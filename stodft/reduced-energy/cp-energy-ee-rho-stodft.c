/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: cp-energy-ee-rho-stodft.c                      */
/*                                                                          */
/* This routine wrapps all functions used within SCF. Nuclei forces are not */
/* calculated.                                                              */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_stodft_local.h"

#include "complex.h"
#define TIME_CP_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* calculate the force on the coeff's (up or down, both) */
/* given an appropriate v_ks (up or down, both).         */
/*==========================================================================*/
void coefForceCalcHybridSCF(CPEWALD *cpewald,int nstate,
                             double *ccreal,double *ccimag,
                             double *fccreal,double  *fccimag,
                             double *cre_scr,double *cim_scr,
                             double *cp_hess_re,double *cp_hess_im,
                             double *zfft,double *zfft_tmp,
                             double *v_ks,double *v_ks_tau,double *ak2_sm,
                             double *eke_ret,double *pvten_cp,
                             int cp_ptens_calc,double *hmati,
                             COMMUNICATE *communicate,
                             int icoef_form,int icoef_orth,int ifcoef_form,
                             int cp_tau_functional,int cp_min_on,
                             PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm,
			     CP *cp, CLASS *class,GENERAL_DATA *general_data)
/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                               */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo = cp->stodftInfo;
  int is,i,iupper;
  int ioff,ncoef1,ioff2;
  int iii,iis,nis;
  int nfft       = cp_sclr_fft_pkg3d_sm->nfft;
  int ncoef      = cp_sclr_fft_pkg3d_sm->ncoef;
  int myid_state = communicate->myid_state;
  int np_states  = communicate->np_states;
  int fftw3dFlag = cpewald->fftw3dFlag;
  //int fftw3dFlag = 1;
  int onebodyMatrixFlag = cpewald->onebodyMatrixFlag;
  int pseudoRealFlag = cp->pseudo.pseudoReal.pseudoRealFlag;

  int  *kastore_sm    =  cpewald->kastr_sm;
  int  *kbstore_sm    =  cpewald->kbstr_sm;
  int  *kcstore_sm    =  cpewald->kcstr_sm;
  MPI_Comm comm_states = communicate->comm_states;

  double tpi;
  double aka,akb,akc,xk,yk,zk,cfact;
  double eke;
  double sum_check,sum_check_tmp;
  double time_st,time_end;
  double *keMatrix = cpewald->keMatrix;
#define DEBUG_OFF
#ifdef DEBUG
  int icount;
  double c_g,g2,anorm,sum,vol,cre_now,cim_now;
  double dx,x_pos,y_pos,z_pos,phase_r,phase_i,arg;
  FILE *fp;
#endif

/* ================================================================= */
/*0) Check the form of the coefficients                              */

  if(icoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in orthogonal form    \n");
    printf("on state processor %d in coef_force_calc_hybrid\n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }//endif

  if(np_states>1)
   if( ifcoef_form==1 || icoef_form == 1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs/coef forces must be in normal (not transposed) \n");
    printf("form on state processor %d in coef_force_calc_hybrid  \n",
           myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }//endif

  cp_sclr_fft_pkg3d_sm->numThreads = communicate->numThreads;

/*=================================================================*/
/*  Find the upper state limit                                     */

  ncoef1 = ncoef - 1;
  iupper = nstate;
  if(nstate%2==1)iupper = nstate-1;

/*=================================================================*/
/*  get the forces on the coefs of each state                      */

  for(is=1;is<=iupper;is+=2 ){
    ioff = (is-1)*ncoef;
    ioff2 = (is)*ncoef;

/*==========================================================================*/
/* 1) get the wave functions in real space two at a time                    */
/*   I) double pack the complex zfft array with two real wavefunctions      */
  
    if(fftw3dFlag==0){
      dble_pack_coef(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
                     zfft,cp_sclr_fft_pkg3d_sm);
    }
    else{
      dble_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
                     zfft,cp_sclr_fft_pkg3d_sm);
    }

/*--------------------------------------------------------------------------*/
/* II) fourier transform the wavefunctions to real space                    */
/*      convention exp(-igr)                                                */

    if(fftw3dFlag==0){
      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
    }
    else{
      para_fft_gen3d_fwd_to_r_fftw3d_threads(zfft,cp_sclr_fft_pkg3d_sm);
    }
    cputime(&time_end);

/*==========================================================================*/
/* 2) get v|psi> in g space and store it in zfft                            */
/*   I) get  v|psi> in real space                                           */
    if(pseudoRealFlag==1){
      memcpy(&zfft_tmp[1],&zfft[1],nfft*sizeof(double));
      cp_vpsi(zfft,v_ks,nfft);
      cp->pseudo.pseudoReal.energyCalcFlag = 1;
      controlEnergyNlppRealThreads(cp,class,general_data,zfft_tmp,zfft,1);
    }
    else cp_vpsi(zfft,v_ks,nfft);
    //printf("v_ks %lg\n",v_ks);

/*--------------------------------------------------------------------------*/
/*  II) fourier transform  to g-space                                       */
/*     convention exp(igr)                                                  */


    if(fftw3dFlag==0){
      para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
    }
    else{
      para_fft_gen3d_bck_to_g_fftw3d_threads(zfft,cp_sclr_fft_pkg3d_sm);
    }

/*==========================================================================*/
/* 3) get forces on coefficients by double unpacking the array zfft         */

    if(fftw3dFlag==0){
      dble_upack_coef_sum(&fccreal[ioff],&fccimag[ioff],
                          &fccreal[ioff2],&fccimag[ioff2],
                          zfft,cp_sclr_fft_pkg3d_sm);
      //printf("fccreal %lg\n",fccreal[ioff+ncoef]);
    }
    else{
      dble_upack_coef_sum_fftw3d(&fccreal[ioff],&fccimag[ioff],
                          &fccreal[ioff2],&fccimag[ioff2],
                          zfft,cp_sclr_fft_pkg3d_sm);
      //printf("fccreal fftw %lg\n",fccreal[ioff+ncoef]);
    }

  }//endfor is

  /*
  if(fftw3dFlag==0){
    for(i=1;i<=ncoef*nstate;i++){
      printf("forceeeee %lg %lg\n",fccreal[i],fccimag[i]);
    }
    exit(0);
  }
  */


/*==========================================================================*/
/*==========================================================================*/
/* 4) if there is an odd number of states, go through                       */
/*      the same procedure using sng_packs                                  */

  if(nstate % 2 != 0){
    is = nstate;
    ioff = (is-1)*ncoef;

/*--------------------------------------------------------------------------*/
/*   I) sngl pack                                                           */


    time_st = omp_get_wtime();
    if(fftw3dFlag==0){
      sngl_pack_coef(&ccreal[ioff],&ccimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);
    }
    else{
      sngl_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);
    }
    //cputime(&time_end);
    time_end = omp_get_wtime();
    stodftInfo->cputime2 += time_end-time_st;

/*--------------------------------------------------------------------------*/
/* II) fourier transform the wavefunctions to real space                    */
/*      convention exp(-igr)                                                */

    if(fftw3dFlag==0){
      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
    }
    else{
      para_fft_gen3d_fwd_to_r_fftw3d_threads(zfft,cp_sclr_fft_pkg3d_sm);
    }

/*==========================================================================*/
/* 5) get v|psi> in g space and store it in zfft                            */
/*   I) get  v|psi> in real space                                           */

    if(pseudoRealFlag==1){
      //cputime(&time_st);
      time_st = omp_get_wtime();
      memcpy(&zfft_tmp[1],&zfft[1],nfft*sizeof(double));
      cp_vpsi(zfft,v_ks,nfft);
      time_end = omp_get_wtime();
      //cputime(&time_end);
      stodftInfo->cputime3 += time_end-time_st;
      cp->pseudo.pseudoReal.energyCalcFlag = 1;
      controlEnergyNlppRealThreads(cp,class,general_data,zfft_tmp,zfft,0);
    }
    else cp_vpsi(zfft,v_ks,nfft);

/*--------------------------------------------------------------------------*/
/*   II) fourier transform the result back to g-space */
/*     convention exp(igr)  */

    if(fftw3dFlag==0){
      para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
    }
    else{
      para_fft_gen3d_bck_to_g_fftw3d_threads(zfft,cp_sclr_fft_pkg3d_sm);
    }

/*==========================================================================*/
/* 6) get forces on coefficients by single unpacking the array zfft         */

    //cputime(&time_st);
    time_st = omp_get_wtime();
    if(fftw3dFlag==0){
      sngl_upack_coef_sum(&fccreal[ioff],&fccimag[ioff],zfft,
			  cp_sclr_fft_pkg3d_sm);
    }
    else{
      sngl_upack_coef_sum_fftw3d(&fccreal[ioff],&fccimag[ioff],zfft,
			  cp_sclr_fft_pkg3d_sm);
    }
    //cputime(&time_end);
    time_end = omp_get_wtime();
    stodftInfo->cputime4 += time_end-time_st;

  }// endif: odd number of states

/*==========================================================================*/
/* 7) If there is an electron KE density dependent functional, calculate    */
/*    this contribution to the force                                        */

/*==========================================================================*/
/* 8) If doing minimization, Fourier transform the KS potential to g-space  */
/*    and unpack it into the diagonal Hessian                               */

/*==========================================================================*/
/* 9) calculate the kinetic energy term and add its contribution to the force*/

  //printf("I'm here kinetic energy!\n");
  //cputime(&time_st);
  time_st = omp_get_wtime();
  tpi = 2.0*M_PI;
  eke = 0.0;
  for(is=1 ; is<= nstate ; is++){
    ioff = (is-1)*ncoef;
    for(i=1; i<= ncoef1 ; i++){
      iis = ioff + i;
      fccreal[iis] -= 2.0*ak2_sm[i]*ccreal[iis];
      fccimag[iis] -= 2.0*ak2_sm[i]*ccimag[iis];
      eke += (2.0*ak2_sm[i]*(ccreal[iis]*ccreal[iis] + ccimag[iis]*ccimag[iis]));
    }/*endfor i*/
   nis = is*ncoef;
   fccimag[nis] = 0.0;
  }/*endfor*/
  //cputime(&time_end);
  time_end = omp_get_wtime();
  stodftInfo->cputime5 += time_end-time_st;


  //debug
  /*
  double sumdebug;
  for(is=0;is<nstate;is++){
    ioff = is*ncoef;
    sumdebug = 0.0;
    for(i=1;i<ncoef;i++){
      sumdebug += fccreal[ioff+i]*fccreal[ioff+i]+fccimag[ioff+i]*fccimag[ioff+i];
    }
    sumdebug *= 2.0;
    sumdebug += fccreal[ioff+ncoef]*fccreal[ioff+ncoef];
    printf("testtttt is %i f norm %lg\n",is,sumdebug);
  }
  printf("fcre %lg\n",fccreal[ncoef]);
  */

  eke *= .50;
  *eke_ret = eke;

/*================================================================================*/
/* 10) If doing minimization, calculat kinetic contribution to diagonal Hessian   */

/*==========================================================================*/
/* 9) calculate kinetic contribution to pressure tensor                     */

/* fine fertig terminado finito */
/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cpGetVksStodft(CPOPTS *cpopts,CPSCR *cpscr,CPEWALD *cpewald,EWALD *ewald,
                COMMUNICATE *communicate,
                CP_COMM_STATE_PKG *cp_comm_state_pkg_up,
                CP_COMM_STATE_PKG *cp_comm_state_pkg_dn,
                STAT_AVG *stat_avg,double *ks_offset,CELL *cell, char *vxc_typ,
                double *pvten_cp,double gc_cut,double alpha_conv_dual,
                int n_interp_pme_dual,
                int cp_gga,int cp_ptens_calc,int nstate_dn,
                int cp_tau_functional,int laplacian_on,
                PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_dens_cp_box,
                PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box,
                PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_lg,
                PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg,
                int cp_dual_grid_opt)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/
#include "../typ_defs/typ_mask.h"

/* assign local pointers */

  double cpu1,cpu2;

  int cp_lda          = cpopts->cp_lda;
  int cp_lsda         = cpopts->cp_lsda;
  int cp_nonint       = cpopts->cp_nonint;
  int cp_para_opt     = cpopts->cp_para_opt;
  int myid_state      = communicate->myid_state;
  int np_states       = communicate->np_states;
  int nstate_up       = cp_comm_state_pkg_up->nstate;
  int cp_lyp          = cpopts->cp_lyp;
  int cp_lypm1        = cpopts->cp_lypm1;
  int nfft            = cp_para_fft_pkg3d_lg->nfft;
  int nfft_proc       = cp_para_fft_pkg3d_lg->nfft_proc;
  int nfft2           = nfft/2;
  int nfft2_proc      = nfft_proc/2;
  int ncoef_l         = cp_para_fft_pkg3d_lg->ncoef_proc;
  int ncoef_l_dens_cp_box = cp_para_fft_pkg3d_dens_cp_box->ncoef_proc;
  int ncoef_l_use     = cp_para_fft_pkg3d_lg->ncoef_use;
  int ncoef_l_use_dens_cp_box = cp_para_fft_pkg3d_dens_cp_box->ncoef_use;
  int icoef_off       = cp_para_fft_pkg3d_lg->icoef_off;
  int iperd           = cell->iperd;

  int *kastore        = ewald->kastr;
  int *kbstore        = ewald->kbstr;
  int *kcstore        = ewald->kcstr;
  int *recv_counts_rho;
  int *displs_rho;

  double *rhocr          =    cpscr->cpscr_rho.rhocr_up;
  double *rhoci          =    cpscr->cpscr_rho.rhoci_up;
  double *rhocr_dens_cp_box =    cpscr->cpscr_rho.rhocr_up_dens_cp_box;
  double *rhoci_dens_cp_box =    cpscr->cpscr_rho.rhoci_up_dens_cp_box;
  double *rho_up         =    cpscr->cpscr_rho.rho_up;
  double *rho_dn         =    cpscr->cpscr_rho.rho_dn;
  double *vextr          =    cpscr->cpscr_loc.vextr;
  double *vexti          =    cpscr->cpscr_loc.vexti;
  double *vextr_dens_cp_box =    cpscr->cpscr_loc.vextr_dens_cp_box;
  double *vexti_dens_cp_box =    cpscr->cpscr_loc.vexti_dens_cp_box;
  double *dvextr         =    cpscr->cpscr_loc.dvextr;
  double *dvexti         =    cpscr->cpscr_loc.dvexti;
  double    *ak2            =    cpewald->ak2;
  double    *ak2_cp_box     =    cpewald->ak2_dens_cp_box;
  double    *v_ks_up        =    cpscr->cpscr_rho.v_ks_up;
  double    *v_ks_dn        =    cpscr->cpscr_rho.v_ks_dn;
  double    *v_ks_tau_up    =    cpscr->cpscr_rho.v_ks_tau_up;
  double    *v_ks_tau_dn    =    cpscr->cpscr_rho.v_ks_tau_dn;
  double    *zfft           =    cpscr->cpscr_wave.zfft;
  double    *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
  double    *bw_r           =    cpscr->cpscr_dual_pme.bw_r;
  double    *bw_i           =    cpscr->cpscr_dual_pme.bw_i;
  double    *clus_corr_r    =    ewald->clus_corr_r;
  double    *hmat_cp        =    cell->hmat_cp;
  double    *hmati          =    cell->hmati;
  double    *hmat           =    cell->hmat;
  double    *eh_ret         =    &(stat_avg->cp_ehart);
  double    *eext_ret       =    &(stat_avg->cp_eext);
  double    *exc_ret        =    &(stat_avg->cp_exc);
  double    *muxc_ret       =    &(stat_avg->cp_muxc);
  MPI_Comm  comm            =    communicate->comm_states;

  int       nfft_dens_cp_box;
  int       nfft_proc_dens_cp_box;

  int       nfft2_proc_dens_cp_box;
  int       nfft2_proc_send;
  int nkf1 = cp_para_fft_pkg3d_lg->nkf1;
  int nkf2 = cp_para_fft_pkg3d_lg->nkf2;
  int nkf3 = cp_para_fft_pkg3d_lg->nkf3;
  int igrid,jgrid;

  /*         Local Variable declarations                                   */
  double ghfact,pi,fpi;
  double aka,akb,akc,xk,yk,zk,tpi;
  double cfact,eh0,vol,vol_cp,eh,eext,exc,muxc;
  double eext_short_dual,eh_short_dual;
  double pre;
  double temp_r,temp_i;
  int i,j,k,iii,igo;

  int fftw3dFlag = cpopts->fftw3dFlag;

/*====================================================================*/
/*  I) calculate hartree potential and add it to the                  */
/*     external potential (sum is over spherically cutoff half space) */

/*--------------------------------------------------------------------*/
/* a) Determine constants and zero energies                           */
  pi = M_PI; fpi = 4.0*pi; tpi = 2.0*pi;
  pre    = (cp_dual_grid_opt == 2 ? 0.25/(alpha_conv_dual*alpha_conv_dual)
                                     : 0.0);

  vol = getdeth(hmat);
  vol_cp = getdeth(hmat_cp);
  eh = 0.0;
  eext = 0.0;
  exc = 0.0;
  muxc = 0.0;

/*--------------------------------------------------------------------*/
/*  c) Get Hartree + external contributions to VKS -- test for CBCs   */
 
  if(iperd==3){
    for(i=1;i<=ncoef_l_use;i++){
      //erf  long range piece on large grid for dual opt
      ghfact = fpi*exp(-ak2[i]*pre)/(ak2[i]*vol);
      eh +=  0.50*ghfact*(rhocr[i]*rhocr[i] + rhoci[i]*rhoci[i]);
      eext +=  vextr[i]*rhocr[i] + vexti[i]*rhoci[i];
      vextr[i] +=  ghfact*rhocr[i];
      vexti[i] +=  ghfact*rhoci[i];
    }//endfor
  }else{
    for(i=1;i<=ncoef_l_use;i++){
      //erf  long range piece on large grid for dual opt
      ghfact = (fpi*exp(-ak2[i]*pre)/ak2[i] + clus_corr_r[i])/vol;
      eh +=  0.50*ghfact*(rhocr[i]*rhocr[i] + rhoci[i]*rhoci[i]);
      eext +=  vextr[i]*rhocr[i] + vexti[i]*rhoci[i];
      vextr[i] += ghfact*rhocr[i];
      vexti[i] += ghfact*rhoci[i];
    }//endfor
  }//endif periodic

  eext *= 2.0;
  if((myid_state+1)==np_states)eext +=  vextr[ncoef_l]*rhocr[ncoef_l];
  eh0 = 2.0*(eh);
  eh  *= 2.0;
/*--------------------------------------------------------------------*/
/*  d) Add in convergent part of long range erf for dualing           */
/*--------------------------------------------------------------------*/
/*  e) Add in cluster correction if necessary */

    if(iperd!=3&&myid_state+1==np_states){
      eh += 0.5*clus_corr_r[ncoef_l]*rhocr[ncoef_l]*rhocr[ncoef_l]/vol;
      vextr[ncoef_l] += clus_corr_r[ncoef_l]*rhocr[ncoef_l]/vol;
    }//endif
/*--------------------------------------------------------------------*/
/*  f) get hartree contribution to pressure tensor                    */
/*--------------------------------------------------------------------*/
/*  i) Increment the energies					      */
  *eext_ret += eext;
  *eh_ret += eh;
/*--------------------------------------------------------------------*/
/* j) Multiply the external potential by the pme g-space weight       */
/*====================================================================*/
/*  II) single pack the ks potential for fourier transform routine    */
  sngl_pack_coef(vextr,vexti,zfft,cp_para_fft_pkg3d_lg);
/*====================================================================*/
/* III) fourier transform ks potential to real space exp(-igr)        */
  para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_lg);
/*====================================================================*/
/* IV) Contract and unpack rho for dual option, otherwise just upack  */
  sngl_upack_rho(zfft,v_ks_up,cp_para_fft_pkg3d_lg);
/*=====================================================================*/
/* V) Coulomb and Hartree erfc calculated on small grid for PME dualing*/
/*====================================================================*/
/* IV) Assign up to down for LSDA because vext is vext        */
  if(cp_lsda==1&&nstate_dn!= 0){
    for(i=1;i<=nfft2_proc;i++){v_ks_dn[i]=v_ks_up[i];}
  }//endif
/*====================================================================*/
/* V) calculate the exchange potential and add to v_ks if necessary   */
/*    on the small grid only for dualing                              */

    igo = 0;
/*--------------------------------------------------------------------*/
/*  a) LDA                                                            */

  if(cp_lda==1){
/*--------------*/
/* Perdew-Zunger*/
/*--------------*/
    if(strcasecmp(vxc_typ,"pz_lda")==0){
      excpot_pz_lda(v_ks_up,rho_up,&exc,&muxc,nfft,nfft_proc,vol_cp,
      	            cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc);
      igo = 1;
    }//endif pz
/*------------------*/
/* Perdew-Wang 1992 */
/*------------------*/
    if(strcasecmp(vxc_typ,"pw_lda")==0 ){
      excpot_pw_lda(v_ks_up,rho_up,&exc,&muxc,nfft,nfft_proc,vol_cp,
		    cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc);
      igo = 1;
    }//endif pw
/*------------------*/
/* Pade */
/*------------------*/
    if(strcasecmp(vxc_typ,"pade_lda")==0 ){
      excpot_pade_lda(v_ks_up,rho_up,&exc,&muxc,nfft,nfft_proc,vol_cp,
		      cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc);
      igo = 1;
    }//endif pade
  }//endif LDA
/*--------------------------------------------------------------------*/
/*  b) LSDA                                                           */

  if(cp_lsda==1){
/*--------------*/
/* Perdew-Zunger*/
/*--------------*/
    if(strcasecmp(vxc_typ,"pz_lsda")==0){
      excpot_pz_lsda(v_ks_up,v_ks_dn,rho_up,rho_dn,&exc,&muxc,nfft,nfft_proc,
                     vol_cp,cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc);
      igo = 1;
    }//endif pz
/*--------------*/
/* Perdew-Wang  */
/*--------------*/
    if(strcasecmp(vxc_typ,"pw_lsda")==0){
      excpot_pw_lsda(v_ks_up,v_ks_dn,rho_up,rho_dn,&exc,&muxc,nfft,nfft_proc,
                     vol_cp,cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc);
      igo = 1;
    }//endif pw
/*---------*/
/* Pade    */
/*---------*/
    if(strcasecmp(vxc_typ,"pade_lsda")==0){
      excpot_pade_lsda(v_ks_up,v_ks_dn,rho_up,rho_dn,&exc,&muxc,nfft,nfft_proc,
                       vol_cp,cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc);
      igo = 1;
    }//endif pade
  }//endif LSDA

/*--------------------------------------------------------------------*/
/*  c) add gradient corrections if necessary                          */

  if(cp_gga==1){
    if(cp_lsda == 0 ){
      grad_corr_lda(cpscr,cpewald,ewald,cell,&exc,&muxc,vol_cp,
                    cpopts,gc_cut,pvten_cp,cp_ptens_calc,
                    cp_tau_functional,laplacian_on,
                    cp_dual_grid_opt,cp_para_fft_pkg3d_lg,nstate_up);
    }else {
      grad_corr_lsda(cpscr,cpewald,ewald,cell,&exc,&muxc,vol_cp,
                     cpopts,gc_cut,pvten_cp,cp_ptens_calc,cp_tau_functional,
                     nstate_up,nstate_dn,laplacian_on,cp_dual_grid_opt,
                     cp_para_fft_pkg3d_lg);
    }// endif cp_lsda
  }//endif cp_gga

/*--------------------------------------------------------------------*/
/*  d Check for implementation of chosen xc functional (barf if not found)*/

  if(igo==0){
    printf("@@@@@@@@@@@@@@@@@@-error-@@@@@@@@@@@@@@@@@@@@@@@ \n");
    printf("user specified correlation functional %s\n",vxc_typ);
    printf("not found in the guts. contact technical support\n");
    printf("@@@@@@@@@@@@@@@@@@-error-@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }//endif
/*--------------------------------------------------------------------*/
/* e) Increment the energies */

  *exc_ret   += exc;
  *muxc_ret  += muxc;
  *ks_offset  = exc-eh-muxc;

/*====================================================================*/
/* Get the KS matrix in form for state level                          */
/*  that is for hybrid parallel, allgather the puppy                  */

  if(np_states>1&&cp_para_opt==0){// hybrid parallel only
    nfft2_proc_send = nfft2_proc;
    recv_counts_rho = cp_para_fft_pkg3d_lg->recv_counts_rho;
    displs_rho      = cp_para_fft_pkg3d_lg->displs_rho;

    for(i=1;i<=nfft2_proc_send;i++)zfft[i] = v_ks_up[i];

    Allgatherv(&(zfft[1]),nfft2_proc_send,MPI_DOUBLE,&(v_ks_up[1]),
                 &recv_counts_rho[1],&(displs_rho[1]),MPI_DOUBLE,0,comm);
    // All the above rho and v_ks using z as leading dimension. Now I 
    // tanspose it to x leading.
    if(fftw3dFlag==1){
      for(i=0;i<nkf3;i++){
        for(j=0;j<nkf2;j++){
          for(k=0;k<nkf1;k++){
            igrid = i*nkf2*nkf1+j*nkf1+k;
            jgrid = k*nkf2*nkf3+j*nkf3+i;
            zfft[jgrid+1] = v_ks_up[igrid+1];
          }//endfor k
        }//endfor j
      }//endfor i
      memcpy(&v_ks_up[1],&(zfft[1]),nfft2*sizeof(double));
    }//endif fftwedFlag
    if(cp_lsda==1){
      for(i=1;i<=nfft2_proc_send;i++)zfft[i]=v_ks_dn[i];
      Allgatherv(&(zfft[1]),nfft2_proc_send,MPI_DOUBLE,&(v_ks_dn[1]),
                 &recv_counts_rho[1],&(displs_rho[1]),MPI_DOUBLE,0,comm);
      if(fftw3dFlag==1){
        for(i=0;i<nkf3;i++){
          for(j=0;j<nkf2;j++){
            for(k=0;k<nkf1;k++){
              igrid = i*nkf2*nkf1+j*nkf1+k;
              jgrid = k*nkf2*nkf3+j*nkf3+i;
              zfft[jgrid+1] = v_ks_dn[igrid+1];
            }//endfor k
          }//endfor j
        }//endfor i
        memcpy(&v_ks_dn[1],&(zfft[1]),nfft2*sizeof(double));
      }//endif fftwedFlag
    }//endif cp_lsda
  }//endif hybrid

/*==========================================================================*/
   }/* end routine*/
/*==========================================================================*/



