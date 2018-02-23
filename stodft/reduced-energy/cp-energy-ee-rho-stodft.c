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

  int is,i,iupper;
  int ioff,ncoef1,ioff2;
  int iii,iis,nis;
  int nfft       = cp_sclr_fft_pkg3d_sm->nfft;
  int ncoef      = cp_sclr_fft_pkg3d_sm->ncoef;
  int myid_state = communicate->myid_state;
  int np_states  = communicate->np_states;
  int fftw3dFlag = cpewald->fftw3dFlag;
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
      para_fft_gen3d_fwd_to_r_fftw3d(zfft,cp_sclr_fft_pkg3d_sm);
    }

/*==========================================================================*/
/* 2) get v|psi> in g space and store it in zfft                            */
/*   I) get  v|psi> in real space                                           */
    if(pseudoRealFlag==1){
      memcpy(&zfft_tmp[1],&zfft[1],nfft*sizeof(double));
      cp_vpsi(zfft,v_ks,nfft);
      cp->pseudo.pseudoReal.energyCalcFlag = 1;
      controlEnergyNlppReal(cp,class,general_data,zfft_tmp,zfft,1);
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
      para_fft_gen3d_bck_to_g_fftw3d(zfft,cp_sclr_fft_pkg3d_sm);
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

    if(fftw3dFlag==0){
      sngl_pack_coef(&ccreal[ioff],&ccimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);
    }
    else{
      sngl_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);
    }

/*--------------------------------------------------------------------------*/
/* II) fourier transform the wavefunctions to real space                    */
/*      convention exp(-igr)                                                */

    if(fftw3dFlag==0){
      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
    }
    else{
      para_fft_gen3d_fwd_to_r_fftw3d(zfft,cp_sclr_fft_pkg3d_sm);
    }

/*==========================================================================*/
/* 5) get v|psi> in g space and store it in zfft                            */
/*   I) get  v|psi> in real space                                           */

    if(pseudoRealFlag==1){
      memcpy(&zfft_tmp[1],&zfft[1],nfft*sizeof(double));
      cp_vpsi(zfft,v_ks,nfft);
      cp->pseudo.pseudoReal.energyCalcFlag = 1;
      controlEnergyNlppReal(cp,class,general_data,zfft_tmp,zfft,0);
    }
    else cp_vpsi(zfft,v_ks,nfft);

/*--------------------------------------------------------------------------*/
/*   II) fourier transform the result back to g-space */
/*     convention exp(igr)  */
    if(fftw3dFlag==0){
      para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
    }
    else{
      para_fft_gen3d_bck_to_g_fftw3d(zfft,cp_sclr_fft_pkg3d_sm);
    }

/*==========================================================================*/
/* 6) get forces on coefficients by single unpacking the array zfft         */

    if(fftw3dFlag==0){
      sngl_upack_coef_sum(&fccreal[ioff],&fccimag[ioff],zfft,
			  cp_sclr_fft_pkg3d_sm);
    }
    else{
      sngl_upack_coef_sum_fftw3d(&fccreal[ioff],&fccimag[ioff],zfft,
			  cp_sclr_fft_pkg3d_sm);
    }

  }// endif: odd number of states

/*==========================================================================*/
/* 7) If there is an electron KE density dependent functional, calculate    */
/*    this contribution to the force                                        */

  if(cp_tau_functional==1){
    coef_force_tau_fun_hybrid(cpewald,nstate,ccreal,ccimag,fccreal,fccimag,
                              cre_scr,cim_scr,zfft,zfft_tmp,v_ks_tau,ak2_sm,
                              pvten_cp,cp_ptens_calc,hmati,communicate,
                              icoef_form,icoef_orth,ifcoef_form,
                              cp_sclr_fft_pkg3d_sm);
  }

/*==========================================================================*/
/* 8) If doing minimization, Fourier transform the KS potential to g-space  */
/*    and unpack it into the diagonal Hessian                               */

/*==========================================================================*/
/* 9) calculate the kinetic energy term and add its contribution to the force*/

  //printf("I'm here kinetic energy!\n");
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




