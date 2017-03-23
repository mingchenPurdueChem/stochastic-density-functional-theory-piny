
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void cp_get_vks(CPOPTS *cpopts,CPSCR *cpscr,CPEWALD *cpewald,EWALD *ewald,
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

   double cpu1,cpu2;
   int i,j,iii,igo;
   int cp_lda          =    cpopts->cp_lda;
   int cp_lsda         =    cpopts->cp_lsda;
   int cp_nonint       =    cpopts->cp_nonint;
   int cp_para_opt     =    cpopts->cp_para_opt;
   int myid_state      =    communicate->myid_state;
   int np_states       =    communicate->np_states;
   int nstate_up       =    cp_comm_state_pkg_up->nstate;
   int cp_lyp          =    cpopts->cp_lyp;
   int cp_lypm1        =    cpopts->cp_lypm1;
   int nfft            =    cp_para_fft_pkg3d_lg->nfft;
   int nfft_proc       =    cp_para_fft_pkg3d_lg->nfft_proc;
   int nfft2           =    nfft/2;
   int nfft2_proc      =    nfft_proc/2;
   int ncoef_l         =    cp_para_fft_pkg3d_lg->ncoef_proc;
   int ncoef_l_dens_cp_box = cp_para_fft_pkg3d_dens_cp_box->ncoef_proc;
   int ncoef_l_use	   = cp_para_fft_pkg3d_lg->ncoef_use;
   int ncoef_l_use_dens_cp_box = cp_para_fft_pkg3d_dens_cp_box->ncoef_use;
   int icoef_off       = cp_para_fft_pkg3d_lg->icoef_off;
   int iperd           = cell->iperd;
   int nfft_dens_cp_box;
   int nfft_proc_dens_cp_box;
   int nfft2_proc_dens_cp_box;
   int nfft2_proc_send;
   int fftw3dFlag = cpopts->fftw3dFlag;

   int *kastore        =    ewald->kastr;
   int *kbstore        =    ewald->kbstr;
   int *kcstore        =    ewald->kcstr;
   int *recv_counts_rho;
   int *displs_rho;

   double ghfact,pi,fpi;
   double aka,akb,akc,xk,yk,zk,tpi;
   double cfact,eh0,vol,vol_cp,eh,eext,exc,muxc;
   double eext_short_dual,eh_short_dual;
   double pre;
   double temp_r,temp_i;


   double *rhocr          =    cpscr->cpscr_rho.rhocr_up;
   double *rhoci          =    cpscr->cpscr_rho.rhoci_up;
   double *rhocr_dens_cp_box =    cpscr->cpscr_rho.rhocr_up_dens_cp_box;
   double *rhoci_dens_cp_box =    cpscr->cpscr_rho.rhoci_up_dens_cp_box;
   double *rho_up         =    cpscr->cpscr_rho.rho_up;
   double *rho_dn         =    cpscr->cpscr_rho.rho_dn;

   double *vextr          =    cpscr->cpscr_loc.vextr;
   double *vexti          =    cpscr->cpscr_loc.vexti;
   double *vextr_loc	  =    cpscr->cpscr_loc.vextr_loc;
   double *vexti_loc      =    cpscr->cpscr_loc.vexti_loc;
   double *vextr_dens_cp_box =    cpscr->cpscr_loc.vextr_dens_cp_box;
   double *vexti_dens_cp_box =    cpscr->cpscr_loc.vexti_dens_cp_box;
   double *vextr_dens_cp_box_loc = cpscr->cpscr_loc.vextr_dens_cp_box_loc;
   double *vexti_dens_cp_box_loc = cpscr->cpscr_loc.vexti_dens_cp_box_loc;
   double *dvextr         =    cpscr->cpscr_loc.dvextr;
   double *dvexti         =    cpscr->cpscr_loc.dvexti;
   double *ak2            =    cpewald->ak2;
   double *ak2_cp_box     =    cpewald->ak2_dens_cp_box;
   double *v_ks_up        =    cpscr->cpscr_rho.v_ks_up;
   double *v_ks_dn        =    cpscr->cpscr_rho.v_ks_dn;
   double *v_ks_tau_up    =    cpscr->cpscr_rho.v_ks_tau_up;
   double *v_ks_tau_dn    =    cpscr->cpscr_rho.v_ks_tau_dn;
   double *zfft           =    cpscr->cpscr_wave.zfft;
   double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
   double *bw_r           =    cpscr->cpscr_dual_pme.bw_r;
   double *bw_i           =    cpscr->cpscr_dual_pme.bw_i;
   double *clus_corr_r    =    ewald->clus_corr_r;

   double *hmat_cp        =    cell->hmat_cp;
   double *hmati          =    cell->hmati;
   double *hmat           =    cell->hmat;
   double *eh_ret         =    &(stat_avg->cp_ehart);
   double *eext_ret       =    &(stat_avg->cp_eext);
   double *exc_ret        =    &(stat_avg->cp_exc);
   double *muxc_ret       =    &(stat_avg->cp_muxc);
   MPI_Comm  comm            =    communicate->comm_states;

   // fftw3d options

   if( cp_dual_grid_opt >= 1){
    nfft_dens_cp_box       = cp_para_fft_pkg3d_dens_cp_box->nfft;
    nfft_proc_dens_cp_box  = cp_para_fft_pkg3d_dens_cp_box->nfft_proc;
    nfft2_proc_dens_cp_box = nfft_proc_dens_cp_box/2;
   }/*endif cp_dual_grid_opt */

/*====================================================================*/
/*  0) Copy the vext back			                      */

  memcpy(&vextr[1],&(vextr_loc[1]),ncoef_l*sizeof(double));
  memcpy(&vexti[1],&(vexti_loc[1]),ncoef_l*sizeof(double));
  if()

/*====================================================================*/
/*  I) calculate hartree potential and add it to the                  */
/*     external potential (sum is over spherically cutoff half space) */

/*--------------------------------------------------------------------*/
/* a) Determine constants and zero energies                           */
   pi     = M_PI; fpi = 4.0*pi; tpi = 2.0*pi;
   pre    = (cp_dual_grid_opt == 2 ? 0.25/(alpha_conv_dual*alpha_conv_dual)
                                     : 0.0);

   vol    = getdeth(hmat);
   vol_cp = getdeth(hmat_cp);

   eh   = 0.0;
   eext = 0.0;
   exc  = 0.0;
   muxc = 0.0;


  if(cp_nonint==0){

/*--------------------------------------------------------------------*/
/*  c) Get Hartree + external contributions to VKS -- test for CBCs   */

    if(iperd == 3 ){
      for(i=1 ; i<= ncoef_l_use; i++){
                  /*erf  long range piece on large grid for dual opt */
        ghfact    = fpi*exp(-ak2[i]*pre)/(ak2[i]*vol);
        eh       +=  0.50*ghfact*(rhocr[i]*rhocr[i] + rhoci[i]*rhoci[i]);
        eext     +=  vextr[i]*rhocr[i] + vexti[i]*rhoci[i];
        vextr[i] +=  ghfact*rhocr[i];
        vexti[i] +=  ghfact*rhoci[i];
      }/*endfor*/
    }else{

      for(i=1 ; i<= ncoef_l_use; i++){
                  /*erf  long range piece on large grid for dual opt */
        ghfact    = (fpi*exp(-ak2[i]*pre)/ak2[i] + clus_corr_r[i])/vol;
        eh       +=  0.50*ghfact*(rhocr[i]*rhocr[i] + rhoci[i]*rhoci[i]);
        eext     +=  vextr[i]*rhocr[i] + vexti[i]*rhoci[i];

        vextr[i] +=  ghfact*rhocr[i];
        vexti[i] +=  ghfact*rhoci[i];
      }/*endfor*/
    }/*endif periodic */

    eext *= 2.0;

    if((myid_state+1) == np_states){eext +=  vextr[ncoef_l]*rhocr[ncoef_l];}

    eh0 = 2.0*(eh);
    eh  *= 2.0;


/*--------------------------------------------------------------------*/
/*  d) Add in convergent part of long range erf for dualing           */

    if( (cp_dual_grid_opt == 2) && ((myid_state+1) == np_states) ){
      eh -= 0.5*fpi*pre*rhocr[ncoef_l]*rhocr[ncoef_l]/vol;
      vextr[ncoef_l] -= fpi*pre*rhocr[ncoef_l]/vol;
    }/*endif*/

/*--------------------------------------------------------------------*/
/*  e) Add in cluster correction if necessary */

    if( (iperd != 3 )&& ((myid_state+1) == np_states) ) {
      eh += 0.5*clus_corr_r[ncoef_l]*rhocr[ncoef_l]*rhocr[ncoef_l]/vol;
      vextr[ncoef_l] += clus_corr_r[ncoef_l]*rhocr[ncoef_l]/vol;
    }/*endif*/

/*--------------------------------------------------------------------*/
/*  g) Only get external contribution to VKS                          */

    for(i=1; i<= ncoef_l_use; i++){
      eext += (vextr[i]*rhocr[i] + vexti[i]*rhoci[i]);
    }/*endfor*/

    eext *= 2.0;

    if((myid_state+1) == np_states) eext += vextr[ncoef_l]*rhocr[ncoef_l];

  }/*endif nonint*/






/*======================================================================*/
    }/*end routine*/
/*======================================================================*/

