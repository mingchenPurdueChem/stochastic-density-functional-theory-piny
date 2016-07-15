/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: density.c                                      */
/*                                                                          */
/* This routine costruct P_N(H)|phi> where P_N is some                      */
/* polynomial.                                                              */
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
/*  rhoCalcRealSto:                                                         */
/*      calculates the real space density (up,down or total)                */
/*      using the stochastic orbitals. Double packing of ffts and other     */
/*      goodies (most studdly coding)                                       */
/*      Even the option to do full g space (cp_rho_calc_full_g)             */
/*      or the studdly hybrid option in (cp_rho_calc_hybrid)                */
/*==========================================================================*/

void rhoCalcRealStoHybrid(CPSCR *cpscr,
                        CPCOEFFS_INFO *cpcoeffs_info,
                        CELL *cell,STODFTINFO *stodftInfo,
			double *creal, double *cimag,double *rho_scr,
                        int icoef_form,int icoef_orth,int nstate,int ncoef,
                        int cp_dual_grid_opt,
                        COMMUNICATE *communicate,
                        PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg,
                        PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_lg,
                        PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box,
                        PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_dens_cp_box,
                        PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"


  int i,ioff,ioff2,is,iupper;
  int myid_state       =    communicate->myid_state;
  int np_states        =    communicate->np_states;
  int nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
  int nfft             =    cp_para_fft_pkg3d_lg->nfft;
  int nfft2            =    nfft/2;
  int nfft2_proc       =    nfft_proc/2;
  int nfft_dens_cp_box,nfft2_dens_cp_box;
  int nfft_proc_dens_cp_box,nfft2_proc_dens_cp_box;
  int ncoef_l_dens_cp_box;

  double vol_cp,rvol_cp;
  double temp_r,temp_i;
  double *zfft           =    cpscr->cpscr_wave.zfft;
  double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
  double *hmati_cp       =    cell->hmati_cp;
  double *hmat_cp        =    cell->hmat_cp;

  MPI_Comm comm_states   =    communicate->comm_states;

  if(cp_dual_grid_opt >= 1){
    nfft_dens_cp_box       =  cp_para_fft_pkg3d_dens_cp_box->nfft;
    nfft2_dens_cp_box      =  nfft_dens_cp_box/2;
    ncoef_l_dens_cp_box    =  cp_para_fft_pkg3d_dens_cp_box->ncoef;
  }/*endif cp_dual_grid_opt*/

  /*
  if(np_states > 1){
    rho_scr = cpscr->cpscr_rho.v_ks_up; 
  }else{
    rho_scr = rho;
  }
  */
/* ================================================================= */
/*0) Check the form of the coefficients                              */

  if(icoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in orthogonal form    \n");
    printf("on state processor %d in cp_rho_calc_hybrid   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1 && icoef_form == 1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in normal (not transposed) \n");
    printf("form on state processor %d in cp_rho_calc_hybrid  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/


/* ================================================================= */
/*1) zero density and gradients if necessary                         */

  if(cp_dual_grid_opt >= 1){
   for(i=1; i<= nfft2_dens_cp_box ; i++){
     rho_scr[i] = 0.0;
   }/*endfor*/
  }else{
   for(i=1; i<= nfft2 ; i++){
     rho_scr[i] = 0.0;
   }/*endfor*/
  }/*endif cp_dual_grid_opt*/

/*==========================================================================*/
/*==========================================================================*/
/*2)sum the density in real space two states at a time                      */
/*  This is done at state level and uses the scalar packages!               */

  iupper = nstate;
  if((nstate % 2) != 0){ 
    iupper = nstate-1; 
  }/* endif */

  for(is = 1; is <= iupper; is = is + 2){
    ioff   = (is-1)*ncoef;
    ioff2 = (is)*ncoef;
/*--------------------------------------------------------------------------*/
/*I) pack the complex zfft array with two wavefunctions (real) 
    the wavefunctions are reperesented in spherically cuttof 
    half g space                                                            */

    dble_pack_coef(&creal[ioff],&cimag[ioff],&creal[ioff2],&cimag[ioff2],
                   zfft,cp_sclr_fft_pkg3d_sm);
/*--------------------------------------------------------------------------*/
/*II) fourier transform the two wavefunctions to real space
     convention exp(-igr)                                                   */

    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
  
/*--------------------------------------------------------------------------*/
/* III) add the square of the two wave functions to the density(real space) */

    if(cp_dual_grid_opt >= 1){
      sum_rho(zfft,rho_scr,cp_sclr_fft_pkg3d_dens_cp_box); 
    }else{
      sum_rho(zfft,rho_scr,cp_sclr_fft_pkg3d_lg);
    }/*endif cp_dual_grid_opt*/

  }/*endfor is*/

/*--------------------------------------------------------------------------*/
/*IV) if there is an odd number of states, use single pack
      for the last state (i caen't be chaengin' the laws of physics
      captn, i've got to be usin the single pack!!!!!!)                     */

  if((nstate%2)!=0){
    ioff = (nstate-1)*ncoef;
    sngl_pack_coef(&creal[ioff],&cimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* V) fourier transform the last wavefunction to real space
       convention exp(-igr)                                                 */

    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*VI) add the square of the last wave function to the density(real space)   */

    if(cp_dual_grid_opt >= 1){ 
      sum_rho(zfft,rho_scr,cp_sclr_fft_pkg3d_dens_cp_box); 
    }else{
      sum_rho(zfft,rho_scr,cp_sclr_fft_pkg3d_lg);
    }/*endif cp_dual_grid_opt*/
  }/*endif*/

/*==============================================================*/
}/*end routine*/
/*==============================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcRhoStoRecipHybrid(CPEWALD *cpewald,CPSCR *cpscr,
                        CPCOEFFS_INFO *cpcoeffs_info,EWALD *ewald,
                        CELL *cell,STODFTINFO *stodftInfo,
                        double *creal, double *cimag,
                        int icoef_form,int icoef_orth,
                        double *rhocr ,double *rhoci,double *rhotemp,double *rho,
                        double *rhocr_dens_cp_box,double *rhoci_dens_cp_box,
                        double *del_rho_x, double *del_rho_y,
                        double *del_rho_z,
                        double *del2_rho,int nstate,int ncoef,int nstate_tot,
                        int cp_gga,int cp_dual_grid_opt,
                        int n_interp_pme_dual,
                        COMMUNICATE *communicate,
                        PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg,
                        PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_lg,
                        PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box,
                        PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_dens_cp_box,
                        PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Now we have real space average density and I want to transfer it to   */
/* reciprocal space and also generate gradient for gga                   */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  #include "../typ_defs/typ_mask.h"

  int iii,ioff,ioff2;
  int is,i,iupper;
  double vol_cp,rvol_cp;
  double temp_r,temp_i;

  /*  Assign local pointers                                           */
  double *zfft           =    cpscr->cpscr_wave.zfft;
  double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
  double *rho_scr;

  double *hmati_cp       =    cell->hmati_cp;
  double *hmat_cp        =    cell->hmat_cp;

  double dbox_rat        =    cpewald->dbox_rat;
  double *bw_r           =    cpscr->cpscr_dual_pme.bw_r;
  double *bw_i           =    cpscr->cpscr_dual_pme.bw_i;

  int cp_elf_calc_frq    =    cpcoeffs_info->cp_elf_calc_frq;

  int   myid_state       =    communicate->myid_state;
  int   np_states        =    communicate->np_states;
  int   laplacian_on     =    cpcoeffs_info->cp_laplacian_on;

  int  *recv_counts_coef =    cp_para_fft_pkg3d_lg->recv_counts_coef;
  int   ncoef_l          =    cp_para_fft_pkg3d_lg->ncoef;
  int   ncoef_l_use      =    cp_para_fft_pkg3d_lg->ncoef_use;
  int   ncoef_l_proc     =    cp_para_fft_pkg3d_lg->ncoef_proc;
  int   nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
  int   nfft             =    cp_para_fft_pkg3d_lg->nfft;
  int   nfft2            =    nfft/2;
  int   nfft2_proc       =    nfft_proc/2;

  int   nfft_dens_cp_box,nfft2_dens_cp_box;
  int   nfft_proc_dens_cp_box,nfft2_proc_dens_cp_box;
  int   ncoef_l_dens_cp_box;

  double integral,int_tmp;
  int    *recv_counts_coef_dens_cp_box;
  int *recvDispls = stodftInfo->recvDispls;
  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls  = stodftInfo->rhoRealDispls;
  int *recvDisplsDens = stodftInfo->recvDisplsDens;

 MPI_Comm comm_states   =    communicate->comm_states;

 if(cp_dual_grid_opt >= 1){
   nfft_dens_cp_box       =  cp_para_fft_pkg3d_dens_cp_box->nfft;
   nfft_proc_dens_cp_box  =  cp_para_fft_pkg3d_dens_cp_box->nfft_proc;
   nfft2_dens_cp_box      =  nfft_dens_cp_box/2;
   nfft2_proc_dens_cp_box =  nfft_proc_dens_cp_box/2;


   ncoef_l_dens_cp_box  =  cp_para_fft_pkg3d_dens_cp_box->ncoef;

   recv_counts_coef_dens_cp_box =
             cp_para_fft_pkg3d_dens_cp_box->recv_counts_coef;

 }/*endif cp_dual_grid_opt*/


/*=========================================================================*/
/*=========================================================================*/
/*  3) get density in g space (on the master process                       */
/*  I)  pack it up                                                         */
/*  rho_scr = rho in scalar                                                */

  if(myid_state==0){
    if(cp_dual_grid_opt >= 1){
      sngl_pack_rho(zfft,rhotemp,cp_sclr_fft_pkg3d_dens_cp_box);
    }else{
      sngl_pack_rho(zfft,rhotemp,cp_sclr_fft_pkg3d_lg);
    }/*endif cp_dual_grid_opt*/


/*--------------------------------------------------------------------------*/
/*  II) back transform to g-space  convention exp(igr)                      */

    if(cp_dual_grid_opt >= 1){
      para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_dens_cp_box);
    }else{
      para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_lg);
    }/*endif cp_dual_grid_opt*/

/*==========================================================================*/
/*  III) unpack the density                                                 */

    if(cp_dual_grid_opt == 0){
      if(np_states == 1){
        sngl_upack_coef(rhocr,rhoci,zfft,cp_sclr_fft_pkg3d_lg);
      }else{
        sngl_upack_coef(zfft_tmp,&zfft_tmp[ncoef_l],zfft,cp_sclr_fft_pkg3d_lg);
    }/*endif*/
    }else{
      if(np_states == 1){
        sngl_upack_coef(rhocr_dens_cp_box,rhoci_dens_cp_box,zfft,
	  	        cp_sclr_fft_pkg3d_dens_cp_box);
      }else{
        sngl_upack_coef(zfft_tmp,&zfft_tmp[ncoef_l_dens_cp_box],zfft,
		        cp_sclr_fft_pkg3d_dens_cp_box);
      }/*endif*/
    }/*endif cp_dual_grid_opt*/
  }//endif

/*==========================================================================*/
/* VII) Scatter rho in g space and real space		                    */

  if(np_states > 1){
    if(cp_dual_grid_opt == 0){
      //Since we already reduce the density in real space, we only need to scatter
      //We need to generate recvDispls 
      Scatterv(&zfft_tmp[1],&recv_counts_coef[1],recvDispls,MPI_DOUBLE,
	       &rhocr[1],recv_counts_coef[myid_state+1],MPI_DOUBLE,0,comm_states);
      Scatterv(&zfft_tmp[ncoef_l+1],&recv_counts_coef[1],recvDispls,MPI_DOUBLE,
	       &rhoci[1],recv_counts_coef[myid_state+1],MPI_DOUBLE,0,comm_states);

      /*
      sngl_pack_coef(rhocr,rhoci,zfft,cp_para_fft_pkg3d_lg);

      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_lg);

      sngl_upack_rho(zfft,rho,cp_para_fft_pkg3d_lg);

      Barrier(comm_states);

      */
      Scatterv(&(rhotemp[1]),rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
               &rho[1],nfft2_proc,MPI_DOUBLE,0,comm_states);

    }else{
      Scatterv(&zfft_tmp[1],&recv_counts_coef_dens_cp_box[1],
	       recvDisplsDens,MPI_DOUBLE,&rhocr[1],
	       recv_counts_coef_dens_cp_box[myid_state+1],
	       MPI_DOUBLE,0,comm_states);
      Scatterv(&zfft_tmp[ncoef_l_dens_cp_box+1],&recv_counts_coef_dens_cp_box[1],
	       recvDisplsDens,MPI_DOUBLE,&rhoci[1],
	       recv_counts_coef_dens_cp_box[myid_state+1],
	       MPI_DOUBLE,0,comm_states);
      /*
      sngl_pack_coef(rhocr_dens_cp_box,rhoci_dens_cp_box,zfft,
                     cp_para_fft_pkg3d_dens_cp_box);

      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_dens_cp_box);

      sngl_upack_rho(zfft,rho,cp_para_fft_pkg3d_dens_cp_box);
      */
      Scatterv(&(rhotemp[1]),rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
               &rho[1],nfft2_proc_dens_cp_box,MPI_DOUBLE,0,comm_states);
    }/*endif cp_dual_grid_opt*/
  }/*endif np_states*/

  /*
  if(np_states>1){
    printf("myid %i rho[1] %lg\n",myid_state,rho[1]);
    if(myid_state==0){
      for(i=1;i<=np_states;i++)printf("i %i recv_counts_coeff %i\n",i,recv_counts_coef[i]);
    }
  }
  */

/*===========================================================================*/
/* IF DUALED put rho real space onto the large grid and fft it to g space    */

 if(cp_dual_grid_opt >= 1){
/* sending density*vol_cp on small grid  */
   control_spread_rho(cpscr,rho,cell,dbox_rat,np_states,
                      n_interp_pme_dual,
                      cp_para_fft_pkg3d_dens_cp_box,
                      cp_para_fft_pkg3d_lg,cp_dual_grid_opt);

  if(np_states == 1){
    para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_lg);
    sngl_upack_coef(rhocr,rhoci,zfft,cp_sclr_fft_pkg3d_lg);
  }else{
    para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_lg);
    sngl_upack_coef(rhocr,rhoci,zfft,cp_para_fft_pkg3d_lg);
  }/*endif np_states*/
 }/*endif cp_dual_grid_opt*/


/*--------------------------------------------------------------------*/
/*  Post-processing for pme grid need to multiply by complex weight factor */

  if((n_interp_pme_dual > 1) && (cp_dual_grid_opt == 2)){
    for(i=1; i<= ncoef_l_use; i++){
     temp_r   =  (rhocr[i]*bw_r[i] - rhoci[i]*bw_i[i]);
     temp_i   =  (rhocr[i]*bw_i[i] + rhoci[i]*bw_r[i]);
     rhocr[i] =  temp_r;
     rhoci[i] =  temp_i;
    }/*endfor*/

    if((myid_state+1) == np_states){rhocr[ncoef_l_proc]*=bw_r[ncoef_l_proc];}

   }/*endif pme grid */

/*===========================================================================*/
/* IV) finish the density in real space by dividing by the volume           */
/* DUALED SYSTEMS only keep the real space density on the cp_grid            */

     vol_cp  = getdeth(hmat_cp);
     rvol_cp = 1.0/vol_cp;

     if(cp_dual_grid_opt >= 1){
      for(i=1 ; i<= nfft2_proc_dens_cp_box;i++){
         rho[i] *= rvol_cp;
      }/*endfor*/
     }else{
      for(i=1 ; i<= nfft2_proc;i++){
         rho[i] *= rvol_cp;
      }/*endfor*/
     }/*endif cp_dual_grid_opt*/

/*==============================================================*/
/* VII) if doing gradient corrections, get gradient of density    */

  if((cp_gga == 1 || cp_elf_calc_frq > 0)) {
   if(cp_dual_grid_opt >= 1){
    control_grad_rho(cpewald,cpscr,ewald,rhocr_dens_cp_box,rhoci_dens_cp_box,
                     del_rho_x,del_rho_y,del_rho_z,del2_rho,
                     hmati_cp,vol_cp,laplacian_on,cp_dual_grid_opt,
                     cp_para_fft_pkg3d_dens_cp_box);
   }else{
    control_grad_rho(cpewald,cpscr,ewald,rhocr,rhoci,
                     del_rho_x,del_rho_y,del_rho_z,del2_rho,
                     hmati_cp,vol_cp,laplacian_on,cp_dual_grid_opt,
                     cp_para_fft_pkg3d_lg);
   }/*endif cp_dual_grid_opt*/
  }/*endif*/

/*==============================================================*/
}/*end routine*/
/*==============================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  rhoCalcRealSto:                                                         */
/*      calculates the real space density (up,down or total)                */
/*      using the stochastic orbitals. Double packing of ffts and other     */
/*      goodies (most studdly coding)                                       */
/*      Even the option to do full g space (cp_rho_calc_full_g)             */
/*      or the studdly hybrid option in (cp_rho_calc_hybrid)                */
/*==========================================================================*/

void rhoCalcRealStoFullg(CPSCR *cpscr,
                        CPCOEFFS_INFO *cpcoeffs_info,
                        CELL *cell,STODFTINFO *stodftInfo,
                        double *creal, double *cimag,double *rho_scr,
                        int icoef_form,int icoef_orth,
                        int nstate,int ncoef,
                        int cp_dual_grid_opt,
                        COMMUNICATE *communicate,
                        PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg,
                        PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box,
                        PARA_FFT_PKG3D *cp_para_fft_pkg3d_sm)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
  #include "../typ_defs/typ_mask.h"

  /* local variables                                                  */

  int ioff,ioff2,i,is,iupper;

  double *zfft           =    cpscr->cpscr_wave.zfft;
  double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;

  int   myid_state       =    communicate->myid_state;
  int   np_states        =    communicate->np_states;

  int   ncoef_l          =    cp_para_fft_pkg3d_lg->ncoef;
  int   ncoef_l_use      =    cp_para_fft_pkg3d_lg->ncoef_use;
  int   ncoef_l_proc     =    cp_para_fft_pkg3d_lg->ncoef_proc;
  int   nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
  int   nfft             =    cp_para_fft_pkg3d_lg->nfft;

  int   nfft2_proc       =    nfft_proc/2;
  int   nfft2            =    nfft/2;

  int   nfft_proc_dens_cp_box,nfft2_proc_dens_cp_box;

  double integral,int_tmp;

  MPI_Comm comm_states   =    communicate->comm_states;

  if(cp_dual_grid_opt >= 1){
   nfft_proc_dens_cp_box   =  cp_para_fft_pkg3d_dens_cp_box->nfft_proc;
   nfft2_proc_dens_cp_box  =  nfft_proc_dens_cp_box/2;
  }/*endif cp_dual_grid_opt*/


/* ================================================================= */
/*0) Check the form of the coefficients                              */

  if(icoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in orthogonal form    \n");
    printf("on state processor %d in cp_rho_calc_full_g   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1 && icoef_form == 0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed (not normal) \n");
    printf("form on state processor %d in cp_rho_calc_full_g  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/* ================================================================= */
/*1) zero density and gradients if necessary                         */

  if(cp_dual_grid_opt>=1){
    for(i=1;i<=nfft2_proc_dens_cp_box;i++)rho_scr[i] = 0.0;
  }else{
    for(i=1;i<=nfft2_proc;i++)rho_scr[i] = 0.0;
  }/*endif cp_dual_grid_opt*/


/*==========================================================================*/
/*==========================================================================*/
/*2)sum the density in real space two states at a time                      */
/*  This is done at state level and uses the scalar packages!               */

  iupper = nstate;
  if((nstate % 2) != 0){ 
     iupper = nstate-1; 
  }/* endif */

  for(is = 1; is <= iupper; is = is + 2){
    ioff   = (is-1)*ncoef;
    ioff2 = (is)*ncoef;

/*--------------------------------------------------------------------------*/
/*I) pack the complex zfft array with two wavefunctions (real) 
  the wavefunctions are reperesented in spherically cuttof 
  half g space                       */

    dble_pack_coef(&creal[ioff],&cimag[ioff],&creal[ioff2],&cimag[ioff2],
  		   zfft,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*II) fourier transform the two wavefunctions to real space
   convention exp(-igr)  */

    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* III) add the square of the two wave functions to the density  */
    if(cp_dual_grid_opt >= 1){
      sum_rho(zfft,rho_scr,cp_para_fft_pkg3d_dens_cp_box);
    }else{
      sum_rho(zfft,rho_scr,cp_para_fft_pkg3d_lg);
    }/*endif cp_dual_grid_opt */
  }/*endfor is*/
/*--------------------------------------------------------------------------*/
/*IV) if there is an odd number of states, use single pack
    for the last state (i caen't be chaengin' the laws of physics
    captn, i've got to be usin the single pack!!!!!!)  */

  if((nstate % 2 ) != 0) {
    ioff = (nstate-1)*ncoef;
    sngl_pack_coef(&creal[ioff],&cimag[ioff],zfft,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* V) fourier transform the last wavefunction to real space
     convention exp(-igr) */

    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*VI) add the square of the last wave function to the density  */

    if(cp_dual_grid_opt >= 1){
      sum_rho(zfft,rho_scr,cp_para_fft_pkg3d_dens_cp_box);
    }else{ 
      sum_rho(zfft,rho_scr,cp_para_fft_pkg3d_lg);
    }/*endif cp_dual_grid_opt*/
  }/*endif*/

/*==============================================================*/
}/*end routine*/
/*==============================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcRhoStoHybrid(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                   CP *cp,int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the routine to calculate the  density from stochastic WF.     */
/* This version is used in scf loop for different chemical potential.    */
/* We shall calculate the real space density for all chem pot and the    */
/* number of electrons. This routine follow the hybrid parallel scheme.  */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  #include "../typ_defs/typ_mask.h"

  EWALD        *ewald        = &(general_data->ewald);
  CELL         *cell         = &(general_data->cell);
  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR        *cpscr        = &(cp->cpscr);
  CPEWALD      *cpewald      = &(cp->cpewald);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  COMMUNICATE   *commCP         = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS  *cpcoeffs_pos   = &(cp->cpcoeffs_pos[ip_now]);
  PSEUDO        *pseudo         = &(cp->pseudo);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

  int cpParaOpt = cpopts->cp_para_opt;
  int cpLsda = cpopts->cp_lsda;
  int cpGGA  = cpopts->cp_gga;
  int numCoeffLarge       = cpcoeffs_info->ncoef_l;
  int numCoeffLargeProc   = cp->cp_para_fft_pkg3d_lg.ncoef_proc;
  int numCoeffLargeProcDensCpBox = cp->cp_para_fft_pkg3d_dens_cp_box.ncoef_proc;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numInterpPmeDual = pseudo->n_interp_pme_dual;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numChemPot     = stodftInfo->numChemPot;
  int numFFTProc        = cp_para_fft_pkg3d_lg->nfft_proc;
  int numFFT            = cp_para_fft_pkg3d_lg->nfft;
  int numFFT2           = numFFT/2;
  int numFFT2Proc       = numFFTProc/2;

  int rhoRealGridNum    = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot    = stodftInfo->rhoRealGridTot;
  int numChemProc       = stodftInfo->numChemProc;
  int numStateStoUp	= stodftInfo->numStateStoUp;
  int numStateStoDn	= stodftInfo->numStateStoDn;
  int occNumber		= stodftInfo->occNumber;
  int myidState		= commCP->myid_state;
  int numProcStates     = commCP->np_states;

  int iCoeff,iChem,iGrid;
  int index;
  int i,j,k;

  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);
  int *densityMap   = stodftInfo->densityMap;
  int *indexChemProc = stodftInfo->indexChemProc;
  int *chemProcIndexInv = stodftInfo->chemProcIndexInv;

  double volCP,rvolCP;
  double numGridTotInv = 1.0/rhoRealGridTot;
  double aveFactUp = occNumber/(double)(numStateStoUp);
  double aveFactDn;

  double *hmatCP    = cell->hmat_cp;
  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;
  double *rhoCoeffReUp   = cpscr->cpscr_rho.rhocr_up;
  double *rhoCoeffImUp   = cpscr->cpscr_rho.rhoci_up;
  double *rhoUp          = cpscr->cpscr_rho.rho_up;
  double *rhoCoeffReUpDensCpBox = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
  double *rhoCoeffImUpDensCpBox = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
  double *divRhoxUp       = cpscr->cpscr_grho.d_rhox_up;
  double *divRhoyUp       = cpscr->cpscr_grho.d_rhoy_up;
  double *divRhozUp       = cpscr->cpscr_grho.d_rhoz_up;
  double *d2RhoUp        = cpscr->cpscr_grho.d2_rho_up;
  double *rhoCoeffReDn   = cpscr->cpscr_rho.rhocr_dn;
  double *rhoCoeffImDn   = cpscr->cpscr_rho.rhoci_dn;
  double *rhoDn          = cpscr->cpscr_rho.rho_dn;
  double *rhoCoeffReDnDensCpBox = cpscr->cpscr_rho.rhocr_dn_dens_cp_box;
  double *rhoCoeffImDnDensCpBox = cpscr->cpscr_rho.rhoci_dn_dens_cp_box;
  double *divRhoxDn       = cpscr->cpscr_grho.d_rhox_dn;
  double *divRhoyDn       = cpscr->cpscr_grho.d_rhoy_dn;
  double *divRhozDn       = cpscr->cpscr_grho.d_rhoz_dn;
  double *d2RhoDn        = cpscr->cpscr_grho.d2_rho_dn;
  double *numElectron = stodftCoefPos->numElectron;
  double *numElectronTemp = (double*)cmalloc(numChemPot*sizeof(double));

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  double **rhoUpChemPot = stodftCoefPos->rhoUpChemPot;
  double **rhoDnChemPot = stodftCoefPos->rhoDnChemPot;

  double *rhoTemp  = (double*)cmalloc(numFFT2*sizeof(double))-1;

  MPI_Comm commStates   =    commCP->comm_states;

/*==========================================================================*/
/* I) Generate density in real space for all chem potentials		    */

  //vol_cp  = getdeth(hmat_cp);
  //rvol_cp = 1.0/vol_cp;
  
  for(iChem=0;iChem<numChemProc;iChem++){
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoUpChemPot[iChem][iGrid] = 0.0;
  }
  if(cpLsda==1&&numStateDnProc!=0){
    for(iChem=0;iChem<numChemProc;iChem++){
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoDnChemPot[iChem][iGrid] = 0.0;
    }
  }

  if(cpLsda==1&&numStateDnProc!=0)aveFactDn = occNumber/(double)(numStateStoDn);

  for(iChem=0;iChem<numChemPot;iChem++){
    rhoCalcRealStoHybrid(cpscr,cpcoeffs_info,
		   cell,stodftInfo,stoWfUpRe[iChem],
		   stoWfUpIm[iChem],rhoTemp,*coefFormUp,*coefOrthUp,
		   numStateUpProc,numCoeff,cpDualGridOptOn,commCP,
		   &(cp->cp_para_fft_pkg3d_lg),
		   &(cp->cp_sclr_fft_pkg3d_lg),
		   &(cp->cp_para_fft_pkg3d_dens_cp_box),
		   &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
		   &(cp->cp_sclr_fft_pkg3d_sm));
    if(numProcStates>1){
      Reduce(&rhoTemp[1],rhoUpChemPot[indexChemProc[iChem]],rhoRealGridNum,MPI_DOUBLE,
  	     MPI_SUM,densityMap[iChem],commStates);
    }
    else memcpy(rhoUpChemPot[iChem],&rhoTemp[1],rhoRealGridNum*sizeof(double));
    //debug
    printf("iChem %i rhoTemp %lg rhoUp %lg\n",iChem,rhoTemp[4001],rhoUpChemPot[iChem][4000]);
    // Calculate the average density, haven't scale by 1/volume
  }
  printf("11111111111111111\n");
  printf("rhoUp %lg\n",rhoUpChemPot[0][4000]);
  if(cpLsda==1&&numStateDnProc!=0){
    for(iChem=0;iChem<numChemPot;iChem++){
      rhoCalcRealStoHybrid(cpscr,cpcoeffs_info,
		   cell,stodftInfo,stoWfDnRe[iChem],
		   stoWfDnIm[iChem],rhoTemp,*coefFormDn,*coefOrthDn,
		   numStateDnProc,numCoeff,cpDualGridOptOn,commCP,
		   &(cp->cp_para_fft_pkg3d_lg),
		   &(cp->cp_sclr_fft_pkg3d_lg),
		   &(cp->cp_para_fft_pkg3d_dens_cp_box),
		   &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
		   &(cp->cp_sclr_fft_pkg3d_sm));
      if(numProcStates>1){
      Reduce(&rhoTemp[1],rhoDnChemPot[indexChemProc[iChem]],rhoRealGridNum,MPI_DOUBLE,
	     MPI_SUM,densityMap[iChem],commStates);
      }
      else memcpy(rhoDnChemPot[iChem],&rhoTemp[1],rhoRealGridNum);
    }  
  }
  // Calculate the average density, haven't scale by 1/volume
  for(i=0;i<10;i++)printf("i %i rhoUp %lg\n",i,rhoUpChemPot[0][i]);
  printf("aveFactUp %lg occNumber %i\n",aveFactUp,occNumber);
  printf("rhoRealGridNum %i rhoRealGridTot %i\n",rhoRealGridNum,rhoRealGridTot);
  for(iChem=0;iChem<numChemProc;iChem++){
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoUpChemPot[iChem][iGrid] *= aveFactUp;
  }
  if(cpLsda==1&&numStateDnProc!=0){
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoDnChemPot[iChem][iGrid] *= aveFactDn;
  }
  FILE *fileDensityTest = fopen("density-test","w");
  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)fprintf(fileDensityTest,"%lg\n",rhoUpChemPot[0][iGrid]);
  fclose(fileDensityTest);

  free(&rhoTemp[1]);

/*==========================================================================*/
/* II) Calculate number of electrons.					    */

  for(iChem=0;iChem<numChemPot;iChem++){
    numElectron[iChem] = 0.0;
    numElectronTemp[iChem] = 0.0;
  }
  
  printf("numChemProc %i\n",numChemProc); 
  for(iChem=0;iChem<numChemProc;iChem++){
    index = chemProcIndexInv[iChem];
    printf("index %i\n",index);
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
      numElectronTemp[index] += rhoUpChemPot[iChem][iGrid];     
    }
    if(cpLsda==1&&numStateDnProc!=0){
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
	numElectronTemp[index] += rhoDnChemPot[iChem][iGrid];
      }
    }
    numElectronTemp[index] *= numGridTotInv;
    printf("numElectronTemp %lg\n",numElectronTemp[index]);
  }
  if(numProcStates>1){
    Allreduce(numElectronTemp,numElectron,numChemPot,MPI_DOUBLE,MPI_SUM,0,commStates);
  }
  else memcpy(numElectron,numElectronTemp,numChemPot*sizeof(double));

  free(numElectronTemp);

//debug
  if(myidState==0){
    double *chemPot = stodftCoefPos->chemPot;
    for(iChem=0;iChem<numChemPot;iChem++)printf("chemPot-numElec %lg %lg\n",chemPot[iChem],numElectron[iChem]);
  }

/*==========================================================================*/
/* III) Interpolate the correct # of electron			            */
  
  
/*==========================================================================*/
/* IV) Interpolate the correct density		                            */
 
/*==========================================================================*/
/* IV) Generate the reciprocal part and all the other things                */

  //We store it in the rhoTemp
  //calcRhoStoRecip();
 
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


