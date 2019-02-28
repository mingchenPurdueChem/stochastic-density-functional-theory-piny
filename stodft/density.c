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
  double sum_sq = 0.0;
  //debug
  //char name[100];
  //FILE *fstowf;

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
    /*
    sprintf(name,"stowf-%i",is-1);
    fstowf = fopen(name,"w");
    for(i=1;i<=nfft2;i++)fprintf(fstowf,"%.16lg\n",-zfft[2*i-1]);
    fclose(fstowf);
    sprintf(name,"stowf-%i",is);
    fstowf = fopen(name,"w");
    for(i=1;i<=nfft2;i++)fprintf(fstowf,"%.16lg\n",-zfft[2*i]);
    fclose(fstowf);
    */
    //printf("1111111111 rho %lg\n",zfft[1985]*zfft[1985]);
    //printf("1111111111 rho %lg\n",zfft[1986]*zfft[1986]);
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
void calcRhoStoRecipFullg(CPEWALD *cpewald,CPSCR *cpscr,
                        CPCOEFFS_INFO *cpcoeffs_info,EWALD *ewald,CELL *cell,
                        double *rhocr ,double *rhoci,double *rho,
                        double *rhocr_dens_cp_box,double *rhoci_dens_cp_box,
                        double *del_rho_x, double *del_rho_y,
                        double *del_rho_z,
                        double *del2_rho,
                        int cp_gga,int cp_dual_grid_opt,
                        int n_interp_pme_dual,
                        COMMUNICATE *communicate,
                        PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg,
                        PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box)
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

/* local variables                                                  */

  int iii,ioff,ioff2;
  int is,i,iupper;
  int cp_elf_calc_frq  = cpcoeffs_info->cp_elf_calc_frq;
  int myid_state       = communicate->myid_state;
  int np_states        = communicate->np_states;
  int laplacian_on     = cpcoeffs_info->cp_laplacian_on;
  int ncoef_l          = cp_para_fft_pkg3d_lg->ncoef;
  int ncoef_l_use      = cp_para_fft_pkg3d_lg->ncoef_use;
  int ncoef_l_proc     = cp_para_fft_pkg3d_lg->ncoef_proc;
  int nfft_proc        = cp_para_fft_pkg3d_lg->nfft_proc;
  int nfft             = cp_para_fft_pkg3d_lg->nfft;
  int nfft2_proc       = nfft_proc/2;
  int nfft2            = nfft/2;
  int nfft_proc_dens_cp_box,nfft2_proc_dens_cp_box;


  double vol_cp,rvol_cp;
  double temp_r,temp_i;

  double *zfft           = cpscr->cpscr_wave.zfft;
  double *zfft_tmp       = cpscr->cpscr_wave.zfft_tmp;
  double *hmati_cp       = cell->hmati_cp;
  double *hmat_cp        = cell->hmat_cp;

  double dbox_rat        = cpewald->dbox_rat;
  double *bw_r           = cpscr->cpscr_dual_pme.bw_r;
  double *bw_i           = cpscr->cpscr_dual_pme.bw_i;

#ifdef WRITE_DENSITY
#include "../proto_defs/proto_friend_lib_entry.h"
  int nkf1,nkf2,nkf3;
  int index,ka,kb,kc;
  FILE *fp_rho;
#endif
#ifdef DEBUG_LSDA
  //double *rhocr_dens_cp_box;
  //double *rhoci_dens_cp_box;
#endif

  MPI_Comm comm_states = communicate->comm_states;

  if(cp_dual_grid_opt >= 1){
    nfft_proc_dens_cp_box   =  cp_para_fft_pkg3d_dens_cp_box->nfft_proc;
    nfft2_proc_dens_cp_box  =  nfft_proc_dens_cp_box/2;
#ifdef DEBUG_LSDA
    if((cp_gga == 1)|| (cp_dual_grid_opt > 1)){
      /*holds the density in g space on cp_box grid*/
      rhocr_dens_cp_box    = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
      rhoci_dens_cp_box    = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
    }/* endif cp_gga */
#endif
  }/*endif cp_dual_grid_opt*/

/*=========================================================================*/
/*=========================================================================*/
/*  3) get density in g space                                              */
/*  I)  pack it up                                                         */

  if(cp_dual_grid_opt >= 1){
    /* sending density*vol_cp on small grid */
    control_spread_rho(cpscr,rho,cell,dbox_rat,np_states,n_interp_pme_dual,
                       cp_para_fft_pkg3d_dens_cp_box,
                       cp_para_fft_pkg3d_lg,cp_dual_grid_opt);
  }else{
    sngl_pack_rho(zfft,rho,cp_para_fft_pkg3d_lg);
  }/*endif cp_dual_grid_opt*/


/*--------------------------------------------------------------------------*/
/*  II) back transform to g-space  convention exp(igr)   */

  para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_lg);

/*--------------------------------------------------------------------------*/
/*  III) unpack the density                                    */

  sngl_upack_coef(rhocr,rhoci,zfft,cp_para_fft_pkg3d_lg);

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
  }

/*=========================================================================*/
/* IF CP_DUAL_GRID_OPT and doing GRADIENT CORRECTIONS you need the density in*/
/*   g-space on the mid size grid corresponding to the CP_BOX              */
/*  3.5) get density in G SPACE for CP_BOX                                 */

  if((cp_dual_grid_opt == 1 && cp_gga == 1) || (cp_dual_grid_opt > 1)){
/*  I)  pack it up                                                         */

    sngl_pack_rho(zfft,rho,cp_para_fft_pkg3d_dens_cp_box);

/*--------------------------------------------------------------------------*/
/*  II) back transform to g-space  convention exp(igr)   */

    para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_dens_cp_box);

/*--------------------------------------------------------------------------*/
/*  III) unpack the density                                                 */

    sngl_upack_coef(rhocr_dens_cp_box,rhoci_dens_cp_box,
                    zfft,cp_para_fft_pkg3d_dens_cp_box);

  }/*endif cp_dual_grid_opt and cp_gga*/

/*--------------------------------------------------------------------------*/
/* IV) finish the density in real space by dividing by the volume */

  vol_cp  = getdeth(hmat_cp);
  rvol_cp = 1.0/vol_cp;
  if(cp_dual_grid_opt >= 1){
    for(i=1 ; i<= nfft2_proc_dens_cp_box;i++)rho[i] *= rvol_cp;
  }else{
    for(i=1 ; i<= nfft2_proc;i++)rho[i] *= rvol_cp;
  }/*endif cp_dual_grid_opt*/

  Barrier(comm_states);
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
void calcRhoStoHybridInterp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
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
  //PARA_FFT_PKG3D *cp_para_fft_pkg3d;

  int cpParaOpt = cpopts->cp_para_opt;
  int cpLsda = cpopts->cp_lsda;
  int cpGGA  = cpopts->cp_gga;
  int realSparseOpt = cpopts->realSparseOpt;
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
  int densityMixFlag	= stodftInfo->densityMixFlag;
  int iScf		= stodftInfo->iScf;
  int myidState		= commCP->myid_state;
  int numProcStates     = commCP->np_states;

  int iCoeff,iChem,iGrid;
  int index;
  int i,j,k;
  int reRunFlag;

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
  double numElecTrue = stodftInfo->numElecTrue;
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
  double *rhoUpCorrect	  = stodftCoefPos->rhoUpCorrect;
  double *divRhoxDn       = cpscr->cpscr_grho.d_rhox_dn;
  double *divRhoyDn       = cpscr->cpscr_grho.d_rhoy_dn;
  double *divRhozDn       = cpscr->cpscr_grho.d_rhoz_dn;
  double *d2RhoDn        = cpscr->cpscr_grho.d2_rho_dn;
  double *numElectron = stodftCoefPos->numElectron;
  double *chemPot = stodftCoefPos->chemPot;
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
  /*
  if(realSparseOpt==0){
    cp_para_fft_pkg3d = &(cp->cp_para_fft_pkg3d_lg);
  }
  else{
    cp_para_fft_pkg3d = &(cp->cp_para_fft_pkg3d_sparse);
  }
  */
  

  if(numProcStates>1)Barrier(commStates);
  
  for(iChem=0;iChem<numChemProc;iChem++){
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoUpChemPot[iChem][iGrid] = 0.0;
  }
  if(cpLsda==1&&numStateDnProc!=0){
    for(iChem=0;iChem<numChemProc;iChem++){
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoDnChemPot[iChem][iGrid] = 0.0;
    }
  }

  if(cpLsda==1&&numStateDnProc!=0)aveFactDn = occNumber/(double)(numStateStoDn);

  //printf("111111111 Finish Initialize Densities\n");
  //fflush(stdout);
  if(numProcStates>1)Barrier(commStates);

  //debug
  //for(iChem=0;iChem<numChemPot;iChem++)printf("myid %i densityMap %i\n",myidState,densityMap[iChem],indexChemProc[iChem]);
  //fflush(stdout);
  //if(numProcStates>1)Barrier(commStates);

  for(iChem=0;iChem<numChemPot;iChem++){
    //debug
    if(myidState==densityMap[iChem]){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoUpChemPot[indexChemProc[iChem]][iGrid] = 0.0;
    }
    //end debug
    rhoCalcRealStoHybrid(cpscr,cpcoeffs_info,
		   cell,stodftInfo,stoWfUpRe[iChem],
		   stoWfUpIm[iChem],rhoTemp,*coefFormUp,*coefOrthUp,
		   numStateUpProc,numCoeff,cpDualGridOptOn,commCP,
		   &(cp->cp_para_fft_pkg3d_lg),
		   &(cp->cp_sclr_fft_pkg3d_lg),
		   &(cp->cp_para_fft_pkg3d_dens_cp_box),
		   &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
		   &(cp->cp_sclr_fft_pkg3d_sm));
    //debug
    /*
    if(checkNanArray(&rhoTemp[1],rhoRealGridTot)==1){
      printf("iChem %i myid %i density is wrong before reduce!\n",iChem,myidState);
    }
    */
    /*
    if(myidState==0){
      printf("iChem %i finish density calc\n",iChem);
      printf("rhoUpChemPot[0] %p\n",rhoUpChemPot[0]);
      fflush(stdout);
    }
    printf("myidState %i iChem %i indexChemProc[iChem] %i densityMap[iChem] %i rhoRealGridTot %i rhoUpChemPot %p rhoUpChemPot again %p numFFT2 %i\n",
	    myidState,iChem,indexChemProc[iChem],densityMap[iChem],rhoRealGridTot,rhoUpChemPot[0],&rhoUpChemPot[0][0],numFFT2);
    */
    //end debug
    if(numProcStates>1){
      Barrier(commStates);
      Reduce(&rhoTemp[1],rhoUpChemPot[indexChemProc[iChem]],rhoRealGridTot,MPI_DOUBLE,
	     MPI_SUM,densityMap[iChem],commStates);
      //MPI_Reduce(&rhoTemp[1],rhoUpChemPot[0],rhoRealGridTot,MPI_DOUBLE,MPI_SUM,densityMap[iChem],commStates);
      Barrier(commStates);
    }
    else memcpy(rhoUpChemPot[iChem],&rhoTemp[1],rhoRealGridTot*sizeof(double));
    //debug
    //printf("iChem %i rhoTemp %lg rhoUp %lg\n",iChem,rhoTemp[4001],rhoUpChemPot[iChem][4000]);
    // Calculate the average density, haven't scale by 1/volume
  }
  //printf("rhoUp %lg\n",rhoUpChemPot[0][4000]);
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
        Reduce(&rhoTemp[1],rhoDnChemPot[indexChemProc[iChem]],rhoRealGridTot,MPI_DOUBLE,
	     MPI_SUM,densityMap[iChem],commStates);
      }
      else memcpy(rhoDnChemPot[iChem],&rhoTemp[1],rhoRealGridTot);
    }  
  }
  if(numProcStates>1)Barrier(commStates);

  // Calculate the average density, haven't scale by 1/volume
  
  //for(i=0;i<10;i++)printf("i %i rhoUp %lg\n",i,rhoUpChemPot[0][i]);
  //printf("aveFactUp %lg occNumber %i\n",aveFactUp,occNumber);
  //printf("rhoRealGridNum %i rhoRealGridTot %i\n",rhoRealGridNum,rhoRealGridTot);
  for(iChem=0;iChem<numChemProc;iChem++){
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoUpChemPot[iChem][iGrid] *= aveFactUp;
  }
  if(cpLsda==1&&numStateDnProc!=0){
    for(iChem=0;iChem<numChemProc;iChem++){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoDnChemPot[iChem][iGrid] *= aveFactDn;
    }
  }

  // debug
  /*
  for(iChem=0;iChem<numChemProc;iChem++){
    if(checkNanArray(rhoUpChemPot[iChem],rhoRealGridTot)==1){
      printf("iChem %i myid %i density is wrong after reduce\n",iChem,myidState);
    }
  }
  */

  // debug
  /*
  FILE *fileDensityTest = fopen("density-test","w");
  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)fprintf(fileDensityTest,"%.16lg\n",rhoUpChemPot[0][iGrid]);
  fclose(fileDensityTest);
  */
  

  free(&rhoTemp[1]);

  if(numProcStates>1)Barrier(commStates);


/*==========================================================================*/
/* II) Calculate number of electrons.					    */

  for(iChem=0;iChem<numChemPot;iChem++){
    numElectron[iChem] = 0.0;
    numElectronTemp[iChem] = 0.0;
  }
  
  //printf("numChemProc %i\n",numChemProc); 
  for(iChem=0;iChem<numChemProc;iChem++){
    index = chemProcIndexInv[iChem];
    //printf("index %i\n",index);
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      numElectronTemp[index] += rhoUpChemPot[iChem][iGrid];     
    }
    if(cpLsda==1&&numStateDnProc!=0){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	numElectronTemp[index] += rhoDnChemPot[iChem][iGrid];
      }
    }
    numElectronTemp[index] *= numGridTotInv;
    //printf("numElectronTemp %lg\n",numElectronTemp[index]);
  }
  if(numProcStates>1){
    Allreduce(numElectronTemp,numElectron,numChemPot,MPI_DOUBLE,MPI_SUM,0,commStates);
  }
  else memcpy(numElectron,numElectronTemp,numChemPot*sizeof(double));

  free(numElectronTemp);
   
  if(myidState==0){
    printf("Min Chem Pot %.16lg # Electron %.16lg\n",chemPot[0],numElectron[0]);
    printf("Max Chem Pot %.16lg # Electron %.16lg\n",chemPot[numChemPot-1],numElectron[numChemPot-1]);

    //for(iChem=0;iChem<numChemPot;iChem++)printf("chemPot-numElec %lg %lg\n",chemPot[iChem],numElectron[iChem]);
  }
  if(numElecTrue<numElectron[0]||numElecTrue>numElectron[numChemPot-1]){
    if(myidState==0){
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      printf("Improper chemical potential range! This step of SCF will\n");
      printf("rerun until chemical potential range is large enough!\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    }
    reRunFlag = 1;
    stodftInfo->reRunFlag = 1;
  }
  else{
    reRunFlag = 0;
    stodftInfo->reRunFlag = 0;
  }
  
  if(numProcStates>1)Barrier(commStates);
  //exit(0);

/*==========================================================================*/
/* III) Interpolate the correct # of electron and the correct density       */

  if(reRunFlag==0){// no rerun
    calcChemPotInterp(cp);
    if(myidState==0){
      printf("Finish Interpolating Chemical Potential\n");
      fflush(stdout);
    }
    /*
    double checkElecNumProc = 0.0;
    double checkElecNum = 0.0;
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
      checkElecNumProc += rhoUpCorrect[iGrid];
    }//endfor iGrid
    if(numProcStates>1){
      Reduce(&checkElecNumProc,&checkElecNum,1,MPI_DOUBLE,MPI_SUM,0,commStates);
    }
    if(myidState==0){
      checkElecNum *= numGridTotInv;
      printf("Test num electron %.16lg\n",checkElecNum);
    }
    */

/*==========================================================================*/
/* IV) Output the density                                                   */

    outputDensity(cp,cell);

/*==========================================================================*/
/* V) Generate the diis density						    */

    if(densityMixFlag>0){
      if(myidState==0){
	printf("Start Mixing Density\n");
      }
      genDensityMix(cp,iScf);
      if(myidState==0){
	printf("Finish Mixing Density\n");
      }
    }
    
/*==========================================================================*/
/* VI) Generate the reciprocal part and all the other things                */
    /*
    FILE *fileRhoReal = fopen("density-init","r");
    FILE *fileRhoRecip = fopen("density-recip-test","w");
    for(iGrid=1;iGrid<=rhoRealGridNum;iGrid++)fscanf(fileRhoReal,"%lg",&rhoUp[iGrid]);
    for(iGrid=1;iGrid<=rhoRealGridNum;iGrid++)rhoUp[iGrid] /= 0.0009250463018013585;
    */
    calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
		       rhoCoeffReUp,rhoCoeffImUp,rhoUp,rhoCoeffReUpDensCpBox,rhoCoeffImUpDensCpBox,
		       divRhoxUp,divRhoyUp,divRhozUp,d2RhoUp,cpGGA,cpDualGridOptOn,numInterpPmeDual,
		       commCP,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));
    if(cpLsda==1&&numStateDnProc!=0){
      calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
			 rhoCoeffReDn,rhoCoeffImDn,rhoDn,rhoCoeffReDnDensCpBox,rhoCoeffImDnDensCpBox,
			 divRhoxDn,divRhoyDn,divRhozDn,d2RhoDn,cpGGA,cpDualGridOptOn,numInterpPmeDual,
			 commCP,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));
      for(iCoeff=1;iCoeff<=numCoeffLargeProc;iCoeff++) {
	rhoCoeffReUp[iCoeff] += rhoCoeffReDn[iCoeff];
	rhoCoeffImUp[iCoeff] += rhoCoeffImDn[iCoeff];
      }/* endfor */
      if(cpDualGridOptOn>=1){
	for(iCoeff=1;iCoeff<=numCoeffLargeProcDensCpBox;iCoeff++){
	  rhoCoeffReUpDensCpBox[iCoeff] += rhoCoeffReDnDensCpBox[iCoeff];
	  rhoCoeffImUpDensCpBox[iCoeff] += rhoCoeffImDnDensCpBox[iCoeff];
	}/* endfor */
      } /* endif */
    }/* endif */

    if(myidState==0){
      printf("Finish Calculating Reciprocal Space Density\n");
    }


    /*
    for(iCoeff=1;iCoeff<=numCoeffLargeProc;iCoeff++){
      fprintf(fileRhoRecip,"%.10lg %.10lg\n",rhoCoeffReUp[iCoeff],rhoCoeffImUp[iCoeff]);
    }
    fclose(fileRhoReal);
    fclose(fileRhoRecip);  
    */
  }//endif reRunFlag
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcRhoStoHybridCheby(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
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
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
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
  int calcFragFlag	= stodftInfo->calcFragFlag;

  int rhoRealGridNum    = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot    = stodftInfo->rhoRealGridTot;
  int numChemProc       = stodftInfo->numChemProc;
  int numStateStoUp	= stodftInfo->numStateStoUp;
  int numStateStoDn	= stodftInfo->numStateStoDn;
  int occNumber		= stodftInfo->occNumber;
  int densityMixFlag	= stodftInfo->densityMixFlag;
  int iScf		= stodftInfo->iScf;
  int myidState		= commCP->myid_state;
  int numProcStates     = commCP->np_states;

  int iCoeff,iChem,iGrid,iProc;
  int index;
  int i,j,k;
  int reRunFlag;
  MPI_Comm comm_states = commCP->comm_states;

  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);
  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;

  double volCP,rvolCP;
  double numGridTotInv = 1.0/rhoRealGridTot;
  double aveFactUp = occNumber/(double)(numStateStoUp);
  double numElecTrue = stodftInfo->numElecTrue;
  double aveFactDn;
  double numElecTest;

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
  double *rhoUpCorrect = stodftCoefPos->rhoUpCorrect;
  double *rhoDnCorrect = stodftCoefPos->rhoDnCorrect;
  double *divRhoxDn       = cpscr->cpscr_grho.d_rhox_dn;
  double *divRhoyDn       = cpscr->cpscr_grho.d_rhoy_dn;
  double *divRhozDn       = cpscr->cpscr_grho.d_rhoz_dn;
  double *d2RhoDn        = cpscr->cpscr_grho.d2_rho_dn;
  double *chemPot = stodftCoefPos->chemPot;
  double *numElectronTemp = (double*)cmalloc(numChemPot*sizeof(double));
  double *rhoUpFragSum;
  double *rhoDnFragSum;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  double **rhoUpChemPot = stodftCoefPos->rhoUpChemPot;
  double **rhoDnChemPot = stodftCoefPos->rhoDnChemPot;

  double *rhoTemp  = (double*)cmalloc(numFFT2*sizeof(double))-1;
  double *rhoReduce;
 
  MPI_Comm commStates   =    commCP->comm_states;

/*==========================================================================*/
/* I) Generate density in real space for all chem potentials		    */


  if(numProcStates>1)Barrier(commStates);
  
  if(myidState==0){
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoUpChemPot[0][iGrid] = 0.0;
    if(cpLsda==1&&numStateDnProc!=0){
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoDnChemPot[0][iGrid] = 0.0;
    }
    rhoReduce = (double*)cmalloc(rhoRealGridTot*sizeof(double));
  }


  if(cpLsda==1&&numStateDnProc!=0)aveFactDn = occNumber/(double)(numStateStoDn);

  if(numProcStates>1)Barrier(commStates);

  rhoCalcRealStoHybrid(cpscr,cpcoeffs_info,
		 cell,stodftInfo,stoWfUpRe[0],
		 stoWfUpIm[0],rhoTemp,*coefFormUp,*coefOrthUp,
		 numStateUpProc,numCoeff,cpDualGridOptOn,commCP,
		 &(cp->cp_para_fft_pkg3d_lg),
		 &(cp->cp_sclr_fft_pkg3d_lg),
		 &(cp->cp_para_fft_pkg3d_dens_cp_box),
		 &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
		 &(cp->cp_sclr_fft_pkg3d_sm));
  if(numProcStates>1){
    Barrier(commStates);
    //printf("%p %i\n",rhoUpChemPot[0],rhoRealGridTot);
    /*
    if(myidState==1){
      for(iGrid=1;iGrid<=rhoRealGridTot;iGrid++){
	printf("%p\n",&rhoTemp[iGrid]);
      }
    }
    */
    /*
    printf("PID %i\n",getpid());
    int idebug = 0;
    while(idebug==0){
      sleep(5);
    }
    */
    Barrier(commStates);
    Reduce(&rhoTemp[1],rhoUpChemPot[0],rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    //Reduce(&rhoTemp[1],rhoReduce,rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    //if(myidState==0) memcpy(rhoUpChemPot[0],rhoReduce,rhoRealGridTot*sizeof(double));
    Barrier(commStates);
  }
  else memcpy(rhoUpChemPot[0],&rhoTemp[1],rhoRealGridTot*sizeof(double));
  if(cpLsda==1&&numStateDnProc!=0){
    rhoCalcRealStoHybrid(cpscr,cpcoeffs_info,
		cell,stodftInfo,stoWfDnRe[0],
		stoWfDnIm[0],rhoTemp,*coefFormDn,*coefOrthDn,
		numStateDnProc,numCoeff,cpDualGridOptOn,commCP,
		&(cp->cp_para_fft_pkg3d_lg),
		&(cp->cp_sclr_fft_pkg3d_lg),
		&(cp->cp_para_fft_pkg3d_dens_cp_box),
		&(cp->cp_sclr_fft_pkg3d_dens_cp_box),
		&(cp->cp_sclr_fft_pkg3d_sm));
    if(numProcStates>1){
      Reduce(&rhoTemp[1],rhoDnChemPot[0],rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    }
    else memcpy(rhoDnChemPot[iChem],&rhoTemp[1],rhoRealGridTot);
  }
  if(numProcStates>1)Barrier(commStates);

  // Calculate the average density, haven't scale by 1/volume
  if(myidState==0){
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoUpChemPot[0][iGrid] *= aveFactUp;
    if(cpLsda==1&&numStateDnProc!=0){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoDnChemPot[0][iGrid] *= aveFactDn;
    }
  }

  // debug
  /*
  FILE *fileDensityTest = fopen("density-test","w");
  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)fprintf(fileDensityTest,"%.16lg\n",rhoUpChemPot[0][iGrid]);
  fclose(fileDensityTest);
  */


  free(&rhoTemp[1]);

  if(numProcStates>1)Barrier(commStates);

/*==========================================================================*/
/* III) Double Check number of electrons.				    */

  if(myidState==0){
    numElecTest = 0.0;
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)numElecTest += rhoUpChemPot[0][iGrid];
    if(cpLsda==1&&numStateDnProc!=0){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)numElecTest += rhoDnChemPot[0][iGrid];
    }
    numElecTest *= numGridTotInv;
    printf("Number of electron is %.16lg.\n",numElecTest);
  }
  
  if(numProcStates>1)Barrier(commStates);


/*==========================================================================*/
/* III) Scatter the density to all nodes			            */

  if(numProcStates>1){
    Scatterv(rhoUpChemPot[0],rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	     rhoUpCorrect,rhoRealGridNum,MPI_DOUBLE,0,comm_states);
  }
  else{
    memcpy(rhoUpCorrect,rhoUpChemPot[0],rhoRealGridNum*sizeof(double));
  }
  if(cpLsda==1&&numStateDnProc!=0){
    if(numProcStates>1){   
      Scatterv(rhoDnChemPot[0],rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	       rhoDnCorrect,rhoRealGridNum,MPI_DOUBLE,0,comm_states);
    }
    else{
      memcpy(rhoDnCorrect,rhoDnChemPot[0],rhoRealGridNum*sizeof(double));
    }
  }
  
  //outputDensity(cp,cell);


/*==========================================================================*/
/* IV) Add the fragment contribution.                                       */

  if(calcFragFlag==1){
    rhoUpFragSum =  fragInfo->rhoUpFragSum;
    //debug
    /*
    for(iProc=0;iProc<numProcStates;iProc++){
      if(myidState==iProc){
        for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)printf("rhostooooo %lg %lg %lg\n",rhoUpCorrect[iGrid],rhoUpFragSum[iGrid],rhoUpCorrect[iGrid]+rhoUpFragSum[iGrid]);
      }
      if(numProcStates>1)Barrier(commStates);
    }
    */
    //printf("rhoUpCorrect %lg %lg rhoUpFragSum %lg %lg\n",rhoUpCorrect[0],rhoUpCorrect[1],rhoUpFragSum[0],rhoUpFragSum[1]);   
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoUpCorrect[iGrid] += rhoUpFragSum[iGrid];
    if(cpLsda==1&&numStateDnProc!=0){
      rhoDnFragSum =  fragInfo->rhoDnFragSum;
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoDnCorrect[iGrid] += rhoDnFragSum[iGrid];
    }
  }

/*==========================================================================*/
/* V) Output the density                                                    */
  //printf("22222222 rho %.16lg %.16lg\n",rhoUpCorrect[1],rhoUpCorrect[2]);

  outputDensity(cp,cell);

/*==========================================================================*/
/* VI) Generate the diis density					    */

  if(densityMixFlag>0){
    if(myidState==0){
      printf("Start Mixing Density\n");
      fflush(stdout);
    }
    genDensityMix(cp,iScf);
    if(myidState==0){
      printf("Finish Mixing Density\n");
      fflush(stdout);
    }
  }
 
    /*
    double sum = 0.0;
    volCP  = getdeth(hmatCP);
    rvolCP = 1.0/volCP;

    FILE *fp_rho = fopen("rho_test","w");
    for(iGrid=1;iGrid<=rhoRealGridTot;iGrid++){
       fprintf(fp_rho,"%.5e\n",rhoUp[iGrid]*rvolCP);
    }
    //printf("rho sum test %lg\n",sum);
    fclose(fp_rho);
    */
  

  //outputDensity(cp,cell);
 
/*==========================================================================*/
/* VII) Generate the reciprocal part and all the other things               */
    /*
    FILE *fileRhoReal = fopen("density-init","r");
    FILE *fileRhoRecip = fopen("density-recip-test","w");
    for(iGrid=1;iGrid<=rhoRealGridNum;iGrid++)fscanf(fileRhoReal,"%lg",&rhoUp[iGrid]);
    for(iGrid=1;iGrid<=rhoRealGridNum;iGrid++)rhoUp[iGrid] /= 0.0009250463018013585;
    */

  //fflush(stdout);
  //exit(0);

  calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
		     rhoCoeffReUp,rhoCoeffImUp,rhoUp,rhoCoeffReUpDensCpBox,rhoCoeffImUpDensCpBox,
		     divRhoxUp,divRhoyUp,divRhozUp,d2RhoUp,cpGGA,cpDualGridOptOn,numInterpPmeDual,
		     commCP,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));
  if(cpLsda==1&&numStateDnProc!=0){
    calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
		       rhoCoeffReDn,rhoCoeffImDn,rhoDn,rhoCoeffReDnDensCpBox,rhoCoeffImDnDensCpBox,
		       divRhoxDn,divRhoyDn,divRhozDn,d2RhoDn,cpGGA,cpDualGridOptOn,numInterpPmeDual,
		       commCP,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));
    for(iCoeff=1;iCoeff<=numCoeffLargeProc;iCoeff++){
      rhoCoeffReUp[iCoeff] += rhoCoeffReDn[iCoeff];
      rhoCoeffImUp[iCoeff] += rhoCoeffImDn[iCoeff];
    }/* endfor */
    if(cpDualGridOptOn>=1){
      for(iCoeff=1;iCoeff<=numCoeffLargeProcDensCpBox;iCoeff++){
	rhoCoeffReUpDensCpBox[iCoeff] += rhoCoeffReDnDensCpBox[iCoeff];
	rhoCoeffImUpDensCpBox[iCoeff] += rhoCoeffImDnDensCpBox[iCoeff];
      }/* endfor */
    } /* endif */
  }/* endif */
  //printf("2222222222222 rhok %.16lg %.16lg\n",rhoCoeffReUp[1],rhoCoeffImUp[1]);

   //exit(0);
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  outputDensity:                                                          */
/*      Gather the density to the master processor and print it to file     */
/*==========================================================================*/
void outputDensity(CP *cp,CELL *cell)
/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"
  COMMUNICATE *commCP		= &(cp->communicate);
  STODFTINFO *stodftInfo	= cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos	= cp->stodftCoefPos;
  CPOPTS *cpopts		= &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);

  int myidState		    = commCP->myid_state;
  int numProcStates	    = commCP->np_states;
  int cpLsda		    = cpopts->cp_lsda;
  int rhoRealGridNum	    = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot	    = stodftInfo->rhoRealGridTot;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int iGrid;
  int chemPotOpt = stodftInfo->chemPotOpt;
  int *rhoRealSendCounts    = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls	    = stodftInfo->rhoRealDispls;
  char *densityFileName	    = stodftInfo->densityFileName;

  double volCP,volCPInv;
  double *hmatCP	    = cell->hmat_cp;
  double *rhoUpCorrect	    = stodftCoefPos->rhoUpCorrect;
  double *rhoDnCorrect	    = stodftCoefPos->rhoDnCorrect;
  double *rhoUpForOutput,*rhoDnForOutput;

  double **rhoUpChemPot = stodftCoefPos->rhoUpChemPot;
  double **rhoDnChemPot = stodftCoefPos->rhoDnChemPot;

  FILE *densityOutputFile;
  MPI_Comm commStates	    = commCP->comm_states;

  volCP  = getdeth(hmatCP);
  volCPInv = 1.0/volCP;
  // test, I'll write a better version later
  FILE *densityFinal;

  //if(chemPotOpt==1){// interpolation way
    if(numProcStates>1){
      if(myidState==0)rhoUpForOutput = (double*)cmalloc(rhoRealGridTot*sizeof(double));
      Barrier(commStates);
      Gatherv(rhoUpCorrect,rhoRealGridNum,MPI_DOUBLE,rhoUpForOutput,
	      rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,0,commStates);
      Barrier(commStates);
      if(cpLsda==1&&numStateDnProc!=0){
	if(myidState==0)rhoDnForOutput = (double*)cmalloc(rhoRealGridTot*sizeof(double));
	Barrier(commStates);
	Gatherv(rhoDnCorrect,rhoRealGridNum,MPI_DOUBLE,rhoDnForOutput,
		rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,0,commStates);
	Barrier(commStates);
      }//endif cpLsda
      if(myidState==0){
	densityOutputFile = cfopen(densityFileName,"a");
        printf("I'm writing the output density\n");
        densityFinal = fopen("density-out","w");
	for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	  //fprintf(densityOutputFile,"%.10lg\n",rhoUpForOutput[iGrid]*volCPInv);
	  //debug
	  fprintf(densityOutputFile,"%.10lg\n",rhoUpForOutput[iGrid]*volCPInv);
          fprintf(densityFinal,"%.16lg\n",rhoUpForOutput[iGrid]);
	}//endfor iGrid
	cfree(rhoUpForOutput);
	if(cpLsda==1&&numStateDnProc!=0){
	  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	    fprintf(densityOutputFile,"%.10lg\n",rhoDnForOutput[iGrid]*volCPInv);
            fprintf(densityFinal,"%.16lg\n",rhoDnForOutput[iGrid]);
	  }//endfor iGrid
	  cfree(rhoDnForOutput);
	}//endif cpLsda
	fclose(densityOutputFile);
        fclose(densityFinal);
      }//endif myState
    }//endif parallel case
    else{
      densityOutputFile = cfopen(densityFileName,"a");
      densityFinal = fopen("density-out","w");
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
	fprintf(densityOutputFile,"%.10lg\n",rhoUpCorrect[iGrid]);
        fprintf(densityFinal,"%.16lg\n",rhoUpCorrect[iGrid]);
      }//endfor iGrid
      if(cpLsda==1&&numStateDnProc!=0){
	for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
	  fprintf(densityOutputFile,"%.10lg\n",rhoDnCorrect[iGrid]);
          fprintf(densityFinal,"%.16lg\n",rhoDnCorrect[iGrid]);
	}//endfor iGrid
      }//endif cpLsda
      fclose(densityOutputFile);
      fclose(densityFinal);
    }//endif sequential case
  //}
  /*
  else if(chemPotOpt==2){
    if(myidState==0){
      densityOutputFile = cfopen(densityFileName,"a");
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	fprintf(densityOutputFile,"%.10lg\n",rhoUpChemPot[0][iGrid]*volCPInv);
      }//endfor iGrid
      if(cpLsda==1&&numStateDnProc!=0){
	for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	  fprintf(densityOutputFile,"%.10lg\n",rhoUpChemPot[0][iGrid]*volCPInv);
	}//endfor iGrid
      }//endif cpLsda
      fclose(densityOutputFile);
    }
  }
  */

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

void calcRhoFilterDiagHybrid(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
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
  int cpGGA  = cpopts->cp_gga;
  int cpLsda                = cpopts->cp_lsda;
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
  int numStateStoUp     = stodftInfo->numStateStoUp;
  int numStateStoDn     = stodftInfo->numStateStoDn;
  int occNumber         = stodftInfo->occNumber;
  int densityMixFlag    = stodftInfo->densityMixFlag;
  int iScf              = stodftInfo->iScf;
  int myidState         = commCP->myid_state;
  int numProcStates     = commCP->np_states;

  int iCoeff,iChem,iGrid;
  int index;
  int i,j,k;
  int reRunFlag;

  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);
  int *densityMap   = stodftInfo->densityMap;
  int *indexChemProc = stodftInfo->indexChemProc;
  int *chemProcIndexInv = stodftInfo->chemProcIndexInv;
  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;

  double volCP,rvolCP;
  double numGridTotInv = 1.0/rhoRealGridTot;
  double aveFactUp = occNumber/(double)(numStateStoUp);
  double numElecTrue = stodftInfo->numElecTrue;
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
  double *chemPot = stodftCoefPos->chemPot;
  double *rhoUpCorrect = stodftCoefPos->rhoUpCorrect;
  double *numElectronTemp = (double*)cmalloc(numChemPot*sizeof(double));

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  double **rhoUpChemPot = stodftCoefPos->rhoUpChemPot;
  double **rhoDnChemPot = stodftCoefPos->rhoDnChemPot;

  double *rhoTemp  = (double*)cmalloc(numFFT2*sizeof(double))-1;
  double *rhoReduce;

  MPI_Comm commStates   =    commCP->comm_states;
  
  if(myidState==0){
    // Just use the first rhoRealGridTot grid points
    rhoReduce = (double*)cmalloc(numFFT2*sizeof(double));
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoUpChemPot[0][iGrid] = 0.0;
  }//endif myidState

/*==========================================================================*/
/* I) Reduce the densities                                             */

  for(iChem=0;iChem<numChemPot;iChem++){
    //debug
    /*
    if(myidState==densityMap[iChem]){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoUpChemPot[indexChemProc[iChem]][iGrid] = 0.0;
    }
    */
    //end debug
    rhoCalcRealStoHybrid(cpscr,cpcoeffs_info,
                   cell,stodftInfo,stoWfUpRe[iChem],
                   stoWfUpIm[iChem],rhoTemp,*coefFormUp,*coefOrthUp,
                   numStateUpProc,numCoeff,cpDualGridOptOn,commCP,
                   &(cp->cp_para_fft_pkg3d_lg),
                   &(cp->cp_sclr_fft_pkg3d_lg),
                   &(cp->cp_para_fft_pkg3d_dens_cp_box),
                   &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
                   &(cp->cp_sclr_fft_pkg3d_sm));
    //debug
    //end debug
    if(numProcStates>1){
      Barrier(commStates);
      Reduce(&rhoTemp[1],rhoReduce,rhoRealGridTot,MPI_DOUBLE,
             MPI_SUM,0,commStates);
      //MPI_Reduce(&rhoTemp[1],rhoUpChemPot[0],rhoRealGridTot,MPI_DOUBLE,MPI_SUM,densityMap[iChem],commStates);
      Barrier(commStates);
      if(myidState==0){
	for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	  rhoUpChemPot[0][iGrid] += rhoReduce[iGrid];
	}
	//printf("iChem %i rho[0] %lg\n",iChem,rhoUpChemPot[0][0]);
      }
      Barrier(commStates);
    }
    else{
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	rhoUpChemPot[0][iGrid] += rhoTemp[iGrid+1];
      }//endfor iGrid
    }
    //else memcpy(rhoUpChemPot[iChem],&rhoTemp[1],rhoRealGridTot*sizeof(double));
    //debug
    //printf("iChem %i rhoTemp %lg rhoUp %lg\n",iChem,rhoTemp[4001],rhoUpChemPot[iChem][4000]);
    // Calculate the average density, haven't scale by 1/volume
  }//endfor iChem

  free(&rhoTemp[1]);
  free(rhoReduce);
  
  if(numProcStates>1)Barrier(commStates);

/*==========================================================================*/
/* II) Scatter the densities                                             */

  if(cpParaOpt==0){
    if(numProcStates>1){
      Scatterv(rhoUpChemPot[0],rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
               rhoUpCorrect,rhoRealGridNum,MPI_DOUBLE,0,commStates);
    }
    else memcpy(rhoUpCorrect,rhoUpChemPot[0],rhoRealGridTot*sizeof(double));
  }

/*==========================================================================*/
/* IV) Output the density                                                   */

  outputDensity(cp,cell);

/*==========================================================================*/
/* V) Generate the diis density                                             */

  if(densityMixFlag>0){
    if(myidState==0){
      printf("Start Mixing Density\n");
      fflush(stdout);
    }
    genDensityMix(cp,iScf);
    if(myidState==0){
      printf("Finish Mixing Density\n");
      fflush(stdout);
    }
  }

/*==========================================================================*/
/* VI) Generate the reciprocal part and all the other things                */
  
  /*FILE *fileRhoReal = fopen("density-init","r");
  FILE *fileRhoRecip = fopen("density-recip-test","w");
  for(iGrid=1;iGrid<=rhoRealGridNum;iGrid++)fscanf(fileRhoReal,"%lg",&rhoUp[iGrid]);
  for(iGrid=1;iGrid<=rhoRealGridNum;iGrid++)rhoUp[iGrid] /= 0.0009250463018013585;
  */

  calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
		     rhoCoeffReUp,rhoCoeffImUp,rhoUp,rhoCoeffReUpDensCpBox,rhoCoeffImUpDensCpBox,
		     divRhoxUp,divRhoyUp,divRhozUp,d2RhoUp,cpGGA,cpDualGridOptOn,numInterpPmeDual,
		     commCP,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));
  if(cpLsda==1&&numStateDnProc!=0){
    calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
		       rhoCoeffReDn,rhoCoeffImDn,rhoDn,rhoCoeffReDnDensCpBox,rhoCoeffImDnDensCpBox,
		       divRhoxDn,divRhoyDn,divRhozDn,d2RhoDn,cpGGA,cpDualGridOptOn,numInterpPmeDual,
		       commCP,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));
    for(iCoeff=1;iCoeff<=numCoeffLargeProc;iCoeff++) {
      rhoCoeffReUp[iCoeff] += rhoCoeffReDn[iCoeff];
      rhoCoeffImUp[iCoeff] += rhoCoeffImDn[iCoeff];
    }/* endfor */
    if(cpDualGridOptOn>=1){
      for(iCoeff=1;iCoeff<=numCoeffLargeProcDensCpBox;iCoeff++){
	rhoCoeffReUpDensCpBox[iCoeff] += rhoCoeffReDnDensCpBox[iCoeff];
	rhoCoeffImUpDensCpBox[iCoeff] += rhoCoeffImDnDensCpBox[iCoeff];
      }/* endfor */
    } /* endif */
  }/* endif */

  /*
  for(iCoeff=1;iCoeff<=numCoeffLargeProc;iCoeff++){
    fprintf(fileRhoRecip,"%.10lg %.10lg\n",rhoCoeffReUp[iCoeff],rhoCoeffImUp[iCoeff]);
  }
  fclose(fileRhoReal);
  fclose(fileRhoRecip);  
  */
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcRhoStoHybridEnergyWindow(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
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
  FRAGINFO *fragInfo            = stodftInfo->fragInfo;
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  //PARA_FFT_PKG3D *cp_para_fft_pkg3d;

  int cpParaOpt = cpopts->cp_para_opt;
  int cpLsda = cpopts->cp_lsda;
  int cpGGA  = cpopts->cp_gga;
  int realSparseOpt = cpopts->realSparseOpt;
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
  int calcFragFlag      = stodftInfo->calcFragFlag;
  int rhoRealGridNum    = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot    = stodftInfo->rhoRealGridTot;
  int numChemProc       = stodftInfo->numChemProc;
  int numStateStoUp	= stodftInfo->numStateStoUp;
  int numStateStoDn	= stodftInfo->numStateStoDn;
  int occNumber		= stodftInfo->occNumber;
  int densityMixFlag	= stodftInfo->densityMixFlag;
  int iScf		= stodftInfo->iScf;
  int myidState		= commCP->myid_state;
  int numProcStates     = commCP->np_states;
  int iCoeff,iChem,iGrid,iState;
  int index;
  int i,j,k;

  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);
  int *densityMap   = stodftInfo->densityMap;
  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;

  double volCP,rvolCP;
  double numGridTotInv = 1.0/rhoRealGridTot;
  double aveFactUp = occNumber/(double)(numStateStoUp);
  double numElecTrue = stodftInfo->numElecTrue;
  double aveFactDn;
  double numElecTest;

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
  double *rhoUpCorrect    = stodftCoefPos->rhoUpCorrect;
  double *divRhoxUp       = cpscr->cpscr_grho.d_rhox_up;
  double *divRhoyUp       = cpscr->cpscr_grho.d_rhoy_up;
  double *divRhozUp       = cpscr->cpscr_grho.d_rhoz_up;
  double *d2RhoUp        = cpscr->cpscr_grho.d2_rho_up;
  double *rhoCoeffReDn   = cpscr->cpscr_rho.rhocr_dn;
  double *rhoCoeffImDn   = cpscr->cpscr_rho.rhoci_dn;
  double *rhoDn          = cpscr->cpscr_rho.rho_dn;
  double *rhoCoeffReDnDensCpBox = cpscr->cpscr_rho.rhocr_dn_dens_cp_box;
  double *rhoCoeffImDnDensCpBox = cpscr->cpscr_rho.rhoci_dn_dens_cp_box;
  double *rhoDnCorrect	  = stodftCoefPos->rhoDnCorrect;
  double *divRhoxDn       = cpscr->cpscr_grho.d_rhox_dn;
  double *divRhoyDn       = cpscr->cpscr_grho.d_rhoy_dn;
  double *divRhozDn       = cpscr->cpscr_grho.d_rhoz_dn;
  double *d2RhoDn        = cpscr->cpscr_grho.d2_rho_dn;
  double *numElectron = stodftCoefPos->numElectron;
  double *chemPot = stodftCoefPos->chemPot;
  double *rhoUpFragSum;
  double *rhoDnFragSum;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  double **rhoUpChemPot = stodftCoefPos->rhoUpChemPot;
  double **rhoDnChemPot = stodftCoefPos->rhoDnChemPot;

  double *rhoTemp  = (double*)cmalloc(numFFT2*sizeof(double))-1;
  double *rhoSum   = (double*)calloc(numFFT2,sizeof(double));

  MPI_Comm commStates   =    commCP->comm_states;
  
  // TEST
  double *wfDetReal = stodftCoefPos->wfDetReal;
  double *rhoLastWindow = (double*)calloc(numFFT2,sizeof(double));
  double *rhoLastWindowDet = (double*)calloc(numFFT2,sizeof(double));

/*==========================================================================*/
/* I) Generate density in real space for all chem potentials		    */

  //vol_cp  = getdeth(hmat_cp);
  //rvol_cp = 1.0/vol_cp;
  /*
  if(realSparseOpt==0){
    cp_para_fft_pkg3d = &(cp->cp_para_fft_pkg3d_lg);
  }
  else{
    cp_para_fft_pkg3d = &(cp->cp_para_fft_pkg3d_sparse);
  }
  */
  
  if(cpLsda==1&&numStateDnProc!=0)aveFactDn = occNumber/(double)(numStateStoDn);

  //printf("111111111 Finish Initialize Densities\n");
  //fflush(stdout);
  if(numProcStates>1)Barrier(commStates);

  //debug
  //for(iChem=0;iChem<numChemPot;iChem++)printf("myid %i densityMap %i\n",myidState,densityMap[iChem],indexChemProc[iChem]);
  //fflush(stdout);

  //if(numProcStates>1)Barrier(commStates);

  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoSum[iGrid] = 0.0;

  double test;
  for(iChem=0;iChem<numChemPot;iChem++){
    printf("!!!!!!!!!!!!! iChem %i\n",iChem);
    rhoCalcRealStoHybrid(cpscr,cpcoeffs_info,
		   cell,stodftInfo,stoWfUpRe[iChem],
		   stoWfUpIm[iChem],rhoTemp,*coefFormUp,*coefOrthUp,
		   numStateUpProc,numCoeff,cpDualGridOptOn,commCP,
		   &(cp->cp_para_fft_pkg3d_lg),
		   &(cp->cp_sclr_fft_pkg3d_lg),
		   &(cp->cp_para_fft_pkg3d_dens_cp_box),
		   &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
		   &(cp->cp_sclr_fft_pkg3d_sm));
    test = 0.0;
    printf("rhoTemp[1] %lg\n",rhoTemp[1]);
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      rhoSum[iGrid] += rhoTemp[iGrid+1];
      test += rhoTemp[iGrid+1];
    }
    //printf("tttttttest %lg\n",test);
    // TEST
    if(iChem==numChemPot-1){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
        rhoLastWindow[iGrid] = rhoTemp[iGrid+1];
      }
    }
  }
  printf("rhoSum[0] %lg\n",rhoSum[0]);
  if(numProcStates>1){
    Barrier(commStates);
    Reduce(&rhoSum[0],rhoUpChemPot[0],rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Barrier(commStates);
  }
  else memcpy(rhoUpChemPot[0],&rhoSum[0],rhoRealGridTot*sizeof(double));
  if(myidState==0)printf("rhoUp %lg\n",rhoUpChemPot[0][0]);
  if(cpLsda==1&&numStateDnProc!=0){
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoSum[iGrid] = 0.0;
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
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoSum[iGrid] += rhoTemp[iGrid+1];
    }
    if(numProcStates>1){
      Barrier(commStates);
      Reduce(&rhoSum[0],rhoDnChemPot[0],rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
      Barrier(commStates);
    }
    else memcpy(rhoDnChemPot[0],&rhoSum[0],rhoRealGridTot*sizeof(double));
  }
  if(numProcStates>1)Barrier(commStates);

  // Calculate the average density, haven't scale by 1/volume
  
  //for(i=0;i<10;i++)printf("i %i rhoUp %lg\n",i,rhoUpChemPot[0][i]);
  //printf("aveFactUp %lg occNumber %i\n",aveFactUp,occNumber);
  //printf("rhoRealGridNum %i rhoRealGridTot %i\n",rhoRealGridNum,rhoRealGridTot);

  if(myidState==0){
    //printf("2222222222222222 aveFactUp %lg rhoUpChemPot[0][0] %lg\n",aveFactUp,rhoUpChemPot[0][0]);
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      rhoUpChemPot[0][iGrid] *= aveFactUp;
      rhoLastWindow[iGrid] *= aveFactUp;
    }
    if(cpLsda==1&&numStateDnProc!=0){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoDnChemPot[0][iGrid] *= aveFactDn;
    }
    //printf("33333333333333 rhoUpChemPot[0][0] %lg\n",rhoUpChemPot[0][0]);
  }
  
  free(&rhoTemp[1]);
  free(&rhoSum[0]);

  // TEST
  // First substract the stochastic orbital from the last energy window
  
  int *ewStateNum = stodftCoefPos->ewStateNum;
  int **ewStateMap = stodftCoefPos->ewStateMap;
  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
    rhoUpChemPot[0][iGrid] -= rhoLastWindow[iGrid];
  }
  for(iState=0;iState<ewStateNum[numChemPot-1];iState++){
    index = ewStateMap[numChemPot-1][iState];
    printf("index %i %lg %lg\n",index,rhoLastWindowDet[0],wfDetReal[index*rhoRealGridTot]);
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      rhoLastWindowDet[iGrid] += wfDetReal[index*rhoRealGridTot+iGrid]*wfDetReal[index*rhoRealGridTot+iGrid];
    }
  }
  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoUpChemPot[0][iGrid] += rhoLastWindowDet[iGrid];
  
  if(numProcStates>1)Barrier(commStates);

/*==========================================================================*/
/* III) Double Check number of electrons.                                   */

  if(myidState==0){
    numElecTest = 0.0;
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
       numElecTest += rhoUpChemPot[0][iGrid];
       //printf("??????? numElecTest %lg %lg\n",rhoUpChemPot[0][iGrid],numElecTest);
    }
    if(cpLsda==1&&numStateDnProc!=0){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)numElecTest += rhoDnChemPot[0][iGrid];
    }
    printf("numElecTest %lg %lg\n",numElecTest,numGridTotInv);
    numElecTest *= numGridTotInv;
    printf("Number of electron is %.16lg.\n",numElecTest);
  }

  if(numProcStates>1)Barrier(commStates);


/*==========================================================================*/
/* III) Scatter the density to all nodes                                    */

  if(numProcStates>1){
    Scatterv(rhoUpChemPot[0],rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
             rhoUpCorrect,rhoRealGridNum,MPI_DOUBLE,0,commStates);
  }
  else{
    memcpy(rhoUpCorrect,rhoUpChemPot[0],rhoRealGridNum*sizeof(double));
  }
  if(cpLsda==1&&numStateDnProc!=0){
    if(numProcStates>1){
      Scatterv(rhoDnChemPot[0],rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
               rhoDnCorrect,rhoRealGridNum,MPI_DOUBLE,0,commStates);
    }
    else{
      memcpy(rhoDnCorrect,rhoDnChemPot[0],rhoRealGridNum*sizeof(double));
    }
  }

/*==========================================================================*/
/* IV) Add the fragment contribution.                                       */

  if(calcFragFlag==1){
    rhoUpFragSum =  fragInfo->rhoUpFragSum;
    //debug
    /*
    for(iProc=0;iProc<numProcStates;iProc++){
      if(myidState==iProc){
        for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)printf("rhostooooo %lg %lg %lg\n",rhoUpCorrect[iGrid],rhoUpFragSum[iGrid],rhoUpCorrect[iGrid]+rhoUpFragSum[iGrid]);
      }
      if(numProcStates>1)Barrier(commStates);
    }
    */
    //printf("rhoUpCorrect %lg %lg rhoUpFragSum %lg %lg\n",rhoUpCorrect[0],rhoUpCorrect[1],rhoUpFragSum[0],rhoUpFragSum[1]);   
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoUpCorrect[iGrid] += rhoUpFragSum[iGrid];
    if(cpLsda==1&&numStateDnProc!=0){
      rhoDnFragSum =  fragInfo->rhoDnFragSum;
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoDnCorrect[iGrid] += rhoDnFragSum[iGrid];
    }
  }

/*==========================================================================*/
/* IV) Output the density                                                   */

    outputDensity(cp,cell);

/*==========================================================================*/
/* V) Generate the diis density						    */

    if(densityMixFlag>0){
      if(myidState==0){
	printf("Start Mixing Density\n");
        fflush(stdout);
      }
      genDensityMix(cp,iScf);
      if(myidState==0){
	printf("Finish Mixing Density\n");
        fflush(stdout);
      }
    }
    
/*==========================================================================*/
/* VI) Generate the reciprocal part and all the other things                */
  /*
  FILE *fileRhoReal = fopen("density-init","r");
  FILE *fileRhoRecip = fopen("density-recip-test","w");
  for(iGrid=1;iGrid<=rhoRealGridNum;iGrid++)fscanf(fileRhoReal,"%lg",&rhoUp[iGrid]);
  for(iGrid=1;iGrid<=rhoRealGridNum;iGrid++)rhoUp[iGrid] /= 0.0009250463018013585;
  */
  calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
		     rhoCoeffReUp,rhoCoeffImUp,rhoUp,rhoCoeffReUpDensCpBox,rhoCoeffImUpDensCpBox,
		     divRhoxUp,divRhoyUp,divRhozUp,d2RhoUp,cpGGA,cpDualGridOptOn,numInterpPmeDual,
		     commCP,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));
  if(cpLsda==1&&numStateDnProc!=0){
    calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
		       rhoCoeffReDn,rhoCoeffImDn,rhoDn,rhoCoeffReDnDensCpBox,rhoCoeffImDnDensCpBox,
		       divRhoxDn,divRhoyDn,divRhozDn,d2RhoDn,cpGGA,cpDualGridOptOn,numInterpPmeDual,
		       commCP,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));
    for(iCoeff=1;iCoeff<=numCoeffLargeProc;iCoeff++) {
      rhoCoeffReUp[iCoeff] += rhoCoeffReDn[iCoeff];
      rhoCoeffImUp[iCoeff] += rhoCoeffImDn[iCoeff];
    }/* endfor */
    if(cpDualGridOptOn>=1){
      for(iCoeff=1;iCoeff<=numCoeffLargeProcDensCpBox;iCoeff++){
	rhoCoeffReUpDensCpBox[iCoeff] += rhoCoeffReDnDensCpBox[iCoeff];
	rhoCoeffImUpDensCpBox[iCoeff] += rhoCoeffImDnDensCpBox[iCoeff];
      }/* endfor */
    } /* endif */
  }/* endif */

  if(myidState==0){
    printf("Finish Calculating Reciprocal Space Density\n");
  }


  FILE *fileRhoRecip = fopen("rho-k","w");
  double *ak2 = cpewald->ak2;
  for(iCoeff=1;iCoeff<=numCoeffLargeProc;iCoeff++){
    fprintf(fileRhoRecip,"%.10lg %.10lg %.10lg\n",ak2[iCoeff],rhoCoeffReUp[iCoeff],rhoCoeffImUp[iCoeff]);
  }
  //fclose(fileRhoReal);
  fclose(fileRhoRecip);  
  

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


