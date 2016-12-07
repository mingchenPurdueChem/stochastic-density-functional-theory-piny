/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: wf-real-frag.c                               */
/*                                                                          */
/* Wrapper routine to FFT wave function to real space                       */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_frag_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void rhoRealCalcDriver(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *classMini,
			CP *cp)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
  CPEWALD *cpewald = &(cpMini->cpewald);
  CPSCR *cpscr = &(cpMini->cpscr);
  CPOPTS *cpopts = &(cpMini->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cpMini->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos = &(cpMini->cpcoeffs_pos[1]);
  COMMUNICATE *communicate = &(cpMini->communicate);
  EWALD *ewald = &(generalDataMini->ewald);
  CELL *cell = &(generalDataMini->cell);

  int cp_norb         = cpopts->cp_norb;
  int cp_lsda         = cpopts.cp_lsda;
  int myid_state      = communicate.myid_state;



  int *icoef_orth_up    = &(cp->cpcoeffs_pos[ip_now].icoef_orth_up);
  int *icoef_form_up    = &(cp->cpcoeffs_pos[ip_now].icoef_form_up);
  int *ifcoef_orth_up   = &(cp->cpcoeffs_pos[ip_now].ifcoef_orth_up);
  int *ifcoef_form_up   = &(cp->cpcoeffs_pos[ip_now].ifcoef_form_up);

  double *ccrealUp        = cp->cpcoeffs_pos[ip_now].cre_up;
  double *ccimagUp        = cp->cpcoeffs_pos[ip_now].cim_up;

/*======================================================================*/
/* I) Check the forms                                                   */

  if(cp_norb>0){
    if((*icoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up Coefs must be in nonorthonormal form under norb \n");
      printf("on state processor %d in cp_elec_energy_ctrl \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
      if((*icoef_orth_dn)!=0){
	printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	printf("Dn Coefs must be in nonorthonormal form under norb \n");
	printf("on state processor %d in cp_elec_energy_ctrl \n",myid_state);
	printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	fflush(stdout);
	exit(1);
      }/*endif*/
    }/*endif*/
  }/*endif*/

/*======================================================================*/
/* III) Initialize Flags, inverse hmat			                */

  (*ifcoef_form_up) = 0;
  (*ifcoef_orth_up) = 1;
  if(cp_lsda==1&&nstate_dn!=0){  
    (*ifcoef_form_dn) = 0;
    (*ifcoef_orth_dn) = 1;
  }
  for(i=1;i<=9;i++){ptens_pvten_tmp[i] = 0.0;}
  gethinv(cell->hmat_cp,cell->hmati_cp,&(cell->vol_cp),iperd);
  gethinv(cell->hmat,cell->hmati,&(cell->vol),iperd);

/*======================================================================*/
/* IV) Allocate memory for big system                                   */


/*======================================================================*/
/* IV) Calculate real space wave functions and densities                */

  rho_scr = rho;

/* ================================================================= */
/*1) zero density and gradients if necessary                         */

  for(i=1;i<= nfft2;i++){rho_scr[i] = 0.0;
 


/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void rhoRealCalcWrapper(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *classMini,
                        CP *cp,double *ccreal, double *ccimag,int icoef_form,int icoef_orth,
			double *rho,double *wfReal,int nstate)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
  CPEWALD *cpewald = &(cpMini->cpewald);
  CPSCR *cpscr = &(cpMini->cpscr);
  CPOPTS *cpopts = &(cpMini->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cpMini->cpcoeffs_info);
  COMMUNICATE *communicate = &(cpMini->communicate);
  EWALD *ewald = &(generalDataMini->ewald);
  CELL *cell = &(generalDataMini->cell);

  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm = &(cpMini->cp_sclr_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_lg = &(cpMini->cp_sclr_fft_pkg3d_lg);

  int iii,ioff,ioff2;
  int is,i,j,k,iupper;
  int igrid;
  int ncoef = cpcoeffs_info->ncoef;
  int   nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
  int   nfft2_proc       =    nfft_proc/2;

  double *zfft           =    cpscr->cpscr_wave.zfft;

/* ================================================================= */
/*1) zero density and gradients if necessary                         */

  for(i=1;i<=nfft2;i++)rho[i] = 0.0;

/*==========================================================================*/
/*==========================================================================*/
/*2)sum the density in real space two states at a time                      */
/*  This is done at state level and uses the scalar packages!               */

  iupper = nstate;
  if(nstate%2!=0){
     iupper = nstate-1;
  }/* endif */

  for(is=1;is<=iupper;is=is+2){
    ioff   = (is-1)*ncoef;
    ioff2 = (is)*ncoef;

/*--------------------------------------------------------------------------*/
/*I) pack the complex zfft array with two wavefunctions (real) 
    the wavefunctions are reperesented in spherically cuttof 
    half g space                                                            */

    dble_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
                      zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*II) fourier transform the two wavefunctions to real space
     convention exp(-igr)                                                   */

    para_fft_gen3d_fwd_to_r_fftw3d(zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* III) Copy the real sapce wave function and add the square of the two     
        wave functions to the density(real space)			    */
 
    for(igrid=0;igrid<nfft2_proc;igrid++){
      wfReal[ioff*nfft2+igrid] = zfft[igrid*2+1];
      wfReal[ioff2*nfft2+igrid] = zfft[igrid*2+2];
      rho[igrid] += zfft[igrid*2+1]*zfft[igrid*2+1]+zfft[igrid*2+2]*zfft[igrid*2+2];
    }
  }/*endfor is*/

/*--------------------------------------------------------------------------*/
/*IV) if there is an odd number of states, use single pack
      for the last state (i caen't be chaengin' the laws of physics
      captn, i've got to be usin the single pack!!!!!!)                     */

  if(nstate%2!=0){
    ioff = (nstate-1)*ncoef;
    sngl_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* V) fourier transform the last wavefunction to real space
       convention exp(-igr)                                                 */
    
    para_fft_gen3d_fwd_to_r_fftw3d(zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*VI) Copy the real sapce wave function and add the square of the last wave 
      function to the density(real space)   */
    
    for(igrid=0;igrid<nfft2_proc;igrid++){
      wfReal[ioff*nfft2+igrid] = zfft[igird*2+1];
      rho[igrid] += zfft[igrid*2+1]*zfft[igrid*2+1];
    }

  }//endif nstat%2

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

