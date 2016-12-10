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
void rhoRealCalcDriverFrag(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *classMini,
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
  int cpLsda         = cpopts.cp_lsda;
  int myidState      = communicate.myid_state;
  int numStateUpProc        = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc        = cpcoeffs_info->nstate_dn_proc;

  int *icoef_orth_up    = &(cpMini->cpcoeffs_pos[1].icoef_orth_up);
  int *icoef_form_up    = &(cpMini->cpcoeffs_pos[1].icoef_form_up);
  int *ifcoef_orth_up   = &(cpMini->cpcoeffs_pos[1].ifcoef_orth_up);
  int *ifcoef_form_up   = &(cpMini->cpcoeffs_pos[1].ifcoef_form_up);
  int *icoef_orth_dn    = &(cpMini->cpcoeffs_pos[1].icoef_orth_dn);
  int *icoef_form_dn    = &(cpMini->cpcoeffs_pos[1].icoef_form_dn);
  int *ifcoef_orth_dn   = &(cpMini->cpcoeffs_pos[1].ifcoef_orth_dn);
  int *ifcoef_form_dn   = &(cpMini->cpcoeffs_pos[1].ifcoef_form_dn);


  double *ccrealUp        = cp->cpcoeffs_pos[ip_now].cre_up;
  double *ccimagUp        = cp->cpcoeffs_pos[ip_now].cim_up;
  double *ccrealDn	  = cp->cpcoeffs_pos[ip_now].cim_dn;
  double *ccimagDn        = cp->cpcoeffs_pos[ip_now].cim_dn;

  double *ccrealUpMini    = cpMini->cpcoeffs_pos[1].cre_up;
  double *ccimagUpMini    = cpMini->cpcoeffs_pos[1].cim_up;
  double *ccrealDnMini    = cpMini->cpcoeffs_pos[1].cre_dn;
  double *ccimagDnMini    = cpMini->cpcoeffs_pos[1].cim_dn;
  double *rhoUpFragProc   = fragInfo->rhoUpFragProc[iFrag];
  double *rhoDnFragProc	  = fragInfo->rhoDnFragProc[iFrag];
  double *coefUpFragProc  = fragInfo->coefUpFragProc[iFrag];
  double *coefDnFragProc  = fragInfo->coefDnFragProc[iFrag];

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
    if(cpLsda==1){
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
  if(cpLsda==1&&numStateDnProc!=0){  
    (*ifcoef_form_dn) = 0;
    (*ifcoef_orth_dn) = 1;
  }
  for(i=1;i<=9;i++){ptens_pvten_tmp[i] = 0.0;}
  gethinv(cell->hmat_cp,cell->hmati_cp,&(cell->vol_cp),iperd);
  gethinv(cell->hmat,cell->hmati,&(cell->vol),iperd);


/*======================================================================*/
/* IV) Calculate real space wave functions and densities for fragments  */

  rhoRealCalcWrapper(generalDataMini,cpMini,classMini,
                     cp,ccrealUpMini,ccimagUpMini,icoef_form_up,icoef_orth_up,
                     rhoUpFragProc,coefUpFragProc,numStateUpFrag);
  if(cpLsda==1&&numStateDnFrag!=0){
    rhoRealCalcWrapper(generalDataMini,cpMini,classMini,
		       cp,ccrealDnMini,ccimagDnMini,icoef_form_dn,icoef_orth_dn,
		       rhoDnFragProc,coefDnFragProc,numStateDnFrag);    
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void rhoRealCalcDriverNoise(GENERAL_DATA *general_data,CP *cp,CLASS *class,int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  CPEWALD *cpewald = &(cp->cpewald);
  CPSCR *cpscr = &(cp->cpscr);
  CPOPTS *cpopts = &(cp->cpopts);
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  CPCOEFFS_INFO *cpcoeffs_info = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos = &(cp->cpcoeffs_pos[ip_now]);
  COMMUNICATE *communicate = &(cp->communicate);
  EWALD *ewald = &(general_data->ewald);
  CELL *cell = &(general_data->cell);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm = &(cp->cp_sclr_fft_pkg3d_sm);

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
  double *noiseWfUpReal  = fragInfo->noiseWfUpReal;
  double *noiseWfDnReal  = fragInfo->noiseWfDnReal;

  double *creal = cp->cpcoeffs_pos[ip_now].cre_up;
  double *cimag = cp->cpcoeffs_pos[ip_now].cim_up;

  MPI_Comm comm_states   =    communicate->comm_states;

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

