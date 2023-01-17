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
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_frag_local.h"
#include "../proto_defs/proto_stodft_local.h"
#include "../proto_defs/proto_vps_params_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void rhoRealCalcDriverFragMol(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *classMini,
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
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;

  int i;
  int cp_norb         = cpopts->cp_norb;
  int cpLsda         = cpopts->cp_lsda;
  int myidState      = communicate->myid_state;
  int numStateUpFrag    = cpcoeffs_info->nstate_up_proc;
  int numStateDnFrag    = cpcoeffs_info->nstate_dn_proc;
  int iFrag		= fragInfo->iFrag;

  int *icoef_orth_up    = &(cpMini->cpcoeffs_pos[1].icoef_orth_up);
  int *icoef_form_up    = &(cpMini->cpcoeffs_pos[1].icoef_form_up);
  int *ifcoef_orth_up   = &(cpMini->cpcoeffs_pos[1].ifcoef_orth_up);
  int *ifcoef_form_up   = &(cpMini->cpcoeffs_pos[1].ifcoef_form_up);
  int *icoef_orth_dn    = &(cpMini->cpcoeffs_pos[1].icoef_orth_dn);
  int *icoef_form_dn    = &(cpMini->cpcoeffs_pos[1].icoef_form_dn);
  int *ifcoef_orth_dn   = &(cpMini->cpcoeffs_pos[1].ifcoef_orth_dn);
  int *ifcoef_form_dn   = &(cpMini->cpcoeffs_pos[1].ifcoef_form_dn);

  double *ccrealUpMini    = cpMini->cpcoeffs_pos[1].cre_up;
  double *ccimagUpMini    = cpMini->cpcoeffs_pos[1].cim_up;
  double *ccrealDnMini    = cpMini->cpcoeffs_pos[1].cre_dn;
  double *ccimagDnMini    = cpMini->cpcoeffs_pos[1].cim_dn;
  double *rhoUpFragProc,*rhoDnFragProc,*coefUpFragProc,*coefDnFragProc;

/*======================================================================*/
/* I) Check the forms                                                   */

  if(cp_norb>0){
    if((*icoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up Coefs must be in nonorthonormal form under norb \n");
      printf("on state processor %d in cp_elec_energy_ctrl \n",myidState);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cpLsda==1){
      if((*icoef_orth_dn)!=0){
	printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	printf("Dn Coefs must be in nonorthonormal form under norb \n");
	printf("on state processor %d in cp_elec_energy_ctrl \n",myidState);
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
  rhoUpFragProc = fragInfo->rhoUpFragProc[iFrag];
  coefUpFragProc = fragInfo->coefUpFragProc[iFrag];
  if(cpLsda==1&&numStateDnFrag!=0){  
    (*ifcoef_form_dn) = 0;
    (*ifcoef_orth_dn) = 1;
    rhoDnFragProc = fragInfo->rhoDnFragProc[iFrag];
    coefDnFragProc = fragInfo->coefDnFragProc[iFrag];
  }
  //gethinv(cell->hmat_cp,cell->hmati_cp,&(cell->vol_cp),iperd);
  //gethinv(cell->hmat,cell->hmati,&(cell->vol),iperd);

  //debug
  /*
  if(stodftInfo->stodftOn!=-1){
    char name[100];
    FILE *fwf;
    int nstate_up = cpcoeffs_info->nstate_up_proc;
    int ncoef = cpcoeffs_info->ncoef;
    int numProcStatesSys = cp->communicate.np_states; 
    int myidStateSys = cp->communicate.myid_state;
    int numFragTot = fragInfo->numFragTot;
    int res = numFragTot%numProcStatesSys;
    int div = numFragTot/numProcStatesSys;
    int indStart,indFrag;
    if(myidStateSys<res)indStart = (div+1)*myidStateSys;
    else indStart = (myidStateSys-res)*div+res*(div+1);
    indFrag = indStart+iFrag;
    //printf("res %i div %i indStart %i indFrag %i\n",res,div,indStart,indFrag);
    sprintf(name,"fragwf-%i",indFrag);
    //fwf = fopen(name,"w");
    fwf = fopen(name,"r");
    for(i=1;i<=nstate_up*ncoef;i++){
      //fprintf(fwf,"%.16lg %.16lg\n",ccrealUpMini[i],ccimagUpMini[i]);
      fscanf(fwf,"%lg",&ccrealUpMini[i]);
      fscanf(fwf,"%lg",&ccimagUpMini[i]);
    }
    fclose(fwf);
  }
  */

/*======================================================================*/
/* IV) Calculate real space wave functions and densities for fragments  */

  rhoRealCalcFragWrapper(generalDataMini,cpMini,classMini,
                     cp,ccrealUpMini,ccimagUpMini,icoef_form_up,icoef_orth_up,
                     rhoUpFragProc,coefUpFragProc,numStateUpFrag);
  if(cpLsda==1&&numStateDnFrag!=0){
    rhoRealCalcFragWrapper(generalDataMini,cpMini,classMini,
		       cp,ccrealDnMini,ccimagDnMini,icoef_form_dn,icoef_orth_dn,
		       rhoDnFragProc,coefDnFragProc,numStateDnFrag);    
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void rhoRealCalcDriverFragUnitCell(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *classMini,
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
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;

  int iState,iGrid;
  int cp_norb         = cpopts->cp_norb;
  int cpLsda         = cpopts->cp_lsda;
  int myidState      = communicate->myid_state;
  int numStateUpFrag    = cpcoeffs_info->nstate_up_proc;
  int numStateDnFrag    = cpcoeffs_info->nstate_dn_proc;
  int iFrag		= fragInfo->iFrag;
  int numGridFragSmall	= fragInfo->numGridFragProcSmall[iFrag];
  int *gridMapSmall	= fragInfo->gridMapProcSmall[iFrag];
  int numGridFrag	= fragInfo->numGridFragProc[iFrag];

  int *icoef_orth_up    = &(cpMini->cpcoeffs_pos[1].icoef_orth_up);
  int *icoef_form_up    = &(cpMini->cpcoeffs_pos[1].icoef_form_up);
  int *ifcoef_orth_up   = &(cpMini->cpcoeffs_pos[1].ifcoef_orth_up);
  int *ifcoef_form_up   = &(cpMini->cpcoeffs_pos[1].ifcoef_form_up);
  int *icoef_orth_dn    = &(cpMini->cpcoeffs_pos[1].icoef_orth_dn);
  int *icoef_form_dn    = &(cpMini->cpcoeffs_pos[1].icoef_form_dn);
  int *ifcoef_orth_dn   = &(cpMini->cpcoeffs_pos[1].ifcoef_orth_dn);
  int *ifcoef_form_dn   = &(cpMini->cpcoeffs_pos[1].ifcoef_form_dn);

  double *ccrealUpMini    = cpMini->cpcoeffs_pos[1].cre_up;
  double *ccimagUpMini    = cpMini->cpcoeffs_pos[1].cim_up;
  double *ccrealDnMini    = cpMini->cpcoeffs_pos[1].cre_dn;
  double *ccimagDnMini    = cpMini->cpcoeffs_pos[1].cim_dn;
  double *rhoUpFragProc,*rhoDnFragProc,*coefUpFragProc,*coefDnFragProc;
  double *wfRealUpSmall,*wfRealDnSmall,*rhoFake;

/*======================================================================*/
/* I) Check the forms                                                   */
  if(cp_norb>0){
    if((*icoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up Coefs must be in nonorthonormal form under norb \n");
      printf("on state processor %d in cp_elec_energy_ctrl \n",myidState);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cpLsda==1){
      if((*icoef_orth_dn)!=0){
	printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	printf("Dn Coefs must be in nonorthonormal form under norb \n");
	printf("on state processor %d in cp_elec_energy_ctrl \n",myidState);
	printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	fflush(stdout);
	exit(1);
      }/*endif*/
    }/*endif*/
  }/*endif*/
  wfRealUpSmall = (double*)cmalloc(numStateUpFrag*numGridFragSmall*sizeof(double));
  if(cpLsda==1&&numStateDnFrag!=0){
    wfRealDnSmall = (double*)cmalloc(numStateUpFrag*numGridFragSmall*sizeof(double));
  }
  //fake density
  rhoFake = (double*)cmalloc(numGridFragSmall*sizeof(double));

/*======================================================================*/
/* III) Initialize Flags, inverse hmat			                */

  (*ifcoef_form_up) = 0;
  (*ifcoef_orth_up) = 1;
  rhoUpFragProc = fragInfo->rhoUpFragProc[iFrag];
  coefUpFragProc = fragInfo->coefUpFragProc[iFrag];
  if(cpLsda==1&&numStateDnFrag!=0){  
    (*ifcoef_form_dn) = 0;
    (*ifcoef_orth_dn) = 1;
    rhoDnFragProc = fragInfo->rhoDnFragProc[iFrag];
    coefDnFragProc = fragInfo->coefDnFragProc[iFrag];
  }
  //gethinv(cell->hmat_cp,cell->hmati_cp,&(cell->vol_cp),iperd);
  //gethinv(cell->hmat,cell->hmati,&(cell->vol),iperd);

  for(iGrid=0;iGrid<numGridFrag;iGrid++)rhoUpFragProc[iGrid] = 0.0;
  if(cpLsda==1&&numStateDnFrag!=0){
    for(iGrid=0;iGrid<numGridFrag;iGrid++)rhoDnFragProc[iGrid] = 0.0;
  }
  for(iState=0;iState<numStateUpFrag;iState++){
    for(iGrid=0;iGrid<numGridFrag;iGrid++)coefUpFragProc[iState*numGridFrag+iGrid] = 0.0;
  }
  if(cpLsda==1&&numStateDnFrag!=0){
    for(iState=0;iState<numStateDnFrag;iState++){
      for(iGrid=0;iGrid<numGridFrag;iGrid++)coefDnFragProc[iState*numGridFrag+iGrid] = 0.0;
    }  
  }

/*======================================================================*/
/* IV) Calculate real space wave functions and densities for fragments  */

  rhoRealCalcFragWrapper(generalDataMini,cpMini,classMini,
                     cp,ccrealUpMini,ccimagUpMini,icoef_form_up,icoef_orth_up,
                     rhoFake,wfRealUpSmall,numStateUpFrag);
  if(cpLsda==1&&numStateDnFrag!=0){
    rhoRealCalcFragWrapper(generalDataMini,cpMini,classMini,
		       cp,ccrealDnMini,ccimagDnMini,icoef_form_dn,icoef_orth_dn,
		       rhoFake,wfRealDnSmall,numStateDnFrag);    
  }

/*======================================================================*/
/* V) Embed the real space wave functions to the bigger cell	        */
  
  //embedWfReal(generalDataMini,cpMini,classMini,cp,wfRealUpSmall,wfRealDnSmall); 

/*======================================================================*/
/* VI) Clean up local memory allocation			                */

  free(wfRealUpSmall);
  if(cpLsda==1&&numStateDnFrag!=0)free(wfRealDnSmall);
  free(rhoFake);


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
  CPCOEFFS_INFO *cpcoeffs_info = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos = &(cp->cpcoeffs_pos[1]);
  COMMUNICATE *communicate = &(cp->communicate);
  EWALD *ewald = &(general_data->ewald);
  CELL *cell = &(general_data->cell);
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;

  int cp_norb         = cpopts->cp_norb;
  int cpLsda          = cpopts->cp_lsda;
  int myidState       = communicate->myid_state;
  int numStateUpProc  = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc  = cpcoeffs_info->nstate_dn_proc;

  double *ccrealUp        = cp->cpcoeffs_pos[ip_now].cre_up;
  double *ccimagUp        = cp->cpcoeffs_pos[ip_now].cim_up;
  double *ccrealDn        = cp->cpcoeffs_pos[ip_now].cim_dn;
  double *ccimagDn        = cp->cpcoeffs_pos[ip_now].cim_dn;

  double *noiseWfUpReal = fragInfo->noiseWfUpReal;
  double *noiseWfDnReal = fragInfo->noiseWfDnReal;

  int *icoef_orth_up    = &(cp->cpcoeffs_pos[1].icoef_orth_up);
  int *icoef_form_up    = &(cp->cpcoeffs_pos[1].icoef_form_up);
  int *ifcoef_orth_up   = &(cp->cpcoeffs_pos[1].ifcoef_orth_up);
  int *ifcoef_form_up   = &(cp->cpcoeffs_pos[1].ifcoef_form_up);
  int *icoef_orth_dn    = &(cp->cpcoeffs_pos[1].icoef_orth_dn);
  int *icoef_form_dn    = &(cp->cpcoeffs_pos[1].icoef_form_dn);
  int *ifcoef_orth_dn   = &(cp->cpcoeffs_pos[1].ifcoef_orth_dn);
  int *ifcoef_form_dn   = &(cp->cpcoeffs_pos[1].ifcoef_form_dn);


/* ================================================================= */
/*0) Check the form of the coefficients                              */

/*======================================================================*/
/* IV) Calculate real space wave functions				*/
  //debug

  rhoRealCalcWrapper(general_data,cp,class,ccrealUp,ccimagUp,icoef_form_up,
		     icoef_orth_up,noiseWfUpReal,numStateUpProc);
  if(cpLsda==1&&numStateDnProc!=0){
    rhoRealCalcWrapper(general_data,cp,class,ccrealUp,ccimagUp,icoef_form_dn,
                       icoef_orth_dn,noiseWfDnReal,numStateDnProc);
  }
    
  /*
  int iState;
  int rhoRealGridTot = stodftInfo->rhoRealGridTot;
  for(iState=0;iState<numStateUpProc;iState++){
    printf("noissss %lg\n",noiseWfUpReal[iState*rhoRealGridTot]);
  }
  */

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void rhoRealCalcWrapper(GENERAL_DATA *general_data,CP *cp,CLASS *class,
			double *ccreal,double *ccimag,int *icoef_form,int *icoef_orth,
			double *wfReal,int nstate)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
  CPEWALD *cpewald = &(cp->cpewald);
  CPSCR *cpscr = &(cp->cpscr);
  CPOPTS *cpopts = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cp->cpcoeffs_info);
  COMMUNICATE *communicate = &(cp->communicate);
  EWALD *ewald = &(general_data->ewald);
  CELL *cell = &(general_data->cell);

  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm = &(cp->cp_sclr_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

  int iii,ioff,ioff2;
  int is,i,j,k,iupper;
  int gridoff,gridoff2;
  int igrid;
  int ncoef = cpcoeffs_info->ncoef;
  int   nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
  int   nfft2_proc       =    nfft_proc/2;
  int nfft = cp_para_fft_pkg3d_lg->nfft;
  int nfft2 = nfft/2;

  double *zfft           =    cpscr->cpscr_wave.zfft;
  double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
  double *hmatCP         =    cell->hmat_cp;
  double volCP		 = getdeth(hmatCP);
  double invVolCP	 = 1.0/invVolCP;


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
    gridoff = (is-1)*nfft2;
    gridoff2 = is*nfft2;

/*--------------------------------------------------------------------------*/
/*I) pack the complex zfft array with two wavefunctions (real) 
    the wavefunctions are reperesented in spherically cuttof 
    half g space                                                            */

    //printf("creal1 %lg creal2 %lg\n",ccreal[ioff+1],ccreal[ioff2+1]);
    dble_pack_coef(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
                      zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*II) fourier transform the two wavefunctions to real space
     convention exp(-igr)                                                   */

    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
    //printf("zfft %lg\n",zfft[1]);

/*--------------------------------------------------------------------------*/
/* III) Copy the real sapce wave function and add the square of the two     
        wave functions to the density(real space)                           */

    for(igrid=0;igrid<nfft2;igrid++){
      wfReal[gridoff+igrid] = zfft[igrid*2+1];
      wfReal[gridoff2+igrid] = zfft[igrid*2+2];
    }
    //printf("wfReal %lg %lg\n",wfReal[(is-1)*nfft2_proc],wfReal[is*nfft2_proc]);
  }/*endfor is*/

/*--------------------------------------------------------------------------*/
/*IV) if there is an odd number of states, use single pack
      for the last state (i caen't be chaengin' the laws of physics
      captn, i've got to be usin the single pack!!!!!!)                     */

  if(nstate%2!=0){
    ioff = (nstate-1)*ncoef;
    gridoff = (nstate-1)*nfft2;
    sngl_pack_coef(&ccreal[ioff],&ccimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* V) fourier transform the last wavefunction to real space
       convention exp(-igr)                                                 */

    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*VI) Copy the real sapce wave function and add the square of the last wave 
      function to the density(real space)   */

    for(igrid=0;igrid<nfft2;igrid++){
      wfReal[gridoff+igrid] = zfft[igrid*2+1];
    }

  }//endif nstat%2
  /*
  printf("nfft2 %i\n",nfft2);
  for(is=0;is<nstate;is++){
    printf("wwwfReal %lg\n",wfReal[is*nfft2_proc]);
  }
  */

  /*
  for(is=0;is<nstate;is++){
    char fname[100]; 
    sprintf(fname,"wf-stodft-%i",is);
    FILE *fout = fopen(fname,"w");
    for(igrid=0;igrid<nfft2;igrid++)fprintf(fout,"%i %.16lg\n",igrid,wfReal[is*nfft2+igrid]);
    fclose(fout);
  }
  */

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void rhoRealCalcFragWrapper(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *classMini,
                        CP *cp,double *ccreal, double *ccimag,int *icoef_form,int *icoef_orth,
			double *rho,double *wfReal,int nstate)
/*========================================================================*/
{/*begin routine*/
/* Fragment calculation uses 3D FFT. So I need to write a sperate wrapper */
/* For it.								  */
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
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cpMini->cp_para_fft_pkg3d_lg);

  int iii,ioff,ioff2;
  int is,i,j,k,iupper;
  int igrid,jgrid;
  int cp_lsda = cpopts->cp_lsda;
  int ncoef = cpcoeffs_info->ncoef;
  int   nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
  int   nfft2_proc       =    nfft_proc/2;
  int nfft = cp_para_fft_pkg3d_lg->nfft;
  int nfft2 = nfft/2;
  int nkf1       = cp_sclr_fft_pkg3d_sm->nkf1;
  int nkf2       = cp_sclr_fft_pkg3d_sm->nkf2;
  int nkf3       = cp_sclr_fft_pkg3d_sm->nkf3;
  int numThreads = cp_para_fft_pkg3d_lg->numThreads;
  int iThread;
  int gridOff1,gridOff2;
  int fragDFTMethod = cp->stodftInfo->fragDFTMethod;

  //printf("nfftttttt2_proc %i\n",nfft2_proc);
  double *zfft           =    cpscr->cpscr_wave.zfft;
  double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
  double *hmatCP         =    cell->hmat_cp;
  double *numOccDetProc;
  double *occPre;
  double volFrag         = getdeth(hmatCP);
  double invVolFrag	 = 1.0/volFrag;
  double prefact;

  double **zfft_threads  =    cpscr->cpscr_wave.zfft_threads;
  double **zfft_tmp_threads = cpscr->cpscr_wave.zfft_tmp_threads;
  double *rho_scr_threads;

  if(cp->stodftInfo->calcLocalTraceOpt == 1){
    prefact = 1.0/sqrt(volFrag);
  }
  else {
    prefact = 1.0/sqrt(2.0*volFrag);
  }
  //prefact = 1.0/sqrt(2.0*volFrag);
  if(cp_lsda==1)prefact = sqrt(invVolFrag);
  rho_scr_threads = (double*)cmalloc((numThreads*nfft2+1)*sizeof(double));

  // DEBUG NOT WORK FOR SPIN DOWN
  if(fragDFTMethod==2){
    numOccDetProc = cpMini->stodftInfo->numOccDetProc;
  }
  occPre = (double*)cmalloc(nstate*sizeof(double));
  switch(fragDFTMethod){
    case 1:
      for(is=0;is<nstate;is++)occPre[is] = 1.0;
      break;
    case 2:
      for(is=0;is<nstate;is++){
        occPre[is] = numOccDetProc[is]/sqrt(2.0);
        //printf("is %i occPre %lg\n",is,occPre[is]);
      }
      break;
  }
  

/* ================================================================= */
/*1) zero density and gradients if necessary                         */

  for(i=0;i<nfft2;i++)rho[i] = 0.0;
  for(i=1;i<=numThreads*nfft2;i++){
    rho_scr_threads[i] = 0.0;
  }

/*==========================================================================*/
/*==========================================================================*/
/*2)sum the density in real space two states at a time                      */
/*  This is done at state level and uses the scalar packages!               */

  iupper = nstate;
  if(nstate%2!=0){
     iupper = nstate-1;
  }/* endif */

  omp_set_num_threads(numThreads);
  #pragma omp parallel private(iThread,is,ioff,ioff2,gridOff1,gridOff2,i,j,k,igrid,jgrid)
  {
    iThread = omp_get_thread_num();
    #pragma omp for
    for(is=1;is<=iupper;is=is+2){
      ioff   = (is-1)*ncoef;
      ioff2 = (is)*ncoef;
      gridOff1 = (is-1)*nfft2_proc;
      gridOff2 = is*nfft2_proc;

/*--------------------------------------------------------------------------*/
/*I) pack the complex zfft array with two wavefunctions (real) 
    the wavefunctions are reperesented in spherically cuttof 
    half g space                                                            */

      //dble_pack_coef(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
      //                  zfft,cp_sclr_fft_pkg3d_sm);
      //printf("is %i %i ccreal %lg %lg\n",is-1,is,ccreal[ioff+2],ccreal[ioff2+2]);
      dble_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
			zfft_threads[iThread],cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*II) fourier transform the two wavefunctions to real space
     convention exp(-igr)                                                   */

      para_fft_gen3d_fwd_to_r_fftw3d(zfft_threads[iThread],cp_sclr_fft_pkg3d_sm,iThread);
      //para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* III) Copy the real sapce wave function and add the square of the two     
        wave functions to the density(real space)			    */
      //transpose zfft_threads xleading
      memcpy(&zfft_tmp_threads[iThread][1],&zfft_threads[iThread][1],nfft*sizeof(double));
      for(i=0;i<nkf1;i++){
	for(j=0;j<nkf2;j++){
	  for(k=0;k<nkf3;k++){
	    igrid = i*nkf2*nkf3+j*nkf3+k; //x leading
	    jgrid = k*nkf2*nkf1+j*nkf1+i; //z leading
	    zfft_threads[iThread][jgrid*2+1] = zfft_tmp_threads[iThread][igrid*2+1];
	    zfft_threads[iThread][jgrid*2+2] = zfft_tmp_threads[iThread][igrid*2+2];
	  }
	}
      }
   
      for(igrid=0;igrid<nfft2_proc;igrid++){
	wfReal[gridOff1+igrid] = zfft_threads[iThread][igrid*2+1]*prefact;
	wfReal[gridOff2+igrid] = zfft_threads[iThread][igrid*2+2]*prefact;
	rho_scr_threads[iThread*nfft2+igrid] 
		+= zfft_threads[iThread][igrid*2+1]*zfft_threads[iThread][igrid*2+1]*occPre[is-1]+
		   zfft_threads[iThread][igrid*2+2]*zfft_threads[iThread][igrid*2+2]*occPre[is];
      }
      //printf("ioff %i ioff2 %i wfReal %lg %lg\n",ioff,ioff2,wfReal[ioff],wfReal[ioff2]);
    }//endfor is
  }//endomp

  for(iThread=0;iThread<numThreads;iThread++){
    #pragma omp parallel for private(igrid)
    for(igrid=0;igrid<nfft2;igrid++){
      rho[igrid] += rho_scr_threads[iThread*nfft2+igrid];
    }
  }

/*--------------------------------------------------------------------------*/
/*IV) if there is an odd number of states, use single pack
      for the last state (i caen't be chaengin' the laws of physics
      captn, i've got to be usin the single pack!!!!!!)                     */

  if(nstate%2!=0){
    ioff = (nstate-1)*ncoef;
    sngl_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],zfft_threads[0],cp_sclr_fft_pkg3d_sm);
    //sngl_pack_coef(&ccreal[ioff],&ccimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* V) fourier transform the last wavefunction to real space
       convention exp(-igr)                                                 */
    
    para_fft_gen3d_fwd_to_r_fftw3d(zfft_threads[0],cp_sclr_fft_pkg3d_sm,0);
    //para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*VI) Copy the real sapce wave function and add the square of the last wave 
      function to the density(real space)   */

    memcpy(&zfft_tmp_threads[0][1],&zfft_threads[0][1],nfft*sizeof(double));
    for(i=0;i<nkf1;i++){
      for(j=0;j<nkf2;j++){
	for(k=0;k<nkf3;k++){
	  igrid = i*nkf2*nkf3+j*nkf3+k; //x leading
	  jgrid = k*nkf2*nkf1+j*nkf1+i; //z leading
	  zfft_threads[0][jgrid*2+1] = zfft_tmp_threads[0][igrid*2+1];
	  zfft_threads[0][jgrid*2+2] = zfft_tmp_threads[0][igrid*2+2];
	}
      }
    }

    
    for(igrid=0;igrid<nfft2_proc;igrid++){
      wfReal[ioff*nfft2+igrid] = zfft_threads[0][igrid*2+1]*prefact;
      rho[igrid] += zfft_threads[0][igrid*2+1]*zfft_threads[0][igrid*2+1]*occPre[nstate-1];
    }

  }//endif nstat%2

  for(igrid=0;igrid<nfft2_proc;igrid++){
    rho[igrid] *= invVolFrag;
    //printf("rhoooooooooooooo %lg\n",rho[igrid]);
  }

  //exit(0);
  /*
  //debug
  double dot;
  int js;
  int icoef;
  for(is=0;is<nstate;is++){
    for(js=is;js<nstate;js++){
      dot = 0.0;
      for(igrid=0;igrid<nfft2_proc;igrid++){
	dot += wfReal[is*nfft2_proc+igrid]*wfReal[js*nfft2_proc+igrid];
      }
      dot *= volFrag/nfft2_proc;
      printf("is %i js %i dot %lg\n",is,js,dot);
    }
  }
  for(is=0;is<nstate;is++){
    for(js=is;js<nstate;js++){
      dot = 0.0;
      for(icoef=1;icoef<ncoef;icoef++){
	dot += ccreal[is*ncoef+icoef]*ccreal[js*ncoef+icoef]+
	       ccimag[is*ncoef+icoef]*ccimag[js*ncoef+icoef];
      }
      dot *= 2.0;
      dot += ccreal[is*ncoef+ncoef]*ccreal[js*ncoef+ncoef];
      printf("is %i js %i dot %lg\n",is,js,dot);
    }
  }
  */

  free(occPre);
  free(rho_scr_threads);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void embedWfReal(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *classMini,
		 CP *cp,double *wfRealUpSmall,double *wfRealDnSmall)
/*========================================================================*/
{/*begin routine*/
/* Loop all grid points in small box and calculate the gaussian conv with */
/* respect to the nbhd in the bigger box.				  */
/*========================================================================*/
/*             Local variable declarations                                */
  CELL *cellMini = &(generalDataMini->cell);
  CPCOEFFS_INFO *cpcoeffsInfoMini = &(cpMini->cpcoeffs_info);
  CPOPTS *cpOptsMini = &(cpMini->cpopts);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cpMini->cp_para_fft_pkg3d_lg);

  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo            = stodftInfo->fragInfo;

  int nkf1    = cp_para_fft_pkg3d_lg->nkf1;
  int nkf2    = cp_para_fft_pkg3d_lg->nkf2;
  int nkf3    = cp_para_fft_pkg3d_lg->nkf3;
  int iState,iGrid,jGrid,kGrid;
  int numGridNbhd = 3;
  int iiGrid,jjGrid,kkGrid,indexa,indexb,indexc;
  int indexaNbhd,indexbNbhd,indexcNbhd;
  int indexaFold,indexbFold,indexcFold,indexFold;
  int indexSmall;
  int cpLsda		= cpOptsMini->cp_lsda;
  int iFrag		= fragInfo->iFrag;
  int numGridSkin	= fragInfo->numGridSkin;
  int numStateUpFrag    = cpcoeffsInfoMini->nstate_up_proc;
  int numStateDnFrag    = cpcoeffsInfoMini->nstate_dn_proc;
  int numGridFragSmall  = fragInfo->numGridFragProcSmall[iFrag];
  int numGridFrag       = fragInfo->numGridFragProc[iFrag];

  int *gridMapSmall     = fragInfo->gridMapProcSmall[iFrag];
  int *numGridFragDim	= fragInfo->numGridFragDim[iFrag];

  double sigma		= fragInfo->gaussianSigma;
  double pre		= -0.5/(sigma*sigma);
  double volMini,volElem,volsqrt;
  double occnumber;

  double *hmatMini = cellMini->hmat;
  double aGrid[3],bGrid[3],cGrid[3];
  double diff[3],distsq;
  double *rhoUpFragProc,*rhoDnFragProc,*coefUpFragProc,*coefDnFragProc;

  //debug
  int i,j,k;
  double sum;
  FILE *fileA;

/*--------------------------------------------------------------------------*/
/*VI) Get the grid vector						    */

  aGrid[0] = hmatMini[1]/nkf1;
  aGrid[1] = hmatMini[2]/nkf1;
  aGrid[2] = hmatMini[3]/nkf1;
  bGrid[0] = hmatMini[4]/nkf2;
  bGrid[1] = hmatMini[5]/nkf2;
  bGrid[2] = hmatMini[6]/nkf2;
  cGrid[0] = hmatMini[7]/nkf3;
  cGrid[1] = hmatMini[8]/nkf3;
  cGrid[2] = hmatMini[9]/nkf3;

/*--------------------------------------------------------------------------*/
/*II) Calculate the convoluted wave functions                               */

  rhoUpFragProc = fragInfo->rhoUpFragProc[iFrag];
  coefUpFragProc = fragInfo->coefUpFragProc[iFrag];
  if(cpLsda==1&&numStateDnFrag!=0){
    rhoDnFragProc = fragInfo->rhoDnFragProc[iFrag];
    coefDnFragProc = fragInfo->coefDnFragProc[iFrag];
  }
  //printf("sigma %lg\n",sigma);
  //printf("alen %lg blen %lg clen %lg\n",aGrid[0],bGrid[1],cGrid[2]);
  /*
  fileA = fopen("test-wf-sm","w");
  for(i=0;i<numGridFragSmall;i++){
    sum = 0.0;
    for(j=0;j<numStateUpFrag;j++){
      fprintf(fileA,"%lg ",wfRealUpSmall[j*numGridFragSmall+i]);
      //sum += wfRealUpSmall[j*numGridFragSmall+i]*wfRealUpSmall[j*numGridFragSmall+i];
    }
    fprintf(fileA,"\n");
    //printf("11111111 rhosm %lg\n",sum);
  }
  fclose(fileA);
  */

  //sigma = 1.0e-3;
  //pre = -0.5/(sigma*sigma);

  for(iState=0;iState<numStateUpFrag;iState++){
    for(iGrid=0;iGrid<nkf3;iGrid++){
      indexc = iGrid+numGridSkin;
      for(jGrid=0;jGrid<nkf2;jGrid++){
	indexb = jGrid+numGridSkin;
	for(kGrid=0;kGrid<nkf1;kGrid++){
	  indexa = kGrid+numGridSkin;
	  indexSmall = iState*numGridFragSmall+iGrid*nkf2*nkf1+jGrid*nkf1+kGrid;
	  for(iiGrid=-numGridNbhd;iiGrid<=numGridNbhd;iiGrid++){
	    indexcFold = indexc+iiGrid;
	    if(indexcFold<0)indexcFold += numGridFragDim[2];
	    if(indexcFold>=numGridFragDim[2]) indexcFold -= numGridFragDim[2];
	    for(jjGrid=-numGridNbhd;jjGrid<=numGridNbhd;jjGrid++){
	      indexbFold = indexb+jjGrid;
	      if(indexbFold<0)indexbFold += numGridFragDim[1];
	      if(indexbFold>=numGridFragDim[1]) indexbFold -= numGridFragDim[1];
	      for(kkGrid=-numGridNbhd;kkGrid<=numGridNbhd;kkGrid++){
	        indexaFold = indexa+kkGrid;
		if(indexaFold<0)indexaFold += numGridFragDim[0];
		if(indexaFold>=numGridFragDim[0]) indexaFold -= numGridFragDim[0];
		diff[0] = iiGrid*cGrid[0]+jjGrid*bGrid[0]+kkGrid*aGrid[0];
                diff[1] = iiGrid*cGrid[1]+jjGrid*bGrid[1]+kkGrid*aGrid[1];
                diff[2] = iiGrid*cGrid[2]+jjGrid*bGrid[2]+kkGrid*aGrid[2];
		distsq = diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2];
		indexFold = indexcFold*numGridFragDim[1]*numGridFragDim[0]+
			    indexbFold*numGridFragDim[0]+indexaFold;
		coefUpFragProc[iState*numGridFrag+indexFold] += exp(distsq*pre)*wfRealUpSmall[indexSmall];
	      }//endfor kkGrid
	    }//endfor jjGrid
	  }//endfor iiGrid
	}//endfor kGrid
      }//endfor jGrid
    }//endfor iGrid
    // Test. Let's recover the wave fun values in small box
    for(iGrid=0;iGrid<numGridFragSmall;iGrid++){
      //coefUpFragProc[iState*numGridFrag+gridMapSmall[iGrid]] = wfRealUpSmall[iState*numGridFragSmall+iGrid];
    }
  }//endfor iState
  if(cpLsda==1&&numStateDnFrag!=0){
    for(iState=0;iState<numStateDnFrag;iState++){
      for(iGrid=0;iGrid<nkf3;iGrid++){
	indexc = iGrid+2;
	for(jGrid=0;jGrid<nkf2;jGrid++){
	  indexb = jGrid+2;
	  for(kGrid=0;kGrid<nkf1;kGrid++){
	    indexa = kGrid+2;
	    indexSmall = iState*numGridFrag+iGrid*nkf2*nkf1+jGrid*nkf1+kGrid;
	    for(iiGrid=-numGridNbhd;iiGrid<=numGridNbhd;iiGrid++){
	      indexcFold = indexc+iiGrid;
	      if(indexcFold<0)indexcFold += numGridFragDim[2];
	      if(indexcFold>=numGridFragDim[2]) indexcFold -= numGridFragDim[2];
	      for(jjGrid=-numGridNbhd;jjGrid<=numGridNbhd;jjGrid++){
		indexbFold = indexb+jjGrid;
		if(indexbFold<0)indexbFold += numGridFragDim[1];
		if(indexbFold>=numGridFragDim[1]) indexbFold -= numGridFragDim[1];
		for(kkGrid=numGridNbhd;kkGrid<=numGridNbhd;kkGrid++){
		  indexaFold = indexa+kkGrid;
		  if(indexaFold<0)indexaFold += numGridFragDim[0];
		  if(indexaFold>=numGridFragDim[0]) indexaFold -= numGridFragDim[0];
		  diff[0] = iiGrid*cGrid[0]+jjGrid*bGrid[0]+kkGrid*aGrid[0];
		  diff[1] = iiGrid*cGrid[1]+jjGrid*bGrid[1]+kkGrid*aGrid[1];
		  diff[2] = iiGrid*cGrid[2]+jjGrid*bGrid[2]+kkGrid*aGrid[2];
		  distsq = diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2];
		  indexFold = indexcFold*numGridFragDim[1]*numGridFragDim[0]+
			      indexbFold*numGridFragDim[0]+indexaFold;
		  coefDnFragProc[iState*numGridFrag+indexFold] += exp(distsq*pre)*wfRealDnSmall[indexSmall];
		}//endfor kkGrid
	      }//endfor jjGrid
	    }//endfor iiGrid
	  }//endfor kGrid
	}//endfor jGrid
      }//endfor iGrid
    }//endfor iState
  }//endif

/*--------------------------------------------------------------------------*/
/* II) Orthonormal by QR decomposition	                                    */
   
  //qrWrapper(coefUpFragProc,numStateUpFrag,numGridFrag);
  /*
  fileA = fopen("test-wf","w");
  for(i=0;i<numGridFrag;i++){
    for(j=0;j<numStateUpFrag;j++)fprintf(fileA,"%lg ",coefUpFragProc[j*numGridFrag+i]);
    fprintf(fileA,"\n");
  }
  fclose(fileA);
  */
  qrWrapper(coefUpFragProc,numStateUpFrag,numGridFrag);

  volMini  = getdeth(hmatMini);
  volElem = volMini/numGridFragSmall;
  volsqrt = 1.0/sqrt(volElem);
  for(iGrid=0;iGrid<numStateUpFrag*numGridFrag;iGrid++)coefUpFragProc[iGrid] *= volsqrt;
  if(cpLsda==1&&numStateDnFrag!=0){
    qrWrapper(coefDnFragProc,numStateDnFrag,numGridFrag);
    for(iGrid=0;iGrid<numStateDnFrag*numGridFrag;iGrid++)coefDnFragProc[iGrid]*volsqrt;
  }
  //test orthonormal
  /*
  int jState;
  double sum;
  for(iState=0;iState<numStateUpFrag;iState++){
    for(jState=iState;jState<numStateUpFrag;jState++){
      sum = 0.0;    
      for(iGrid=0;iGrid<numGridFrag;iGrid++){
        sum += coefUpFragProc[iState*numGridFrag+iGrid]*coefUpFragProc[jState*numGridFrag+iGrid];
      }
      sum *= volElem;
      printf("iState %i jState %i sum %lg\n",iState,jState,sum);
    }
  }
  */

/*--------------------------------------------------------------------------*/
/* III) Calculate density	                                            */

  occnumber = 1;
  if(cpLsda==0)occnumber = 2;
  for(iState=0;iState<numStateUpFrag;iState++){
    for(iGrid=0;iGrid<numGridFrag;iGrid++){
      rhoUpFragProc[iGrid] += occnumber*coefUpFragProc[iState*numGridFrag+iGrid]*
				coefUpFragProc[iState*numGridFrag+iGrid];
    }//endfor iGrid
  }//endfor iState
  /*
  sum = 0.0;
  for(iGrid=0;iGrid<numGridFrag;iGrid++){
    printf("111111 rhoFrag %lg\n",rhoUpFragProc[iGrid]);
    sum += rhoUpFragProc[iGrid];
  }
  printf("sum %lg\n",sum*volElem);
  fflush(stdout);
  exit(0);  
  */
  if(cpLsda==1&&numStateDnFrag!=0){
    for(iState=0;iState<numStateDnFrag;iState++){
      for(iGrid=0;iGrid<numGridFrag;iGrid++){
	rhoDnFragProc[iGrid] += occnumber*coefDnFragProc[iGrid]*coefDnFragProc[iGrid];
      }//endfor iGrid
    }//endfor iState
  }//endif  

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void qrWrapper(double *waveFun, int nstate,int ngrid) 
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
  int m = ngrid;
  int n = nstate;
  int lda = m;
  int lwork = 64*n;
  int info;
  double *work = (double*)cmalloc(lwork*sizeof(double));
  double *A = (double*)cmalloc(m*n*sizeof(double));
  double *tau = (double*)cmalloc(n*sizeof(double));

  memcpy(A,waveFun,m*n*sizeof(double));
  
  /*
  FILE *fileA = fopen("test-wf","w");
  int i,j,k;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++)fprintf(fileA,"%lg ",A[j*m+i]);
    fprintf(fileA,"\n");
  }
  fclose(fileA);
  fflush(stdout);
  exit(0);
  */

  dgeqrf_(&m,&n,A,&lda,tau,work,&lwork,&info);
  if(info<0){
    printf("In qr decomposition, the %i'th element has illegal value!\n",-info);
    fflush(stdout);
    exit(0);
  }
  
  /*
  FILE *fileA = fopen("test-wf","w");
  int i,j,k;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++)fprintf(fileA,"%lg ",A[j*m+i]);
    fprintf(fileA,"\n");
  }
  fclose(fileA);
  fflush(stdout);
  exit(0);
  */
  
  //int i;
  //for(i=0;i<n;i++)printf("tau %lg\n",tau[i]);
  // generate q matrix


  int lwork2 = m*32;
  double *work2 = (double*)cmalloc(2*lwork2*sizeof(double));
  dorgqr_(&m,&n,&n,A,&lda,tau,&work2[0],&lwork2,&info);
  if(info<0){
    printf("In q calculation, the %i'th element has illegal value!\n",-info);
    fflush(stdout);
    exit(0);
  }

  
  /*
  FILE *fileA = fopen("test-wf","w");
  int i,j,k;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++)fprintf(fileA,"%lg ",A[j*m+i]);
    fprintf(fileA,"\n");
  }
  fclose(fileA);
  fflush(stdout);
  exit(0);
  */
  memcpy(waveFun,A,m*n*sizeof(double));
  free(work);
  free(A);
  free(tau);
  free(work2);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void noiseRealReGen(GENERAL_DATA *general_data,CP *cp,CLASS *class,int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR        *cpscr        = &(cp->cpscr);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  FRAGINFO	*fragInfo    = stodftInfo->fragInfo;	
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  COMMUNICATE   *communicate      = &(cp->communicate);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm = &(cp->cp_sclr_fft_pkg3d_sm);

  int iStat,iGrid,iOff,iOff2,iProc;
  int cpLsda = cpopts->cp_lsda;
  int numRandNum;
  int numStatUpProc = cpcoeffs_info->nstate_up_proc;
  int numStatDnProc = cpcoeffs_info->nstate_dn_proc;
  int nfft          = cp_para_fft_pkg3d_lg->nfft;
  int nfft2         = nfft/2;
  int numCoeff = cpcoeffs_info->ncoef;
  int numStatUpTot = numStatUpProc*numCoeff;
  int numStatDnTot = numStatDnProc*numCoeff;
  int numProcStates             = communicate->np_states;
  int myidState                 = communicate->myid_state;
  MPI_Comm comm_states   =    communicate->comm_states;

  double ranValue = sqrt(nfft2);
  double *randNumSeedTot = stodftInfo->randSeedTot;
  double *randNum;
  double *noiseWfUpReal = fragInfo->noiseWfUpReal;
  double *noiseWfDnReal = fragInfo->noiseWfDnReal;


  numRandNum = numStatUpProc*nfft2;
  if(cpLsda==1)numRandNum += numStatDnProc*nfft2;
  randNum = (double*)cmalloc(numRandNum*sizeof(double));

  // Generate the random number seeds from the given random number seed
  if(myidState==0){
#ifdef MKL_RANDOM
    VSLStreamStatePtr stream;
    int errcode;
    int seed = (int)(seed = stodftInfo->randSeed);
    errcode = vslNewStream(&stream,VSL_BRNG_MCG31,seed);
    errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,
                            numProcStates,randNumSeedTot,0.0,100000.0);
#endif
#ifndef MKL_RANDOM
    double seed = stodftInfo->randSeed;
    printf("seed!!!!!!!! %lg\n",seed);
    //whatever random number is good, I'm using Gaussian in this case
    //double seed = 8.3;
    //double seed = 2.5;
    int iseed;
    //printf("numRandTot %i numStateUpTot %i numCoeff %i\n",numRandTot,numStatUpTot,numCoeff);
    //fflush(stdout);
    double x;
    gaussran2(numProcStates,&iseed,&iseed,&seed,randNumSeedTot);
    for(iProc=0;iProc<numProcStates;iProc++){
      x = randNumSeedTot[iProc]*randNumSeedTot[iProc];
      if(x>=1.0)randNumSeedTot[iProc] = x*100.0;
      else randNumSeedTot[iProc] = 100.0/x;
    }
#endif
  }
  // Bcast the random number seeds
  if(numProcStates>1){
    Bcast(randNumSeedTot,numProcStates,MPI_DOUBLE,0,comm_states);
    Barrier(comm_states);
  }

#ifdef MKL_RANDOM
  VSLStreamStatePtr streamNew;
  int errcodeNew;
  int seedNew = (int)randNumSeedTot[myidState];
  errcode = vslNewStream(&streamNew,VSL_BRNG_MCG31,seedNew);
  errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,streamNew,
                          numRandNum,randNum,-1.0,1.0);
#endif
#ifndef MKL_RANDOM
  double seedNew = randNumSeedTot[myidState];
  int iseedNew;
  gaussran2(numRandNum,&iseedNew,&iseedNew,&seedNew,randNum);
  //debug
#endif
  
  /*
  char fileNameRand[100];
  FILE *fileRand;
  
  sprintf(fileNameRand,"rand-%i",myidState);
  fileRand = fopen(fileNameRand,"w");
  for(iGrid=0;iGrid<numStatUpProc*nfft2;iGrid++){
    if(randNum[iGrid]<0.0)fprintf(fileRand,"-1.0\n");
    else fprintf(fileRand,"1.0\n");
  }
  fclose(fileRand);
  */
  
 
  /*
  //char fileNameRand[100];
  //FILE *fileRand;
  printf("I'm reading noise orbital\n");
  sprintf(fileNameRand,"rand-%i",myidState);
  fileRand = fopen(fileNameRand,"r");
  for(iStat=0;iStat<numStatUpProc;iStat++){
    for(iGrid=0;iGrid<nfft2;iGrid++){
      fscanf(fileRand,"%lg",&randNum[iStat*nfft2+iGrid]);
      //if(randNum[iStat*nfft2+iGrid]<0.0)fprintf(fileRand,"-1.0\n");
      //else fprintf(fileRand,"1.0\n");
    }
  }
  fclose(fileRand);
  */
  
  //printf("randNum[1] %lg\n",randNum[1]);
  for(iGrid=0;iGrid<numStatUpProc*nfft2;iGrid++){
    if(randNum[iGrid]<0.0)noiseWfUpReal[iGrid] = -ranValue;
    else noiseWfUpReal[iGrid] = ranValue;
  }
  if(cpLsda==1){
    for(iGrid=0;iGrid<numStatDnProc*nfft2;iGrid++){
      if(randNum[iGrid]<0.0)noiseWfDnReal[iGrid] = -ranValue;
      else noiseWfDnReal[iGrid] = ranValue;
    }
  }
  /* 
  char fileNameRand[100];
  FILE *fileRand;
  double test;
  sprintf(fileNameRand,"rand-%i",myidState);
  //fileRand = fopen(fileNameRand,"r");
  fileRand = fopen("rand-all","r");
  for(iGrid=0;iGrid<numStatUpProc*nfft2;iGrid++){
    fscanf(fileRand,"%lg",&test);
    if(test<0.0)noiseWfUpReal[iGrid] = -ranValue;
    else noiseWfUpReal[iGrid] = ranValue;
  }
  */

  //DEBUG
  /*
  int *gridMapProc = fragInfo->gridMapProc[0];
  double *wfTemp = (double*)cmalloc(2*nfft2*sizeof(double))-1;
  double *wfTemp_temp = (double*)cmalloc(2*nfft2*sizeof(double))-1;

  double *cre_temp = (double*)cmalloc(numCoeff*sizeof(double))-1;
  double *cim_temp = (double*)cmalloc(numCoeff*sizeof(double))-1;

  FILE *ftest = fopen("noise-test-k","w");
  for(iStat=0;iStat<numStatUpProc;iStat++){
    for(iGrid=0;iGrid<nfft2;iGrid++){
      wfTemp[2*iGrid+1] = noiseWfUpReal[iStat*nfft2+gridMapProc[iGrid]];
      //wfTemp[2*iGrid+1] = noiseWfUpReal[iStat*nfft2+iGrid];
      wfTemp[2*iGrid+2] = 0.0;
    }
    for(iGrid=0;iGrid<numCoeff;iGrid++){
      cre_temp[iGrid+1] = 0.0;
      cim_temp[iGrid+1] = 0.0;
    }
    para_fft_gen3d_bck_to_g(wfTemp,wfTemp_temp,cp_sclr_fft_pkg3d_sm);
    sngl_upack_coef_sum(cre_temp,cim_temp,wfTemp,
                            cp_sclr_fft_pkg3d_sm);
    
    for(iGrid=0;iGrid<numCoeff-1;iGrid++){
      cre_temp[iGrid+1] *= 0.25;
      cim_temp[iGrid+1] *= 0.25;
    }
    cre_temp[numCoeff] *= 0.5;
    cim_temp[numCoeff] = 0.0;
    
    for(iGrid=0;iGrid<numCoeff;iGrid++){
      fprintf(ftest,"%.16lg %.16lg\n",cre_temp[iGrid+1],cim_temp[iGrid+1]);
    }
  }
  fclose(ftest);
  */

  free(randNum);
  if(numProcStates>1)Barrier(comm_states);
  //fflush(stdout);
  //exit(0);
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void noiseFilterGen(GENERAL_DATA *general_data,CP *cp,CLASS *class,
                    int ip_now)
/*========================================================================*/
{/*begin routine*/
/*************************************************************************/
/* Instead of using real space stochastic orbitals, we decompose these   */
/* orbitals by filters of energy windows corresponding to a h_KS.        */
/* The h_KS is constructed by fragment density.                          */
/*************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  CELL *cell             = &(general_data->cell);
  PTENS *ptens           = &(general_data->ptens);
  STAT_AVG *stat_avg     = &(general_data->stat_avg);
  EWALD *ewald           = &(general_data->ewald); 

  STODFTINFO *stodftInfo         = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos   = cp->stodftCoefPos;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CPEWALD *cpewald              = &(cp->cpewald);
  CPSCR *cpscr                  = &(cp->cpscr);
  PSEUDO *pseudo                = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal       = &(pseudo->pseudoReal);
  COMMUNICATE *communicate      = &(cp->communicate);

  FOR_SCR      *for_scr         = &(class->for_scr);
  EWD_SCR      *ewd_scr         = &(class->ewd_scr);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  CLATOMS_POS *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  ATOMMAPS *atommaps        = &(class->atommaps);

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
  int iperd             = cell->iperd;
  int checkPerdSize             = cpopts->icheck_perd_size;
  int checkDualSize             = cpopts->icheck_dual_size;

  int rhoRealGridNum    = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot    = stodftInfo->rhoRealGridTot;
  int numChemProc       = stodftInfo->numChemProc;
  int numStateStoUp     = stodftInfo->numStateStoUp;
  int numStateStoDn     = stodftInfo->numStateStoDn;
  int occNumber         = stodftInfo->occNumber;
  int densityMixFlag    = stodftInfo->densityMixFlag;
  int readCoeffFlag     = stodftInfo->readCoeffFlag;
  int myidState         = communicate->myid_state;
  int numProcStates     = communicate->np_states;
  int pseudoRealFlag    = pseudoReal->pseudoRealFlag;

  int iCoeff,iChem,iGrid;
  int iCell;
  int index;
  int i,j,k;
  int reRunFlag;
  MPI_Comm comm_states = communicate->comm_states;

  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);
  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;

  double volCP,rvolCP;
  double tolEdgeDist        = cpopts->tol_edge_dist;

  double *rhoUpRead,*rhoDnRead;
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
  double *rhoUpFragSumCpy = fragInfo->rhoUpFragSumCpy;
  double *rhoDnFragSumCpy = fragInfo->rhoDnFragSumCpy;
  double *ptensPvtenTmp    = ptens->pvten_tmp;

  double *rhoUpRealBackup;
  double *rhoDnRealBackup;
/*======================================================================*/
/* I) Copy the fragment density into rho and calculate k-space density  */

  if(readCoeffFlag!=-1){
    rhoUpRealBackup = (double*)cmalloc(rhoRealGridNum*sizeof(double));
    memcpy(rhoUpRealBackup,rhoUpCorrect,rhoRealGridNum*sizeof(double));
    if(cpLsda==1&&numStateDnProc!=0){
      rhoDnRealBackup = (double*)cmalloc(rhoRealGridNum*sizeof(double));
      memcpy(rhoDnRealBackup,rhoDnCorrect,rhoRealGridNum*sizeof(double));
    }
  }

  // I don't need to do anything if the fragment density is the initial guess
  memcpy(rhoUpCorrect,rhoUpFragSumCpy,rhoRealGridNum*sizeof(double));
  memcpy(&rhoUp[1],rhoUpCorrect,rhoRealGridNum*sizeof(double));
  if(cpLsda==1&&numStateDnProc!=0){
    memcpy(rhoDnCorrect,rhoDnFragSumCpy,rhoRealGridNum*sizeof(double));
    memcpy(&rhoDn[1],rhoDnCorrect,rhoRealGridNum*sizeof(double));
  }

  if(numProcStates>1)Barrier(comm_states);
  printf("rhoUpCorrect %lg\n",rhoUpCorrect[0]);

  calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
                     rhoCoeffReUp,rhoCoeffImUp,rhoUp,rhoCoeffReUpDensCpBox,rhoCoeffImUpDensCpBox,
                     divRhoxUp,divRhoyUp,divRhozUp,d2RhoUp,cpGGA,cpDualGridOptOn,numInterpPmeDual,
                     communicate,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));
  if(cpLsda==1&&numStateDnProc!=0){
    calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
                       rhoCoeffReDn,rhoCoeffImDn,rhoDn,rhoCoeffReDnDensCpBox,rhoCoeffImDnDensCpBox,
                       divRhoxDn,divRhoyDn,divRhozDn,d2RhoDn,cpGGA,cpDualGridOptOn,numInterpPmeDual,
                       communicate,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));
    for(iCoeff=1;iCoeff<=numCoeffLargeProc;iCoeff++){
      rhoCoeffReUp[iCoeff] += rhoCoeffReDn[iCoeff];
      rhoCoeffImUp[iCoeff] += rhoCoeffImDn[iCoeff];
    }// endfor 
    if(cpDualGridOptOn>=1){
      for(iCoeff=1;iCoeff<=numCoeffLargeProcDensCpBox;iCoeff++){
        rhoCoeffReUpDensCpBox[iCoeff] += rhoCoeffReDnDensCpBox[iCoeff];
        rhoCoeffImUpDensCpBox[iCoeff] += rhoCoeffImDnDensCpBox[iCoeff];
      }// endfor
    } // endif
  }// endif 
  //printf("rhoUpCorrect %lg\n",rhoUpCorrect[1]);
  //printf("rhoCoeffReUp[1] %lg %lg\n",rhoCoeffReUp[1],rhoCoeffImUp[1]);

/*======================================================================*/
/* II) Filter the stochastic orbital                                    */

  if((iperd<3||iperd==4)&&checkPerdSize==1){
    cp_boundary_check(cell,clatoms_info,clatoms_pos,tolEdgeDist);
  }//endif
  if(cpDualGridOptOn>=1&&checkDualSize==1){
    cp_dual_check(cell,clatoms_info,clatoms_pos,
                  atommaps->cp_atm_lst,tolEdgeDist);
  }//endif

  cpcoeffs_pos->ifcoef_form_up = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1;
  if(cpLsda==1&&numStateDnProc>0){
    cpcoeffs_pos->ifcoef_form_dn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }//endif

  for(iCell=1;iCell<=9;iCell++){ptensPvtenTmp[iCell] = 0.0;}
  gethinv(cell->hmat_cp,cell->hmati_cp,&(cell->vol_cp),iperd);
  gethinv(cell->hmat,cell->hmati,&(cell->vol),iperd);

  calcLocalPseudoScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

  calcKSPot(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  /*
  for(iGrid=1;iGrid<numFFT2;iGrid++){
    if(isnan(cpscr->cpscr_rho.v_ks_up[1])==1){
      printf("v_ks_up %lg\n",cpscr->cpscr_rho.v_ks_up[1]);
    }
  }
  */

  if(pseudoRealFlag==1){
    if(myidState==0)printf("**Calculating Real Space Non-local Pseudopotential...\n");
    pseudoReal->forceCalcFlag = 1;
    initRealNlppWf(cp,class,general_data);
    allocRealNl(cp,class);
    pseudoReal->forceCalcFlag = 0;
  }
  //printf("vnlPhiAtomGridRe[1] %lg\n",pseudoReal->vnlPhiAtomGridRe[1]);  

  //genStoOrbitalEnergyWindow(class,general_data,cp,ip_now);

  /*
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iCoeff=1;iCoeff<numCoeff*numStateUpProc;iCoeff++){
      printf("filllllllllter %i %i %lg %lg\n",iChem,iCoeff,stodftCoefPos->stoWfUpRe[iChem][iCoeff],stodftCoefPos->stoWfUpIm[iChem][iCoeff]);
    }
  }
  fflush(stdout);
  exit(0);
  */

/*======================================================================*/
/* III) Restore the correct density                                     */

  if(readCoeffFlag!=-1){  
    memcpy(rhoUpCorrect,rhoUpRealBackup,rhoRealGridNum*sizeof(double));
    memcpy(&rhoUp[1],rhoUpCorrect,rhoRealGridNum*sizeof(double));
    if(cpLsda==1&&numStateDnProc!=0){
      memcpy(rhoDnCorrect,rhoDnRealBackup,rhoRealGridNum*sizeof(double));
      memcpy(&rhoDn[1],rhoDnCorrect,rhoRealGridNum*sizeof(double));
    }
    if(numProcStates>1)Barrier(comm_states);
    calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
                       rhoCoeffReUp,rhoCoeffImUp,rhoUp,rhoCoeffReUpDensCpBox,rhoCoeffImUpDensCpBox,
                       divRhoxUp,divRhoyUp,divRhozUp,d2RhoUp,cpGGA,cpDualGridOptOn,numInterpPmeDual,
                       communicate,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));
    if(cpLsda==1&&numStateDnProc!=0){
      calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
                         rhoCoeffReDn,rhoCoeffImDn,rhoDn,rhoCoeffReDnDensCpBox,rhoCoeffImDnDensCpBox,
                         divRhoxDn,divRhoyDn,divRhozDn,d2RhoDn,cpGGA,cpDualGridOptOn,numInterpPmeDual,
                         communicate,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));
      for(iCoeff=1;iCoeff<=numCoeffLargeProc;iCoeff++){
        rhoCoeffReUp[iCoeff] += rhoCoeffReDn[iCoeff];
        rhoCoeffImUp[iCoeff] += rhoCoeffImDn[iCoeff];
      }// endfor 
      if(cpDualGridOptOn>=1){
        for(iCoeff=1;iCoeff<=numCoeffLargeProcDensCpBox;iCoeff++){
          rhoCoeffReUpDensCpBox[iCoeff] += rhoCoeffReDnDensCpBox[iCoeff];
          rhoCoeffImUpDensCpBox[iCoeff] += rhoCoeffImDnDensCpBox[iCoeff];
        }// endfor
      } // endif
    }// endif 
  }//endif
  /*
  else{
    free(rhoUpFragSumCpy);
    if(cpLsda==1&&numStateDnProc!=0){
      free(rhoDnFragSumCpy);
    }
  }
  */

  if(readCoeffFlag!=-1){
    free(rhoUpRealBackup);
    if(cpLsda==1&&numStateDnProc!=0){
      free(rhoDnRealBackup);
    }
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void noiseFilterGenFake(GENERAL_DATA *general_data,CP *cp,CLASS *class,
                    int ip_now)
/*========================================================================*/
{/*begin routine*/
/*************************************************************************/
/* Fake version of noiseFilterGen. The filter is constructed with        */
/* deterministic orbitals.                                               */
/*************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS  *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR        *cpscr        = &(cp->cpscr);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  NEWTONINFO   *newtonInfo      = stodftInfo->newtonInfo;
  COMMUNICATE   *commCP         = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS  *cpcoeffs_pos   = &(cp->cpcoeffs_pos[ip_now]);
  MPI_Comm commStates = commCP->comm_states;
  CELL *cell = &(general_data->cell);

  int iPoly,iState,jState,iCoeff,iChem,iOff,iOff2,iGrid;
  int iProc;
  int numChemPot = stodftInfo->numChemPot;
  int polynormLength = stodftInfo->polynormLength;
  int expanType = stodftInfo->expanType;
  int numProcStates = commCP->np_states;
  int myidState = commCP->myid_state;
  int cpLsda = cpopts->cp_lsda;
  int cpGGA  = cpopts->cp_gga;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int totalPoly;
  int numCoeffUpTot   = numStateUpProc*numCoeff;
  int numCoeffDnTot   = numStateDnProc*numCoeff;
  int numStatesDet = stodftInfo->numStatesDet;
  int rhoRealGridTot = stodftInfo->rhoRealGridTot;

  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *forceFormUp  = &(cpcoeffs_pos->ifcoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *forceFormDn  = &(cpcoeffs_pos->ifcoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *forceOrthUp  = &(cpcoeffs_pos->ifcoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);
  int *forceOrthDn  = &(cpcoeffs_pos->ifcoef_orth_dn);

  double dot;

  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;
  double *coeffForceReUp = cpcoeffs_pos->fcre_up;
  double *coeffForceImUp = cpcoeffs_pos->fcim_up;

  double *wfDetBackupUpRe = stodftCoefPos->wfDetBackupUpRe;
  double *wfDetBackupUpIm = stodftCoefPos->wfDetBackupUpIm;
  double *wfReTemp,*wfImTemp;
  double *wfReProjTemp,*wfImProjTemp;
  double *ksEigv = (double*)cmalloc(numStatesDet*sizeof(double));
  double *hmatCP = cell->hmat_cp;
  double *chemPot = stodftCoefPos->chemPot;


  FILE *feigv = fopen("orbital-e","r");

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

 
/*======================================================================*/
/* I) Set flags                     */

  *forceFormUp = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1; // temp keep this flag for debug filter
  if(cpLsda==1&&numStateDnProc!=0){
    *forceFormDn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }


/*======================================================================*/
/* IV) Generate random orbital                                          */

  genNoiseOrbitalReal(cp,cpcoeffs_pos);

/*======================================================================*/
/* II) Filtering by deterministic orbitals                     */

  wfReTemp = (double*)cmalloc(numCoeff*sizeof(double));
  wfImTemp = (double*)cmalloc(numCoeff*sizeof(double));
  wfReProjTemp = (double*)cmalloc(numCoeff*sizeof(double));
  wfImProjTemp = (double*)cmalloc(numCoeff*sizeof(double));

  for(iChem=0;iChem<numChemPot;iChem++){
    for(iCoeff=1;iCoeff<=numCoeffUpTot;iCoeff++){
      stoWfUpRe[iChem][iCoeff] = 0.0;
      stoWfUpIm[iChem][iCoeff] = 0.0;
    }
  }

  //numChemPot = 1000; // TEST
  int junk;
  int windowIndex;
  int index;
  stodftCoefPos->ewStateNum = (int*)cmalloc(numChemPot*sizeof(int));
  stodftCoefPos->ewStateMap = (int**)cmalloc(numChemPot*sizeof(int*));
  int *ewStateNum = stodftCoefPos->ewStateNum;
  int **ewStateMap = stodftCoefPos->ewStateMap;
  for(iChem=0;iChem<numChemPot;iChem++){
    ewStateMap[iChem] = (int*)cmalloc(numStatesDet*sizeof(int));
  }
  for(iState=0;iState<numStatesDet;iState++){
    fscanf(feigv,"%i",&junk);
    fscanf(feigv,"%lg",&ksEigv[iState]);
  }

  double energyMin = stodftInfo->energyMin;
  // test only
  energyMin = ksEigv[0]-0.1;
  double dE = (chemPot[numChemPot-2]-energyMin)/(numChemPot-1.0);
  /*
  for(iState=0;iState<numStatesDet;iState++){
    if(ksEigv[iState]<chemPot[numChemPot-2]){
      printf("iState %i ksEigv %lg chemPot %lg energyMin %lg\n",iState,ksEigv[iState],chemPot[numChemPot-2],energyMin);
      windowIndex = (int)((ksEigv[iState]-energyMin)/dE);
      ewStateMap[windowIndex][ewStateNum[windowIndex]] = iState;
      ewStateNum[windowIndex] += 1;
    }
    else{
      ewStateMap[numChemPot-1][ewStateNum[numChemPot-1]] = iState;
      ewStateNum[numChemPot-1] += 1;
    }
  }
  */
  for(iChem=0;iChem<numChemPot;iChem++){
    ewStateNum[iChem] = 0;
  }
  for(iState=0;iState<numStatesDet;iState++){
    if(ksEigv[iState]<chemPot[0]){
      ewStateMap[0][ewStateNum[0]] = iState;
      ewStateNum[0] += 1;
      if(iState==numStatesDet-1)stodftInfo->homoIndex = 0;
    }
  }
  for(iChem=1;iChem<numChemPot-1;iChem++){
    for(iState=0;iState<numStatesDet;iState++){
      if(ksEigv[iState]>=chemPot[iChem-1]&&ksEigv[iState]<chemPot[iChem]){
        ewStateMap[iChem][ewStateNum[iChem]] = iState;
        ewStateNum[iChem] += 1;
        if(iState==numStatesDet-1)stodftInfo->homoIndex = iChem;
      }
    }
  }
  for(iState=0;iState<numStatesDet;iState++){
    if(ksEigv[iState]>chemPot[numChemPot-2]){
      ewStateMap[numChemPot-1][ewStateNum[numChemPot-1]] = iState;
      ewStateNum[numChemPot-1] += 1;
      if(iState==numStatesDet-1)stodftInfo->homoIndex = numStatesDet-1;
    }
  }

  for(iChem=0;iChem<numChemPot;iChem++){
    printf("iChem %i ewStateNum %i\n",iChem,ewStateNum[iChem]);
  }
 

  for(iState=0;iState<numStateUpProc;iState++){
    iOff = iState*numCoeff;
    memcpy(wfReTemp,&coeffReUp[iOff+1],numCoeff*sizeof(double));
    memcpy(wfImTemp,&coeffImUp[iOff+1],numCoeff*sizeof(double));
    for(iChem=0;iChem<numChemPot-1;iChem++){
      for(jState=0;jState<ewStateNum[iChem];jState++){
        dot = 0.0;
        index = ewStateMap[iChem][jState];
        iOff2 = index*numCoeff;
        #pragma omp parallel for reduction(+:dot) private(iCoeff)
        for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
          dot += wfDetBackupUpRe[iOff2+iCoeff]*wfReTemp[iCoeff]+
                 wfDetBackupUpIm[iOff2+iCoeff]*wfImTemp[iCoeff];
        }
        dot = (dot*2.0+wfDetBackupUpRe[iOff2+numCoeff-1]*wfReTemp[numCoeff-1])*0.5;
        #pragma omp parallel for private(iCoeff)
        for(iCoeff=0;iCoeff<numCoeff;iCoeff++){
          stoWfUpRe[iChem][iOff+iCoeff+1] += wfDetBackupUpRe[iOff2+iCoeff]*dot;
          stoWfUpIm[iChem][iOff+iCoeff+1] += wfDetBackupUpIm[iOff2+iCoeff]*dot;
        }
      }
    }
    memcpy(&stoWfUpRe[numChemPot-1][iOff+1],&coeffReUp[iOff+1],numCoeff*sizeof(double));
    memcpy(&stoWfUpIm[numChemPot-1][iOff+1],&coeffImUp[iOff+1],numCoeff*sizeof(double));
    for(iChem=0;iChem<numChemPot-1;iChem++){
      #pragma omp parallel for private(iCoeff)
      for(iCoeff=0;iCoeff<numCoeff;iCoeff++){
        stoWfUpRe[numChemPot-1][iOff+iCoeff+1] -= stoWfUpRe[iChem][iOff+iCoeff+1];
        stoWfUpIm[numChemPot-1][iOff+iCoeff+1] -= stoWfUpIm[iChem][iOff+iCoeff+1];
      }
    }
  }

  //debug
  /*
  double testNumElec;
  testNumElec = 0.0;
  for(iChem=0;iChem<numChemPot-1;iChem++){
    for(iState=0;iState<numStateUpProc;iState++){
      iOff = iState*numCoeff;
      dot = 0.0;
      for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
        dot += stoWfUpRe[iChem][iOff+iCoeff+1]*stoWfUpRe[iChem][iOff+iCoeff+1]+
               stoWfUpIm[iChem][iOff+iCoeff+1]*stoWfUpIm[iChem][iOff+iCoeff+1];
      }
      dot = dot*2.0+stoWfUpRe[iChem][iOff+numCoeff]*stoWfUpRe[iChem][iOff+numCoeff];
      testNumElec += dot;
    }
  }
  testNumElec /= numStateUpProc;
  printf("testNumElec %lg\n",testNumElec*2.0);
  */

  free(wfReTemp);
  free(wfImTemp);
  free(wfReProjTemp);
  free(wfImProjTemp);
  free(ksEigv);


/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void noiseFilterRealReGen(CLASS *class,GENERAL_DATA *general_data,
                          CP *cp,double *creal, double *cimag,double *noiseWfReal,
                          int nstate)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  STODFTINFO    *stodftInfo      = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos   = cp->stodftCoefPos;
  FRAGINFO      *fragInfo        = stodftInfo->fragInfo;
  COMMUNICATE   *commCP          = &(cp->communicate);
  CPSCR         *cpscr           = &(cp->cpscr);
  CELL          *cell            = &(general_data->cell);
  COMMUNICATE   *communicate     = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
    
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm = &(cp->cp_sclr_fft_pkg3d_sm);
  
  int i,ioff,ioff2,is,iupper;
  int myid_state       =    communicate->myid_state;
  int np_states        =    communicate->np_states;
  int nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
  int nfft             =    cp_para_fft_pkg3d_lg->nfft;
  int nfft2            =    nfft/2;
  int nfft2_proc       =    nfft_proc/2;
  int ncoef            =    cpcoeffs_info->ncoef;


  double vol_cp,rvol_cp;
  double temp_r,temp_i;
  double *zfft           =    cpscr->cpscr_wave.zfft;
  double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
  double *hmati_cp       =    cell->hmati_cp;
  double *hmat_cp        =    cell->hmat_cp;
  double sum_sq = 0.0;

/*======================================================================*/
/* I) Check the form of the coefficients */

  MPI_Comm comm_states   =    communicate->comm_states;
     
/*==========================================================================*/
/*==========================================================================*/
/*2)sum the density in real space two states at a time                      */
/*  This is done at state level and uses the scalar packages!               */

  iupper = nstate;
  if((nstate % 2) != 0){
    iupper = nstate-1;
  }// endif

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

    for(i=0;i<nfft2;i++){
      noiseWfReal[(is-1)*nfft2+i] = zfft[i*2+1];
      noiseWfReal[is*nfft2+i] = zfft[i*2+2];
    }
  }//endfor is

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
    for(i=0;i<nfft2;i++){
      noiseWfReal[(nstate-1)*nfft2+i] = zfft[i*2+1];
    }
  }//endif

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

