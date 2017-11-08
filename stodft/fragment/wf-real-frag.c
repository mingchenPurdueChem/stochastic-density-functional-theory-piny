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
  int igrid;
  int cp_lsda = cpopts->cp_lsda;
  int ncoef = cpcoeffs_info->ncoef;
  int   nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
  int   nfft2_proc       =    nfft_proc/2;
  int nfft = cp_para_fft_pkg3d_lg->nfft;
  int nfft2 = nfft/2;
  int gridOff1,gridOff2;

  //printf("nfftttttt2_proc %i\n",nfft2_proc);

  double *zfft           =    cpscr->cpscr_wave.zfft;
  double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
  double *hmatCP         =    cell->hmat_cp;
  double volFrag         = getdeth(hmatCP);
  double invVolFrag	 = 1.0/volFrag;
  double prefact;

  prefact = 1.0/sqrt(2.0*volFrag);
  if(cp_lsda==1)prefact = sqrt(invVolFrag);


/* ================================================================= */
/*1) zero density and gradients if necessary                         */

  for(i=0;i<nfft2;i++)rho[i] = 0.0;

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
    gridOff1 = (is-1)*nfft2_proc;
    gridOff2 = is*nfft2_proc;

/*--------------------------------------------------------------------------*/
/*I) pack the complex zfft array with two wavefunctions (real) 
    the wavefunctions are reperesented in spherically cuttof 
    half g space                                                            */

    dble_pack_coef(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
                      zfft,cp_sclr_fft_pkg3d_sm);

    //dble_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
    //                  zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*II) fourier transform the two wavefunctions to real space
     convention exp(-igr)                                                   */

    //para_fft_gen3d_fwd_to_r_fftw3d(zfft,cp_sclr_fft_pkg3d_sm);
    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* III) Copy the real sapce wave function and add the square of the two     
        wave functions to the density(real space)			    */
 
    for(igrid=0;igrid<nfft2_proc;igrid++){
      wfReal[gridOff1+igrid] = zfft[igrid*2+1]*prefact;
      wfReal[gridOff2+igrid] = zfft[igrid*2+2]*prefact;
      rho[igrid] += zfft[igrid*2+1]*zfft[igrid*2+1]+zfft[igrid*2+2]*zfft[igrid*2+2];
    }
    //printf("ioff %i ioff2 %i wfReal %lg %lg\n",ioff,ioff2,wfReal[ioff],wfReal[ioff2]);
  }/*endfor is*/

/*--------------------------------------------------------------------------*/
/*IV) if there is an odd number of states, use single pack
      for the last state (i caen't be chaengin' the laws of physics
      captn, i've got to be usin the single pack!!!!!!)                     */

  if(nstate%2!=0){
    ioff = (nstate-1)*ncoef;
    //sngl_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);
    sngl_pack_coef(&ccreal[ioff],&ccimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* V) fourier transform the last wavefunction to real space
       convention exp(-igr)                                                 */
    
    //para_fft_gen3d_fwd_to_r_fftw3d(zfft,cp_sclr_fft_pkg3d_sm);
    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*VI) Copy the real sapce wave function and add the square of the last wave 
      function to the density(real space)   */
    
    for(igrid=0;igrid<nfft2_proc;igrid++){
      wfReal[ioff*nfft2+igrid] = zfft[igrid*2+1]*prefact;
      rho[igrid] += zfft[igrid*2+1]*zfft[igrid*2+1];
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
#endif
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
  free(randNum);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

