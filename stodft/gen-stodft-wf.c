/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: gen-stodft-wf.c                              */
/*                                                                          */
/* Subroutine to filter noise orbitals                                      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_stodft_local.h"
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genStoOrbitalInterp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
	    CP *cp,int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS  *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR	       *cpscr	     = &(cp->cpscr);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  NEWTONINFO   *newtonInfo      = stodftInfo->newtonInfo;
  COMMUNICATE   *commCP	        = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS  *cpcoeffs_pos   = &(cp->cpcoeffs_pos[ip_now]);
  MPI_Comm commStates = commCP->comm_states;
  

  int iPoly,iState,iState2,iCoeff,iChem;
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

  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *forceFormUp  = &(cpcoeffs_pos->ifcoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *forceFormDn  = &(cpcoeffs_pos->ifcoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *forceOrthUp  = &(cpcoeffs_pos->ifcoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);
  int *forceOrthDn  = &(cpcoeffs_pos->ifcoef_orth_dn);
  
  double energyMin,energyMax;
  double Smin,Smax; 
  double energyDiff;
  double scale;
  double length;

  double *chemPot = stodftCoefPos->chemPot;
  double *sampPoint = newtonInfo->sampPoint;
  double *sampPointUnscale = newtonInfo->sampPointUnscale;
  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;
  double *scrCoeffReUp   = cpscr->cpscr_wave.cre_up;
  double *scrCoeffImUp   = cpscr->cpscr_wave.cim_up;
  double *rhoCoeffReUp   = cpscr->cpscr_rho.rhocr_up;
  double *rhoCoeffImUp   = cpscr->cpscr_rho.rhoci_up;
  double *rhoUp          = cpscr->cpscr_rho.rho_up;
  double *rhoCoeffReUpDensCpBox = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
  double *rhoCoeffImUpDensCpBox = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
  double *divRhoxUp       = cpscr->cpscr_grho.d_rhox_up;
  double *divRhoyUp       = cpscr->cpscr_grho.d_rhoy_up;
  double *divRhozUp       = cpscr->cpscr_grho.d_rhoz_up;
  double *d2RhoUp        = cpscr->cpscr_grho.d2_rho_up;
  double *scrCoeffReDn   = cpscr->cpscr_wave.cre_dn;
  double *scrCoeffImDn   = cpscr->cpscr_wave.cim_dn;
  double *rhoCoeffReDn   = cpscr->cpscr_rho.rhocr_dn;
  double *rhoCoeffImDn   = cpscr->cpscr_rho.rhoci_dn;
  double *rhoDn          = cpscr->cpscr_rho.rho_dn;
  double *rhoCoeffReDnDensCpBox = cpscr->cpscr_rho.rhocr_dn_dens_cp_box;
  double *rhoCoeffImDnDensCpBox = cpscr->cpscr_rho.rhoci_dn_dens_cp_box;
  double *divRhoxDn       = cpscr->cpscr_grho.d_rhox_dn;
  double *divRhoyDn       = cpscr->cpscr_grho.d_rhoy_dn;
  double *divRhozDn       = cpscr->cpscr_grho.d_rhoz_dn;
  double *d2RhoDn         = cpscr->cpscr_grho.d2_rho_dn;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  double *coeffReUpBackup = stodftCoefPos->coeffReUpBackup;
  double *coeffImUpBackup = stodftCoefPos->coeffImUpBackup;
  
/*==========================================================================*/
/* 0) Check the forms				    */

/*======================================================================*/
/* 0.05) Check the approximations in the methods                        */

/*======================================================================*/
/* II) In parallel, transpose coefs back to normal form                 */

/*======================================================================*/
/* III) Set flags			    */

  *forceFormUp = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1; // temp keep this flag for debug filter
  if(cpLsda==1&&numStateDnProc!=0){
    *forceFormDn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }

/*======================================================================*/
/* IV) Calculate Emax and Emin                                          */

  //if(myidState==0){
  genEnergyMax(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
  genEnergyMin(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
  //}
  if(numProcStates>1){
    Barrier(commStates);
    Bcast(&(stodftInfo->energyMin),1,MPI_DOUBLE,0,commStates);
    Bcast(&(stodftInfo->energyMax),1,MPI_DOUBLE,0,commStates);
  }
  energyMin = stodftInfo->energyMin;
  energyMax = stodftInfo->energyMax;
  //printf("energyMax %lg energyMin %lg\n",energyMax,energyMin);
  Barrier(commStates);
  stodftInfo->energyDiff = energyMax-energyMin;
  energyDiff = stodftInfo->energyDiff;
  stodftInfo->energyMean = 0.5*(energyMin+energyMax);

/*======================================================================*/
/* V) Generate Coeffcients for Polynormial interpolation                */
  
  if(expanType==2){
    Smin = newtonInfo->Smin;
    Smax = newtonInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    newtonInfo->scale = scale;
    if(stodftInfo->filterDiagFlag==1){
      // regenerate chemical potentials
      if(numProcStates>1){
	Barrier(commStates);
	Bcast(&(stodftInfo->eigValMin),1,MPI_DOUBLE,0,commStates);
      }
      int iScf = stodftInfo->iScf;
      if(iScf>1){ //Not the first SCF step
	double eigValMin = stodftInfo->eigValMin;
	double eigValLUMO = stodftInfo->eigValLUMO;
	if(myidState==0){
	  printf("eigValMin %lg\n",eigValMin);
	  printf("eigValLUMO %lg\n",eigValLUMO);  
	}
	stodftInfo->gapInit = eigValLUMO-eigValMin;
	stodftInfo->chemPotInit = 0.5*(eigValLUMO+eigValMin);
	genChemPotInterpPoints(stodftInfo,stodftCoefPos);
      }
    }
    if(myidState==0)genNewtonHermit(stodftInfo,stodftCoefPos);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      polynormLength = stodftInfo->polynormLength;
      totalPoly = polynormLength*numChemPot;
      if(myidState!=0){
	stodftCoefPos->expanCoeff = (double*)crealloc(stodftCoefPos->expanCoeff,
				totalPoly*sizeof(double));
	newtonInfo->sampPoint = (double*)crealloc(newtonInfo->sampPoint,
				polynormLength*sizeof(double));
	newtonInfo->sampPointUnscale = (double*)crealloc(newtonInfo->sampPointUnscale,
				polynormLength*sizeof(double));
      }
      Barrier(commStates);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPoint,polynormLength,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPointUnscale,polynormLength,MPI_DOUBLE,0,commStates);
    }
    /*
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      sampPointUnscale[iPoly] = (sampPoint[iPoly]-Smin)/scale+energyMin;
      //printf("sampunscale %lg\n",sampPointUnscale[iPoly]);
    }
    genCoeffNewtonHermit(stodftInfo,stodftCoefPos);
    */
  }

/*======================================================================*/
/* IV) Generate random orbital                                          */

  genNoiseOrbital(cp,cpcoeffs_pos);
  //debug
  /*
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    coeffReUp[iCoeff] = stodftCoefPos->creTest[iCoeff];
    coeffImUp[iCoeff] = stodftCoefPos->cimTest[iCoeff];
  }
  */
  
  
  
/*======================================================================*/
/* V) Filter the stochastic orbitals			*/

  switch(expanType){
    case 2:
      filterNewtonPolyHerm(cp,class,general_data,ip_now);
      break;
  }

//debug print wave function
  
  //Barrier(commStates);
  /*
  char wfname[100];
  //sprintf(wfname,"/scratch/mingchen/tmp/sto-wf-save-%i",myidState);
  FILE *filePrintWF = NULL;
  //filePrintWF = fopen(wfname,"w");
  filePrintWF = fopen("sto-wf-save","w");
  if(filePrintWF!=NULL){
    for(iChem=0;iChem<numChemPot;iChem++){
      for(iState=0;iState<numStateUpProc;iState++){
	for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
	  fprintf(filePrintWF,"%.16lg %.16lg\n",
          stoWfUpRe[iChem][iState*numCoeff+iCoeff],stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
        }//endfor iCoeff
      }//endfor iState
    }//endfor iChem
    fclose(filePrintWF);
  }
  else{
    printf("myid %i I can't open files to write stochastic orbitals!\n",myidState);
    fflush(stdout);
  }
  */
  
  
  double testfilter = 0.0;
  for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
    testfilter += stoWfUpRe[0][iCoeff]*stoWfUpRe[0][iCoeff]+
		stoWfUpIm[0][iCoeff]*stoWfUpIm[0][iCoeff];
  }
  testfilter += stoWfUpRe[0][numCoeff]*stoWfUpRe[0][numCoeff];
  printf("testfilter %lg\n",testfilter);
  
  //Barrier(commStates);
    
  

/*======================================================================*/
/* VI) Calculate Energy	                                            */

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genStoOrbitalCheby(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
	    CP *cp,int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS  *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR	       *cpscr	     = &(cp->cpscr);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  NEWTONINFO   *newtonInfo      = stodftInfo->newtonInfo;
  COMMUNICATE   *commCP	        = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS  *cpcoeffs_pos   = &(cp->cpcoeffs_pos[ip_now]);
  MPI_Comm commStates = commCP->comm_states;
  

  int iPoly,iState,iState2,iCoeff,iChem;
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

  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *forceFormUp  = &(cpcoeffs_pos->ifcoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *forceFormDn  = &(cpcoeffs_pos->ifcoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *forceOrthUp  = &(cpcoeffs_pos->ifcoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);
  int *forceOrthDn  = &(cpcoeffs_pos->ifcoef_orth_dn);
  
  double energyMin,energyMax;
  double chemPotInit = stodftInfo->chemPotInit;
  double gapInit = stodftInfo->gapInit;
  double Smin,Smax; 
  double energyDiff;
  double scale;
  double length;

  double *sampPoint = newtonInfo->sampPoint;
  double *sampPointUnscale = newtonInfo->sampPointUnscale;
  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;
  double *scrCoeffReUp   = cpscr->cpscr_wave.cre_up;
  double *scrCoeffImUp   = cpscr->cpscr_wave.cim_up;
  double *rhoCoeffReUp   = cpscr->cpscr_rho.rhocr_up;
  double *rhoCoeffImUp   = cpscr->cpscr_rho.rhoci_up;
  double *rhoUp          = cpscr->cpscr_rho.rho_up;
  double *rhoCoeffReUpDensCpBox = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
  double *rhoCoeffImUpDensCpBox = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
  double *divRhoxUp       = cpscr->cpscr_grho.d_rhox_up;
  double *divRhoyUp       = cpscr->cpscr_grho.d_rhoy_up;
  double *divRhozUp       = cpscr->cpscr_grho.d_rhoz_up;
  double *d2RhoUp        = cpscr->cpscr_grho.d2_rho_up;
  double *scrCoeffReDn   = cpscr->cpscr_wave.cre_dn;
  double *scrCoeffImDn   = cpscr->cpscr_wave.cim_dn;
  double *rhoCoeffReDn   = cpscr->cpscr_rho.rhocr_dn;
  double *rhoCoeffImDn   = cpscr->cpscr_rho.rhoci_dn;
  double *rhoDn          = cpscr->cpscr_rho.rho_dn;
  double *rhoCoeffReDnDensCpBox = cpscr->cpscr_rho.rhocr_dn_dens_cp_box;
  double *rhoCoeffImDnDensCpBox = cpscr->cpscr_rho.rhoci_dn_dens_cp_box;
  double *divRhoxDn       = cpscr->cpscr_grho.d_rhox_dn;
  double *divRhoyDn       = cpscr->cpscr_grho.d_rhoy_dn;
  double *divRhozDn       = cpscr->cpscr_grho.d_rhoz_dn;
  double *d2RhoDn         = cpscr->cpscr_grho.d2_rho_dn;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  double *coeffReUpBackup = stodftCoefPos->coeffReUpBackup;
  double *coeffImUpBackup = stodftCoefPos->coeffImUpBackup;
  
/*======================================================================*/
/* I) Set flags			    */

  *forceFormUp = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1; // temp keep this flag for debug filter
  if(cpLsda==1&&numStateDnProc!=0){
    *forceFormDn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }
  stodftInfo->filterFlag = 1;

/*======================================================================*/
/* II) Calculate Emax and Emin                                          */

  //if(myidState==0){
  genEnergyMax(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
  genEnergyMin(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
  
  //}
  if(numProcStates>1){
    Barrier(commStates);
    Bcast(&(stodftInfo->energyMin),1,MPI_DOUBLE,0,commStates);
    Bcast(&(stodftInfo->energyMax),1,MPI_DOUBLE,0,commStates);
  }
  energyMin = stodftInfo->energyMin;
  energyMax = stodftInfo->energyMax;
  //printf("energyMax %lg energyMin %lg\n",energyMax,energyMin);
  Barrier(commStates);
  stodftInfo->energyDiff = energyMax-energyMin;
  energyDiff = stodftInfo->energyDiff;
  stodftInfo->energyMean = 0.5*(energyMin+energyMax);


/*======================================================================*/
/* III) Generate Length of Polynomial Chain		                */
  
  if(expanType==2){
    Smin = newtonInfo->Smin;
    Smax = newtonInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    newtonInfo->scale = scale;
    //cheat my code
    stodftInfo->numChemPot = 2;
    free(&(stodftCoefPos->chemPot[0]));
    stodftCoefPos->chemPot = (double*)cmalloc(2*sizeof(double));
    stodftCoefPos->chemPot[0] = chemPotInit-gapInit*0.5;
    stodftCoefPos->chemPot[1] = chemPotInit+gapInit*0.5;
    if(myidState==0)genNewtonHermit(stodftInfo,stodftCoefPos);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      Barrier(commStates);
    }
    stodftInfo->numChemPot = 1;
    free(&(stodftCoefPos->chemPot[0]));
    stodftCoefPos->chemPot = (double*)cmalloc(sizeof(double));    
    //finish cheating my code
  }

/*======================================================================*/
/* IV) Calculate the True Chemical Potential                            */

  calcChemPotCheby(cp,class,general_data,ip_now);  

/*======================================================================*/
/* V) Generate Coeffcients for Polynormial interpolation with Correct   */
/*    Chemical Potential.						*/
  
  
  if(expanType==2){
    Smin = newtonInfo->Smin;
    Smax = newtonInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    newtonInfo->scale = scale;
    if(myidState==0)genNewtonHermitTrueChemPot(stodftInfo,stodftCoefPos);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      polynormLength = stodftInfo->polynormLength;
      totalPoly = polynormLength*numChemPot;
      if(myidState!=0){
	stodftCoefPos->expanCoeff = (double*)crealloc(stodftCoefPos->expanCoeff,
			    totalPoly*sizeof(double));
	newtonInfo->sampPoint = (double*)crealloc(newtonInfo->sampPoint,
			    polynormLength*sizeof(double));
	newtonInfo->sampPointUnscale = (double*)crealloc(newtonInfo->sampPointUnscale,
			    polynormLength*sizeof(double));
      }
      Barrier(commStates);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPoint,polynormLength,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPointUnscale,polynormLength,MPI_DOUBLE,0,commStates);
    }
    
    /*
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      sampPointUnscale[iPoly] = (sampPoint[iPoly]-Smin)/scale+energyMin;
      //printf("sampunscale %lg\n",sampPointUnscale[iPoly]);
    }
    genCoeffNewtonHermit(stodftInfo,stodftCoefPos);
    */
  }

/*======================================================================*/
/* IV) Generate random orbital                                          */

  genNoiseOrbitalReal(cp,cpcoeffs_pos);

/*======================================================================*/
/* V) Filter the stochastic orbitals			*/

  switch(expanType){
    case 2:
      filterNewtonPolyHerm(cp,class,general_data,ip_now);
      break;
  }
  stodftInfo->filterFlag = 0;
//debug print wave function
  //Barrier(commStates);
  /*
  char wfname[100];
  //sprintf(wfname,"/scratch/mingchen/tmp/sto-wf-save-%i",myidState);
  sprintf(wfname,"sto-wf-save-%i",myidState);
  FILE *filePrintWF = NULL;
  filePrintWF = fopen(wfname,"w");
  //filePrintWF = fopen("sto-wf-save","w");
  numChemPot = stodftInfo->numChemPot;
  if(filePrintWF!=NULL){
    for(iChem=0;iChem<numChemPot;iChem++){
      for(iState=0;iState<numStateUpProc;iState++){
	for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
	  fprintf(filePrintWF,"%.16lg %.16lg\n",
          stoWfUpRe[iChem][iState*numCoeff+iCoeff],stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
	}//endfor iCoeff
      }//endfor iState
    }//endfor iChem
    fclose(filePrintWF);
  }
  else{
    printf("myid %i I can't open files to write stochastic orbitals!\n",myidState);
    fflush(stdout);
  }
  */
  //Barrier(commStates);
/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genStoOrbitalFake(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
            CP *cp,int ip_now)
/*========================================================================*/
{/*begin routine*/
/* When we do the filtering for test cases, we can construct the filter   */
/* by deterministic orbitals. DEBUG ONLY NO PARALLEL			  */
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
  
  int iPoly,iState,jState,iCoeff,iChem,iOff,iOff2;
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
 
  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *forceFormUp  = &(cpcoeffs_pos->ifcoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *forceFormDn  = &(cpcoeffs_pos->ifcoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *forceOrthUp  = &(cpcoeffs_pos->ifcoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);
  int *forceOrthDn  = &(cpcoeffs_pos->ifcoef_orth_dn);
  int *numStateUpAllProc;

  double dot;

  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;
  double *wfDetBackupUpRe = stodftCoefPos->wfDetBackupUpRe;
  double *wfDetBackupUpIm = stodftCoefPos->wfDetBackupUpIm;
  double *wfReTemp,*wfImTemp;
  double *wfReProjTemp,*wfImProjTemp;
 
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

  numStateUpAllProc = (int*)cmalloc(numProcStates*sizeof(int));
  if(numProcStates>1){
    Allgather(&numStateUpProc,1,MPI_INT,numStateUpAllProc,1,MPI_INT,0,commStates);
  }
  else numStateUpAllProc[0] = numStateUpProc;
  wfReTemp = (double*)cmalloc(numCoeff*sizeof(double));
  wfImTemp = (double*)cmalloc(numCoeff*sizeof(double));
  wfReProjTemp = (double*)cmalloc(numCoeff*sizeof(double));
  wfImProjTemp = (double*)cmalloc(numCoeff*sizeof(double));

  for(iCoeff=1;iCoeff<=numCoeffUpTot;iCoeff++){
    stoWfUpRe[0][iCoeff] = 0.0;
    stoWfUpIm[0][iCoeff] = 0.0;
  }

  //printf("11111111 numStatesDet %i\n",numStatesDet);
  //double *wfNow = (double*)calloc(numCoeff*sizeof(double));
  //printf("Test projjjjj %i %lg %lg\n",numStatesDet,wfDetBackupUpRe[0],wfDetBackupUpIm[0]);
  //printf("Test proj 2 %lg %lg\n",coeffReUp[1],coeffImUp[1]);
  printf("wfDetBackupUpReeeeeeeee %.16lg %.16lg\n",wfDetBackupUpRe[0],wfDetBackupUpIm[0]);
  for(iProc=0;iProc<numProcStates;iProc++){
    for(iState=0;iState<numStateUpAllProc[iProc];iState++){
      iOff = iState*numCoeff;
      if(myidState==iProc){
        memcpy(wfReTemp,&coeffReUp[iOff+1],numCoeff*sizeof(double));
        memcpy(wfImTemp,&coeffImUp[iOff+1],numCoeff*sizeof(double));      
	//printf("iState %i iProc %i wfReTemp %lg\n",iState,iProc,wfReTemp[0]);
      }
      if(numProcStates>1){
	Bcast(wfReTemp,numCoeff,MPI_DOUBLE,iProc,commStates);
        Bcast(wfImTemp,numCoeff,MPI_DOUBLE,iProc,commStates);
      }
      //printf("wfReTemppppp %lg\n",wfReTemp[0]);
      for(iCoeff=0;iCoeff<numCoeff;iCoeff++){
        wfReProjTemp[iCoeff] = 0.0;
        wfImProjTemp[iCoeff] = 0.0;
      }
      for(jState=0;jState<numStatesDet;jState++){
	dot = 0.0;
	iOff2 = jState*numCoeff;
	for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
	  dot += wfDetBackupUpRe[iOff2+iCoeff]*wfReTemp[iCoeff]+
	  	 wfDetBackupUpIm[iOff2+iCoeff]*wfImTemp[iCoeff];
	}
	dot = (dot*2.0+wfDetBackupUpRe[iOff2+numCoeff-1]*wfReTemp[numCoeff-1])*0.5;
        for(iCoeff=0;iCoeff<numCoeff;iCoeff++){
          wfReProjTemp[iCoeff] += wfDetBackupUpRe[iOff2+iCoeff]*dot;
          wfImProjTemp[iCoeff] += wfDetBackupUpIm[iOff2+iCoeff]*dot;
        }
	/*
	for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
	  stoWfUpRe[0][iOff+iCoeff] += wfDetBackupUpRe[iOff2+iCoeff]*dot;
	  stoWfUpIm[0][iOff+iCoeff] += wfDetBackupUpIm[iOff2+iCoeff]*dot;
	}
	*/
      }//endfor jState
      if(numProcStates>1){
	Reduce(wfReProjTemp,&stoWfUpRe[0][iOff+1],numCoeff,MPI_DOUBLE,
	       MPI_SUM,iProc,commStates);
	Reduce(wfImProjTemp,&stoWfUpIm[0][iOff+1],numCoeff,MPI_DOUBLE,
	       MPI_SUM,iProc,commStates);
      }
      else{
	memcpy(&stoWfUpRe[0][iOff+1],wfReProjTemp,numCoeff*sizeof(double));
	memcpy(&stoWfUpIm[0][iOff+1],wfImProjTemp,numCoeff*sizeof(double));
      }
    }//endfor iState
  }//endfor iProc
  free(wfReTemp);
  free(wfImTemp);
  free(wfReProjTemp);
  free(wfImProjTemp);
  /*
  FILE *ftest = fopen("wf-test","w");
  for(iCoeff=1;iCoeff<=numStateUpProc*numCoeff;iCoeff++){
    fprintf(ftest,"%.16lg %.16lg\n",stoWfUpRe[0][iCoeff],stoWfUpIm[0][iCoeff]);
  }
  fclose(ftest);
  for(iState=0;iState<numStateUpProc;iState++){
    for(jState=0;jState<numStatesDet;jState++){
      dot = 0.0;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	dot += stoWfUpRe[0][iState*numCoeff+iCoeff]*wfDetBackupUpRe[jState*numCoeff+iCoeff-1]+
	       stoWfUpIm[0][iState*numCoeff+iCoeff]*wfDetBackupUpIm[jState*numCoeff+iCoeff-1];
      }
      dot *= 2.0;
      dot += stoWfUpRe[0][iState*numCoeff+numCoeff]*wfDetBackupUpRe[jState*numCoeff+numCoeff-1];
      printf("11111111 iState %i jState %i dot %lg\n",iState,jState,dot*sqrt(0.5));
    }
  }
  */

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

