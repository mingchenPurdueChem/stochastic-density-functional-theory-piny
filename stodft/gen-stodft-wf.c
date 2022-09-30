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
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_stodft_local.h"
#include "../proto_defs/proto_frag_local.h"
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genStoOrbitalInterp(CLASS *class,GENERAL_DATA *general_data,
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
  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;
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
  int smearOpt        = stodftInfo->smearOpt;
  int calcLocalTraceOpt = stodftInfo->calcLocalTraceOpt;

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

  if(myidState==0){
    genEnergyMax(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
    genEnergyMin(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
  }
  if(numProcStates>1){
    Barrier(commStates);
    Bcast(&(stodftInfo->energyMin),1,MPI_DOUBLE,0,commStates);
    Bcast(&(stodftInfo->energyMax),1,MPI_DOUBLE,0,commStates);
  }
  energyMin = stodftInfo->energyMin;
  energyMax = stodftInfo->energyMax;
  //printf("energyMax %lg energyMin %lg\n",energyMax,energyMin);
  if(numProcStates>1)Barrier(commStates);
  stodftInfo->energyDiff = energyMax-energyMin;
  energyDiff = stodftInfo->energyDiff;
  stodftInfo->energyMean = 0.5*(energyMin+energyMax);

/*======================================================================*/
/* V) Generate Coeffcients for Polynormial interpolation                */
 
  if(expanType==1){
    Smin = chebyshevInfo->Smin;
    Smax = chebyshevInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    chebyshevInfo->scale = scale;
    if(stodftInfo->filterDiagFlag==1){
      // regenerate chemical potentials
      if(numProcStates>1){
	Barrier(commStates);
	Bcast(&(stodftInfo->eigValMin),1,MPI_DOUBLE,0,commStates);
      }
      int iScfTrue = stodftInfo->iScfTrue;
      if(iScfTrue>1){ //Not the first SCF step
        int numStatePrintUp = stodftInfo->numStatePrintUp;
	double eigValMin = stodftCoefPos->energyLevel[0];
	double eigValMax = stodftCoefPos->energyLevel[numStatePrintUp-1];
	if(myidState==0){
	  printf("eigValMin %lg\n",eigValMin);
	  printf("eigValMax %lg\n",eigValMax);  
	}
	//stodftInfo->gapInit = eigValMax-eigValMin;
	stodftInfo->gapInit = eigValMax-energyMin;
	//stodftInfo->chemPotInit = 0.5*(eigValMax+eigValMin);
	stodftInfo->chemPotInit = 0.5*(eigValMax+energyMin);
	genChemPotInterpPoints(stodftInfo,stodftCoefPos);
      }
    }
    if(myidState==0)genChebyHermit(stodftInfo,stodftCoefPos,0);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      polynormLength = stodftInfo->polynormLength;
      totalPoly = polynormLength*numChemPot;
      if(myidState!=0){
	stodftCoefPos->expanCoeff = (double*)crealloc(stodftCoefPos->expanCoeff,
				totalPoly*sizeof(double));
        if(smearOpt>0&&stodftInfo->filterDiagFlag==0){
          stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));
        }
      }
      Barrier(commStates);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      if(smearOpt>0&&stodftInfo->filterDiagFlag==0){
        Bcast(stodftCoefPos->entropyExpanCoeff,polynormLength,MPI_DOUBLE,0,commStates);
      }
      Barrier(commStates);
    }
  }
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
      int iScfTrue = stodftInfo->iScfTrue;
      if(iScfTrue>1){ //Not the first SCF step
        int numStatePrintUp = stodftInfo->numStatePrintUp;
	double eigValMin = stodftCoefPos->energyLevel[0];
	double eigValMax = stodftCoefPos->energyLevel[numStatePrintUp-1];
	if(myidState==0){
	  printf("eigValMin %lg\n",eigValMin);
	  printf("eigValMax %lg\n",eigValMax);  
	}
	//stodftInfo->gapInit = eigValMax-eigValMin;
	stodftInfo->gapInit = eigValMax-energyMin;
	//stodftInfo->chemPotInit = 0.5*(eigValMax+eigValMin);
	stodftInfo->chemPotInit = 0.5*(eigValMax+energyMin);	
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
        if(smearOpt>0&&stodftInfo->filterDiagFlag==0){
          stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));
        }
      }
      Barrier(commStates);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPoint,polynormLength,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPointUnscale,polynormLength,MPI_DOUBLE,0,commStates);
      if(smearOpt>0&&stodftInfo->filterDiagFlag==0){
        Bcast(stodftCoefPos->entropyExpanCoeff,polynormLength,MPI_DOUBLE,0,commStates);
      }
      Barrier(commStates);
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
    case 1:
      stodftInfo->storeChebyMomentsFlag = 0;
      filterChebyPolyHerm(cp,class,general_data,ip_now);
      break;
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
  
  
  /*
  double testfilter = 0.0;
  for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
    testfilter += stoWfUpRe[0][iCoeff]*stoWfUpRe[0][iCoeff]+
		stoWfUpIm[0][iCoeff]*stoWfUpIm[0][iCoeff];
  }
  testfilter += stoWfUpRe[0][numCoeff]*stoWfUpRe[0][numCoeff];
  printf("testfilter %lg\n",testfilter);
  */
  
  //Barrier(commStates);
    
  

/*======================================================================*/
/* VI) Calculate Energy	                                            */

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genStoOrbitalCheby(CLASS *class,GENERAL_DATA *general_data,
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
  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;
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
  int smearOpt        = stodftInfo->smearOpt;

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

  //double *sampPoint;
  //double *sampPointUnscale;
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
  //double *scrCoeffReDn   = cpscr->cpscr_wave.cre_dn;
  //double *scrCoeffImDn   = cpscr->cpscr_wave.cim_dn;
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

  //timing
  double timeStart1,timeEnd1;
  double timeStart2,timeEnd2;
  double timeStart3,timeEnd3;
  double timeStart4,timeEnd4;
  double timeStart5,timeEnd5;
  double timeStart6,timeEnd6;
  double diffTime1 = 0.0;
  double diffTime2 = 0.0;
  double diffTime3 = 0.0;
  double diffTime4 = 0.0;
  double diffTime5 = 0.0;
  double diffTime6 = 0.0;


  
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

  timeStart1 = omp_get_wtime();
  //if(myidState==0){
  genEnergyMax(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
  //fflush(stdout);
  //exit(0);
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
  if(numProcStates>1)Barrier(commStates);
  stodftInfo->energyDiff = energyMax-energyMin;
  energyDiff = stodftInfo->energyDiff;
  stodftInfo->energyMean = 0.5*(energyMin+energyMax);
  timeEnd1 = omp_get_wtime();
  diffTime1 = timeEnd1-timeStart1;


/*======================================================================*/
/* III) Generate Length of Polynomial Chain		                */
  
  timeStart2 = omp_get_wtime();

  if(expanType==1){
    Smin = chebyshevInfo->Smin;
    Smax = chebyshevInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    chebyshevInfo->scale = scale;
    //cheat my code
    stodftInfo->numChemPot = 2;
    free(&(stodftCoefPos->chemPot[0]));
    stodftCoefPos->chemPot = (double*)cmalloc(2*sizeof(double));
    stodftCoefPos->chemPot[0] = chemPotInit-gapInit*0.5;
    stodftCoefPos->chemPot[1] = chemPotInit+gapInit*0.5;
    if(myidState==0)genChebyHermit(stodftInfo,stodftCoefPos,0);
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
  if(expanType==2){
    //sampPoint = newtonInfo->sampPoint;
    //sampPointUnscale = newtonInfo->sampPointUnscale;
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
  timeEnd2 = omp_get_wtime();
  diffTime2 = timeEnd2-timeStart2;

/*======================================================================*/
/* IV) Calculate the True Chemical Potential                            */

  timeStart3 = omp_get_wtime();
  calcChemPotCheby(cp,class,general_data,ip_now);
  timeEnd3 = omp_get_wtime();
  diffTime3 = timeEnd3-timeStart3;

/*======================================================================*/
/* V) Generate Coeffcients for Polynormial interpolation with Correct   */
/*    Chemical Potential.						*/
  
  timeStart4 = omp_get_wtime();
  if(expanType==1){
    Smin = chebyshevInfo->Smin;
    Smax = chebyshevInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    chebyshevInfo->scale = scale;
    if(myidState==0)genChebyHermitTrueChemPot(stodftInfo,stodftCoefPos,0);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      polynormLength = stodftInfo->polynormLength;
      totalPoly = polynormLength*numChemPot;
      if(myidState!=0){
	stodftCoefPos->expanCoeff = (double*)crealloc(stodftCoefPos->expanCoeff,
			    totalPoly*sizeof(double));
        if(smearOpt>0){
          stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));          
        }
      }
      Barrier(commStates);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      if(smearOpt>0){
        Bcast(stodftCoefPos->entropyExpanCoeff,polynormLength,MPI_DOUBLE,0,commStates);
      }
      Barrier(commStates);
    }
  }
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
        if(smearOpt>0){
          stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));          
        }
      }
      Barrier(commStates);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPoint,polynormLength,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPointUnscale,polynormLength,MPI_DOUBLE,0,commStates);
      if(smearOpt>0){
        Bcast(stodftCoefPos->entropyExpanCoeff,polynormLength,MPI_DOUBLE,0,commStates);
      }
      Barrier(commStates);
    }
    
    /*
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      sampPointUnscale[iPoly] = (sampPoint[iPoly]-Smin)/scale+energyMin;
      //printf("sampunscale %lg\n",sampPointUnscale[iPoly]);
    }
    genCoeffNewtonHermit(stodftInfo,stodftCoefPos);
    */
  }
  timeEnd4 = omp_get_wtime();
  diffTime4 = timeEnd4-timeStart4;


/*======================================================================*/
/* IV) Generate random orbital                                          */

  timeStart5 = omp_get_wtime();
  genNoiseOrbitalReal(cp,cpcoeffs_pos);
  timeEnd5 = omp_get_wtime();
  diffTime5 = timeEnd5-timeStart5;

/*======================================================================*/
/* V) Filter the stochastic orbitals			*/

  timeStart6 = omp_get_wtime();
  switch(expanType){
    case 1:
      stodftInfo->storeChebyMomentsFlag = 0;
      filterChebyPolyHerm(cp,class,general_data,ip_now);
      break;
    case 2:
      filterNewtonPolyHerm(cp,class,general_data,ip_now);
      break;
  }
  stodftInfo->filterFlag = 0;
  timeEnd6 = omp_get_wtime();
  diffTime6 = timeEnd6-timeStart6;

  printf("Gen-stowf time myid %i spec-range %.8lg gen-poly-length %.8lg gen-chempot %.8lg gen-poly-coeff %.8lg gen-rand %.8lg filter %.8lg\n",myidState,diffTime1,diffTime2,diffTime3,diffTime4,diffTime5,diffTime6);

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
void genStoOrbitalFake(CLASS *class,GENERAL_DATA *general_data,
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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genStoOrbitalEnergyWindow(CLASS *class,GENERAL_DATA *general_data,
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
  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;
  COMMUNICATE   *commCP	        = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS  *cpcoeffs_pos   = &(cp->cpcoeffs_pos[ip_now]);
  MPI_Comm commStates = commCP->comm_states;
  

  int iPoly,iState,iState2,iCoeff,iChem;
  int numChemPot = stodftInfo->numChemPot;
  int numChemPotTemp;
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
  int numCoeffUpTotal   = numStateUpProc*numCoeff;
  int numCoeffDnTotal   = numStateDnProc*numCoeff;
  int smearOpt        = stodftInfo->smearOpt;

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

  double *chemPot = stodftCoefPos->chemPot;
  //double *sampPoint = newtonInfo->sampPoint;
  //double *sampPointUnscale = newtonInfo->sampPointUnscale;
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

  //timing
  double timeStart1,timeEnd1;
  double timeStart2,timeEnd2;
  double timeStart3,timeEnd3;
  double timeStart4,timeEnd4;
  double timeStart5,timeEnd5;
  double timeStart6,timeEnd6;
  double diffTime1 = 0.0;
  double diffTime2 = 0.0;
  double diffTime3 = 0.0;
  double diffTime4 = 0.0;
  double diffTime5 = 0.0;
  double diffTime6 = 0.0;

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
  if(numProcStates>1)Barrier(commStates);
  stodftInfo->energyDiff = energyMax-energyMin;
  energyDiff = stodftInfo->energyDiff;
  stodftInfo->energyMean = 0.5*(energyMin+energyMax);

/*======================================================================*/
/* III) Generate Length of Polynomial Chain                             */

  timeStart2 = omp_get_wtime();
  if(expanType==1){
    Smin = chebyshevInfo->Smin;
    Smax = chebyshevInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    chebyshevInfo->scale = scale;
    //cheat my code
    numChemPotTemp = stodftInfo->numChemPot;
    stodftInfo->numChemPot = 2;
    free(&(stodftCoefPos->chemPot[0]));
    stodftCoefPos->chemPot = (double*)cmalloc(2*sizeof(double));
    stodftCoefPos->chemPot[0] = chemPotInit-gapInit*0.5;
    stodftCoefPos->chemPot[1] = chemPotInit+gapInit*0.5;
    if(myidState==0)genChebyHermit(stodftInfo,stodftCoefPos,0);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      Barrier(commStates);
    }
    //printf("numChemPotTemp %i\n",numChemPotTemp);
    stodftInfo->numChemPot = numChemPotTemp;
    free(&(stodftCoefPos->chemPot[0]));
    stodftCoefPos->chemPot = (double*)cmalloc(numChemPotTemp*sizeof(double));
    //finish cheating my code
  }
  if(expanType==2){
    Smin = newtonInfo->Smin;
    Smax = newtonInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    newtonInfo->scale = scale;
    //cheat my code
    numChemPotTemp = stodftInfo->numChemPot;
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
    //printf("numChemPotTemp %i\n",numChemPotTemp);
    stodftInfo->numChemPot = numChemPotTemp;
    free(&(stodftCoefPos->chemPot[0]));
    stodftCoefPos->chemPot = (double*)cmalloc(numChemPotTemp*sizeof(double));
    //finish cheating my code
  }
  timeEnd2 = omp_get_wtime();
  diffTime2 = timeEnd2-timeStart2;

/*======================================================================*/
/* IV) Calculate the True Chemical Potential                            */

  timeStart3 = omp_get_wtime();
  calcChemPotCheby(cp,class,general_data,ip_now);
  timeEnd3 = omp_get_wtime();
  diffTime3 = timeEnd3-timeStart3;

/*======================================================================*/
/* V) Generate list of chemical potential for energy windows.           */

  genChemPotInterpPoints(stodftInfo,stodftCoefPos);

/*======================================================================*/
/* V) Generate Coeffcients for Polynormial interpolation                */

  if(expanType==1){
    Smin = chebyshevInfo->Smin;
    Smax = chebyshevInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    chebyshevInfo->scale = scale;
    if(myidState==0)genChebyHermitTrueChemPot(stodftInfo,stodftCoefPos,1);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      polynormLength = stodftInfo->polynormLength;
      totalPoly = polynormLength*numChemPot;
      if(myidState!=0){
        stodftCoefPos->expanCoeff = (double*)crealloc(stodftCoefPos->expanCoeff,
                            totalPoly*sizeof(double));
        if(smearOpt>0){
          stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));
        }
      }
      Barrier(commStates);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      if(smearOpt>0){
        Bcast(stodftCoefPos->entropyExpanCoeff,polynormLength,MPI_DOUBLE,0,commStates);
      }
      Barrier(commStates);
    }
  }
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
        if(smearOpt>0){
            stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));        
        }
      }
      Barrier(commStates);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPoint,polynormLength,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPointUnscale,polynormLength,MPI_DOUBLE,0,commStates);
      if(smearOpt>0){
        Bcast(stodftCoefPos->entropyExpanCoeff,polynormLength,MPI_DOUBLE,0,commStates);
      }
      Barrier(commStates);
    }
  }

/*======================================================================*/
/* IV) Generate random orbital                                          */

  genNoiseOrbitalReal(cp,cpcoeffs_pos);
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
    case 1:
      stodftInfo->storeChebyMomentsFlag = 0;
      filterChebyPolyHerm(cp,class,general_data,ip_now);
      break;
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
  printf("1111111111111111111111111111\n");
  //filePrintWF = fopen(wfname,"w");
  //filePrintWF = fopen("sto-wf-save","w");
  //if(filePrintWF!=NULL){
    for(iChem=0;iChem<numChemPot;iChem++){
      sprintf(wfname,"sto-wf-save-%i",iChem);
      filePrintWF = fopen(wfname,"w");
      for(iState=0;iState<numStateUpProc;iState++){
	for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
	  fprintf(filePrintWF,"%.16lg %.16lg\n",
          stoWfUpRe[iChem][iState*numCoeff+iCoeff],stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
        }//endfor iCoeff
      }//endfor iState
      fclose(filePrintWF);
    }//endfor iChem
    //fclose(filePrintWF);
  //}
  //else{
  //  printf("myid %i I can't open files to write stochastic orbitals!\n",myidState);
  //  fflush(stdout);
  //}
  fflush(stdout);
  exit(0);
  */
  
/*======================================================================*/
/* VI) Calculate filtered states w.r.t. energy windows                  */
    
  //Barrier(commStates);
  
  /*  
  if(myidState==0)printf("Filter by window...\n");
  for(iChem=numChemPot-1;iChem>0;iChem--){
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      stoWfUpRe[iChem][iCoeff] -= stoWfUpRe[iChem-1][iCoeff];
      stoWfUpIm[iChem][iCoeff] -= stoWfUpIm[iChem-1][iCoeff];
    }//endfor iCoeff
  }//endfor iChem
  if(cpLsda==1){
    for(iChem=numChemPot-1;iChem>0;iChem--){
      for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
	stoWfDnRe[iChem][iCoeff] -= stoWfDnRe[iChem-1][iCoeff];
	stoWfDnIm[iChem][iCoeff] -= stoWfDnIm[iChem-1][iCoeff];
      }//endfor iCoeff
    }//endfor iChem
  }
  */

/*======================================================================*/
/* VI) Calculate Energy	                                            */

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genStoOrbitalEnergyWindowFake(CLASS *class,GENERAL_DATA *general_data,
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
  int fragWindowFlag = stodftInfo->fragWindowFlag;
  int homoIndex;

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
  double *chemPot = stodftCoefPos->chemPot;

  double *wfDetBackupUpRe = stodftCoefPos->wfDetBackupUpRe;
  double *wfDetBackupUpIm = stodftCoefPos->wfDetBackupUpIm;
  double *wfReTemp,*wfImTemp;
  double *wfReProjTemp,*wfImProjTemp;
  double *ksEigv = (double*)cmalloc(numStatesDet*sizeof(double));
  double *hmatCP = cell->hmat_cp;
 
  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  FILE *feigv = fopen("orbital-e","r");
  // TEST. Let's calculate the theoretical estimation
  //FILE *frealwf = fopen("wf-real","r");
  double *wfDetReal = stodftCoefPos->wfDetReal;

  /*
  double volCPRev  = 1.0/getdeth(hmatCP);
  for(iState=0;iState<numStatesDet;iState++){
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      fscanf(frealwf,"%lg",&wfDetReal[iState*rhoRealGridTot+iGrid]);
    }
  }
  fclose(frealwf);
  */

  
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
  double energyOrbMin = ksEigv[0]-0.0001;
  double energyOrbMax = ksEigv[numStatesDet-1]+0.0001;
  double dE = (energyOrbMax-energyOrbMin)/numChemPot;
  printf("11111111111111 dE %lg\n",dE);
  for(iState=0;iState<numStatesDet;iState++){
    windowIndex = (int)((ksEigv[iState]-energyOrbMin)/dE);
    ewStateMap[windowIndex][ewStateNum[windowIndex]] = iState;
    ewStateNum[windowIndex] += 1;
    printf("windowIndex %i\n",windowIndex);
  }

  // TEST:
  /*
  double *rhoTemp = (double*)calloc(rhoRealGridTot,sizeof(double));
  double *rhoStd = (double*)calloc(rhoRealGridTot,sizeof(double));
  double pre = sqrt(2.0);
  FILE *fileRhoStd = fopen("rho-std-theory","w");
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoTemp[iGrid] = 0.0;
    for(iState=0;iState<ewStateNum[iChem];iState++){
      index = ewStateMap[iChem][iState];
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
        rhoTemp[iGrid] += wfDetReal[index*rhoRealGridTot+iGrid]*wfDetReal[index*rhoRealGridTot+iGrid];
      }
    }
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoStd[iGrid] += rhoTemp[iGrid]*rhoTemp[iGrid];
  }
  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
    fprintf(fileRhoStd,"%.10lg\n",sqrt(rhoStd[iGrid])*pre);
  }
  fclose(fileRhoStd);
  fflush(stdout);
  exit(0);
  */

  //TEST Kinetic energy matrix
  /*
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    coeffForceReUp[iCoeff] = 0.0;
    coeffForceImUp[iCoeff] = 0.0;
  }
  double *keMatrix = (double*)calloc(numStatesDet*numStatesDet,sizeof(double));
  double *ak2Kinetic = cp->cpewald.ak2Kinetic;
  double temp;
  FILE *fkeMat = fopen("ke-mat","w");
  FILE *fwfDetTemp = fopen("wf-k","w");
  for(iState=0;iState<numStatesDet;iState++){
  //printf("ak2_sm[1] %lg\n",ak2_sm[1]);
    iOff = iState*numCoeff;
    for(iCoeff=0;iCoeff<numCoeff;iCoeff++){
      fprintf(fwfDetTemp,"%.16lg %.16lg\n",wfDetBackupUpRe[iOff+iCoeff],wfDetBackupUpIm[iOff+iCoeff]);
    }
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      coeffForceReUp[iCoeff] = 0.5*ak2Kinetic[iCoeff]*wfDetBackupUpRe[iOff+iCoeff-1];
      coeffForceImUp[iCoeff] = 0.5*ak2Kinetic[iCoeff]*wfDetBackupUpIm[iOff+iCoeff-1];
      //eke += (2.0*ak2Kinetic[i]*(ccreal[iis]*ccreal[iis] + ccimag[iis]*ccimag[iis]));
    }//endfor iCoeff
    for(jState=iState;jState<numStatesDet;jState++){
      temp = 0.0;
      iOff2 = jState*numCoeff;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
        temp += 2.0*(coeffForceReUp[iCoeff]*wfDetBackupUpRe[iOff2+iCoeff-1]+
                coeffForceImUp[iCoeff]*wfDetBackupUpIm[iOff2+iCoeff-1]);
      }//endfor iCoeff
      keMatrix[iState*numStatesDet+jState] = temp;
      keMatrix[jState*numStatesDet+iState] = temp;
    }//endfor jState
  }//endfor iState
  for(iState=0;iState<numStatesDet;iState++){
    for(jState=0;jState<numStatesDet;jState++){
      fprintf(fkeMat,"%.8lg ",keMatrix[iState*numStatesDet+jState]);
    }
    fprintf(fkeMat,"\n");
  }
  fclose(fkeMat);
  fclose(fwfDetTemp);
  fflush(stdout);
  exit(0);
  */
  
  //printf("numStateUpProc %i\n",numStateUpProc);
  for(iState=0;iState<numStateUpProc;iState++){
    iOff = iState*numCoeff;
    memcpy(wfReTemp,&coeffReUp[iOff+1],numCoeff*sizeof(double));
    memcpy(wfImTemp,&coeffImUp[iOff+1],numCoeff*sizeof(double));
    for(iChem=0;iChem<numChemPot;iChem++){
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
  }

  free(wfReTemp);
  free(wfImTemp);
  free(wfReProjTemp);
  free(wfImProjTemp);
  //for(iChem=0;iChem<numChemPot;iChem++){
  //  free(ewStateMap[iChem]);
  //}
  //free(ewStateMap);
  //free(ewStateNum);
  free(ksEigv);
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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genStoOrbitalEnergyWindowFrag(CLASS *class,GENERAL_DATA *general_data,
            CP *cp,GENERAL_DATA *generalDataMini,
            CP *cpMini,CLASS *classMini,int ip_now)
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
  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;
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
  int fragWindowFlag = stodftInfo->fragWindowFlag;
  int homoIndex;
  int smearOpt        = stodftInfo->smearOpt;

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

  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;
  double *coeffForceReUp = cpcoeffs_pos->fcre_up;
  double *coeffForceImUp = cpcoeffs_pos->fcim_up;
  double *chemPot = stodftCoefPos->chemPot;

  double *wfDetBackupUpRe = stodftCoefPos->wfDetBackupUpRe;
  double *wfDetBackupUpIm = stodftCoefPos->wfDetBackupUpIm;
  double *wfReTemp,*wfImTemp;
  double *wfReProjTemp,*wfImProjTemp;
  double *ksEigv = (double*)cmalloc(numStatesDet*sizeof(double));
  double *hmatCP = cell->hmat_cp;
 
  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;


  double timeStart1,timeEnd1;
  double timeStart2,timeEnd2;
  double timeStart3,timeEnd3;
  double timeStart4,timeEnd4;
  double timeStart5,timeEnd5;
  double timeStart6,timeEnd6;
  double diffTime1 = 0.0;
  double diffTime2 = 0.0;
  double diffTime3 = 0.0;
  double diffTime4 = 0.0;
  double diffTime5 = 0.0;
  double diffTime6 = 0.0;

  /*
  double volCPRev  = 1.0/getdeth(hmatCP);
  */
/*======================================================================*/
/* I) Set flags                                                         */

  *forceFormUp = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1; // temp keep this flag for debug filter
  if(cpLsda==1&&numStateDnProc!=0){
    *forceFormDn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }
  
  //genChemPotInterpPoints(stodftInfo,stodftCoefPos);
  //combineStoUCEnergyWindow(cp,general_data,class,cpMini,
  //             generalDataMini,classMini,ip_now);
  //energyCorrect(cpMini,generalDataMini,classMini,cp,class,ip_now);    
  //numChemPot -= 1;

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
  if(numProcStates>1)Barrier(commStates);
  stodftInfo->energyDiff = energyMax-energyMin;
  energyDiff = stodftInfo->energyDiff;
  stodftInfo->energyMean = 0.5*(energyMin+energyMax);

/*======================================================================*/
/* III) Calculate chemPot Values for Projectors                         */

  genChemPotInterpPoints(stodftInfo,stodftCoefPos);

/*======================================================================*/
/* IV) Generate Length of Polynomial Chain                              */

  timeStart2 = omp_get_wtime();
  if(expanType==1){
    Smin = chebyshevInfo->Smin;
    Smax = chebyshevInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    chebyshevInfo->scale = scale;
    //cheat my code
    if(myidState==0)genChebyHermit(stodftInfo,stodftCoefPos,2);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      polynormLength = stodftInfo->polynormLength;
      totalPoly = polynormLength*numChemPot;
      if(myidState!=0){
        stodftCoefPos->expanCoeff = (double*)crealloc(stodftCoefPos->expanCoeff,
                                totalPoly*sizeof(double));
        if(smearOpt>0){
          stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));
        }
      }
      Barrier(commStates);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      if(smearOpt>0){
        Bcast(stodftCoefPos->entropyExpanCoeff,polynormLength,MPI_DOUBLE,0,commStates);
      }
      Barrier(commStates);
    }
    //printf("numChemPotTemp %i\n",numChemPotTemp);
  }
  if(expanType==2){
    /*
    Smin = newtonInfo->Smin;
    Smax = newtonInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    newtonInfo->scale = scale;
    //cheat my code
    if(myidState==0)genNewtonHermit(stodftInfo,stodftCoefPos);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      Barrier(commStates);
    }
    */
    printf("We don't support Newton's Interpolation in EW+Frag scheme.\n");
    fflush(stdout);
    exit(0);
    //printf("numChemPotTemp %i\n",numChemPotTemp);
  }
  timeEnd2 = omp_get_wtime();
  diffTime2 = timeEnd2-timeStart2;

/*======================================================================*/
/* V) Project Noise Orbitals with Projectors                            */
/*    Reminder: Chebyshev Moments are also generated in this step.      */


  // 1. Generate noise orbitals
  genNoiseOrbitalReal(cp,cpcoeffs_pos);

  // 2. Filtering
 
  switch(expanType){
    case 1:
      stodftInfo->storeChebyMomentsFlag = 1;
      filterChebyPolyHerm(cp,class,general_data,ip_now);
      stodftInfo->storeChebyMomentsFlag = 0;
      break;
    default:
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("Only Chebyshev polynormial can be used for ew+frag!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }



/*======================================================================*/
/* VI) Project Stochastic Orbitals with Fragment Density Matrix         */
 
  combineStoUCEnergyWindow(cp,general_data,class,cpMini,
               generalDataMini,classMini,ip_now);

/*======================================================================*/
/* VII) Generate Energy and Force Corrections                           */

  energyCorrect(cpMini,generalDataMini,classMini,cp,class,ip_now);    

/*======================================================================*/
/* IX) Generate Correct Chemical Potential                              */

  calcChemPotChebyEWFrag(cp,class,general_data,ip_now);

/*======================================================================*/
/* X) Redo the Polynormial Fitting with Correct Chemical Potential      */

  if(expanType==1){
    if(myidState==0)genChebyHermitTrueChemPot(stodftInfo,stodftCoefPos,3);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      polynormLength = stodftInfo->polynormLength;
      totalPoly = polynormLength*numChemPot;
      if(myidState!=0){
        stodftCoefPos->expanCoeff = (double*)crealloc(stodftCoefPos->expanCoeff,
                            totalPoly*sizeof(double));
      if(smearOpt>0){
          stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));
       }
      }
      Barrier(commStates);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      if(smearOpt>0){
        Bcast(stodftCoefPos->entropyExpanCoeff,polynormLength,MPI_DOUBLE,0,commStates);
      }
      Barrier(commStates);
    }
  }   
        
/*======================================================================*/
/* XI) Redo filtering with \sqrt{F*P_i}                                 */

  // 1. Generate noise orbitals
  genNoiseOrbitalReal(cp,cpcoeffs_pos);
  
  // 2. Filtering
  switch(expanType){
    case 1:
      filterChebyPolyHerm(cp,class,general_data,ip_now);
      break;
  }

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genStoOrbitalEnergyWindowFragFake(CLASS *class,GENERAL_DATA *general_data,
            CP *cp,GENERAL_DATA *generalDataMini,
            CP *cpMini,CLASS *classMini,int ip_now)
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
  int fragWindowFlag = stodftInfo->fragWindowFlag;
  int homoIndex;

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
  double *chemPot = stodftCoefPos->chemPot;

  double *wfDetBackupUpRe = stodftCoefPos->wfDetBackupUpRe;
  double *wfDetBackupUpIm = stodftCoefPos->wfDetBackupUpIm;
  double *wfReTemp,*wfImTemp;
  double *wfReProjTemp,*wfImProjTemp;
  double *ksEigv = (double*)cmalloc(numStatesDet*sizeof(double));
  double *hmatCP = cell->hmat_cp;
 
  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  FILE *feigv = fopen("orbital-e","r");
  // TEST. Let's calculate the theoretical estimation
  //FILE *frealwf = fopen("wf-real","r");
  double *wfDetReal = stodftCoefPos->wfDetReal;

  /*
  double volCPRev  = 1.0/getdeth(hmatCP);
  for(iState=0;iState<numStatesDet;iState++){
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      fscanf(frealwf,"%lg",&wfDetReal[iState*rhoRealGridTot+iGrid]);
    }
  }
  fclose(frealwf);
  */

  
  int junk;
  for(iState=0;iState<numStatesDet;iState++){
    fscanf(feigv,"%i",&junk);
    fscanf(feigv,"%lg",&ksEigv[iState]);
  }
  fclose(feigv);
  stodftInfo->energyMin = ksEigv[0]-0.1;
  
/*======================================================================*/
/* III) Set flags                           */

  *forceFormUp = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1; // temp keep this flag for debug filter
  if(cpLsda==1&&numStateDnProc!=0){
    *forceFormDn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }
  
  genChemPotInterpPoints(stodftInfo,stodftCoefPos);
  combineStoUCEnergyWindow(cp,general_data,class,cpMini,
               generalDataMini,classMini,ip_now);
  energyCorrect(cpMini,generalDataMini,classMini,cp,class,ip_now);    
  //numChemPot -= 1;

/*======================================================================*/
/* IV) Generate random orbital                                          */

  genNoiseOrbitalReal(cp,cpcoeffs_pos);


/*======================================================================*/
/* II) Filtering by deterministic orbitals                     */

  wfReTemp = (double*)cmalloc(numCoeff*sizeof(double));
  wfImTemp = (double*)cmalloc(numCoeff*sizeof(double));
  wfReProjTemp = (double*)cmalloc(numCoeff*sizeof(double));
  wfImProjTemp = (double*)cmalloc(numCoeff*sizeof(double));

  /*
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iCoeff=1;iCoeff<=numCoeffUpTot;iCoeff++){
      stoWfUpRe[iChem][iCoeff] = 0.0;
      stoWfUpIm[iChem][iCoeff] = 0.0;
    }
  }
  */
   
  //numChemPot = 1000; // TEST
  int windowIndex;
  int index;
  stodftCoefPos->ewStateNum = (int*)cmalloc(numChemPot*sizeof(int));
  stodftCoefPos->ewStateMap = (int**)cmalloc(numChemPot*sizeof(int*));
  int *ewStateNum = stodftCoefPos->ewStateNum;
  int **ewStateMap = stodftCoefPos->ewStateMap;
  for(iChem=0;iChem<numChemPot;iChem++){
    ewStateMap[iChem] = (int*)cmalloc(numStatesDet*sizeof(int));
  }
  double energyMin = stodftInfo->energyMin;
  energyMin = ksEigv[0]-0.1;
  printf("energyMin %lg\n",energyMin);
  double dE = (chemPot[numChemPot-2]-energyMin)/(numChemPot-1.0);
  /*
  for(iState=0;iState<numStatesDet;iState++){
    if(ksEigv[iState]<chemPot[numChemPot-2]){
      windowIndex = (int)((ksEigv[iState]-energyMin)/dE);
      ewStateMap[windowIndex][ewStateNum[windowIndex]] = iState;
      ewStateNum[windowIndex] += 1;
      if(iState==numStatesDet-1)homoIndex = windowIndex;
    }
    else{
      ewStateMap[numChemPot-1][ewStateNum[numChemPot-1]] = iState;
      ewStateNum[numChemPot-1] += 1;
      if(iState==numStatesDet-1)homoIndex = numChemPot-1;
    }
  }
  */
  homoIndex = stodftInfo->homoIndex;
  
  printf("homoIndex %i\n",homoIndex);
  printf("chemPot window %lg HOMO %lg HOMO window index %i\n",chemPot[numChemPot-2],ksEigv[numStatesDet-1],homoIndex);
  printf("test stowf %lg %lg %lg %lg\n",stoWfUpRe[homoIndex][1],stoWfUpRe[homoIndex][100],stoWfUpIm[homoIndex][1],stoWfUpIm[homoIndex][100]);
  for(iState=0;iState<numStateUpProc;iState++){
    iOff = iState*numCoeff;
    memcpy(wfReTemp,&stoWfUpRe[homoIndex][iOff+1],numCoeff*sizeof(double));
    memcpy(wfImTemp,&stoWfUpIm[homoIndex][iOff+1],numCoeff*sizeof(double));
    #pragma omp parallel for private(iCoeff)
    for(iCoeff=0;iCoeff<numCoeff;iCoeff++){
      stoWfUpRe[homoIndex][iOff+iCoeff+1] = 0.0;
      stoWfUpIm[homoIndex][iOff+iCoeff+1] = 0.0;
    }
    for(jState=0;jState<numStatesDet;jState++){
      dot = 0.0;
      iOff2 = jState*numCoeff;
      #pragma omp parallel for reduction(+:dot) private(iCoeff)
      for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
        dot += wfDetBackupUpRe[iOff2+iCoeff]*wfReTemp[iCoeff]+
               wfDetBackupUpIm[iOff2+iCoeff]*wfImTemp[iCoeff];
      }
      dot = (dot*2.0+wfDetBackupUpRe[iOff2+numCoeff-1]*wfReTemp[numCoeff-1])*0.5;
      #pragma omp parallel for private(iCoeff)
      for(iCoeff=0;iCoeff<numCoeff;iCoeff++){
        stoWfUpRe[homoIndex][iOff+iCoeff+1] += wfDetBackupUpRe[iOff2+iCoeff]*dot;
        stoWfUpIm[homoIndex][iOff+iCoeff+1] += wfDetBackupUpIm[iOff2+iCoeff]*dot;
      }
    }//endfor jState
  }//endfor iState
  
  for(iChem=homoIndex+1;iChem<numChemPot;iChem++){
    #pragma omp parallel for private(iCoeff)
    for(iCoeff=0;iCoeff<numCoeffUpTot;iCoeff++){
      stoWfUpRe[iChem][iCoeff+1] = 0.0;
      stoWfUpIm[iChem][iCoeff+1] = 0.0;
    }
  }

  printf("test stowf %lg %lg %lg %lg\n",stoWfUpRe[homoIndex][1],stoWfUpRe[homoIndex][100],stoWfUpIm[homoIndex][1],stoWfUpIm[homoIndex][100]);


  free(wfReTemp);
  free(wfImTemp);
  free(wfReProjTemp);
  free(wfImProjTemp);
  //for(iChem=0;iChem<numChemPot;iChem++){
  //  free(ewStateMap[iChem]);
  //}
  //free(ewStateMap);
  //free(ewStateNum);
  free(ksEigv);
  
  FILE *ftest = fopen("wf-test","w");
  for(iCoeff=1;iCoeff<=numStateUpProc*numCoeff;iCoeff++){
    fprintf(ftest,"%.16lg %.16lg\n",stoWfUpRe[homoIndex][iCoeff],stoWfUpIm[homoIndex][iCoeff]);
  }
  /*
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

#ifdef FAST_FILTER   
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genStoOrbitalChebyTest(CLASS *class,GENERAL_DATA *general_data,
	                       CP *cp,
                               CLASS *class2,GENERAL_DATA *general_data2,
                               CP *cp2,
                               int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS  *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR	       *cpscr	     = &(cp->cpscr);
  CPSCR        *cpscr2       = &(cp2->cpscr);

  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTINFO   *stodftInfo2   = cp2->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  STODFTCOEFPOS *stodftCoefPos2 = cp2->stodftCoefPos;

  NEWTONINFO   *newtonInfo      = stodftInfo->newtonInfo;
  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;

  COMMUNICATE   *commCP	        = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_INFO *cpcoeffs_info2  = &(cp2->cpcoeffs_info);

  CPCOEFFS_POS  *cpcoeffs_pos   = &(cp->cpcoeffs_pos[ip_now]);
  CPCOEFFS_POS  *cpcoeffs_pos2   = &(cp2->cpcoeffs_pos[ip_now]);

  MPI_Comm commStates = commCP->comm_states;
  

  int iPoly,iState,iState2,iCoeff,iChem;
  int numChemPot = stodftInfo->numChemPot;
  int numChemPot2 = stodftInfo2->numChemPot;
  int polynormLength = stodftInfo->polynormLength;
  int expanType = stodftInfo->expanType;
  int numProcStates = commCP->np_states;
  int myidState = commCP->myid_state;
  int cpLsda = cpopts->cp_lsda;
  int cpGGA  = cpopts->cp_gga;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numStateUpProc2 = cpcoeffs_info2->nstate_up_proc;
  int numStateDnProc2 = cpcoeffs_info2->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int totalPoly;
  int numCoeffUpTot   = numStateUpProc*numCoeff;
  int numCoeffDnTot   = numStateDnProc*numCoeff;
  int rhoRealGridTot    = stodftInfo->rhoRealGridTot;
  int smearOpt = stodftInfo->smearOpt;

  int numStatePrintUp = stodftInfo2->numStatePrintUp;
  int numStatePrintDn = stodftInfo2->numStatePrintDn;

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

  //double *sampPoint = newtonInfo->sampPoint;
  //double *sampPointUnscale = newtonInfo->sampPointUnscale;
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
  //double *scrCoeffReDn   = cpscr->cpscr_wave.cre_dn;
  //double *scrCoeffImDn   = cpscr->cpscr_wave.cim_dn;
  double *rhoCoeffReDn   = cpscr->cpscr_rho.rhocr_dn;
  double *rhoCoeffImDn   = cpscr->cpscr_rho.rhoci_dn;
  double *rhoDn          = cpscr->cpscr_rho.rho_dn;
  double *rhoCoeffReDnDensCpBox = cpscr->cpscr_rho.rhocr_dn_dens_cp_box;
  double *rhoCoeffImDnDensCpBox = cpscr->cpscr_rho.rhoci_dn_dens_cp_box;
  double *divRhoxDn       = cpscr->cpscr_grho.d_rhox_dn;
  double *divRhoyDn       = cpscr->cpscr_grho.d_rhoy_dn;
  double *divRhozDn       = cpscr->cpscr_grho.d_rhoz_dn;
  double *d2RhoDn         = cpscr->cpscr_grho.d2_rho_dn;

  double *vksUp           = cpscr->cpscr_rho.v_ks_up;
  double *vksUp2          = cpscr2->cpscr_rho.v_ks_up;
  double *vksDn           = cpscr->cpscr_rho.v_ks_dn;
  double *vksDn2          = cpscr2->cpscr_rho.v_ks_dn;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  double *coeffReUpBackup = stodftCoefPos->coeffReUpBackup;
  double *coeffImUpBackup = stodftCoefPos->coeffImUpBackup;

  //timing
  double timeStart1,timeEnd1;
  double timeStart2,timeEnd2;
  double timeStart3,timeEnd3;
  double timeStart4,timeEnd4;
  double timeStart5,timeEnd5;
  double timeStart6,timeEnd6;
  double timeStart01,timeEnd01;
  double timeStart02,timeEnd02;
  double timeStart03,timeEnd03;
  double timeStart04,timeEnd04;
  double diffTime1 = 0.0;
  double diffTime2 = 0.0;
  double diffTime3 = 0.0;
  double diffTime4 = 0.0;
  double diffTime5 = 0.0;
  double diffTime6 = 0.0;
  double diffTime01 = 0.0;
  double diffTime02 = 0.0;
  double diffTime03 = 0.0;
  double diffTime04 = 0.0;

  double timeStartAll,timeEndAll;
  double diffTimeAll;


/*======================================================================*/
/* I) Copy the KS potential */

  memcpy(&vksUp2[1],&vksUp[1],rhoRealGridTot*sizeof(double));
  if(cpLsda==1&&numStateDnProc>0){
    memcpy(&vksDn2[1],&vksDn[1],rhoRealGridTot*sizeof(double));
  }

/*======================================================================*/
/* I) Filter Diag */

  // free stodft part 
  for(iChem=0;iChem<numChemPot;iChem++){
    cfree(&(stoWfUpRe[iChem][1]));
    cfree(&(stoWfUpIm[iChem][1]));
  }
  cfree(&coeffReUp[1]);
  cfree(&coeffImUp[1]);
  cfree(&(cpcoeffs_pos->fcre_up[1]));
  cfree(&(cpcoeffs_pos->fcim_up[1]));

  // allocate filter diag part
  for(iChem=0;iChem<numChemPot2;iChem++){
    stodftCoefPos2->stoWfUpRe[iChem] = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
    stodftCoefPos2->stoWfUpIm[iChem] = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  }
  cpcoeffs_pos2->cre_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  cpcoeffs_pos2->cim_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  cpcoeffs_pos2->fcre_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  cpcoeffs_pos2->fcim_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;

  stodftInfo2->iScf = stodftInfo->iScf;
  stodftInfo2->iScfTrue = stodftInfo->iScfTrue;
  timeStart01 = omp_get_wtime();
  genStoOrbitalInterp(class2,general_data2,cp2,ip_now);
  timeEnd01 = omp_get_wtime();
  diffTime01 = timeEnd01-timeStart01;

  timeStart02 = omp_get_wtime();
  orthDiagDriver(cp2,class2,general_data2,ip_now);
  timeEnd02 = omp_get_wtime();
  diffTime02 = timeEnd02-timeStart02;

  timeStart03 = omp_get_wtime();
  broadcastWfDet(cp2,class2,general_data2,cp); 
  timeEnd03 = omp_get_wtime();
  diffTime03 = timeEnd03-timeStart03;

  // free filter diag part
  // Done in filter diag
 
  // allocate stodft part
  for(iChem=0;iChem<numChemPot;iChem++){
    stoWfUpRe[iChem] = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
    stoWfUpIm[iChem] = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  }
  cpcoeffs_pos->cre_up = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  cpcoeffs_pos->cim_up = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  cpcoeffs_pos->fcre_up = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  cpcoeffs_pos->fcim_up = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  coeffReUp = cpcoeffs_pos->cre_up;
  coeffImUp = cpcoeffs_pos->cim_up;
  
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

  timeStart1 = omp_get_wtime();
  //if(myidState==0){
  genEnergyMax(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
  //fflush(stdout);
  //exit(0);
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
  if(numProcStates>1)Barrier(commStates);
  stodftInfo->energyDiff = energyMax-energyMin;
  energyDiff = stodftInfo->energyDiff;
  stodftInfo->energyMean = 0.5*(energyMin+energyMax);
  timeEnd1 = omp_get_wtime();
  diffTime1 = timeEnd1-timeStart1;

/*======================================================================*/
/* III) Generate Length of Polynomial Chain		                */
  
  timeStart2 = omp_get_wtime();
  if(expanType==1){
    Smin = chebyshevInfo->Smin;
    Smax = chebyshevInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    chebyshevInfo->scale = scale;
    //cheat my code
    stodftInfo->numChemPot = 2;
    free(&(stodftCoefPos->chemPot[0]));
    stodftCoefPos->chemPot = (double*)cmalloc(2*sizeof(double));
    stodftCoefPos->chemPot[0] = chemPotInit-gapInit*0.5;
    stodftCoefPos->chemPot[1] = chemPotInit+gapInit*0.5;
    if(myidState==0)genChebyHermit(stodftInfo,stodftCoefPos,0);
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
  timeEnd2 = omp_get_wtime();
  diffTime2 = timeEnd2-timeStart2;

/*======================================================================*/
/* IV) Calculate the True Chemical Potential                            */
  
  timeStart3 = omp_get_wtime();
  calcChemPotCheby(cp,class,general_data,ip_now);
  timeEnd3 = omp_get_wtime();
  diffTime3 = timeEnd3-timeStart3;

/*======================================================================*/
/* V) Generate Coeffcients for Polynormial interpolation with Correct   */
/*    Chemical Potential.						*/
  
  timeStart4 = omp_get_wtime();
  if(expanType==1){
    Smin = chebyshevInfo->Smin;
    Smax = chebyshevInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    chebyshevInfo->scale = scale;
    if(myidState==0)genChebyHermitTrueChemPot(stodftInfo,stodftCoefPos,0);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      polynormLength = stodftInfo->polynormLength;
      totalPoly = polynormLength*numChemPot;
      if(myidState!=0){
        stodftCoefPos->expanCoeff = (double*)crealloc(stodftCoefPos->expanCoeff,
                            totalPoly*sizeof(double));
        if(smearOpt>0){
          stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));
        }
      }
      Barrier(commStates);
      //printf("%p\n",stodftCoefPos->expanCoeff);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      //printf("smearOpt %i\n",smearOpt);
      if(smearOpt>0){
        Bcast(stodftCoefPos->entropyExpanCoeff,polynormLength,MPI_DOUBLE,0,commStates);
      }
    }
    //printf("11111111111111111\n");
  }
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
        if(smearOpt>0){
          stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));
        }
      }
      Barrier(commStates);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPoint,polynormLength,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPointUnscale,polynormLength,MPI_DOUBLE,0,commStates);
      if(smearOpt>0){
        Bcast(stodftCoefPos->entropyExpanCoeff,polynormLength,MPI_DOUBLE,0,commStates);
      }
      Barrier(commStates);
    }
    
    /*
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      sampPointUnscale[iPoly] = (sampPoint[iPoly]-Smin)/scale+energyMin;
      //printf("sampunscale %lg\n",sampPointUnscale[iPoly]);
    }
    genCoeffNewtonHermit(stodftInfo,stodftCoefPos);
    */
  }
  timeEnd4 = omp_get_wtime();
  diffTime4 = timeEnd4-timeStart4;


/*======================================================================*/
/* IV) Generate random orbital                                          */

  timeStart5 = omp_get_wtime();
  genNoiseOrbitalReal(cp,cpcoeffs_pos);
  timeEnd5 = omp_get_wtime();
  diffTime5 = timeEnd5-timeStart5;

/*======================================================================*/
/* V) Filter the stochastic orbitals			*/

  timeStart6 = omp_get_wtime();
  switch(expanType){
    case 1:
      filterChebyPolyHermFake(cp,class,general_data,ip_now,0);
      break;
    case 2:
      filterNewtonPolyHermFake(cp,class,general_data,ip_now);
      break;
  }
  stodftInfo->filterFlag = 0;
  timeEnd6 = omp_get_wtime();
  diffTime6 = timeEnd6-timeStart6;

  printf("Gen-stowf time myid %i spec-range %.8lg gen-poly-length %.8lg gen-chempot %.8lg gen-poly-coeff %.8lg gen-rand %.8lg filter %.8lg fd_filter %.8lg fd_diag %.8lg fd_bcast %.8lg\n",myidState,diffTime1,diffTime2,diffTime3,diffTime4,diffTime5,diffTime6,diffTime01,diffTime02,diffTime03);


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
void genStoOrbitalEnergyWindowTest(CLASS *class,GENERAL_DATA *general_data,
                               CP *cp,
                               CLASS *class2,GENERAL_DATA *general_data2,
                               CP *cp2,
                               int ip_now)
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
  STODFTCOEFPOS *stodftCoefPos2 = cp2->stodftCoefPos;

  NEWTONINFO   *newtonInfo      = stodftInfo->newtonInfo;
  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;
  COMMUNICATE   *commCP	        = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_INFO *cpcoeffs_info2 = &(cp2->cpcoeffs_info);

  CPCOEFFS_POS  *cpcoeffs_pos   = &(cp->cpcoeffs_pos[ip_now]);
  CPCOEFFS_POS  *cpcoeffs_pos2  = &(cp2->cpcoeffs_pos[ip_now]);

  MPI_Comm commStates = commCP->comm_states;

  CPSCR        *cpscr2       = &(cp2->cpscr);
  STODFTINFO   *stodftInfo2   = cp2->stodftInfo;

  int iPoly,iState,iState2,iCoeff,iChem;
  int numChemPot = stodftInfo->numChemPot;
  int numChemPotTemp;
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
  int numCoeffUpTotal   = numStateUpProc*numCoeff;
  int numCoeffDnTotal   = numStateDnProc*numCoeff;
  int rhoRealGridTot    = stodftInfo->rhoRealGridTot;
  int smearOpt = stodftInfo->smearOpt;

  int numChemPot2 = stodftInfo2->numChemPot;
  int numStateUpProc2 = cpcoeffs_info2->nstate_up_proc;
  int numStateDnProc2 = cpcoeffs_info2->nstate_dn_proc;

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

  double *vksUp           = cpscr->cpscr_rho.v_ks_up;
  double *vksUp2          = cpscr2->cpscr_rho.v_ks_up;
  double *vksDn           = cpscr->cpscr_rho.v_ks_dn;
  double *vksDn2          = cpscr2->cpscr_rho.v_ks_dn;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  double *coeffReUpBackup = stodftCoefPos->coeffReUpBackup;
  double *coeffImUpBackup = stodftCoefPos->coeffImUpBackup;

  //timing
  double timeStart1,timeEnd1;
  double timeStart2,timeEnd2;
  double timeStart3,timeEnd3;
  double timeStart4,timeEnd4;
  double timeStart5,timeEnd5;
  double timeStart6,timeEnd6;
  double diffTime1 = 0.0;
  double diffTime2 = 0.0;
  double diffTime3 = 0.0;
  double diffTime4 = 0.0;
  double diffTime5 = 0.0;
  double diffTime6 = 0.0;


/*======================================================================*/
/* I) Copy the KS potential */

  memcpy(&vksUp2[1],&vksUp[1],rhoRealGridTot*sizeof(double));
  if(cpLsda==1&&numStateDnProc>0){
    memcpy(&vksDn2[1],&vksDn[1],rhoRealGridTot*sizeof(double));
  }

/*======================================================================*/
/* I) Filter Diag */
  // free stodft part 
  for(iChem=0;iChem<numChemPot;iChem++){
    cfree(&(stoWfUpRe[iChem][1]));
    cfree(&(stoWfUpIm[iChem][1]));
  }
  cfree(&coeffReUp[1]);
  cfree(&coeffImUp[1]);
  cfree(&(cpcoeffs_pos->fcre_up[1]));
  cfree(&(cpcoeffs_pos->fcim_up[1]));

  // allocate filter diag part
  for(iChem=0;iChem<numChemPot2;iChem++){
    stodftCoefPos2->stoWfUpRe[iChem] = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
    stodftCoefPos2->stoWfUpIm[iChem] = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  }
  cpcoeffs_pos2->cre_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  cpcoeffs_pos2->cim_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  cpcoeffs_pos2->fcre_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  cpcoeffs_pos2->fcim_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;

  stodftInfo2->iScf = stodftInfo->iScf;
  stodftInfo2->iScfTrue = stodftInfo->iScfTrue;
  genStoOrbitalInterp(class2,general_data2,cp2,ip_now);
  orthDiagDriver(cp2,class2,general_data2,ip_now);

  broadcastWfDet(cp2,class2,general_data2,cp);

  // free filter diag part
  // Done in filter diag

  // allocate stodft part
  for(iChem=0;iChem<numChemPot;iChem++){
    stoWfUpRe[iChem] = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
    stoWfUpIm[iChem] = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  }
  cpcoeffs_pos->cre_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  cpcoeffs_pos->cim_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  cpcoeffs_pos->fcre_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  cpcoeffs_pos->fcim_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  coeffReUp = cpcoeffs_pos->cre_up;
  coeffImUp = cpcoeffs_pos->cim_up;

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
  if(numProcStates>1)Barrier(commStates);
  stodftInfo->energyDiff = energyMax-energyMin;
  energyDiff = stodftInfo->energyDiff;
  stodftInfo->energyMean = 0.5*(energyMin+energyMax);

/*======================================================================*/
/* III) Generate Length of Polynomial Chain                             */

  timeStart2 = omp_get_wtime();
  if(expanType==1){
    Smin = chebyshevInfo->Smin;
    Smax = chebyshevInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    chebyshevInfo->scale = scale;
    //cheat my code
    numChemPotTemp = stodftInfo->numChemPot;
    stodftInfo->numChemPot = 2;
    free(&(stodftCoefPos->chemPot[0]));
    stodftCoefPos->chemPot = (double*)cmalloc(2*sizeof(double));
    stodftCoefPos->chemPot[0] = chemPotInit-gapInit*0.5;
    stodftCoefPos->chemPot[1] = chemPotInit+gapInit*0.5;
    if(myidState==0)genChebyHermit(stodftInfo,stodftCoefPos,0);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      Barrier(commStates);
    }
    stodftInfo->numChemPot = numChemPotTemp;
    free(&(stodftCoefPos->chemPot[0]));
    stodftCoefPos->chemPot = (double*)cmalloc(numChemPotTemp*sizeof(double));
    //finish cheating my code
  }
  if(expanType==2){
    Smin = newtonInfo->Smin;
    Smax = newtonInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    newtonInfo->scale = scale;
    //cheat my code
    numChemPotTemp = stodftInfo->numChemPot;
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
    //printf("numChemPotTemp %i\n",numChemPotTemp);
    stodftInfo->numChemPot = numChemPotTemp;
    free(&(stodftCoefPos->chemPot[0]));
    stodftCoefPos->chemPot = (double*)cmalloc(numChemPotTemp*sizeof(double));
    //finish cheating my code
  }
  timeEnd2 = omp_get_wtime();
  diffTime2 = timeEnd2-timeStart2;

/*======================================================================*/
/* IV) Calculate the True Chemical Potential                            */

  timeStart3 = omp_get_wtime();
  calcChemPotCheby(cp,class,general_data,ip_now);
  timeEnd3 = omp_get_wtime();
  diffTime3 = timeEnd3-timeStart3;

/*======================================================================*/
/* V) Generate list of chemical potential for energy windows.           */

  genChemPotInterpPoints(stodftInfo,stodftCoefPos);

/*======================================================================*/
/* V) Generate Coeffcients for Polynormial interpolation                */
  
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
        if(smearOpt>0){
          stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));
        }
      }
      Barrier(commStates);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPoint,polynormLength,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPointUnscale,polynormLength,MPI_DOUBLE,0,commStates);
      if(smearOpt>0){
        Bcast(stodftCoefPos->entropyExpanCoeff,polynormLength,MPI_DOUBLE,0,commStates);
      }
      Barrier(commStates);
    }
  }

/*======================================================================*/
/* IV) Generate random orbital                                          */

  genNoiseOrbitalReal(cp,cpcoeffs_pos);
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
      //filterNewtonPolyHerm(cp,class,general_data,ip_now);
      filterNewtonPolyHermFake(cp,class,general_data,ip_now);
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
  
/*======================================================================*/
/* VI) Calculate filtered states w.r.t. energy windows                  */
    
  //Barrier(commStates);
  
  /*  
  if(myidState==0)printf("Filter by window...\n");
  for(iChem=numChemPot-1;iChem>0;iChem--){
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      stoWfUpRe[iChem][iCoeff] -= stoWfUpRe[iChem-1][iCoeff];
      stoWfUpIm[iChem][iCoeff] -= stoWfUpIm[iChem-1][iCoeff];
    }//endfor iCoeff
  }//endfor iChem
  if(cpLsda==1){
    for(iChem=numChemPot-1;iChem>0;iChem--){
      for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
	stoWfDnRe[iChem][iCoeff] -= stoWfDnRe[iChem-1][iCoeff];
	stoWfDnIm[iChem][iCoeff] -= stoWfDnIm[iChem-1][iCoeff];
      }//endfor iCoeff
    }//endfor iChem
  }
  */

/*======================================================================*/
/* VI) Calculate Energy	                                            */

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genStoOrbitalInterpTest(CLASS *class,GENERAL_DATA *general_data,CP *cp,
            CLASS *class2,GENERAL_DATA *general_data2,CP *cp2,
            int ip_now)
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
  STODFTCOEFPOS *stodftCoefPos2 = cp2->stodftCoefPos;

  NEWTONINFO   *newtonInfo      = stodftInfo->newtonInfo;
  COMMUNICATE   *commCP	        = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_INFO *cpcoeffs_info2 = &(cp2->cpcoeffs_info);

  CPCOEFFS_POS  *cpcoeffs_pos   = &(cp->cpcoeffs_pos[ip_now]);
  CPCOEFFS_POS  *cpcoeffs_pos2  = &(cp2->cpcoeffs_pos[ip_now]);

  MPI_Comm commStates = commCP->comm_states;
  
  CPSCR        *cpscr2       = &(cp2->cpscr);
  STODFTINFO   *stodftInfo2   = cp2->stodftInfo;

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
  int rhoRealGridTot    = stodftInfo->rhoRealGridTot;

  int numChemPot2 = stodftInfo2->numChemPot;
  int numStateUpProc2 = cpcoeffs_info2->nstate_up_proc;
  int numStateDnProc2 = cpcoeffs_info2->nstate_dn_proc;
  int smearOpt        = stodftInfo->smearOpt;

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

  double *vksUp           = cpscr->cpscr_rho.v_ks_up;
  double *vksUp2          = cpscr2->cpscr_rho.v_ks_up;
  double *vksDn           = cpscr->cpscr_rho.v_ks_dn;
  double *vksDn2          = cpscr2->cpscr_rho.v_ks_dn;

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
/* I) Copy the KS potential */

  memcpy(&vksUp2[1],&vksUp[1],rhoRealGridTot*sizeof(double));
  if(cpLsda==1&&numStateDnProc>0){
    memcpy(&vksDn2[1],&vksDn[1],rhoRealGridTot*sizeof(double));
  }

/*======================================================================*/
/* I) Filter Diag */
  // free stodft part 
  for(iChem=0;iChem<numChemPot;iChem++){
    cfree(&(stoWfUpRe[iChem][1]));
    cfree(&(stoWfUpIm[iChem][1]));
  }
  cfree(&coeffReUp[1]);
  cfree(&coeffImUp[1]);
  cfree(&(cpcoeffs_pos->fcre_up[1]));
  cfree(&(cpcoeffs_pos->fcim_up[1]));

  // allocate filter diag part
  for(iChem=0;iChem<numChemPot2;iChem++){
    stodftCoefPos2->stoWfUpRe[iChem] = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
    stodftCoefPos2->stoWfUpIm[iChem] = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  }
  cpcoeffs_pos2->cre_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  cpcoeffs_pos2->cim_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  cpcoeffs_pos2->fcre_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  cpcoeffs_pos2->fcim_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;

  stodftInfo2->iScf = stodftInfo->iScf;
  stodftInfo2->iScfTrue = stodftInfo->iScfTrue;
  genStoOrbitalInterp(class2,general_data2,cp2,ip_now);
  orthDiagDriver(cp2,class2,general_data2,ip_now);

  broadcastWfDet(cp2,class2,general_data2,cp);

  // free filter diag part
  // Done in filter diag

  // allocate stodft part
  for(iChem=0;iChem<numChemPot;iChem++){
    stoWfUpRe[iChem] = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
    stoWfUpIm[iChem] = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  }
  cpcoeffs_pos->cre_up = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  cpcoeffs_pos->cim_up = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  cpcoeffs_pos->fcre_up = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  cpcoeffs_pos->fcim_up = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  coeffReUp = cpcoeffs_pos->cre_up;
  coeffImUp = cpcoeffs_pos->cim_up;

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

  if(myidState==0){
    genEnergyMax(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
    genEnergyMin(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
  }
  if(numProcStates>1){
    Barrier(commStates);
    Bcast(&(stodftInfo->energyMin),1,MPI_DOUBLE,0,commStates);
    Bcast(&(stodftInfo->energyMax),1,MPI_DOUBLE,0,commStates);
  }
  energyMin = stodftInfo->energyMin;
  energyMax = stodftInfo->energyMax;
  //printf("energyMax %lg energyMin %lg\n",energyMax,energyMin);
  if(numProcStates>1)Barrier(commStates);
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
        int numStatePrintUp = stodftInfo->numStatePrintUp;
	double eigValMin = stodftCoefPos->energyLevel[0];
	double eigValMax = stodftCoefPos->energyLevel[numStatePrintUp-1];
	if(myidState==0){
	  printf("eigValMin %lg\n",eigValMin);
	  printf("eigValMax %lg\n",eigValMax);  
	}
	//stodftInfo->gapInit = eigValMax-eigValMin;
	stodftInfo->gapInit = eigValMax-energyMin;
	//stodftInfo->chemPotInit = 0.5*(eigValMax+eigValMin);
	stodftInfo->chemPotInit = 0.5*(eigValMax+energyMin);
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
        if(smearOpt>0){
          stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));
        }
      }
      Barrier(commStates);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPoint,polynormLength,MPI_DOUBLE,0,commStates);
      Bcast(newtonInfo->sampPointUnscale,polynormLength,MPI_DOUBLE,0,commStates);
      if(smearOpt>0){
          stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));
      }
      Barrier(commStates);
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
      filterNewtonPolyHermFake(cp,class,general_data,ip_now);
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
  
  
  /*
  double testfilter = 0.0;
  for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
    testfilter += stoWfUpRe[0][iCoeff]*stoWfUpRe[0][iCoeff]+
		stoWfUpIm[0][iCoeff]*stoWfUpIm[0][iCoeff];
  }
  testfilter += stoWfUpRe[0][numCoeff]*stoWfUpRe[0][numCoeff];
  printf("testfilter %lg\n",testfilter);
  */
  
  //Barrier(commStates);
    
  

/*======================================================================*/
/* VI) Calculate Energy	                                            */

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genStoOrbitalEnergyWindowFragTest(CLASS *class,GENERAL_DATA *general_data,
            CP *cp,GENERAL_DATA *generalDataMini,
            CP *cpMini,CLASS *classMini,
            CLASS *class2,GENERAL_DATA *general_data2,CP *cp2,int ip_now)
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
  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;
  COMMUNICATE   *commCP         = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_INFO *cpcoeffs_info2 = &(cp2->cpcoeffs_info);

  CPCOEFFS_POS  *cpcoeffs_pos   = &(cp->cpcoeffs_pos[ip_now]);
  CPCOEFFS_POS  *cpcoeffs_pos2  = &(cp2->cpcoeffs_pos[ip_now]);

  MPI_Comm commStates = commCP->comm_states;
  CELL *cell = &(general_data->cell);

  CPSCR        *cpscr2       = &(cp2->cpscr);
  STODFTINFO   *stodftInfo2   = cp2->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos2 = cp2->stodftCoefPos;

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
  int fragWindowFlag = stodftInfo->fragWindowFlag;
  int homoIndex;
  int smearOpt = stodftInfo->smearOpt;
  int smearOptBackup;

  int numChemPot2 = stodftInfo2->numChemPot;
  int numStateUpProc2 = cpcoeffs_info2->nstate_up_proc;
  int numStateDnProc2 = cpcoeffs_info2->nstate_dn_proc;

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

  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;
  double *coeffForceReUp = cpcoeffs_pos->fcre_up;
  double *coeffForceImUp = cpcoeffs_pos->fcim_up;
  double *chemPot = stodftCoefPos->chemPot;

  double *wfDetBackupUpRe = stodftCoefPos->wfDetBackupUpRe;
  double *wfDetBackupUpIm = stodftCoefPos->wfDetBackupUpIm;
  double *wfReTemp,*wfImTemp;
  double *wfReProjTemp,*wfImProjTemp;
  double *ksEigv = (double*)cmalloc(numStatesDet*sizeof(double));
  double *hmatCP = cell->hmat_cp;
 
  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  double *vksUp           = cpscr->cpscr_rho.v_ks_up;
  double *vksUp2          = cpscr2->cpscr_rho.v_ks_up;
  double *vksDn           = cpscr->cpscr_rho.v_ks_dn;
  double *vksDn2          = cpscr2->cpscr_rho.v_ks_dn;

  double timeStart1,timeEnd1;
  double timeStart2,timeEnd2;
  double timeStart3,timeEnd3;
  double timeStart4,timeEnd4;
  double timeStart5,timeEnd5;
  double timeStart6,timeEnd6;
  double diffTime1 = 0.0;
  double diffTime2 = 0.0;
  double diffTime3 = 0.0;
  double diffTime4 = 0.0;
  double diffTime5 = 0.0;
  double diffTime6 = 0.0;

  /*
  double volCPRev  = 1.0/getdeth(hmatCP);
  */
/*======================================================================*/
/* I) Copy the KS potential */

  memcpy(&vksUp2[1],&vksUp[1],rhoRealGridTot*sizeof(double));
  if(cpLsda==1&&numStateDnProc>0){
    memcpy(&vksDn2[1],&vksDn[1],rhoRealGridTot*sizeof(double));
  }

/*======================================================================*/
/* I) Filter Diag */

  // free stodft part 
  for(iChem=0;iChem<numChemPot;iChem++){
    cfree(&(stoWfUpRe[iChem][1]));
    cfree(&(stoWfUpIm[iChem][1]));
  }
  cfree(&coeffReUp[1]);
  cfree(&coeffImUp[1]);
  cfree(&(cpcoeffs_pos->fcre_up[1]));
  cfree(&(cpcoeffs_pos->fcim_up[1]));

  // allocate filter diag part
  for(iChem=0;iChem<numChemPot2;iChem++){
    stodftCoefPos2->stoWfUpRe[iChem] = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
    stodftCoefPos2->stoWfUpIm[iChem] = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  }
  cpcoeffs_pos2->cre_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  cpcoeffs_pos2->cim_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  cpcoeffs_pos2->fcre_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;
  cpcoeffs_pos2->fcim_up = (double*)cmalloc(numStateUpProc2*numCoeff*sizeof(double))-1;

  stodftInfo2->iScf = stodftInfo->iScf;
  stodftInfo2->iScfTrue = stodftInfo->iScfTrue;
  genStoOrbitalInterp(class2,general_data2,cp2,ip_now);
  orthDiagDriver(cp2,class2,general_data2,ip_now);

  broadcastWfDet(cp2,class2,general_data2,cp);

  // free filter diag part
  // Done in filter diag

  // allocate stodft part
  for(iChem=0;iChem<numChemPot;iChem++){
    stoWfUpRe[iChem] = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
    stoWfUpIm[iChem] = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  }
  cpcoeffs_pos->cre_up = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  cpcoeffs_pos->cim_up = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  cpcoeffs_pos->fcre_up = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  cpcoeffs_pos->fcim_up = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  coeffReUp = cpcoeffs_pos->cre_up;
  coeffImUp = cpcoeffs_pos->cim_up;

/*======================================================================*/
/* I) Set flags                                                         */

  *forceFormUp = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1; // temp keep this flag for debug filter
  if(cpLsda==1&&numStateDnProc!=0){
    *forceFormDn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }
  
  //genChemPotInterpPoints(stodftInfo,stodftCoefPos);
  //combineStoUCEnergyWindow(cp,general_data,class,cpMini,
  //             generalDataMini,classMini,ip_now);
  //energyCorrect(cpMini,generalDataMini,classMini,cp,class,ip_now);    
  //numChemPot -= 1;

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
  if(numProcStates>1)Barrier(commStates);
  stodftInfo->energyDiff = energyMax-energyMin;
  energyDiff = stodftInfo->energyDiff;
  stodftInfo->energyMean = 0.5*(energyMin+energyMax);

/*======================================================================*/
/* III) Calculate chemPot Values for Projectors                         */

  genChemPotInterpPoints(stodftInfo,stodftCoefPos);

/*======================================================================*/
/* IV) Generate Length of Polynomial Chain                              */

  timeStart2 = omp_get_wtime();
  if(expanType==1){
    Smin = chebyshevInfo->Smin;
    Smax = chebyshevInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    chebyshevInfo->scale = scale;
    //cheat my code
    if(myidState==0)genChebyHermit(stodftInfo,stodftCoefPos,2);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      polynormLength = stodftInfo->polynormLength;
      totalPoly = polynormLength*numChemPot;
      if(myidState!=0){
        stodftCoefPos->expanCoeff = (double*)crealloc(stodftCoefPos->expanCoeff,
                                totalPoly*sizeof(double));
        if(smearOpt>0){
            stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));
        }
      }
      Barrier(commStates);     
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      if(smearOpt>0){
        Bcast(stodftCoefPos->entropyExpanCoeff,polynormLength,MPI_DOUBLE,0,commStates);
      }
      Barrier(commStates);
    }
    //printf("numChemPotTemp %i\n",numChemPotTemp);
  }
  if(expanType==2){
    /*
    Smin = newtonInfo->Smin;
    Smax = newtonInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    newtonInfo->scale = scale;
    //cheat my code
    if(myidState==0)genNewtonHermit(stodftInfo,stodftCoefPos);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      Barrier(commStates);
    }
    */
    printf("We don't support Newton's Interpolation in EW+Frag scheme.\n");
    fflush(stdout);
    exit(0);
    //printf("numChemPotTemp %i\n",numChemPotTemp);
  }
  timeEnd2 = omp_get_wtime();
  diffTime2 = timeEnd2-timeStart2;

/*======================================================================*/
/* V) Project Noise Orbitals with Projectors                            */
/*    Reminder: Chebyshev Moments are also generated in this step.      */


  // 1. Generate noise orbitals
  genNoiseOrbitalReal(cp,cpcoeffs_pos);
  /*
  char fileName[100];
  sprintf(fileName,"noise-%i",myidState);
  FILE *fnoise = fopen(fileName,"r");
  for(iCoeff=1;iCoeff<=numCoeffUpTot;iCoeff++){
    fscanf(fnoise,"%lg",&coeffReUp[iCoeff]);
    fscanf(fnoise,"%lg",&coeffImUp[iCoeff]);
  }
  fclose(fnoise);
  */

  // 2. Filtering
 
  switch(expanType){
    case 1:
      stodftInfo->storeChebyMomentsFlag = 1;
      // At this step we need to turn off smearing 
      // since we have not initialize entropy coefficient
      smearOptBackup = stodftInfo->smearOpt;
      stodftInfo->smearOpt = 0;
      filterChebyPolyHermFake(cp,class,general_data,ip_now,1);
      stodftInfo->smearOpt = smearOptBackup;
      //filterChebyPolyHerm(cp,class,general_data,ip_now);
      stodftInfo->storeChebyMomentsFlag = 0;
      break;
    default:
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("Only Chebyshev polynormial can be used for ew+frag!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }



/*======================================================================*/
/* VI) Project Stochastic Orbitals with Fragment Density Matrix         */
 
  combineStoUCEnergyWindow(cp,general_data,class,cpMini,
               generalDataMini,classMini,ip_now);

/*======================================================================*/
/* VII) Generate Energy and Force Corrections                           */

  energyCorrect(cpMini,generalDataMini,classMini,cp,class,ip_now);    

/*======================================================================*/
/* IX) Generate Correct Chemical Potential                              */

  calcChemPotChebyEWFrag(cp,class,general_data,ip_now);

/*======================================================================*/
/* X) Redo the Polynormial Fitting with Correct Chemical Potential      */

  if(expanType==1){
    if(myidState==0)genChebyHermitTrueChemPot(stodftInfo,stodftCoefPos,3);
    if(numProcStates>1){
      Barrier(commStates);
      Bcast(&(stodftInfo->polynormLength),1,MPI_INT,0,commStates);
      polynormLength = stodftInfo->polynormLength;
      totalPoly = polynormLength*numChemPot;
      if(myidState!=0){
        stodftCoefPos->expanCoeff = (double*)crealloc(stodftCoefPos->expanCoeff,
                            totalPoly*sizeof(double));
        if(smearOpt>0){
          stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));
        }
      }
      Barrier(commStates);
      Bcast(stodftCoefPos->expanCoeff,totalPoly,MPI_DOUBLE,0,commStates);
      if(smearOpt>0){
        Bcast(stodftCoefPos->entropyExpanCoeff,polynormLength,MPI_DOUBLE,0,commStates);
      }
      Barrier(commStates);
    }
  }   
        
/*======================================================================*/
/* XI) Redo filtering with \sqrt{F*P_i}                                 */

  // 1. Generate noise orbitals
  genNoiseOrbitalReal(cp,cpcoeffs_pos);
  /*
  sprintf(fileName,"noise-%i",myidState);
  fnoise = fopen(fileName,"r");
  for(iCoeff=1;iCoeff<=numCoeffUpTot;iCoeff++){
    fscanf(fnoise,"%lg",&coeffReUp[iCoeff]);
    fscanf(fnoise,"%lg",&coeffImUp[iCoeff]);
  }
  fclose(fnoise);
  */

  
  // 2. Filtering
  switch(expanType){
    case 1:
      filterChebyPolyHermFake(cp,class,general_data,ip_now,0);
      break;
  }

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



#endif

