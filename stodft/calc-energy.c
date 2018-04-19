/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: calc-energy.c                                  */
/*                                                                          */
/* This routine constructs filters of Fermi function.                       */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_stodft_local.h"

#define TIME_CP_OFF
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcEnergyChemPot(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                       CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate average energy for each chemical potential, then            */
/* Interpolate the energy for the correct chemical potential             */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  COMMUNICATE *communicate      = &(cp->communicate);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  CPEWALD *cpewald              = &(cp->cpewald);
  CELL *cell			= &(general_data->cell);
  CLATOMS_INFO *clatoms_info	= &(class->clatoms_info);
  PSEUDO *pseudo		= &(cp->pseudo);
  PSEUDO_REAL *pseudoReal	= &(pseudo->pseudoReal);

  int cpLsda         = cpopts->cp_lsda;
  int numStateStoUp  = stodftInfo->numStateStoUp;
  int atomForceFlag  = stodftInfo->atomForceFlag;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numChemPot = stodftInfo->numChemPot;
  int occNumber = stodftInfo->occNumber;
  int myidState         = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int numAtomTot = clatoms_info->natm_tot;
  int pseudoRealFlag = pseudoReal->pseudoRealFlag;
  int iState,iCoeff,iChem,iAtom;
  int ioff,iis;

  double tpi = 2.0*M_PI;
  double eke,ekeDn;
  double chemPotTrue = stodftInfo->chemPotTrue;
  double energyKineticTemp,energyNLTemp;

  double *energyKe  = stodftInfo->energyKe;
  double *energyPNL = stodftInfo->energyPNL;
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  double *ak2_sm  =  cpewald->ak2_sm;
  double *chemPot = stodftCoefPos->chemPot;
  double *fx = clatoms_pos->fx;
  double *fy = clatoms_pos->fy;
  double *fz = clatoms_pos->fz;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;
  double **fxNl	     = stodftCoefPos->fxNl;
  double **fyNl	     = stodftCoefPos->fxNl;
  double **fzNl	     = stodftCoefPos->fxNl;

  MPI_Comm commStates = communicate->comm_states;


/*--------------------------------------------------------------------------*/
/* I) Generate kinetic energy and nonlocal pseudopotential energy for       */
/*    each chemical potential.                                              */

  get_ak2_sm(cpewald,cell);

  for(iChem=0;iChem<numChemPot;iChem++){
    stat_avg->vrecip = 0.0;
    stat_avg->cp_enl = 0.0;
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      cre_up[iCoeff] = stoWfUpRe[iChem][iCoeff];
      cim_up[iCoeff] = stoWfUpIm[iChem][iCoeff];
      fcre_up[iCoeff] = 0.0;
      fcim_up[iCoeff] = 0.0;
    }//endfor iCoeff
    if(cpLsda==1&&numStateDnProc!=0){
      for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
        cre_dn[iCoeff] = stoWfDnRe[iChem][iCoeff];
        cim_dn[iCoeff] = stoWfDnIm[iChem][iCoeff];
        fcre_dn[iCoeff] = 0.0;
        fcim_dn[iCoeff] = 0.0;
      }//endfor iCoeff
    }//endif cpLsda
    
    
    //printf("stowf %lg %lg\n",stoWfUpRe[iChem][219],stoWfUpIm[iChem][219]);
    eke = 0.0;
    for(iState=0;iState<numStateUpProc;iState++){
      ioff = iState*numCoeff;
      for(iCoeff=1;iCoeff<=numCoeff-1;iCoeff++){
        iis = ioff+iCoeff;
        eke += 2.0*ak2_sm[iCoeff]*(cre_up[iis]*cre_up[iis]+cim_up[iis]*cim_up[iis]);
      }//endfor iCoeff
    }//endfor iState
    eke *= occNumber*0.5;
    if(cpLsda==1&&numStateDnProc!=0){
      ekeDn = 0.0;
      for(iState=0;iState<numStateDnProc;iState++){
        ioff = iState*numCoeff;
        for(iCoeff=1;iCoeff<=numCoeff-1;iCoeff++){
          iis = ioff+iCoeff;
          ekeDn += 2.0*ak2_sm[iCoeff]*(cre_dn[iis]*cre_dn[iis]+cim_dn[iis]*cim_dn[iis]);
        }//endfor iCoeff
      }//endfor iState
      eke += ekeDn*occNumber*0.5;
    }//endif cpLsda
    //stat_avg->cp_eke = eke;
    
    // test ke
    // read in stowf-frag
    
    //pp 
    //calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    if(pseudoRealFlag==0){
      calcNonLocalPseudoScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    }
    else{
      for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
	fcre_up[iCoeff] = 0.0;
	fcim_up[iCoeff] = 0.0;
      }//endfor iCoeff
      if(cpLsda==1&&numStateDnProc!=0){
	for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
	  fcre_dn[iCoeff] = 0.0;
	  fcim_dn[iCoeff] = 0.0;
	}//endfor iCoeff
      }//endif cpLsda     
      calcCoefForceScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    }
    //calcCoefForceExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    stat_avg->cp_eke = eke;
    stat_avg->cp_enl *= occNumber;

    if(numProcStates>1){
      Reduce(&(stat_avg->cp_enl),&energyNLTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);
      Reduce(&(stat_avg->cp_eke),&energyKineticTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);
    }
    else{
      energyNLTemp = stat_avg->cp_enl;
      energyKineticTemp = stat_avg->cp_eke;
    }

    if(myidState==0){
      energyKe[iChem] = energyKineticTemp/numStateStoUp;
      energyPNL[iChem] = energyNLTemp/numStateStoUp;
      //printf("iChem %i chemPot %lg K %lg NL %lg\n",iChem,chemPot[iChem],energyKineticTemp,energyNLTemp);
    }
  }//endfor iChem
  /*
  // debug
  printf("Hartree Energy: %.6lg\n",stat_avg->cp_ehart);
  printf("External Potential Energy(Local Pseudopotential): %.6lg\n",stat_avg->cp_eext);
  printf("Exchange-Correlation Energy: %.6lg\n",stat_avg->cp_exc);
  */

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcTotEnergy(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                  CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate average energy for each chemical potential, then            */
/* Interpolate the energy for the correct chemical potential	     */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  COMMUNICATE *communicate      = &(cp->communicate);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  FRAGINFO *fragInfo		= stodftInfo->fragInfo;
  CLATOMS_INFO *clatoms_info	= &(class->clatoms_info);

  int numChemPot = stodftInfo->numChemPot;
  int chemPotOpt = stodftInfo->chemPotOpt;
  int atomForceFlag  = stodftInfo->atomForceFlag;
  int myidState         = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int calcFragFlag = stodftInfo->calcFragFlag;
  int numAtomTot = clatoms_info->natm_tot;
  int iState,iCoeff,iChem,iAtom;

  double chemPotTrue = stodftInfo->chemPotTrue;
  double energyKineticTemp,energyNLTemp;
  double energyHartTemp,energyExtTemp,energyExcTemp;
  double energyKeTrue,energyPNLTrue,energyTotElec;
  double energyKeNoCor,energyPNLNoCor;

  double *chemPot = stodftCoefPos->chemPot;
  double *energyKe = stodftInfo->energyKe;
  double *energyPNL = stodftInfo->energyPNL;
  double *lagFunValue = (double*)cmalloc(numChemPot*sizeof(double));
  double *fxTemp,*fyTemp,*fzTemp;
  double **fxNl      = stodftCoefPos->fxNl;
  double **fyNl      = stodftCoefPos->fxNl;
  double **fzNl      = stodftCoefPos->fxNl;

  MPI_Comm commStates = communicate->comm_states; 

/*--------------------------------------------------------------------------*/
/* II) Interpolate the correct energy		                */
    
  if(myidState==0){
    if(chemPotOpt==1){
      //debug
      /*
      printf("chemPotTrue %lg\n",chemPotTrue);
      for(iChem=0;iChem<numChemPot;iChem++){
	printf("%lg %lg\n",chemPot[iChem],energyKNL[iChem]);
      }
      */
      energyKeTrue = calcLagrangeInterpFun(numChemPot,chemPotTrue,chemPot,energyKe,lagFunValue);
      energyPNLTrue = calcLagrangeInterpFun(numChemPot,chemPotTrue,chemPot,energyPNL,lagFunValue);
    }
    if(chemPotOpt==2){
      energyKeTrue = energyKe[0];
      energyPNLTrue = energyPNL[0];
    }
  }

/*--------------------------------------------------------------------------*/
/* III) Add fragmentation correction		                            */

  if(calcFragFlag==1&&myidState==0){
    energyKeNoCor = energyKeTrue;
    energyKeTrue += fragInfo->keCor;
    energyPNLNoCor = energyPNLTrue;
    energyPNLTrue += fragInfo->vnlCor;
  }
  if(myidState==0){ // do this for either using frag or not
    general_data->stat_avg.cp_enl = energyPNLTrue;
    general_data->stat_avg.cp_eke = energyKeTrue;
  }

/*--------------------------------------------------------------------------*/
/* IV) Reduce all the other energy terms calculated from density            */
 
  if(numProcStates>1){
    Reduce(&(stat_avg->cp_ehart),&energyHartTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);   
    Reduce(&(stat_avg->cp_eext),&energyExtTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&(stat_avg->cp_exc),&energyExcTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);
    // Let's copy it back to stat_avg since I need to use these energies after SCF
    stat_avg->cp_ehart = energyHartTemp;
    stat_avg->cp_eext  = energyExtTemp;
    stat_avg->cp_exc   = energyExcTemp;
  }
  else{
    energyHartTemp = stat_avg->cp_ehart;
    energyExtTemp = stat_avg->cp_eext;
    energyExcTemp = stat_avg->cp_exc;
  }

/*--------------------------------------------------------------------------*/
/* V) Output the energy Term						    */

  if(myidState==0){
    energyTotElec = energyKeTrue+energyPNLTrue+energyHartTemp+energyExtTemp+energyExcTemp;
    stodftInfo->energyElecTot = energyTotElec;
    printf("==============================================\n");
    printf("Output Energy\n");
    printf("==============================================\n");
    printf("Kinetic Energy:	 %.16lg %.16lg\n",energyKeNoCor,energyKeTrue);
    printf("NL Pseudopotential:  %.16lg %.16lg %.16lg\n",energyPNLNoCor,energyPNLTrue,fragInfo->vnlCor);
    printf("Hartree Energy:      %.16lg\n",energyHartTemp);
    printf("Ext Energy:          %.16lg\n",energyExtTemp); 
    printf("Ex-Cor Energy:       %.16lg\n",energyExcTemp); 
    printf("Total Elec Energy    %.16lg\n",energyTotElec);
    printf("Elec Energy Diff	 %.16lg\n",
	    stodftInfo->energyElecTot-stodftInfo->energyElecTotOld);
    printf("==============================================\n");
  }
  if(numProcStates>1) Bcast(&stodftInfo->energyElecTot,1,MPI_DOUBLE,0,commStates);
  free(lagFunValue);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcKNEEnergyFilterDiag(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                       CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate average energy for each chemical potential, then            */
/* Interpolate the energy for the correct chemical potential             */
/*************************************************************************/
/*=======================================================================*/
/*=======================================================================*/
/*         Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  COMMUNICATE *communicate      = &(cp->communicate);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  CPEWALD *cpewald              = &(cp->cpewald);

  int cpLsda         = cpopts->cp_lsda;
  int numStateStoUp  = stodftInfo->numStateStoUp;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numChemPot = stodftInfo->numChemPot;
  int occNumber = stodftInfo->occNumber;
  int myidState         = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int iState,iCoeff,iChem;
  int ioff,iis;

  double tpi = 2.0*M_PI;
  double eke,ekeDn;
  double chemPotTrue = stodftInfo->chemPotTrue;
  double energyKineticTemp,energyNLTemp;

  double *energyKe  = stodftInfo->energyKe;
  double *energyPNL = stodftInfo->energyPNL;
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  double *ak2_sm  =  cpewald->ak2_sm;
  double *chemPot = stodftCoefPos->chemPot;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  MPI_Comm commStates = communicate->comm_states;

/*--------------------------------------------------------------------------*/
/* I) Generate kinetic energy and nonlocal pseudopotential energy for       */
/*    each chemical potential.                                              */
  if(myidState==0){
    energyKe[0] = 0.0;
    energyPNL[0] = 0.0;
  }
  Barrier(commStates);
  
  for(iChem=0;iChem<numChemPot;iChem++){
    stat_avg->vrecip = 0.0;
    stat_avg->cp_enl = 0.0;
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      cre_up[iCoeff] = stoWfUpRe[iChem][iCoeff];
      cim_up[iCoeff] = stoWfUpIm[iChem][iCoeff];
      fcre_up[iCoeff] = 0.0;
      fcim_up[iCoeff] = 0.0;
    }//endfor iCoeff
    if(cpLsda==1&&numStateDnProc!=0){
      for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
        cre_dn[iCoeff] = stoWfDnRe[iChem][iCoeff];
        cim_dn[iCoeff] = stoWfDnIm[iChem][iCoeff];
        fcre_dn[iCoeff] = 0.0;
        fcim_dn[iCoeff] = 0.0;
      }//endfor iCoeff
    }//endif cpLsda

    
    eke = 0.0;
    for(iState=0;iState<numStateUpProc;iState++){
      ioff = iState*numCoeff;
      for(iCoeff=1;iCoeff<=numCoeff-1;iCoeff++){
        iis = ioff+iCoeff;
        eke += 2.0*ak2_sm[iCoeff]*(cre_up[iis]*cre_up[iis]+cim_up[iis]*cim_up[iis]);
      }//endfor iCoeff
    }//endfor iState
    eke *= 0.5;
    if(cpLsda==1&&numStateDnProc!=0){
      ekeDn = 0.0;
      for(iState=0;iState<numStateDnProc;iState++){
        ioff = iState*numCoeff;
        for(iCoeff=1;iCoeff<=numCoeff-1;iCoeff++){
          iis = ioff+iCoeff;
          ekeDn += 2.0*ak2_sm[iCoeff]*(cre_dn[iis]*cre_dn[iis]+cim_dn[iis]*cim_dn[iis]);
        }//endfor iCoeff
      }//endfor iState
      eke += ekeDn*0.5;
    }//endif cpLsda
    printf("eke %lg\n",eke);
    stat_avg->cp_eke = eke;
    
    //printf("before cp_eke %lg\n",stat_avg->cp_eke);
    
    calcCoefForceWrapSCF(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    //calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    //calcCoefForceExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    
    if(numProcStates>1){
      Reduce(&(stat_avg->cp_enl),&energyNLTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);
      Reduce(&eke,&energyKineticTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);
    }
    else{
      energyNLTemp = stat_avg->cp_enl;
      energyKineticTemp = eke;
    }

    if(myidState==0){
      //energyNLTemp /= numStateStoUp;
      //energyKineticTemp /= numStateStoUp;
      energyKe[0] += energyKineticTemp;
      energyPNL[0] += energyNLTemp;
      //printf("iChem %i chemPot %lg K %lg NL %lg\n",iChem,chemPot[iChem],energyKineticTemp,energyNLTemp);
    }
  }//endfor iChem
  /*
  // debug
  printf("Hartree Energy: %.6lg\n",stat_avg->cp_ehart);
  printf("External Potential Energy(Local Pseudopotential): %.6lg\n",stat_avg->cp_eext);
  printf("Exchange-Correlation Energy: %.6lg\n",stat_avg->cp_exc);
  */


/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcTotEnergyFilterDiag(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                  CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate average energy for each chemical potential, then            */
/* Interpolate the energy for the correct chemical potential             */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  COMMUNICATE *communicate      = &(cp->communicate);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);

  int numChemPot = stodftInfo->numChemPot;
  int myidState         = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int iState,iCoeff,iChem;

  double chemPotTrue = stodftInfo->chemPotTrue;
  double energyKineticTemp,energyNLTemp;
  double energyHartTemp,energyExtTemp,energyExcTemp;
  double energyTrue,energyTotElec;

  double *chemPot = stodftCoefPos->chemPot;
  double *energyKe = stodftInfo->energyKe;
  double *energyPNL = stodftInfo->energyPNL;
  double *lagFunValue = (double*)cmalloc(numChemPot*sizeof(double));

  MPI_Comm commStates = communicate->comm_states;

/*--------------------------------------------------------------------------*/
/* III) Reduce all the other energy terms calculated from density           */

  if(numProcStates>1){
    Reduce(&(stat_avg->cp_ehart),&energyHartTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&(stat_avg->cp_eext),&energyExtTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&(stat_avg->cp_exc),&energyExcTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);
  }
  else{
    energyHartTemp = stat_avg->cp_ehart;
    energyExtTemp = stat_avg->cp_eext;
    energyExcTemp = stat_avg->cp_exc;
  }

/*--------------------------------------------------------------------------*/
/* IV) Output the energy Term                                               */

  if(myidState==0){
    energyTotElec = energyKe[0]+energyPNL[0]+energyHartTemp
	    +energyExtTemp+energyExcTemp;
    printf("==============================================\n");
    printf("Output Energy\n");
    printf("==============================================\n");
    printf("Kinetic Energy:      %.20lg\n",energyKe[0]);
    printf("NLPP:	         %.20lg\n",energyPNL[0]);
    printf("Hartree Energy:      %.20lg\n",energyHartTemp);
    printf("Ext Energy:          %.20lg\n",energyExtTemp);
    printf("Ex-Cor Energy:       %.20lg\n",energyExcTemp);
    printf("Total Elec Energy    %.20lg\n",energyTotElec);
    printf("==============================================\n");
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/



