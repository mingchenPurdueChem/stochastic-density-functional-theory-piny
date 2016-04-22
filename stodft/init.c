/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: init.c                                         */
/*                                                                          */
/* This routine initialize the stochastic dft calculation.                  */
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

#define TIME_CP_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initStodft(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate expanCoeff for different chemical potentials                */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  CPOPTS *cpopts = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cp->cpcoeffs_info);
  STODFTINFO *stodftInfo       = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos = cp->stodftCoefPos;
  NEWTONINFO *newtonInfo;

  int cpLsda         = cpopts->cp_lsda;
  int expanType      = stodftInfo->expanType;
  int numOrbital     = stodftInfo->numOrbital;
  int polynormLength = stodftInfo->polynormLength;
  int numChemPot     = stodftInfo->numChemPot;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numStateUpTot  = numStateUpProc*numCoeff;
  int numStateDnTot  = numStateDnProc*numCoeff;
  int totalPoly	     = polynormLength*numChemPot;
  int fermiFunType   = stodftInfo->fermiFunType;
  
  int iChem,iSamp;
  
  double energyMax = stodftInfo->energyMax;
  double energyMin = stodftInfo->energyMin;
  double Smin = -2.0;
  double Smax = 2.0;
  double energyDiff;

/*==========================================================================*/
/* I) General parameters and malloc					    */
  
  stodftInfo->energyDiff = energyMax-energyMin;
  energyDiff = stodftInfo->energyDiff;
  stodftInfo->energyMean = (energyMax+energyMin)*0.5;
  stodftInfo->vpsAtomListFlag = 0;

  stodftCoefPos->wfInUp = (double complex*)cmalloc(numStateUpTot*sizeof(double complex));
  stodftCoefPos->wfOutUp = (double complex*)cmalloc(numStateUpTot*sizeof(double complex));
  stodftCoefPos->stoWfUp = (double complex**)cmalloc(numChempot*sizeof(double complex*));
  for(iChem=0;iChem<numChemPot;iChem++){
    stodftCoefPos->stoWfUp[iChem] = (double complex*)cmalloc(numStateUpTot*sizeof(double complex))
  }
  if(cpLsda==1&&numStateDnProc!=0){
    stodftCoefPos->wfInDn = (double complex*)cmalloc(numStateDnTot*sizeof(double complex));
    stodftCoefPos->wfOutDn = (double complex*)cmalloc(numStateDnTot*sizeof(double complex));
    stodftCoefPos->stoWfDn = (double complex**)cmalloc(numChempot*sizeof(double complex*));
    for(iChem=0;iChem<numChemPot;iChem++){
      stodftCoefPos->stoWfDn[iChem] = (double complex*)cmalloc(numStateDnTot*sizeof(double complex))
    }//endfor iChem
  }//endif

  if(expanType==2&&fermiFunType==1)stodftInfo->fermiFunction = &fermiExpReal;
  if(expanType==2&&fermiFunType==2)stodftInfo->fermiFunction = &fermiErfcReal;
  if(expanType==3&&fermiFunType==1)stodftInfo->fermiFunction = &fermiExpComplex;
  if(expanType==3&&fermiFunType==1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("We haven't implement erfc type of Fermi function \n");
    printf("for non-Hermitian KS Hamiltonian. Please use the \n");
    printf("exponential type in this case!\n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(0);
  }


/*==========================================================================*/
/* II) Malloc by expension type						    */

  switch(expanType){
    case 2:
      stodftCoefPos->expanCoeff = (double *)cmalloc(totalPoly*sizeof(double));
      stodftInfo->newtonInfo = (NEWTONINFO *)cmalloc(sizeof(NEWTONINFO));
      newtonInfo = stodftInfo->newtonInfo;
      newtonInfo->samplePoint = (double *)cmalloc(polynormLength*sizeof(double));
      newtonInfo->sampPointUnscale = (double *)cmalloc(polynormLength*sizeof(double));
      newtonInfo->Smin = Smin;
      newtonInfo->Smax = Smax;
      newtonInfo->scale = (Smax-Smin)/energyDiff;      
      
      break;
    case 3:
      stodftCoefPos->expanCoeff = (double complex*)cmalloc(totalPoly*sizeof(double));
      stodftInfo->newtonInfo = (NEWTONINFO *)cmalloc(sizeof(NEWTONINFO));
      newtonInfo = stodftInfo->newtonInfo;
      newtonInfo->samplePoint = (double complex*)cmalloc(polynormLength*sizeof(double));
      newtonInfo->sampPointUnscale = (double complex*)cmalloc(polynormLength*sizeof(double));
      newtonInfo->Smin = Smin;
      newtonInfo->Smax = Smax;
      newtonInfo->scale = (Smax-Smin)/energyDiff;
      break;
  }

/*==========================================================================*/
/* III) Initialize utility data						    */
  
  switch(expanType){
    case 2:
      genSampNewtonHermit(stodftInfo,stodftCoefPos);
      break;
    // I'll do chebyshev and non-Hermitain Newtonian later    
  }

/*==========================================================================*/
/* IV) Initialize Flags							    */

/*==========================================================================*/
/* V) Calculate the non-local pseudopotential list                          */

  if(stodftInfo->vpsAtomListFlag==0||cp_dual_grid_opt_on >= 1){
    control_vps_atm_list(pseudo,cell,clatoms_pos,clatoms_info,
                         atommaps,ewd_scr,for_scr,cp_dual_grid_opt_on,
                         stodftInfo->vpsAtomListFlag);
    stodftInfo->vpsAtomListFlag = 1;
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/







