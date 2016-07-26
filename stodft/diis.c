/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: diis.c	                                    */
/*                                                                          */
/* This routine calculate diis density update.                              */
/*	                                                                    */
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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genDensityDiis(CP *cp,int iScf)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the routine to calculate the diis density                     */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPSCR *cpscr		        = &(cp->cpscr);
  COMMUNICATE *communicate      = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  
  int iGrid;
  int rhoRealGridNum = stodftInfo->rhoRealGridNum;
  int numDiis = stodftInfo->numDiis;
  int numStepMix = stodftInfo->numStepMix;
  int cpLsda = cpopts->cp_lsda;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;

  double mixRatio1 = stodftInfo->mixRatio;
  double mixRatio2 = 1.0-mixRatio1;

  double *rhoUpCorrect = stodftCoefPos->rhoUpCorrect;
  double *rhoDnCorrect = stodftCoefPos->rhoDnCorrect;
  double *rhoUp = cpscr->cpscr_rho.rho_up;
  double *rhoDn = cpscr->cpscr_rho.rho_dn;

  double **rhoUpBank = stodftCoefPos->rhoUpBank;
  double **rhoDnBank = stodftCoefPos->rhoDnBank;
  double **rhoUpErr  = stodftCoefPos->rhoUpErr;
  double **rhoDnErr  = stodftCoefPos->rhoDnErr;

  updateBank(rhoUpCorrect,rhoUpBank,iScf);
  updateErr(rhoUpBank,rhoUpErr,iScf);

  if(iScf==0){//Initial Step
  }
  else if(iScf<=numStepMix){//Mixing Step
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
      rhoUp[iGrid+1] = rhoUpBank[0][iGrid]*mixRatio1+rhoUpBank[1][iGrid]*mixRatio2;
    }
  }
  else{//diis
    
  }
  if(cpLsda==1&&numStateDnProc>0){
    updateBank(rhoDnCorrect,rhoDnBank,iScf);
    updateErr(rhoDnBank,rhoDnErr,iScf);

    if(iScf==0){//Initial Step
    }
    else if(iScf<=numStepMix){//Mixing Step
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
	rhoDn[iGrid+1] = rhoDnBank[0][iGrid]*mixRatio1+rhoDnBank[1][iGrid]*mixRatio2;
      }
    }
    else{//diis
    }

  }


/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

