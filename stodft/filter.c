/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: filters.c                                      */
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
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_stodft_local.h"

#define TIME_CP_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void filterNewtonPolyHerm(CP *cp,CLASS *class,GENERAL_DATA *general_data,
			  int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CLATOMS_POS*  clatoms_pos     = &(class->clatoms_pos[ip_now]);

  NEWTONINFO *newtonInfo = stodftInfo->newtonInfo;

  int expanType	     = stodftInfo->expanType;
  int polynormLength = stodftInfo->polynormLength;
  int numChemPot     = stodftInfo->numChemPot;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int cpLsda	     = cpopts->cp_lsda;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int imu,iCoeff,iPoly,indexStart;
  int startIndex;

  double energyDiff  = stodftInfo->energyDiff;
  double polyCoeff;
  double *sampPoint = (double*)newtonInfo->sampPoint;

  double *expanCoeff = (double*)stodftCoefPos->expanCoeff;
  double complex *wfInUp   = stodftCoefPos->wfInUp;
  double complex *wfInDn   = stodftCoefPos->wfInDn;
  double complex *wfOutUp  = stodftCoefPos->wfOutUp;
  double complex *wfOutDn  = stodftCoefPos->wfOutDn;
  double complex **stoWfUp = stodftCoefPos->stoWfUp;
  double complex **stoWfDn = stodftCoefPos->stoWfDn;
  
/*==========================================================================*/
/* 0) Copy the initial stochastic orbital */
  
  for(imu=0;imu<numChemPot;imu++){
    for(iCoeff=0;iCoeff<numCoeffUpTotal;iCoeff++){
      stoWfUp[imu][iCoeff] = expanCoeff[imu]*wfInUp[iCoeff];
    }//endfor iCoeff
    if(cpLsda==1&&numStateDnProc!=0){
      for(iCoeff=0;iCoeff<numCoeffDnTotal;iCoeff++){
        stoWfDn[imu][iCoeff] = expanCoeff[imu]*wfInDn[iCoeff];
      }//endfor iCoeff      
    }//endif 
  }//endfor imu

/*==========================================================================*/
/* 1) Loop over all polynomial terms (iPoly=0<=>polynomial order 1) */
  
  for(iPoly=1;iPoly<polynormLength;iPoly++){
    normHNewtonNoHerm(cp,class,general_data,
                 cpcoeffs_pos,clatoms_pos,sampPoint[iPoly-1]);
    for(imu=0;imu<numChemPot;imu++){
      polyCoeff = expanCoeff[iPoly*numChemPot+imu];
      for(iCoeff=0;iCoeff<numCoeffUpTotal;iCoeff++){
	stoWfUp[imu][iCoeff] += polyCoeff*wfOutUp[iCoeff];	
      }//endfor iCoeff
      if(cpLsda==1&&numStateDnProc!=0){
        for(iCoeff=0;iCoeff<numCoeffUpTotal;iCoeff++){
          stoWfDn[imu][iCoeff] += polyCoeff*wfOutDn[iCoeff];
        }//endfor iCoeff        
      }//endif 
    }//endfor imu
  }//endfor iPoly

/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void filterNewtonPolyNoHerm(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                          int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CLATOMS_POS*  clatoms_pos     = &(class->clatoms_pos[ip_now]);

  NEWTONINFO *newtonInfo = stodftInfo->newtonInfo;

  int expanType	     = stodftInfo->expanType;
  int polynormLength = stodftInfo->polynormLength;
  int numChemPot     = stodftInfo->numChemPot;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int cpLsda	     = cpopts->cp_lsda;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int imu,iCoeff,iPoly,indexStart;
  int startIndex;

  double energyDiff  = stodftInfo->energyDiff;
  double polyCoeff;
  double complex *sampPoint = (double complex*)newtonInfo->sampPoint;

  double complex *expanCoeff = (double complex*)stodftCoefPos->expanCoeff;
  double complex *wfInUp   = stodftCoefPos->wfInUp;
  double complex *wfInDn   = stodftCoefPos->wfInDn;
  double complex *wfOutUp  = stodftCoefPos->wfOutUp;
  double complex *wfOutDn  = stodftCoefPos->wfOutDn;
  double complex **stoWfUp = stodftCoefPos->stoWfUp;
  double complex **stoWfDn = stodftCoefPos->stoWfDn;
  
/*==========================================================================*/
/* 0) Copy the initial stochastic orbital */
  
  for(imu=0;imu<numChemPot;imu++){
    for(iCoeff=0;iCoeff<numCoeffUpTotal;iCoeff++){
      stoWfUp[imu][iCoeff] = expanCoeff[imu]*wfInUp[iCoeff];
    }//endfor iCoeff
    if(cpLsda==1&&numStateDnProc!=0){
      for(iCoeff=0;iCoeff<numCoeffDnTotal;iCoeff++){
        stoWfDn[imu][iCoeff] = expanCoeff[imu]*wfInDn[iCoeff];
      }//endfor iCoeff
    }
  }//endfor imu

/*==========================================================================*/
/* 1) Loop over all polynomial terms (iPoly=0<=>polynomial order 1) */
  
  for(iPoly=1;iPoly<polynormLength;iPoly++){
    normHNewtonNoHerm(cp,class,general_data,
                 cpcoeffs_pos,clatoms_pos,sampPoint[iPoly-1]);
    for(imu=0;imu<numChemPot;imu++){
      polyCoeff = expanCoeff[iPoly*numChemPot+imu];
      for(iCoeff=0;iCoeff<numCoeffUpTotal;iCoeff++){
        stoWfUp[imu][iCoeff] += polyCoeff*wfOutUp[iCoeff];
      }//endfor iCoeff
      if(cpLsda==1&&numStateDnProc!=0){
        for(iCoeff=0;iCoeff<numCoeffUpTotal;iCoeff++){
          stoWfDn[imu][iCoeff] += polyCoeff*wfOutDn[iCoeff];
        }//endfor iCoeff
      }//endif
    }//endfor imu
  }//endfor iPoly

/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/



