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
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define TIME_CP_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void filterNewtonPolyHerm(CP *cp,int ip_now,EWALD *ewald,EWD_SCR *ewd_scr,
			    CELL *cell,CLATOMS_INFO *clatoms_info,
			    CLATOMS_POS *clatoms_pos,ATOMMAPS *atommaps,
			    STAT_AVG *stat_avg,PTENS *ptens,SIMOPTS *simopts,
			    FOR_SCR *for_scr)
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
  CPEWALD *cpewald		= &(cp->cpewald);
  CPSCR *cpscr			= &(cp->cpscr);
  CPOPTS *cpopts		= &(cp->cpopts);
  PSEUDO *pseudo		= &(cp->pseudo);
  COMMUNICATE *communicate	= &(cp->communicate);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg); 

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
  double *sampPointRe = newtonInfo->sampPointRe;
  double *sampPointIm = newtonInfo->sampPointIm;

  double *expanCoeff = stodftCoefPos->expanCoeffRe;
  double *wfInReUp   = stodftCoefPos->wfInReUp;
  double *wfInImUp   = stodftCoefPos->wfInImUp;
  double *wfInReDn   = stodftCoefPos->wfInReDn;
  double *wfInImDn   = stodftCoefPos->wfInImDn;
  double *wfOutReUp  = stodftCoefPos->wfOutReUp;
  double *wfOutImUp  = stodftCoefPos->wfOutImUp;
  double *wfOutReDn  = stodftCoefPos->wfOutReDn;
  double *wfOutImDn  = stodftCoefPos->wfOutImDn;
  double **stoWfReUp = stodftCoefPos->stoWfReUp;
  double **stoWfImUp = stodftCoefPos->stoWfImUp;
  double **stoWfReDn = stodftCoefPos->stoWfReDn;
  double **stoWfImDn = stodftCoefPos->stoWfImDn;
  
/*==========================================================================*/
/* 0) Copy the initial stochastic orbital */
  
  for(imu=0;imu<numChemPot;imu++){
    for(iCoeff=0;iCoeff<numCoeffUpTotal;iCoeff++){
      stoWfReUp[imu][iCoeff] = expanCoeff[imu]*wfInReUp[iCoeff];
      stoWfImUp[imu][iCoeff] = expanCoeff[imu]*wfInImUp[iCoeff];
    }//endfor iCoeff
  }//endfor imu
  if(cpLsda==1&&numStateDnProc!=0){
    for(imu=0;imu<numChemPot;imu++){
      for(iCoeff=0;iCoeff<numCoeffDnTotal;iCoeff++){
	stoWfReDn[imu][iCoeff] = expanCoeff[imu]*wfInReDn[iCoeff];
        stoWfImDn[imu][iCoeff] = expanCoeff[imu]*wfInImDn[iCoeff];
      }//endfor iCoeff
    }//endfor imu
  }//endif

/*==========================================================================*/
/* 1) Loop over all polynomial terms (iPoly=0<=>polynomial order 1) */
  
  for(iPoly=1;iPoly<polynormLength;iPoly++){
    normHNewton(cp,cpcoeffs_pos,cpcoeffs_info,cell,clatoms_info,clatoms_pos,
	    ewald,ewd_scr,atommaps,for_scr,stat_avg,ptens,sampPointRe[iPoly-1]);
    for(imu=0;imu<numChemPot;imu++){
      polyCoeff = expanCoeff[iPoly*numChemPot+imu];
      for(iCoeff=0;iCoeff<numCoeffUpTotal;iCoeff++){
	stoWfReUp[imu][iCoeff] += polyCoeff*wfOutReUp[iCoeff];
        stoWfImUp[imu][iCoeff] += polyCoeff*wfOutImUp[iCoeff];
      }//endfor iCoeff
    }//endfor imu
  }//endfor iPoly
  if(cpLsda==1&&numStateDnProc!=0){
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      normHNewton(cp,cpcoeffs_pos,cpcoeffs_info,cell,clatoms_info,clatoms_pos,
		      ewald,ewd_scr,atommaps,for_scr,stat_avg,ptens,);
      for(imu=0;imu<numChemPot;imu++){
	polyCoeff = expanCoeff[iPoly*numChemPot+imu];
	for(iCoeff=0;iCoeff<numCoeffUpTotal;iCoeff++){
	  stoWfReDn[imu][iCoeff] += polyCoeff*wfOutReDn[iCoeff];
	  stoWfImDn[imu][iCoeff] += polyCoeff*wfOutImDn[iCoeff];
	}//endfor iCoeff
      }//endfor imu
    }//endfor iPoly
  }//endif

  return;
  
/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void filterNewtonPolyNoHerm(CP *cp,int ip_now,EWALD *ewald,EWD_SCR *ewd_scr,
			    CELL *cell,CLATOMS_INFO *clatoms_info,
			    CLATOMS_POS *clatoms_pos,ATOMMAPS *atommaps,
			    STAT_AVG *stat_avg,PTENS *ptens,SIMOPTS *simopts,
			    FOR_SCR *for_scr)
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
  CPEWALD *cpewald	= &(cp->cpewald);
  CPSCR *cpscr		= &(cp->cpscr);
  CPOPTS *cpopts	= &(cp->cpopts);
  PSEUDO *pseudo	= &(cp->pseudo);
  COMMUNICATE *communicate  = &(cp->communicate);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg); 

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
  double *sampPointRe = newtonInfo->sampPointRe;
  double *sampPointIm = newtonInfo->sampPointIm;

  double *expanCoeffRe = stodftCoefPos->expanCoeffRe;
  double *expanCoeffIm = stodftCoefPos->expanCoeffIm;
  double *wfInReUp   = stodftCoefPos->wfInReUp;
  double *wfInImUp   = stodftCoefPos->wfInImUp;
  double *wfInReDn   = stodftCoefPos->wfInReDn;
  double *wfInImDn   = stodftCoefPos->wfInImDn;
  double *wfOutReUp  = stodftCoefPos->wfOutReUp;
  double *wfOutImUp  = stodftCoefPos->wfOutImUp;
  double *wfOutReDn  = stodftCoefPos->wfOutReDn;
  double *wfOutImDn  = stodftCoefPos->wfOutImDn;
  double **stoWfReUp = stodftCoefPos->stoWfReUp;
  double **stoWfImUp = stodftCoefPos->stoWfImUp;
  double **stoWfReDn = stodftCoefPos->stoWfReDn;
  double **stoWfImDn = stodftCoefPos->stoWfImDn;
  
/*==========================================================================*/
/* 0) Copy the initial stochastic orbital */
  
  for(imu=0;imu<numChemPot;imu++){
    for(iCoeff=0;iCoeff<numCoeffUpTotal;iCoeff++){
      stoWfReUp[imu][iCoeff] = expanCoeffRe[imu]*wfInReUp[iCoeff]-expanCoeffIm[imu]*wfInImUp[iCoeff];
      stoWfImUp[imu][iCoeff] = expanCoeffRe[imu]*wfInImUp[iCoeff]+expanCoeffIm[imu]*wfInReUp[iCoeff];
    }//endfor iCoeff
  }//endfor imu
  if(cpLsda==1&&numStateDnProc!=0){
    for(imu=0;imu<numChemPot;imu++){
      for(iCoeff=0;iCoeff<numCoeffDnTotal;iCoeff++){
    stoWfReDn[imu][iCoeff] = expanCoeff[imu]*wfInReDn[iCoeff];
        stoWfImDn[imu][iCoeff] = expanCoeff[imu]*wfInImDn[iCoeff];
      }//endfor iCoeff
    }//endfor imu
  }//endif

/*==========================================================================*/
/* 1) Loop over all polynomial terms (iPoly=0<=>polynomial order 1) */
  
  for(iPoly=1;iPoly<polynormLength;iPoly++){
    switch(expanType){
      case 1:
    //normHChebyshev(); 
    break;
      case 2:
    normHNewton(cp,cpcoeffs_pos,cpcoeffs_info,cell,clatoms_info,clatoms_pos,
	    ewald,ewd_scr,atommaps,for_scr,stat_avg,ptens,sampPointRe[iPoly-1]);
    break;
    }//endswitch expanType    
    for(imu=0;imu<numChemPot;imu++){
      polyCoeff = expanCoeff[iPoly*numChemPot+imu];
      for(iCoeff=0;iCoeff<numCoeffUpTotal;iCoeff++){
    stoWfReUp[imu][iCoeff] += polyCoeff*wfOutReUp[iCoeff];
        stoWfImUp[imu][iCoeff] += polyCoeff*wfOutImUp[iCoeff];
      }//endfor iCoeff
    }//endfor imu
  }//endfor iPoly
  if(cpLsda==1&&numStateDnProc!=0){
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      switch(expanType){
    case 1:
      //normHChebyshev();
      break;
    case 2:
      normHNewton(cp,cpcoeffs_pos,cpcoeffs_info,cell,clatoms_info,clatoms_pos,
	      ewald,ewd_scr,atommaps,for_scr,stat_avg,ptens,);
      break;
      }//endswitch expanType    
      for(imu=0;imu<numChemPot;imu++){
    polyCoeff = expanCoeff[iPoly*numChemPot+imu];
    for(iCoeff=0;iCoeff<numCoeffUpTotal;iCoeff++){
      stoWfReDn[imu][iCoeff] += polyCoeff*wfOutReDn[iCoeff];
      stoWfImDn[imu][iCoeff] += polyCoeff*wfOutImDn[iCoeff];
    }//endfor iCoeff
      }//endfor imu
    }//endfor iPoly
  }//endif

  return;
  
/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/



