/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: normh.c                                        */
/*                                                                          */
/* This routine costruct P_N(H)|phi> where P_N is some                      */
/* polynomial.                                                              */
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
#define TIME_CP_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void normHNewtonHerm(CP *cp,CLASS *class,GENERAL_DATA *general_data,
		 CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos,double zn)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Let's first build a simple version. We will first apply this on       */
/* Hermitian operator. So all the sample points will locate on the real  */
/* axis. zn is the interpolation point. For details please read		 */
/* Ashkenazi et.al. JChemPhys 103, 10005(1995).				 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  COMMUNICATE *communicate      = &(cp->communicate);

  NEWTONINFO *newtonInfo = stodftInfo->newtonInfo;
  double Smin		= newtonInfo->Smin;
  double Smax		= newtonInfo->Smax;
  double scale		= newtonInfo->scale;
  double energyMean	= stodftInfo->energyMean;
  double energyDiff	= stodftInfo->energyDiff;
  double prefact	= -scale*energyMean-zn;
  double scale1		= -scale*0.25;
  double scale2		= -scale*0.5;

  int cpLsda         = cpopts->cp_lsda;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int cpDualGridOptOn  = cpopts->cp_dual_grid_opt;
  int numCoeffM1     = numCoeff-1;
  int incx = 1;
  int incy = 1;
  int iState,iCoeff,iCoeffStart,index1,index2;
  
  double *expanCoeff = (double*)stodftCoefPos->expanCoeff;
  
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  double dot;

  
/*==========================================================================*/
/* 1) Calculate the H/sigma|phi> */
  //control_vps_atm_list will be done somewhere else (perhaps in density calculation?)

  //calcCoefForceWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  /*
  calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  calcCoefForceWrapReduce(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  */
  calcCoefForceWrapSCF(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

  for(iState=0;iState<numStateUpProc;iState++){
    iCoeffStart = iState*numCoeff;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      index1 = iCoeffStart+iCoeff;
      fcre_up[index1] *= scale1;
      fcim_up[index1] *= scale1;	
    }//endfor iCoeff
    index1 = iCoeffStart+numCoeff;
    fcre_up[index1] *= scale2;
  }//endfor iState
  if(cpLsda==1&&numStateDnProc!=0){
    for(iState=0;iState<numStateDnProc;iState++){
      iCoeffStart = iState*numCoeff;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	index1 = iCoeffStart+iCoeff;
	fcre_dn[index1] *= scale1;
	fcim_dn[index1] *= scale1;
      }//endfor iCoeff
      index1 = iCoeffStart+numCoeff;
      fcim_dn[index1] *= scale2;
    }//endfor iState
  }

/*==========================================================================*/
/* 2) Calculate P_(n+1)(H)|phi>. Everything store in fcre(im)_up(dn) */

  DAXPY(&numCoeffUpTotal,&prefact,&cre_up[1],&incx,&fcre_up[1],&incy);
  DAXPY(&numCoeffUpTotal,&prefact,&cim_up[1],&incx,&fcim_up[1],&incy);
  if(cpLsda==1&&numStateDnProc!=0){
    DAXPY(&numCoeffDnTotal,&prefact,&cre_dn[1],&incx,&fcre_dn[1],&incy);
    DAXPY(&numCoeffDnTotal,&prefact,&cim_dn[1],&incx,&fcim_dn[1],&incy);
  }

/*==========================================================================*/
/* 3) Copy the force back to the coefficients */
  //debug
  /*
  for(iState=0;iState<numStateUpProc;iState++){
    dot = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      dot += cre_up[iState*numCoeff+iCoeff]*fcre_up[iState*numCoeff+iCoeff]+
	    cim_up[iState*numCoeff+iCoeff]*fcim_up[iState*numCoeff+iCoeff];
    }
    dot *= 2.0;
    dot += cre_up[iState*numCoeff+numCoeff]*fcre_up[iState*numCoeff+numCoeff];
    dot *= 0.5;
    printf("zn %lg dot %lg\n",zn,dot);
  } 
  */
  for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
    cre_up[iCoeff] = fcre_up[iCoeff];
    cim_up[iCoeff] = fcim_up[iCoeff];
  }
  if(cpLsda==1&&numStateDnProc!=0){
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      cre_dn[iCoeff] = fcre_dn[iCoeff];
      cim_dn[iCoeff] = fcim_dn[iCoeff];
    }
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void normHCheby(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                 CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos,int iPoly)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Different from the Newton case, the first step is diffierent from all */
/* the others. Also, to save memory, we shall return force vector for    */
/* T_(n+1).								 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  COMMUNICATE *communicate      = &(cp->communicate);

  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;
  double Smin           = chebyshevInfo->Smin;
  double Smax           = chebyshevInfo->Smax;
  double scale          = chebyshevInfo->scale;
  double energyMean     = stodftInfo->energyMean;
  double energyDiff     = stodftInfo->energyDiff;
  double prefact        = -scale*energyMean;
  double scale1         = -scale*0.25;
  double scale2         = -scale*0.5;

  int cpLsda         = cpopts->cp_lsda;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int cpDualGridOptOn  = cpopts->cp_dual_grid_opt;
  int numCoeffM1     = numCoeff-1;
  int incx = 1;
  int incy = 1;
  int iState,iCoeff,iCoeffStart,index1,index2;

  double *expanCoeff = (double*)stodftCoefPos->expanCoeff;

  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  double *wfUpRe1 = stodftCoefPos->wfUpRe1;
  double *wfUpIm1 = stodftCoefPos->wfUpIm1;
  double *wfDnRe1 = stodftCoefPos->wfDnRe1;
  double *wfDnIm1 = stodftCoefPos->wfDnIm1;

  double dot;

/*==========================================================================*/
/* 1) Calculate the H(norm)|phi> */
  //control_vps_atm_list will be done somewhere else (perhaps in density calculation?)

  //calcCoefForceWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  /*
  calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  calcCoefForceWrapReduce(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  */
  calcCoefForceWrapSCF(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

/*==========================================================================*/
/* 2) Calculate T_(n+1)|phi> */


  if(iPoly==1){//iPoly=1
    for(iState=0;iState<numStateUpProc;iState++){
      iCoeffStart = iState*numCoeff;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	index1 = iCoeffStart+iCoeff;
	fcre_up[index1] = scale1*fcre_up[index1]+prefact*cre_up[index1];
	fcim_up[index1] = scale1*fcim_up[index1]+prefact*cim_up[index1];
      }//endfor iCoeff
      index1 = iCoeffStart+numCoeff;
      fcre_up[index1] = scale2*fcre_up[index1]+prefact*cre_up[index1];
    }//endfor iState
    if(cpLsda==1&&numStateDnProc!=0){
      for(iState=0;iState<numStateDnProc;iState++){
	iCoeffStart = iState*numCoeff;
	for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	  index1 = iCoeffStart+iCoeff;
	  fcre_dn[index1] = scale1*fcre_dn[index1]+prefact*cre_dn[index1];
	  fcim_dn[index1] = scale1*fcim_dn[index1]+prefact*cim_dn[index1];
	}//endfor iCoeff
	index1 = iCoeffStart+numCoeff;
	fcim_dn[index1] = scale2*fcre_dn[index1]+prefact*cre_dn[index1];
      }//endfor iState
    }
  }
  else{//iPoly>1
    scale1 *= 2.0;
    scale2 *= 2.0;
    prefact *= 2.0;
    for(iState=0;iState<numStateUpProc;iState++){
      iCoeffStart = iState*numCoeff;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
        index1 = iCoeffStart+iCoeff;
        fcre_up[index1] = scale1*fcre_up[index1]+prefact*cre_up[index1]-wfUpRe1[index1];
        fcim_up[index1] = scale1*fcim_up[index1]+prefact*cim_up[index1]-wfUpIm1[index1];
      }//endfor iCoeff
      index1 = iCoeffStart+numCoeff;
      fcre_up[index1] = scale2*fcre_up[index1]+prefact*cre_up[index1]-wfUpRe1[index1];
    }//endfor iState
    if(cpLsda==1&&numStateDnProc!=0){
      for(iState=0;iState<numStateDnProc;iState++){
        iCoeffStart = iState*numCoeff;
        for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
          index1 = iCoeffStart+iCoeff;
          fcre_dn[index1] = scale1*fcre_dn[index1]+prefact*cre_dn[index1]-wfDnRe1[index1];
          fcim_dn[index1] = scale1*fcim_dn[index1]+prefact*cim_dn[index1]-wfDnIm1[index1];
        }//endfor iCoeff
        index1 = iCoeffStart+numCoeff;
        fcim_dn[index1] = scale2*fcre_dn[index1]+prefact*cre_dn[index1]-wfDnRe1[index1];
      }//endfor iState
    }
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


