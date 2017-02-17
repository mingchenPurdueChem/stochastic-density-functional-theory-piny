/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: energy-cor.c                                 */
/*                                                                          */
/* This routine calculate the kinectic energy, non-local pseudo-potential   */
/* energy and nuclei force correction.					    */
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
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_frag_local.h"
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void energyCorrect(CP *cpMini,GENERAL_DATA *generalDataMini,CLASS *classMini,
		   CP *cp,int ip_now)
/*========================================================================*/
{/*begin routine*/
/**************************************************************************/
/* This function calculate the energy correction, including kinetic	  */
/* energy, non-local pseudo potential energy and nuclei forces comes from */
/* non-local pseudo potential.						  */
/**************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  CPOPTS *cpOpts = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cp->cpcoeffs_info);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  COMMUNICATE *commCP = &(cp->communicate);

  int myidState             = commCP->myid_state;
  int numProcStates         = commCP->np_states;
  int numFragProc	    = fragInfo->numFragProc;
  int iFrag;
  double keCorProc = 0.0;
  MPI_Comm commStates   =    commCP->comm_states;


  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInfo->iFrag = iFrag;

/*======================================================================*/
/* I) Kinetic energy	                                                */
    
    calcKECor(cpMini,generalDataMini,cp,&keCorProc);

/*======================================================================*/
/* II) Non-local pseudo potential energy                                */


/*======================================================================*/
/* III) Non-local pseudo potential force                                */
  }

/*======================================================================*/
/* I) Reduce everything                                                 */

  if(numProcStates>1){
    Reduce(&keCorProc,&(fragInfo->keCor),1,MPI_DOUBLE,MPI_SUM,0,commStates);
  }


/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcKECor(CP *cpMini,GENERAL_DATA *generalDataMini,CP *cp,double *keCorProc)
/*========================================================================*/
{/*begin routine*/
/**************************************************************************/
/* This function calculates the kinetic energy correction. \sum_f ke_f    */
/* -\sum_f a_f^T B_f a_f where ke_f = \sum_i <psi_i^f|K|psi_i^f> , a_f(i) */
/* = <kai|psi_f^i> , and B_f(i,j)=<psi_f^i|K|psi_f^j>.			  */
/**************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  CPOPTS *cpOpts = &(cpMini->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cpMini->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos = &(cpMini->cpcoeffs_pos[1]);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cpMini->cp_para_fft_pkg3d_lg);
  COMMUNICATE *commCP = &(cpMini->communicate);
  STAT_AVG *statAvg = &(generalDataMini->stat_avg);
  

  int iState,jState,iCoeff;
  int iFrag = fragInfo->iFrag;
  int cpLsda = cpOpts->cp_lsda;
  int numFragProc           = fragInfo->numFragProc;
  int numFragTot            = fragInfo->numFragTot;
  int numStateUp = cpcoeffs_info->nstate_up_proc;
  int numStateDn = cpcoeffs_info->nstate_dn_proc;
  double *keMatrixUp,*keMatrixDn;
  double *wfProjUp,*wfProjDn;
  double *temp;
  double keCor;
  double ke = statAvg->kinet_cp;

/*======================================================================*/
/* I) Calculate the matrix                                              */

  calcKEMatrix(cpMini,cp);

/*======================================================================*/
/* I) Allocate Local Memory                                             */

  
  keMatrixUp = fragInfo->keMatrixUp[iFrag];
  wfProjUp = fragInfo->wfProjUp[iFrag];
  temp = (double*)cmalloc(numStateUp*sizeof(double));
  dsymvWrapper('U',numStateUp,1.0,keMatrixUp,wfProjUp,1,0.0,temp,1);
  keCor = ddotBlasWrapper(numStateUp,temp,1,wfProjUp,1);
  free(temp);
  if(cpLsda==1&&numStateDn!=0){
    keMatrixDn = fragInfo->keMatrixDn[iFrag];
    wfProjDn = fragInfo->wfProjDn[iFrag];
    temp = (double*)cmalloc(numStateDn*sizeof(double));
    dsymvWrapper('U',numStateDn,1.0,keMatrixDn,wfProjDn,1,0.0,temp,1);
    keCor += ddotBlasWrapper(numStateDn,temp,1,wfProjDn,1);
    free(temp);
  }
  *keCorProc += ke-keCor;

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcKEMatrix(CP *cpMini,CP *cp)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;  
  CPOPTS *cpOpts = &(cpMini->cpopts);
  CPEWALD *cpEwald = &(cpMini->cpewald);
  CPCOEFFS_INFO *cpcoeffs_info = &(cpMini->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos = &(cpMini->cpcoeffs_pos[1]);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cpMini->cp_para_fft_pkg3d_lg);
  COMMUNICATE *commCP = &(cpMini->communicate);

  int iState,jState,iCoeff;
  int numStateUp = cpcoeffs_info->nstate_up_proc;
  int numStateDn = cpcoeffs_info->nstate_dn_proc;
  int numCoeff = cpcoeffs_info->ncoef;
  int cpLsda = cpOpts->cp_lsda;
  int numCoeffUpTot = numStateUp*numCoeff;
  int numCoeffDnTot = numStateDn*numCoeff;
  int index,index1,index2,index3;
  int iFrag = fragInfo->iFrag;
  
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *ak2Small = cpEwald->ak2_sm;
  double *coefForceRe,*coefForceIm;
  double *keMatrixUp;
  double *keMatrixDn;
  double tpi = 2.0*M_PI;

/*======================================================================*/
/* I) Allocate Local Memory                                             */


  coefForceRe = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  coefForceIm = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;


/*======================================================================*/
/* II) Calculate Spin up matrix                                         */

  keMatrixUp = fragInfo->keMatrixUp[iFrag];

  for(iState=0;iState<numStateUp;iState++){
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      index = iState*numCoeff+iCoeff;
      coefForceRe[index] = 0.5*ak2Small[iCoeff]*cre_up[index];
      coefForceIm[index] = 0.5*ak2Small[iCoeff]*cim_up[index];
    }//endfor iCoeff
    coefForceRe[iState*numCoeff+numCoeff] = 0.0;
    coefForceIm[iState*numCoeff+numCoeff] = 0.0;
  }//endfor iState

  for(iState=0;iState<numStateUp;iState++){
    for(jState=iState;jState<numStateUp;jState++){
      index = iState*numStateUp+jState;
      index1 = jState*numStateUp+iState;
      keMatrixUp[index] = 0.0;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	index2 = iState*numCoeff+iCoeff;
	index3 = jState*numCoeff+iCoeff;
        // We should have *2 here since we have CC* and C*C but the wave functions
	// are normalized to 2, so there is another 0.5 to bring it back to normal 
	// But don't forget to scale everything by occupied number at the end of the day
	keMatrixUp[index] += coefForceRe[index2]*cre_up[index3]+coefForceIm[index2]*cim_up[index3];
      }//endfor iCoeff
      keMatrixUp[index1] = keMatrixUp[index];
    }//endfor jState
  }//endfor iState
  
/*======================================================================*/
/* III) Calculate Spin down matrix					*/


  if(cpLsda==1&&numStateDn!=0){// spin down
    keMatrixDn = fragInfo->keMatrixDn[iFrag];
    for(iState=0;iState<numStateDn;iState++){
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	index = iState*numCoeff+iCoeff;
	coefForceRe[index] = 0.5*ak2Small[iCoeff]*cre_dn[index];
	coefForceIm[index] = 0.5*ak2Small[iCoeff]*cim_dn[index];
      }
      coefForceRe[iState*numCoeff+numCoeff] = 0.0;
      coefForceIm[iState*numCoeff+numCoeff] = 0.0;
    }

    for(iState=0;iState<numStateDn;iState++){
      for(jState=iState;jState<numStateDn;jState++){
	index = iState*numStateDn+jState;
	index1 = jState*numStateDn+iState;
	keMatrixDn[index] = 0.0;
	for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	  index2 = iState*numCoeff+iCoeff;
	  index3 = jState*numCoeff+iCoeff;
	  // We should have *2 here since we have CC* and C*C but the wave functions
	  // are normalized to 2, so there is another 0.5 to bring it back to normal 
	  // But don't forget to scale everything by occupied number at the end of the day
	  keMatrixDn[index] += coefForceRe[index2]*cre_dn[index3]+coefForceIm[index2]*cim_dn[index3];
	}//endfor iCoeff
	keMatrixDn[index1] = keMatrixDn[index];
      }//endfor jState
    }//endfor iState   
  }//endif 


  /*
  //eke = 0.0;
  for(is=1 ; is<= nstate ; is++){
    ioff = (is-1)*ncoef;
    for(i=1; i<= ncoef1 ; i++){
      iis = ioff + i;
      fccreal[iis] -= 2.0*ak2_sm[i]*ccreal[iis];
      fccimag[iis] -= 2.0*ak2_sm[i]*ccimag[iis];
      //eke += (2.0*ak2_sm[i]*(ccreal[iis]*ccreal[iis] + ccimag[iis]*ccimag[iis]));
    }//endfor i
   nis = is*ncoef;
   fccimag[nis] = 0.0;
  }//endfor
  */

/*======================================================================*/
/* IV) free local memories						*/

  free(coefForceRe);
  free(coefForceIm);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

