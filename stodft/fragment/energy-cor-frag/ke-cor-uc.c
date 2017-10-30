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
/* energy and nuclei force correction.                                      */
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
void calcKECorUC(CP *cpMini,GENERAL_DATA *generalDataMini,CLASS *classMini, 
		 CP *cp,double *keCorProc)
/*========================================================================*/
{/*begin routine*/
/**************************************************************************/
/* This function calculates the kinetic energy correction. \sum_f ke_f    */
/* -\sum_f a_f^T B_f a_f where ke_f = \sum_i \int_A dr<psi_i^f|K|psi_i^f>,*/
/* a_f(i) = <kai|psi_f^i> , and B_f(i,j)=\int_A dr<psi_f^i|r><r|K|psi_f^j>*/
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


  int iState,jState,iCoeff,iStoc;
  int iFrag = fragInfo->iFrag;
  int cpLsda = cpOpts->cp_lsda;
  int numFragProc           = fragInfo->numFragProc;
  int numFragTot            = fragInfo->numFragTot;
  int numStateUp = cpcoeffs_info->nstate_up_proc;
  int numStateDn = cpcoeffs_info->nstate_dn_proc;
  int numStateStoUp = stodftInfo->numStateStoUp;
  int numStateStoDn = stodftInfo->numStateStoDn;
  int occNumber = stodftInfo->occNumber;
  double *keMatrixUp,*keMatrixDn;
  double *wfProjUp,*wfProjDn;
  double *temp;
  double keCor,keCorUp,keCorDn;
  double ke;

/*======================================================================*/
/* I) Calculate the matrix                                              */

  calcKEMatrixUC(generalDataMini,cpMini,classMini,cp,&ke);

/*======================================================================*/
/* I) Allocate Local Memory                                             */

  //printf("ke %lg\n",ke);
  keMatrixUp = fragInfo->keMatrixUp[iFrag];
  wfProjUp = fragInfo->wfProjUp[iFrag];
  //debug
  /*
  for(iState=0;iState<numStateUp;iState++){
    for(jState=iState;jState<numStateUp;jState++){
      printf("iState %i jState %i Matrix %lg\n",iState,jState,keMatrixUp[iState*numStateUp+jState]);
    }
  }
  for(iStoc=0;iStoc<numStateStoUp;iStoc++){
    for(iState=0;iState<numStateUp;iState++){
      printf("iStoch %i iState %i wfProj %lg\n",iStoc,iState,wfProjUp[iStoc*numStateUp+iState]);
    }
  }
  */
  temp = (double*)cmalloc(numStateUp*sizeof(double));
  keCorUp = 0.0;
  for(iStoc=0;iStoc<numStateStoUp;iStoc++){
    dsymvWrapper('U',numStateUp,1.0,keMatrixUp,numStateUp,&wfProjUp[iStoc*numStateUp],1,0.0,temp,1);
    keCorUp += ddotBlasWrapper(numStateUp,temp,1,&wfProjUp[iStoc*numStateUp],1);
  }
  keCorUp /= numStateStoUp;
  keCor = keCorUp;
  free(temp);
  if(cpLsda==1&&numStateDn!=0){
    keMatrixDn = fragInfo->keMatrixDn[iFrag];
    wfProjDn = fragInfo->wfProjDn[iFrag];
    temp = (double*)cmalloc(numStateDn*sizeof(double));
    for(iStoc=0;iStoc<numStateStoDn;iStoc++){
      dsymvWrapper('U',numStateDn,1.0,keMatrixDn,numStateDn,&wfProjDn[iStoc*numStateDn],1,0.0,temp,1);
      keCorDn += ddotBlasWrapper(numStateDn,temp,1,&wfProjDn[iStoc*numStateDn],1);
    }
    keCorDn /= numStateStoDn;
    keCor += keCorDn;
    free(temp);
  }
  printf("ke %lg keCor %lg\n",ke,keCor);
  *keCorProc += ke-occNumber*keCor;

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcKEMatrixUC(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *classMini, 
		    CP *cp, double *ke)
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

  int iState,jState,iCoeff,iGrid;
  int numStateUp = cpcoeffs_info->nstate_up_proc;
  int numStateDn = cpcoeffs_info->nstate_dn_proc;
  int numCoeff = cpcoeffs_info->ncoef;
  int cpLsda = cpOpts->cp_lsda;
  int numCoeffUpTot = numStateUp*numCoeff;
  int numCoeffDnTot = numStateDn*numCoeff;
  int numAlloc = MAX(numCoeffUpTot,numCoeffDnTot);
  int index,index1,index2,index3;
  int iFrag = fragInfo->iFrag;
  int numGrid = fragInfo->numGridFragProc[iFrag];
  int numGridSmall = fragInfo->numGridFragProcSmall[iFrag];
  int occNumber = stodftInfo->occNumber;

  int *gridMapProcSmall    = fragInfo->gridMapProcSmall[iFrag];

  double tpi = 2.0*M_PI;
  double volMini;
  double keLocal;
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  double *ak2Small = cpEwald->ak2_sm;
  double *coefForceRe,*coefForceIm;
  double *coefTemp,*wfTemp;
  double *keMatrixUp;
  double *keMatrixDn;
  double *coefUpFragCoreProc;
  double *coefDnFragCoreProc;
  double *coefUpFragProc,*coefDnFragProc;
  double *coefUpTemp,*coefDnTemp;
  double *hmatCpMini = generalDataMini->cell.hmat_cp;
  

/*======================================================================*/
/* I) Allocate Local Memory                                             */

  fragInfo->rhoUpFragProc[iFrag] = (double*)cmalloc(numGrid*sizeof(double));
  fragInfo->coefUpFragProc[iFrag] = (double*)cmalloc(numStateUp*numGrid*sizeof(double));
  coefUpFragProc = fragInfo->coefUpFragProc[iFrag];
  coefUpFragCoreProc = fragInfo->coefUpFragCoreProc[iFrag];  
  if(cpLsda==1&&numStateDn!=0){
    fragInfo->rhoDnFragProc[iFrag] = (double*)cmalloc(numGrid*sizeof(double));
    fragInfo->coefDnFragProc[iFrag] = (double*)cmalloc(numStateDn*numGrid*sizeof(double));
    coefDnFragProc = fragInfo->coefDnFragProc[iFrag];
    coefDnFragCoreProc = fragInfo->coefDnFragCoreProc[iFrag];
  }

  coefTemp = (double*)cmalloc(numCoeff*sizeof(double));
  wfTemp = (double*)cmalloc(numGridSmall*sizeof(double));
  //coefForceRe = (double*)cmalloc((numAlloc+1)*sizeof(double));
  //coefForceIm = (double*)cmalloc((numAlloc+1)*sizeof(double));
  //debug
  /*
  coefUpTemp = (double*)cmalloc(numStateUp*numGrid*sizeof(double));
  if(cpLsda==1&&numStateDn!=0)coefDnTemp = (double*)cmalloc(numStateUp*numGrid*sizeof(double));

  rhoRealCalcDriverFragMol(generalDataMini,cpMini,classMini,cp);
  
  memcpy(coefUpTemp,coefUpFragProc,numStateUp*numGrid*sizeof(double));
  if(cpLsda==1&&numStateDn!=0)memcpy(coefDnTemp,coefDnFragProc,numStateUp*numGrid*sizeof(double));
  */

/*======================================================================*/
/* II) Calculate kinetic energy density                                 */

  keMatrixUp = fragInfo->keMatrixUp[iFrag];
  keLocal = 0;

  //printf("cre %lg\n",cre_up[1]);
  //printf("ak2Small %lg\n",ak2Small[1]);
  for(iState=0;iState<numStateUp;iState++){
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      index = iState*numCoeff+iCoeff;
      fcre_up[index] = 0.5*ak2Small[iCoeff]*cre_up[index];
      fcim_up[index] = 0.5*ak2Small[iCoeff]*cim_up[index];
    }//endfor iCoeff
    fcre_up[iState*numCoeff+numCoeff] = 0.0;
    fcim_up[iState*numCoeff+numCoeff] = 0.0;
    memcpy(coefTemp,&fcre_up[iState*numCoeff+1],numCoeff*sizeof(double));
    memcpy(&fcre_up[iState*numCoeff+1],&cre_up[iState*numCoeff+1],numCoeff*sizeof(double));
    memcpy(&cre_up[iState*numCoeff+1],coefTemp,numCoeff*sizeof(double));    
    memcpy(coefTemp,&fcim_up[iState*numCoeff+1],numCoeff*sizeof(double));
    memcpy(&fcim_up[iState*numCoeff+1],&cim_up[iState*numCoeff+1],numCoeff*sizeof(double));
    memcpy(&cim_up[iState*numCoeff+1],coefTemp,numCoeff*sizeof(double));
  }//endfor iState
  //printf("cre %lg fcre %lg\n",cim_up[2],fcim_up[1]);
  if(cpLsda==1&&numStateDn!=0){
    for(iState=0;iState<numStateDn;iState++){
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	index = iState*numCoeff+iCoeff;
	fcre_dn[index] = 0.5*ak2Small[iCoeff]*cre_dn[index];
	fcim_dn[index] = 0.5*ak2Small[iCoeff]*cim_dn[index];
      }//endfor iCoeff
      fcre_dn[iState*numCoeff+numCoeff] = 0.0;
      fcim_dn[iState*numCoeff+numCoeff] = 0.0;
      memcpy(coefTemp,&fcre_dn[iState*numCoeff+1],numCoeff*sizeof(double));
      memcpy(&fcre_dn[iState*numCoeff+1],&cre_dn[iState*numCoeff+1],numCoeff*sizeof(double));
      memcpy(&cre_dn[iState*numCoeff+1],coefTemp,numCoeff*sizeof(double));
      memcpy(coefTemp,&fcim_dn[iState*numCoeff+1],numCoeff*sizeof(double));
      memcpy(&fcim_dn[iState*numCoeff+1],&cim_dn[iState*numCoeff+1],numCoeff*sizeof(double));
      memcpy(&cim_dn[iState*numCoeff+1],coefTemp,numCoeff*sizeof(double));
    }//endfor iState
  }

  int stodftOnTemp = stodftInfo->stodftOn;
  stodftInfo->stodftOn = -1;
  rhoRealCalcDriverFragMol(generalDataMini,cpMini,classMini,cp);
  stodftInfo->stodftOn = stodftOnTemp;

  memcpy(&cre_up[1],&fcre_up[1],numStateUp*numCoeff*sizeof(double));
  memcpy(&cim_up[1],&fcim_up[1],numStateUp*numCoeff*sizeof(double));
  if(cpLsda==1&&numStateDn!=0){
    memcpy(&cre_dn[1],&fcre_dn[1],numStateUp*numCoeff*sizeof(double));
    memcpy(&cim_dn[1],&fcim_dn[1],numStateUp*numCoeff*sizeof(double));
  }

/*======================================================================*/
/* III) Calculate Spin up matrix                                        */

  volMini = getdeth(hmatCpMini)/numGrid;

  /*
  printf("coefUpFragProc %lg\n",coefUpFragProc[0]);

  double test = 0.0;
  for(iState=0;iState<numStateUp;iState++){
    test += ddotBlasWrapper(numGrid,&coefUpFragProc[iState*numGrid],1,&coefUpTemp[iState*numGrid],1)*volMini;
  }
  printf("testke %lg\n",test*2.0);
  */

  for(iState=0;iState<numStateUp;iState++){
    for(iGrid=0;iGrid<numGridSmall;iGrid++){
      wfTemp[iGrid] = coefUpFragProc[iState*numGrid+gridMapProcSmall[iGrid]];
    }
    for(jState=iState;jState<numStateUp;jState++){
      index = iState*numStateUp+jState;
      index1 = jState*numStateUp+iState;
      keMatrixUp[index] = ddotBlasWrapper(numGridSmall,&wfTemp[0],1,&coefUpFragCoreProc[jState*numGridSmall],1)*volMini;
      if(jState==iState){
	keLocal += keMatrixUp[index];
	//printf("iState %i keLocal %lg\n",iState,keMatrixUp[index]);
      }
      keMatrixUp[index1] = keMatrixUp[index];
    }//endfor jState
  }//endfor iState
  printf("iFrag %i keLocal %lg\n",iFrag,keLocal);  
  //fflush(stdout);
  //exit(0);

/*======================================================================*/
/* III) Calculate Spin down matrix                                      */


  if(cpLsda==1&&numStateDn!=0){// spin down
    keMatrixDn = fragInfo->keMatrixDn[iFrag];
    for(iState=0;iState<numStateDn;iState++){
      for(iGrid=0;iGrid<numGridSmall;iGrid++){
        wfTemp[iGrid] = coefDnFragProc[iState*numGrid+gridMapProcSmall[iGrid]];
      }
      for(jState=iState;jState<numStateDn;jState++){
        index = iState*numStateDn+jState;
        index1 = jState*numStateDn+iState;
        keMatrixDn[index] = ddotBlasWrapper(numGridSmall,&wfTemp[0],1,&coefDnFragCoreProc[jState*numGridSmall],1)*volMini;
        if(jState==iState)keLocal += keMatrixDn[index];
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
  *ke = keLocal*occNumber;

/*======================================================================*/
/* IV) free local memories                                              */

  //free(coefForceRe);
  //free(coefForceIm);
  free(fragInfo->rhoUpFragProc[iFrag]);
  free(fragInfo->coefUpFragProc[iFrag]);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


