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
void calcKECorUC(CP *cpMini,GENERAL_DATA *generalDataMini,CP *cp,double *keCorProc)
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

  calcKEMatrixUC(cpMini,cp,&ke);

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
  //printf("ke %lg keCor %lg\n",ke,keCor);
  *keCorProc += ke-occNumber*keCor;

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcKEMatrix(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *classMini,CP *cp)
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
  int numGrid = fragInfo->numGridFragProc[iFrag];
  int numGridSmall = fragInfo->numGridFragProcSmall[iFrag];
  int numAlloc = MAX(numCoeffUpTot,numCoeffDnTot);
  int index,index1,index2,index3;
  int iFrag = fragInfo->iFrag;
  int *gridMapProcSmall    = fragInfo->gridMapProcSmall[iFrag];

  double tpi = 2.0*M_PI;
  double volMini;
  double occNum;
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
  double *coefTemp,wfTemp;
  double *keMatrixUp;
  double *keMatrixDn;
  double *coefUpFragCoreProc;
  double *coefDnFragCoreProc;
  double *coefUpFragProc,coefDnFragProc;
  double *hmatCpMini = generalDataMini->cell.hmat_cp;
  

/*======================================================================*/
/* I) Allocate Local Memory                                             */

  fragInfo->rhoUpFragProc[iFrag] = (double*)cmalloc(numGrid*sizeof(double));
  fragInfo->coefUpFragProc[iFrag] = (double*)cmalloc(numStateUp*numGrid*sizeof(double));
  coefUpFragProc = fragInfo->coefUpFragProc[iFrag];
  if(cpLsda==1&&numStateDn!=0){
    fragInfo->rhoDnFragProc[iFrag] = (double*)cmalloc(numGrid*sizeof(double));
    fragInfo->coefDnFragProc[iFrag] = (double*)cmalloc(numStateDn*numGrid*sizeof(double));
    coefDnFragProc = fragInfo->coefDnFragProc[iFrag];
  }

  coefTemp = (double*)cmalloc(numCoeff*sizeof(double));
  wfTemp = (double*)cmalloc(numGridSmall*sizeof(double));
  //coefForceRe = (double*)cmalloc((numAlloc+1)*sizeof(double));
  //coefForceIm = (double*)cmalloc((numAlloc+1)*sizeof(double));

/*======================================================================*/
/* II) Calculate kinetic energy density                                 */

  keMatrixUp = fragInfo->keMatrixUp[iFrag];
  occNum = 2.0;
  if(cpLsda==1&&numStateDn!=0) occNum = 1.0;

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
  }//endfor iState
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
    }//endfor iState
  }

  rhoRealCalcDriverFragMol(generalDataMini,cpMini,classMini,cp);

  memcpy(&cre_up[1],&fcre_up[1],numStateUp*numCoeff*sizeof(double));
  if(cpLsda==1&&numStateDn!=0){
    memcpy(&cre_dn[1],&fcre_dn[1],numStateUp*numCoeff*sizeof(double));
  }

/*======================================================================*/
/* III) Calculate Spin up matrix                                        */

  volMini = getdeth(hmatCpMini)/numGrid;

  for(iState=0;iState<numStateUp;iState++){
    for(iGrid=0;iGrid<numGridSmall;iGrid++){
      wfTemp[iGrid] = coefUpFragProc[iState*numGrid+gridMapProcSmall[iGrid]];
    }
    for(jState=iState;jState<numStateUp;jState++){
      index = iState*numStateUp+jState;
      index1 = jState*numStateUp+iState;
      keMatrixUp[index] = ddotBlasWrapper(numGridSmall,&wfTemp[0],1,&coefUpFragCoreProc[jState*numGridSmall],1)*occNum*volMini;
      keMatrixUp[index1] = keMatrixUp[index];
    }//endfor jState
  }//endfor iState

/*======================================================================*/
/* III) Calculate Spin down matrix                                      */


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
/* IV) free local memories                                              */

  free(coefForceRe);
  free(coefForceIm);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genKeDensity(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *classMini,
		  CP *cpMini)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
  CPEWALD *cpewald = &(cpMini->cpewald);
  CPSCR *cpscr = &(cpMini->cpscr);
  CPOPTS *cpopts = &(cpMini->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cpMini->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos = &(cpMini->cpcoeffs_pos[1]);
  COMMUNICATE *communicate = &(cpMini->communicate);
  EWALD *ewald = &(generalDataMini->ewald);
  CELL *cell = &(generalDataMini->cell);
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;

  int i;
  int cp_norb         = cpopts->cp_norb;
  int cpLsda         = cpopts->cp_lsda;
  int myidState      = communicate->myid_state;
  int numStateUpFrag    = cpcoeffs_info->nstate_up_proc;
  int numStateDnFrag    = cpcoeffs_info->nstate_dn_proc;
  int iFrag             = fragInfo->iFrag;

  int *icoef_orth_up    = &(cpMini->cpcoeffs_pos[1].icoef_orth_up);
  int *icoef_form_up    = &(cpMini->cpcoeffs_pos[1].icoef_form_up);
  int *ifcoef_orth_up   = &(cpMini->cpcoeffs_pos[1].ifcoef_orth_up);
  int *ifcoef_form_up   = &(cpMini->cpcoeffs_pos[1].ifcoef_form_up);
  int *icoef_orth_dn    = &(cpMini->cpcoeffs_pos[1].icoef_orth_dn);
  int *icoef_form_dn    = &(cpMini->cpcoeffs_pos[1].icoef_form_dn);
  int *ifcoef_orth_dn   = &(cpMini->cpcoeffs_pos[1].ifcoef_orth_dn);
  int *ifcoef_form_dn   = &(cpMini->cpcoeffs_pos[1].ifcoef_form_dn);

  double *ccrealUpMini    = cpMini->cpcoeffs_pos[1].cre_up;
  double *ccimagUpMini    = cpMini->cpcoeffs_pos[1].cim_up;
  double *ccrealDnMini    = cpMini->cpcoeffs_pos[1].cre_dn;
  double *ccimagDnMini    = cpMini->cpcoeffs_pos[1].cim_dn;
  double *rhoUpFragProc,*rhoDnFragProc,*coefUpFragProc,*coefDnFragProc;

/*======================================================================*/
/* I) Check the forms                                                   */

  if(cp_norb>0){
    if((*icoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up Coefs must be in nonorthonormal form under norb \n");
      printf("on state processor %d in cp_elec_energy_ctrl \n",myidState);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cpLsda==1){
      if((*icoef_orth_dn)!=0){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Dn Coefs must be in nonorthonormal form under norb \n");
        printf("on state processor %d in cp_elec_energy_ctrl \n",myidState);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
    }/*endif*/
  }/*endif*/

/*======================================================================*/
/* III) Initialize Flags, inverse hmat                                  */

  (*ifcoef_form_up) = 0;
  (*ifcoef_orth_up) = 1;
  rhoUpFragProc = fragInfo->rhoUpFragProc[iFrag];
  coefUpFragProc = fragInfo->coefUpFragProc[iFrag];
  if(cpLsda==1&&numStateDnFrag!=0){
    (*ifcoef_form_dn) = 0;
    (*ifcoef_orth_dn) = 1;
    rhoDnFragProc = fragInfo->rhoDnFragProc[iFrag];
    coefDnFragProc = fragInfo->coefDnFragProc[iFrag];
  }

/*======================================================================*/
/* IV) Calculate real space wave functions and densities for fragments  */

  rhoRealCalcFragWrapper(generalDataMini,cpMini,classMini,
                     cp,ccrealUpMini,ccimagUpMini,icoef_form_up,icoef_orth_up,
                     rhoUpFragProc,coefUpFragProc,numStateUpFrag);
  if(cpLsda==1&&numStateDnFrag!=0){
    rhoRealCalcFragWrapper(generalDataMini,cpMini,classMini,
                       cp,ccrealDnMini,ccimagDnMini,icoef_form_dn,icoef_orth_dn,
                       rhoDnFragProc,coefDnFragProc,numStateDnFrag);
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

