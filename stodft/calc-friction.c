/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: calc-nuclei-force.c                            */
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
#include "../proto_defs/proto_stodft_local.h"

#define TIME_CP_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcElectronFricDet(CLASS *class,GENERAL_DATA *general_data,CP *cp,BONDED *bonded,
                        CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is a routine to calculate electron friction tensor after SCF     */
/* loop finishes. (No fragment correction in the current implementation) */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPSCR *cpscr                  = &(cp->cpscr);
  PSEUDO *pseudo                = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal       = &(pseudo->pseudoReal);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  COMMUNICATE *communicate      = &(cp->communicate);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  CPEWALD *cpewald              = &(cp->cpewald);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  ATOMMAPS *atommaps            = &(class->atommaps);
  EWALD *ewald                  = &(general_data->ewald);
  PTENS *ptens                  = &(general_data->ptens);
  FRAGINFO *fragInfo            = stodftInfo->fragInfo;
  METALLIC *metallic            = stodftInfo->metallic;
  PARA_FFT_PKG3D *cp_para_fft_pkg3d;

  int pseudoRealFlag = pseudoReal->pseudoRealFlag;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int myidState         = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int numAtomTot = clatoms_info->natm_tot;
  int iState,jState,iCoeff,iChem,iAtom,iGrid;
  int iDim,iProc,jAtom,jDim,iMat;
  int ioff,iis;
  int smearOpt = stodftInfo->smearOpt;
  int numStateFric = metallic->numStateFric;
  int numAtomFricProc = metallic->numAtomFricProc;
  int numAtomFric = metallic->numAtomFric;
  int atomType,atomIndex;
  int hDevSendCount;

  MPI_Comm commStates   =    communicate->comm_states;
  
  int *atomFricIndProc = metallic->atomFricIndProc;
  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *numNlppAtom = pseudoReal->numNlppAtom;
  int *hDevRecvCounts;
  int *hDevRecvDispls;
  int *numAtomFricProcAllProc;

  double tpi = 2.0*M_PI;
  double chemPotTrue = stodftInfo->chemPotUpMetallic;
  double energyTotElec,energyTot;
  double energyExtTemp,energyExcTemp,energyHartTemp;
  double entropy = stodftInfo->entropy;
  double smearTemperature = stodftInfo->smearTemperature;
  double sigma = metallic->sigma;

  //double *energyKe  = stodftInfo->energyKe;
  //double *energyPNL = stodftInfo->energyPNL;
  double *ksStateChemPotRe = metallic->ksStateChemPotRe;
  double *ksStateChemPotIm = metallic->ksStateChemPotIm;
  double *ksEnergyFric = metallic->ksEnergyFric;
  double *gauValue;
  double *hDevMatTotal;

  double *vlocDevMat;
  double *hDevMat;
  double *fricTensor;
  FILE *ftensor;

  if(myidState==0){
    printf("==============================================\n");
    printf("Calculate Electron Friction\n");
    printf("==============================================\n");
  }
/*======================================================================*/
/* 0) Allocation                                                        */

  hDevMat = (double*)cmalloc(numAtomFricProc*3*numStateFric*numStateFric*
                             sizeof(double));
  vlocDevMat = (double*)cmalloc(numAtomFricProc*3*numStateFric*numStateFric*
                             sizeof(double));

  for(iMat=0;iMat<numAtomFricProc*3*numStateFric*numStateFric;iMat++){
    hDevMat[iMat] = 0.0;
    vlocDevMat[iMat] = 0.0;
  }

/*======================================================================*/
/* I) Calculate nlpp of <m|dVnl/dR|n>                                   */

  if(pseudoRealFlag==1){
    calcNlppRealFriction(class,general_data,cp,hDevMat);
  }
  // else: FUTURE DEVELOPMENT


/*======================================================================*/
/* II) Calculate local pp of <m|dVloc/dR|n>                             */

  calcLocalPotFriction(class,general_data,cp,vlocDevMat);

  for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
    //atomIndex = atomFricIndProc[iAtom];
    //atomType = iAtomAtomType[atomIndex+1]-1;

/*======================================================================*/
/* III) Calculate nlpp of <m|dVnl/dR|n>                                */

    //if(numNlppAtom[atomType]>=1){
    //  calcNlppFriction(class,general_data,cp,atomIndex,iAtom);
    //}

/*======================================================================*/
/* IV) Calculate <m|dh/dR|n>                                            */
  
    for(iState=0;iState<numStateFric;iState++){
      for(jState=0;jState<numStateFric;jState++){
        for(iDim=0;iDim<3;iDim++){
          //hDevMat[(iAtom*3+iDim)*numStateFric*numStateFric+iState*numStateFric+jState] 
          //       = vlocDevMat[iDim*numStateFric*numStateFric+iState*numStateFric+jState]+
          //         vnlDevMat[iDim*numStateFric*numStateFric+iState*numStateFric+jState];
          printf("iAtom %i iDim %i iState %i jState %i nlmat %lg locmat %lg\n",iAtom, iDim,iState,jState,hDevMat[(iAtom*3+iDim)*numStateFric*numStateFric+iState*numStateFric+jState],vlocDevMat[(iAtom*3+iDim)*numStateFric*numStateFric+iState*numStateFric+jState]);
          hDevMat[(iAtom*3+iDim)*numStateFric*numStateFric+iState*numStateFric+jState] += 
                     vlocDevMat[(iAtom*3+iDim)*numStateFric*numStateFric+iState*numStateFric+jState];
        }//endfor iDim
      }//endfor jState
    }//endfor iState
  }//endfor iAtom
/*======================================================================*/
/* V) Calculate friction                                                */
  
  //Reduce hDevMat to the master process
  if(numProcStates>1){
    hDevSendCount = numAtomFricProc*3*numStateFric*numStateFric;
    
    if(myidState==0){
      hDevMatTotal = (double*)cmalloc(numAtomFric*3*numStateFric*numStateFric*
                                     sizeof(double));
      numAtomFricProcAllProc = (int*)cmalloc(numProcStates*sizeof(int));
        
    }
    Gather(&numAtomFricProc,1,MPI_INT,numAtomFricProcAllProc,1,MPI_INT,0,commStates);
    if(myidState==0){
      hDevRecvCounts = (int*)cmalloc(numProcStates*sizeof(int));
      hDevRecvDispls = (int*)cmalloc(numProcStates*sizeof(int));
      for(iProc=0;iProc<numProcStates;iProc++){
        hDevRecvCounts[iProc] = numAtomFricProcAllProc[iProc]*3*numStateFric*numStateFric;
      }
      hDevRecvDispls[0] = 0;
      for(iProc=1;iProc<numProcStates;iProc++){
        hDevRecvDispls[iProc] = hDevRecvDispls[iProc-1]+hDevRecvCounts[iProc-1];
      }
    }
    Gatherv(hDevMat,hDevSendCount,MPI_DOUBLE,hDevMatTotal,hDevRecvCounts,
            hDevRecvDispls,MPI_DOUBLE,0,commStates);
  }
  else{
    hDevMatTotal = hDevMat;
  }
  
  //Calculate friction
  if(myidState==0){
    gauValue = (double*)cmalloc(numStateFric*sizeof(double));
    fricTensor = (double*)cmalloc(numAtomFric*numAtomFric*9);
    for(iState=0;iState<numStateFric;iState++){
      gauValue[iState] = gaussianReal(ksEnergyFric[iState],chemPotTrue,1.0/sigma);
      printf("ksEnergyFric %lg chemPotTrue %lg 1.0/sigma %lg gauValue %lg\n",ksEnergyFric[iState],chemPotTrue,1.0/sigma,gauValue[iState]);
    }

    for(iAtom=0;iAtom<numAtomFric;iAtom++){
      for(iDim=0;iDim<3;iDim++){
        for(jAtom=0;jAtom<numAtomFric;jAtom++){
          for(jDim=0;jDim<3;jDim++){
            for(iState=0;iState<numStateFric;iState++){
              for(jState=0;jState<numStateFric;jState++){
                fricTensor[(iAtom*3+iDim)*numAtomFric*3+(jAtom*3+jDim)] += 
                           gauValue[iState]*gauValue[jState]*
                           hDevMatTotal[(iAtom*3+iDim)*numStateFric*numStateFric+iState*numStateFric+jState]*
                           hDevMatTotal[(jAtom*3+jDim)*numStateFric*numStateFric+jState*numStateFric+iState];
              }//endfor jState
            }//endfor iState
            fricTensor[(iAtom*3+iDim)*numAtomFric*3+(jAtom*3+jDim)] *= -M_PI;
          }//endfor jDim
        }//endfor jAtom
      }//endfor iDim
    }//endfor iAtom
    // Print tensor
    ftensor = fopen("friction-tensor","w");
    for(iAtom=0;iAtom<numAtomFric;iAtom++){
      for(iDim=0;iDim<3;iDim++){
        for(jAtom=0;jAtom<numAtomFric;jAtom++){
          for(jDim=0;jDim<3;jDim++){
            fprintf(ftensor,"%.16lg\n",fricTensor[(iAtom*3+iDim)*numAtomFric*3+(jAtom*3+jDim)]);
          }
        }
      }
    }
    fclose(ftensor);   
  }//endif myidState

  if(myidState==0){
    free(hDevMatTotal);
  }
  free(vlocDevMat);
  free(hDevMat);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcElectronFric(CLASS *class,GENERAL_DATA *general_data,CP *cp,BONDED *bonded,
                        CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is a routine to calculate electron friction tensor after SCF     */
/* loop finishes. (No fragment correction in the current implementation) */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR        *cpscr        = &(cp->cpscr);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  COMMUNICATE   *commCP         = &(cp->communicate);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  METALLIC *metallic            = stodftInfo->metallic;

  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm = &(cp->cp_sclr_fft_pkg3d_sm);

  int numStateStoUp     = stodftInfo->numStateStoUp;
  int numStateStoDn     = stodftInfo->numStateStoDn;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;  
  int cpLsda = cpopts->cp_lsda;
  int numAtomFric = metallic->numAtomFric;
  int iAtom,jAtom,iDim,jDim,index1,index2;
  int iState,iCoeff;
  int myidState         = commCP->myid_state;
  int numProcStates     = commCP->np_states;

  FILE *ftensor;
  MPI_Comm commStates   =    commCP->comm_states;


  double dot;
  double aveFactUp = -M_PI/(double)(numStateStoUp);

  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;
  double *daHChiReUp,*daHChiImUp,*daHChiReDn,*daHChiImDn;
  double *daHXiReUp,*daHXiImUp,*daHXiReDn,*daHXiImDn;
  double *fricTensor;
  double *fricTensorReduce;


/*======================================================================*/
/* I) Memory allocation                                                 */
  fricTensor = (double*)cmalloc(numAtomFric*numAtomFric*9);

  daHChiReUp = (double*)cmalloc(numAtomFric*3*numCoeffUpTotal*sizeof(double));
  daHChiImUp = (double*)cmalloc(numAtomFric*3*numCoeffUpTotal*sizeof(double));
  daHXiReUp = (double*)cmalloc(numAtomFric*3*numCoeffUpTotal*sizeof(double));
  daHXiImUp = (double*)cmalloc(numAtomFric*3*numCoeffUpTotal*sizeof(double));


  if(cpLsda==1){
    daHChiReDn = (double*)cmalloc(numAtomFric*3*numCoeffDnTotal*sizeof(double));
    daHChiImDn = (double*)cmalloc(numAtomFric*3*numCoeffDnTotal*sizeof(double));
    daHXiReDn = (double*)cmalloc(numAtomFric*3*numCoeffDnTotal*sizeof(double));
    daHXiImDn = (double*)cmalloc(numAtomFric*3*numCoeffDnTotal*sizeof(double));
  }

/*======================================================================*/
/* II) Generate coefficients                                            */

  calcChebyCoeffWrapper(stodftInfo,stodftCoefPos,4);

/*======================================================================*/
/* III) Generate noise orbitals |X>                                     */

  genNoiseOrbitalReal(cp,cpcoeffs_pos);
  
/*======================================================================*/
/* IV) Calculate |Y>=(D_a H)|X>                                         */

  genDHPhi(cp,class,general_data,
           coeffReUp,coeffImUp,coeffReDn,coeffImDn,
           daHChiReUp,daHChiImUp,daHChiReDn,daHChiImDn); 

/*======================================================================*/
/* V) Calculate |Z>=P(u)|Y>                                             */

  filterChebyPolyFric(cp,class,general_data,1,numAtomFric*3,
                      daHChiReUp,daHChiImUp,daHChiReDn,daHChiImDn);


/*======================================================================*/
/* VI) Generate noise orbitals |X> again                                */

  genNoiseOrbitalReal(cp,cpcoeffs_pos);


/*======================================================================*/
/* VII) Calculate |A>=P(u)|X>                                           */


  filterChebyPolyFric(cp,class,general_data,1,1,
                      coeffReUp,coeffImUp,coeffReDn,coeffImDn);
  
/*======================================================================*/
/* VIII) Calculate |B>=(D_b H)|A>                                       */

  genDHPhi(cp,class,general_data,
           coeffReUp,coeffImUp,coeffReDn,coeffImDn,
           daHXiReUp,daHXiImUp,daHXiReDn,daHXiImDn);
  
/*======================================================================*/
/* IX) Calculate <Z|B>                                                  */

  for(iAtom=0;iAtom<numAtomFric;iAtom++){
    for(iDim=0;iDim<3;iDim++){
      for(jAtom=0;jAtom<numAtomFric;iAtom++){
        for(jDim=0;jDim<3;jDim++){
          index1 = iAtom*3*numCoeffUpTotal+iDim*numCoeffUpTotal;
          index2 = jAtom*3*numCoeffUpTotal+jDim*numCoeffUpTotal;
          fricTensor[(iAtom*3+iDim)*numAtomFric*3+jAtom*3+jDim] = 0;
          for(iState=0;iState<numStateUpProc;iState++){
            dot = 0.0;
            for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
              dot += daHChiReUp[index1+iState*numCoeff+iCoeff]*daHXiReUp[index2+iState*numCoeff+iCoeff]+
                     daHChiImUp[index1+iState*numCoeff+iCoeff]*daHXiImUp[index2+iState*numCoeff+iCoeff];
            }
            dot *= 2.0;
            dot += daHChiReUp[index1+iState*numCoeff+numCoeff-1]*daHXiReUp[index2+iState*numCoeff+numCoeff-1];
            fricTensor[(iAtom*3+iDim)*numAtomFric*3+jAtom*3+jDim] += dot;
          }
        }
      }
    }
  }

  if(numProcStates>1){
    if(myidState==0)fricTensorReduce = (double*)cmalloc(numAtomFric*numAtomFric*9);
    Reduce(fricTensor,fricTensorReduce,numAtomFric*numAtomFric*9,MPI_DOUBLE,MPI_SUM,0,commStates);
    ftensor = fopen("friction-tensor","w");
    for(iAtom=0;iAtom<numAtomFric;iAtom++){
      for(iDim=0;iDim<3;iDim++){
        for(jAtom=0;jAtom<numAtomFric;jAtom++){
          for(jDim=0;jDim<3;jDim++){
            fprintf(ftensor,"%.16lg\n",fricTensorReduce[(iAtom*3+iDim)*numAtomFric*3+(jAtom*3+jDim)]*aveFactUp);
          }
        }
      }
    }
    fclose(ftensor);
  }


/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

