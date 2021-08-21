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
  CELL *cell                    = &(general_data->cell);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  EWALD *ewald                  = &(general_data->ewald);
  PTENS *ptens                  = &(general_data->ptens);
  FRAGINFO *fragInfo            = stodftInfo->fragInfo;
  METALLIC *metallic            = stodftInfo->metallic;
  PARA_FFT_PKG3D *cp_para_fft_pkg3d;

  int pseudoRealFlag = pseudoReal->pseudoRealFlag;
  int cpLsda         = cpopts->cp_lsda;
  int realSparseOpt  = cpopts->realSparseOpt;
  int cpGGA  = cpopts->cp_gga;
  int numStateStoUp  = stodftInfo->numStateStoUp;
  int atomForceFlag  = stodftInfo->atomForceFlag;
  int chemPotOpt     = stodftInfo->chemPotOpt;
  int calcFragFlag   = stodftInfo->calcFragFlag;
  int energyWindowOn = stodftInfo->energyWindowOn;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int numCoeffLarge       = cpcoeffs_info->ncoef_l;
  int numCoeffLargeProc;
  //int numCoeffLargeProc   = cp->cp_para_fft_pkg3d_lg.ncoef_proc;
  int numCoeffLargeProcDensCpBox = cp->cp_para_fft_pkg3d_dens_cp_box.ncoef_proc;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int rhoRealGridNum    = stodftInfo->rhoRealGridNum;
  int numInterpPmeDual = pseudo->n_interp_pme_dual;
  int iperd           = cell->iperd;
  int numChemPot = stodftInfo->numChemPot;
  int occNumber = stodftInfo->occNumber;
  int myidState         = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int numAtomTot = clatoms_info->natm_tot;
  int iState,jState,iCoeff,iChem,iAtom,iGrid;
  int ioff,iis;
  int smearOpt = stodftInfo->smearOpt;
  int numStateFric = metallic->numStateFric;
  int numAtomFricProc = metallic->numAtomFricProc;
  int atomType,atomIndex;
  int hDevSendCount;
  
  int *atomFricIndProc = metallic->atomFricIndProc;
  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *numNlppAtom = pseudoReal->numNlppAtom;
  int *hDevRecvCounts;
  int *hDevRecvDispls;
  int *numAtomFricProcAllProc;

  double tpi = 2.0*M_PI;
  double eke,ekeDn;
  double chemPotTrue = stodftInfo->chemPotTrue;
  double energyKe   = stat_avg->cp_eke;
  double energyPnl  = stat_avg->cp_enl;
  double energyHart = stat_avg->cp_ehart;
  double energyEext = stat_avg->cp_eext;
  double energyExc  = stat_avg->cp_exc;
  double vol        = cell->vol;
  double volInv     = 1.0/vol;
  double energyTotElec,energyTot;
  double energyExtTemp,energyExcTemp,energyHartTemp;
  double vInter;
  double vrecip;
  double vself,vbgr;
  double vrecipLocal;
  double entropy = stodftInfo->entropy;
  double smearTemperature = stodftInfo->smearTemperature;
  double sigma = metallic->sigma;

  //double *energyKe  = stodftInfo->energyKe;
  //double *energyPNL = stodftInfo->energyPNL;
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  double *rhoCoeffReUp   = cpscr->cpscr_rho.rhocr_up;
  double *rhoCoeffImUp   = cpscr->cpscr_rho.rhoci_up;
  double *rhoUp          = cpscr->cpscr_rho.rho_up;
  double *rhoCoeffReUpDensCpBox = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
  double *rhoCoeffImUpDensCpBox = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
  double *divRhoxUp       = cpscr->cpscr_grho.d_rhox_up;
  double *divRhoyUp       = cpscr->cpscr_grho.d_rhoy_up;
  double *divRhozUp       = cpscr->cpscr_grho.d_rhoz_up;
  double *d2RhoUp        = cpscr->cpscr_grho.d2_rho_up;
  double *rhoCoeffReDn   = cpscr->cpscr_rho.rhocr_dn;
  double *rhoCoeffImDn   = cpscr->cpscr_rho.rhoci_dn;
  double *rhoDn          = cpscr->cpscr_rho.rho_dn;
  double *rhoCoeffReDnDensCpBox = cpscr->cpscr_rho.rhocr_dn_dens_cp_box;
  double *rhoCoeffImDnDensCpBox = cpscr->cpscr_rho.rhoci_dn_dens_cp_box;
  double *divRhoxDn       = cpscr->cpscr_grho.d_rhox_dn;
  double *divRhoyDn       = cpscr->cpscr_grho.d_rhoy_dn;
  double *divRhozDn       = cpscr->cpscr_grho.d_rhoz_dn;
  double *d2RhoDn        = cpscr->cpscr_grho.d2_rho_dn;
  double *ak2_sm  =  cpewald->ak2_sm;
  double *chemPot = stodftCoefPos->chemPot;
  double *fx = clatoms_pos->fx;
  double *fy = clatoms_pos->fy;
  double *fz = clatoms_pos->fz;
  double *vnlFxCor;
  double *vnlFyCor;
  double *vnlFzCor;
  double *lagFunValue = (double*)cmalloc(numChemPot*sizeof(double));
  double *hmatCP    = cell->hmat_cp;
  double *ksStateChemPotRe = metallic->ksStateChemPotRe;
  double *ksStateChemPotIm = metallic->ksStateChemPotIm;
  double *ksEnergyFric = metallic->ksEnergyFric;
  double *gauValue;
  double *hDevMatTotal;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;
  double **fxNl;
  double **fyNl;
  double **fzNl;

  double *fxTemp;
  double *fyTemp;
  double *fzTemp;
  double *fxNuclei,*fyNuclei,*fzNuclei;
  double *fxNlTrue,*fyNlTrue,*fzNlTrue;
  double *fxUnCor,*fyUnCor,*fzUnCor;
  double *fxLoc,*fyLoc,*fzLoc;

  double *vlocDevMat = metallic->vlocDevMat;
  double *vnlDevMat = metallic->vnlDevMat;
  double *hDevMat = metallic->hDevMat;

/*======================================================================*/
/* 0) Allocation                                                        */

  hDevMat = (double*)cmalloc(numAtomFricProc*);


/*======================================================================*/
/* I) Calculate nlpp of <m|dVnl/dR|n>                                   */

  if(pseudoRealFlag==1){
    calcNlppRealFriction(class,general_data,cp,hDevMat);
  }
  // else: FUTURE DEVELOPMENT


/*======================================================================*/
/* II) Calculate local pp of <m|dVloc/dR|n>                             */

  calcLocalPotFriction(class,general_data,cp,atomIndex);


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
      hDevMatTotal = (double)cmalloc(numAtomFric*3*numStateFric*numStateFric*
                                     sizeof(double));
      numAtomFricProcAllProc = (int*)cmalloc(numProcStates*sizeof(int));
        
    }
    Gather(&numAtomFricProc,1,MPI_INT,numAtomFricProcAllProc,1,MPI_INT,0);
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
            hDevRecvDispls,MPI_DOUBLE);
  }
  else{
    hDevMatAll = hDevMat;
  }
  
  //Calculate friction
  if(myidState==0){
    gauValue = (dobule*)cmalloc(numStateFric*sizeof(double));
    for(iState=0;iState<numStateFric;iState++){
      gauValue[iState] = gaussianReal(ksEnergyFric[iState],chemPotTrue,1.0/sigma);
    }

    for(iAtom=0;iAtom<numAtomFric;iAtom++){
      for(iDim=0;iDim<3;iDim++){
        for(jAtom=0;jAtom<numAtomFric;jAtom++){
          for(jDim=0;jDim<3;jDim++){
            for(iState=0;iState<numStateFric;iState++){
              for(jState=0;jState<numStateFric;jState++){
                fricTensor[(iAtom*3+iDim)*numAtomFric*3+(jAtom*3+jDim)] += 
                           gauValue(iState)*gauValue(jState)*
                           hDevMatAll[(iAtom*3+iDim)*numStateFric*numStateFric+iState*numStateFric+jState]*
                           hDevMatAll[(jAtom*3+jDim)*numStateFric*numStateFric+jState*numStateFric+iState];
              }//endfor jState
            }//endfor iState
            fricTensor[(iAtom*3+iDim)*numAtomFric*3+(jAtom*3+jDim)] *= -M_PI;
          }//endfor jDim
        }//endfor jAtom
      }//endfor iDim
    }//endfor iAtom
    // Print tensor
  }//endif myidState

  if(myidState==0){
    free(hDevMatTotal);
  }

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


/*======================================================================*/
/* I) Generate noise orbitals |X>                                       */
  

/*======================================================================*/
/* II) Calculate |Y>=(D_a H)|X>                                         */


/*======================================================================*/
/* III) Calculate |Z>=P(u)|Y>                                           */


/*======================================================================*/
/* IV) Calculate |A>=P(u)|X>                                            */


/*======================================================================*/
/* V) Calculate |B>=(D_b H)|A>                                          */


/*======================================================================*/
/* VI) Calculate <Z|B>                                                  */


/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

