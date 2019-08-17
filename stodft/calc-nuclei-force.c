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
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_stodft_local.h"

#define TIME_CP_OFF
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcEnergyForce(CLASS *class,GENERAL_DATA *general_data,CP *cp,BONDED *bonded,
			CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is a routine to calculate energy and nuclei forces after SCF     */
/* loop finishes. The fragment force correction is also included here.   */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPSCR *cpscr		        = &(cp->cpscr);
  PSEUDO *pseudo	        = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal	= &(pseudo->pseudoReal);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  COMMUNICATE *communicate      = &(cp->communicate);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  CPEWALD *cpewald              = &(cp->cpewald);
  CELL *cell                    = &(general_data->cell);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  EWALD *ewald			= &(general_data->ewald);
  PTENS *ptens			= &(general_data->ptens);
  FRAGINFO *fragInfo            = stodftInfo->fragInfo;
  PARA_FFT_PKG3D *cp_para_fft_pkg3d;
  
  int pseudoRealFlag = pseudoReal->pseudoRealFlag;
  int cpLsda         = cpopts->cp_lsda;
  int realSparseOpt  = cpopts->realSparseOpt;
  int cpGGA  = cpopts->cp_gga;
  int numStateStoUp  = stodftInfo->numStateStoUp;
  int atomForceFlag  = stodftInfo->atomForceFlag;
  int chemPotOpt     = stodftInfo->chemPotOpt;
  int calcFragFlag   = stodftInfo->calcFragFlag;
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
  int iperd	      = cell->iperd;
  int numChemPot = stodftInfo->numChemPot;
  int occNumber = stodftInfo->occNumber;
  int myidState         = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int numAtomTot = clatoms_info->natm_tot;
  int iState,jState,iCoeff,iChem,iAtom,iGrid;
  int ioff,iis;

  double tpi = 2.0*M_PI;
  double eke,ekeDn;
  double chemPotTrue = stodftInfo->chemPotTrue;
  double energyKe   = stat_avg->cp_eke;
  double energyPnl  = stat_avg->cp_enl;
  double energyHart = stat_avg->cp_ehart;
  double energyEext = stat_avg->cp_eext;
  double energyExc  = stat_avg->cp_exc;
  double vol	    = cell->vol;
  double volInv	    = 1.0/vol;
  double energyTotElec,energyTot;
  double energyExtTemp,energyExcTemp,energyHartTemp;
  double vInter;
  double vrecip;
  double vself,vbgr;
  double vrecipLocal;

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

  MPI_Comm commStates = communicate->comm_states;

/*======================================================================*/
/* 0) Initialize force calculation                                      */
  //if(realSparseOpt==0)cp_para_fft_pkg3d = &(cp->cp_para_fft_pkg3d_lg);
  //else cp_para_fft_pkg3d = &(cp->cp_para_fft_pkg3d_sparse);
  cp_para_fft_pkg3d = &(cp->cp_para_fft_pkg3d_lg);
  numCoeffLargeProc = cp_para_fft_pkg3d->ncoef_proc;

  if(myidState==0){
    energyHartTemp = stat_avg->cp_ehart;
    energyExtTemp = stat_avg->cp_eext;
    energyExcTemp = stat_avg->cp_exc;
  }
  class->energy_ctrl.iget_full_inter = 1;
  class->energy_ctrl.iget_res_inter = 0;
  energy_control_initial(class,bonded,general_data);
  if(calcFragFlag==1){
    vnlFxCor = fragInfo->vnlFxCor;
    vnlFyCor = fragInfo->vnlFyCor;
    vnlFzCor = fragInfo->vnlFzCor;
  }

/*======================================================================*/
/* I) Recalculate k space density                                       */

  for(iGrid=1;iGrid<=rhoRealGridNum;iGrid++){
    rhoUp[iGrid] *= vol;
  }
  if(cpLsda==1&&numStateDnProc!=0){
    for(iGrid=1;iGrid<=rhoRealGridNum;iGrid++){
      rhoDn[iGrid] *= vol;
    }    
  }

  calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
                     rhoCoeffReUp,rhoCoeffImUp,rhoUp,rhoCoeffReUpDensCpBox,rhoCoeffImUpDensCpBox,
                     divRhoxUp,divRhoyUp,divRhozUp,d2RhoUp,cpGGA,cpDualGridOptOn,numInterpPmeDual,
                     communicate,cp_para_fft_pkg3d,&(cp->cp_para_fft_pkg3d_dens_cp_box));
  if(cpLsda==1&&numStateDnProc!=0){
    calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
                       rhoCoeffReDn,rhoCoeffImDn,rhoDn,rhoCoeffReDnDensCpBox,rhoCoeffImDnDensCpBox,
                       divRhoxDn,divRhoyDn,divRhozDn,d2RhoDn,cpGGA,cpDualGridOptOn,numInterpPmeDual,
                       communicate,cp_para_fft_pkg3d,&(cp->cp_para_fft_pkg3d_dens_cp_box));
    for(iCoeff=1;iCoeff<=numCoeffLargeProc;iCoeff++){
      rhoCoeffReUp[iCoeff] += rhoCoeffReDn[iCoeff];
      rhoCoeffImUp[iCoeff] += rhoCoeffImDn[iCoeff];
    }/* endfor */
    if(cpDualGridOptOn>=1){
      for(iCoeff=1;iCoeff<=numCoeffLargeProcDensCpBox;iCoeff++){
        rhoCoeffReUpDensCpBox[iCoeff] += rhoCoeffReDnDensCpBox[iCoeff];
        rhoCoeffImUpDensCpBox[iCoeff] += rhoCoeffImDnDensCpBox[iCoeff];
      }/* endfor */
    } /* endif */
  }/* endif */


/*======================================================================*/
/* I) Calculate Local pp	                                        */

  fxNuclei = (double *)cmalloc(numAtomTot*sizeof(double));
  fyNuclei = (double *)cmalloc(numAtomTot*sizeof(double));
  fzNuclei = (double *)cmalloc(numAtomTot*sizeof(double));
  fxLoc = (double *)cmalloc(numAtomTot*sizeof(double));
  fyLoc = (double *)cmalloc(numAtomTot*sizeof(double));
  fzLoc = (double *)cmalloc(numAtomTot*sizeof(double));
  fxUnCor = (double *)cmalloc(numAtomTot*sizeof(double));
  fyUnCor = (double *)cmalloc(numAtomTot*sizeof(double));
  fzUnCor = (double *)cmalloc(numAtomTot*sizeof(double));

  for(iAtom=1;iAtom<=numAtomTot;iAtom++){
    fx[iAtom] = 0.0;
    fy[iAtom] = 0.0;
    fz[iAtom] = 0.0;
  }

  //debug
  /*
  FILE *fp_rhok = fopen("rho_bm_k","r");
  int ncoef_l = cp->cp_para_fft_pkg3d_lg.ncoef;
  for(iCoeff=1;iCoeff<=ncoef_l;iCoeff++){
    fscanf(fp_rhok,"%lg",&(cp->cpscr.cpscr_rho.rhocr_up[iCoeff]));
    fscanf(fp_rhok,"%lg",&(cp->cpscr.cpscr_rho.rhoci_up[iCoeff]));
    //printf("rho k %lg %lg\n",cp->cpscr.cpscr_rho.rhocr_up[1],cp->cpscr.cpscr_rho.rhoci_up[1]);
  }
  fclose(fp_rhok);
  */

  calcLocExtPostScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  //printf("fx[1] %lg fy[1] %lg fz[1] %lg\n",fx[1],fy[1],fz[1]);
  vrecipLocal = stat_avg->vrecip;
  //printf("vrecipLocal %.16lg\n",vrecipLocal);

  if(numProcStates==1){
    memcpy(&fxLoc[0],&fx[1],numAtomTot*sizeof(double));
    memcpy(&fyLoc[0],&fy[1],numAtomTot*sizeof(double));
    memcpy(&fzLoc[0],&fz[1],numAtomTot*sizeof(double));
  }
  else{
    Reduce(&fx[1],&fxLoc[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&fy[1],&fyLoc[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&fz[1],&fzLoc[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&(stat_avg->vrecip),&vrecipLocal,1,MPI_DOUBLE,MPI_SUM,0,commStates);
  }
  //debug
  //printf("ffffff myidState %i fxLoc %.16lg fyLoc %.16lg fzLoc %.16lg\n",
  //       myidState,fx[1],fy[1],fz[1]);
  /*
  if(myidState==0){
    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      printf("fxloc %lg fyloc %lg fzloc %lg\n",
	      fxLoc[iAtom],fyLoc[iAtom],fzLoc[iAtom]);
    }
  }
  */

/*======================================================================*/
/* II) Calculate nl pp force+energy                                     */
/*--------------------------------------------------------------------------*/
/* 0) Allocate temp arrays				                    */

  fxNl = (double**)cmalloc(numChemPot*sizeof(double*));
  fyNl = (double**)cmalloc(numChemPot*sizeof(double*));
  fzNl = (double**)cmalloc(numChemPot*sizeof(double*));
  for(iChem=0;iChem<numChemPot;iChem++){
    fxNl[iChem] = (double*)cmalloc(numAtomTot*sizeof(double));
    fyNl[iChem] = (double*)cmalloc(numAtomTot*sizeof(double));
    fzNl[iChem] = (double*)cmalloc(numAtomTot*sizeof(double));
  }
  fxNlTrue = (double*)cmalloc(numAtomTot*sizeof(double));
  fyNlTrue = (double*)cmalloc(numAtomTot*sizeof(double));
  fzNlTrue = (double*)cmalloc(numAtomTot*sizeof(double));

/*--------------------------------------------------------------------------*/
/* i) Copy the stochastic wave function and reset the force		    */
  for(iChem=0;iChem<numChemPot;iChem++){
    stat_avg->vrecip = 0.0;
    stat_avg->cp_enl = 0.0;
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      cre_up[iCoeff] = stoWfUpRe[iChem][iCoeff];
      cim_up[iCoeff] = stoWfUpIm[iChem][iCoeff];
      //printf("1111111 Re %lg Im %lg\n",stoWfUpRe[iChem][iCoeff],stoWfUpIm[iChem][iCoeff]);
      fcre_up[iCoeff] = 0.0;
      fcim_up[iCoeff] = 0.0;
    }//endfor iCoeff
    if(cpLsda==1&&numStateDnProc!=0){
      for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
        cre_dn[iCoeff] = stoWfDnRe[iChem][iCoeff];
        cim_dn[iCoeff] = stoWfDnIm[iChem][iCoeff];
        fcre_dn[iCoeff] = 0.0;
        fcim_dn[iCoeff] = 0.0;
      }//endfor iCoeff
    }//endif cpLsda

    //debug
    /*
    int numGrid = fragInfo->numGridFragProc[0];
    int numGridBig = (cp->cp_para_fft_pkg3d_lg.nfft_proc)/2;
    int *gridMapProc = fragInfo->gridMapProc[0];
    double *projRealWF = fragInfo->projRealWF;
    double norm1,norm2,dot;
    double *noiseWfUpReal = fragInfo->noiseWfUpReal;
    double *noiseWfDnReal = fragInfo->noiseWfDnReal;

    rhoRealCalcDriverNoise(general_data,cp,class,1);
    for(iState=0;iState<numStateUpProc;iState++){
      norm1 = ddotBlasWrapper(numGrid,&projRealWF[iState*numGrid],1,&projRealWF[iState*numGrid],1);
      norm2 = ddotBlasWrapper(numGridBig,&noiseWfUpReal[iState*numGridBig],1,&noiseWfUpReal[iState*numGridBig],1);
      norm1 = sqrt(norm1);
      norm2 = sqrt(norm2);
      dot = 0.0;
      for(iGrid=0;iGrid<numGrid;iGrid++){
	dot += noiseWfUpReal[iState*numGridBig+gridMapProc[iGrid]]*projRealWF[iState*numGrid+iGrid];
      }
      printf("iState %i dot %lg\n",iState,dot/norm1/norm2);
    }
    */
    
    //debug
    /*
    int jState;
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      cre_up[iCoeff] = 0.0;
      cim_up[iCoeff] = 0.0;
    }
    for(iState=0;iState<4;iState++){
      FILE *fwfread = fopen("wf-det","r");
      for(jState=0;jState<4;jState++){
	for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
	  fscanf(fwfread,"%lg",&(stoWfUpRe[0][iState*numCoeff*4+jState*numCoeff+iCoeff]));
	  fscanf(fwfread,"%lg",&(stoWfUpIm[0][iState*numCoeff*4+jState*numCoeff+iCoeff]));
	}
      }
      fclose(fwfread);
    }
    double presq = sqrt(2.0);
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      stoWfUpRe[0][iCoeff] *= presq;
      stoWfUpIm[0][iCoeff] *= presq;
      cre_up[iCoeff] = stoWfUpRe[0][iCoeff];
      cim_up[iCoeff] = stoWfUpIm[0][iCoeff];
    }
    */
    
    //pp 
    for(iAtom=1;iAtom<=numAtomTot;iAtom++){
      fx[iAtom] = 0.0;
      fy[iAtom] = 0.0;
      fz[iAtom] = 0.0;
    }

/*--------------------------------------------------------------------------*/
/* ii) Calculate nl pp force						    */
    //debug
    
    /*
    for(iState=0;iState<numStateUpProc;iState++){
      cpcoeffs_info->nstate_up_proc = 1;
      fy[1] = 0.0;
      for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
	cre_up[iCoeff] = 0.0;
	cim_up[iCoeff] = 0.0;
        fcre_up[iCoeff] = 0.0;
        fcim_up[iCoeff] = 0.0;
      }
      for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
	cre_up[iCoeff] = stoWfUpRe[0][iState*numCoeff+iCoeff];
        cim_up[iCoeff] = stoWfUpIm[0][iState*numCoeff+iCoeff];
      }
      calcNlPseudoPostScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
      printf("fyyyyyy %lg\n",fy[1]);
    }
    exit(0);
    */
    
    
    //calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    if(pseudoRealFlag==0){
      calcNlPseudoPostScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    }
    else{
      pseudoReal->forceCalcFlag = 1;
      pseudoReal->nlppForceOnly = 1;
      calcCoefForcePosScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);   
      pseudoReal->forceCalcFlag = 0;
      pseudoReal->nlppForceOnly = 0;
      if(numProcStates>1){
	Barrier(commStates);
	Reduce(&fx[1],&fxNl[iChem][0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
	Reduce(&fy[1],&fyNl[iChem][0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
	Reduce(&fz[1],&fzNl[iChem][0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
	Barrier(commStates);
	if(myidState==0){
	  memcpy(&fx[1],&fxNl[iChem][0],numAtomTot*sizeof(double));
	  memcpy(&fy[1],&fyNl[iChem][0],numAtomTot*sizeof(double));
	  memcpy(&fz[1],&fzNl[iChem][0],numAtomTot*sizeof(double));
	}
      }
    }
    //calcCoefForceExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    // force already reduce
    /*
    for(iAtom=1;iAtom<=numAtomTot;iAtom++){
      fx[iAtom] *= occNumber;
      fy[iAtom] *= occNumber;
      fz[iAtom] *= occNumber;
      //printf("fx %lg fy %lg fz %lg\n",fx[iAtom],fy[iAtom],fz[iAtom]);
    }
    */

/*--------------------------------------------------------------------------*/
/* iii) Reduce forces to the master proc                                    */

    // force already reduce in calcNlPseudoPostScf
    /*
    if(numProcStates>1){
      printf("iChem %i\n",iChem);
      Barrier(commStates);
      Reduce(&fx[1],&fxNl[iChem][0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
      Reduce(&fy[1],&fyNl[iChem][0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
      Reduce(&fz[1],&fzNl[iChem][0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
      Barrier(commStates);
      printf("111111 myid %i %lg %lg\n",myidState,fxNl[0][0],fx[1]);
    }
    else{
      memcpy(fxNl[iChem],&fx[1],numAtomTot*sizeof(double));
      memcpy(fyNl[iChem],&fy[1],numAtomTot*sizeof(double));
      memcpy(fzNl[iChem],&fz[1],numAtomTot*sizeof(double));
      printf("fx[1] %lg fy[1] %lg fz[1] %lg\n",fx[1],fy[1],fz[1]);
    }
    */

/*--------------------------------------------------------------------------*/
/* iv) Calculate the average values                                         */

    if(myidState==0){
      for(iAtom=1;iAtom<=numAtomTot;iAtom++){
        fx[iAtom] *= occNumber;
        fy[iAtom] *= occNumber;
        fz[iAtom] *= occNumber;
        //printf("fx %lg fy %lg fz %lg\n",fx[iAtom],fy[iAtom],fz[iAtom]);
      }
      memcpy(fxNl[iChem],&fx[1],numAtomTot*sizeof(double));
      memcpy(fyNl[iChem],&fy[1],numAtomTot*sizeof(double));
      memcpy(fzNl[iChem],&fz[1],numAtomTot*sizeof(double));      

      for(iAtom=0;iAtom<numAtomTot;iAtom++){
	fxNl[iChem][iAtom] /= numStateStoUp;
	fyNl[iChem][iAtom] /= numStateStoUp;
	fzNl[iChem][iAtom] /= numStateStoUp;
	//printf("fxNl %lg fyNl %lg fzNl %lg\n",fxNl[iChem][iAtom],fyNl[iChem][iAtom],fzNl[iChem][iAtom]);
      }
      //printf("iChem %i chemPot %lg K %lg NL %lg\n",iChem,chemPot[iChem],energyKineticTemp,energyNLTemp);
      //debug
      /*
      for(iAtom=0;iAtom<numAtomTot;iAtom++){
	fxNl[iChem][iAtom] *= numStateStoUp/(double)(occNumber);
        fyNl[iChem][iAtom] *= numStateStoUp/(double)(occNumber);
        fzNl[iChem][iAtom] *= numStateStoUp/(double)(occNumber);
      }
      */
    }
  }//endfor iChem

/*======================================================================*/
/* III) Interpolat nl pp and force                                      */


  if(myidState==0){
    if(chemPotOpt==1){
      //debug
      /*
      printf("chemPotTrue %lg\n",chemPotTrue);
      for(iChem=0;iChem<numChemPot;iChem++){
        printf("%lg %lg\n",chemPot[iChem],energyKNL[iChem]);
      }
      */
      // Force from non-local pp
      // Transpose first
      fxTemp = (double *)cmalloc(numChemPot*sizeof(double));
      fyTemp = (double *)cmalloc(numChemPot*sizeof(double));
      fzTemp = (double *)cmalloc(numChemPot*sizeof(double));
      for(iAtom=0;iAtom<numAtomTot;iAtom++){
	for(iChem=0;iChem<numChemPot;iChem++){
	  fxTemp[iChem] = fxNl[iChem][iAtom];
	  fyTemp[iChem] = fyNl[iChem][iAtom];
	  fzTemp[iChem] = fzNl[iChem][iAtom];
	}
	fxNlTrue[iAtom] = calcLagrangeInterpFun(numChemPot,chemPotTrue,chemPot,fxTemp,lagFunValue);
	fyNlTrue[iAtom] = calcLagrangeInterpFun(numChemPot,chemPotTrue,chemPot,fyTemp,lagFunValue);
	fzNlTrue[iAtom] = calcLagrangeInterpFun(numChemPot,chemPotTrue,chemPot,fzTemp,lagFunValue);
      }
      free(fxTemp);
      free(fyTemp);
      free(fzTemp);
    }
    if(chemPotOpt==2){
      for(iAtom=0;iAtom<numAtomTot;iAtom++){
	fxNlTrue[iAtom] = fxNl[0][iAtom];
	fyNlTrue[iAtom] = fyNl[0][iAtom];
	fzNlTrue[iAtom] = fzNl[0][iAtom];
	//printf("fxNlTrue %.16lg fyNlTrue %.16lg fzNlTrue %.16lg vnlFxCor %.16lg vnlFyCor %.16lg vnlFzCor %.16lg\n",
	//	fxNlTrue[iAtom],fyNlTrue[iAtom],fzNlTrue[iAtom],vnlFxCor[iAtom],vnlFyCor[iAtom],vnlFzCor[iAtom]);
      }//endfor iAtom
    }//endif chemPotOpt
    //Correct the non local force by fragment       
  }//endif myidState
 
/*======================================================================*/
/* IV) Add fragmentation correction and local contribution              */

  if(calcFragFlag==1&&myidState==0){
    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      // Add contribution from local pp
      /*
      fxNlTrue[iAtom] += fxBackup[iAtom];
      fyNlTrue[iAtom] += fyBackup[iAtom];
      fzNlTrue[iAtom] += fzBackup[iAtom];
      */
      // backup the uncorrected	version
      fxUnCor[iAtom] = fxNlTrue[iAtom];
      fyUnCor[iAtom] = fyNlTrue[iAtom];
      fzUnCor[iAtom] = fzNlTrue[iAtom];
      // fragmentation correction
      fxNlTrue[iAtom] += vnlFxCor[iAtom];
      fyNlTrue[iAtom] += vnlFyCor[iAtom];
      fzNlTrue[iAtom] += vnlFzCor[iAtom];
      //printf("NL cor force %.16lg %.16lg %.16lg\n",fxNlTrue[iAtom],fyNlTrue[iAtom],fzNlTrue[iAtom]);
    }
  }

  /*
  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    printf("NL cor force %.16lg %.16lg %.16lg\n",fxNlTrue[iAtom],fyNlTrue[iAtom],fzNlTrue[iAtom]);
  }
  */

/*======================================================================*/
/* V)   Ewald self and background terms                                */

  vself       = 0.0;
  vbgr        = 0.0;
  stat_avg->vrecip = vrecipLocal;

  for(iAtom=1;iAtom<=numAtomTot;iAtom++){
    fx[iAtom] = 0.0;
    fy[iAtom] = 0.0;
    fz[iAtom] = 0.0;
  }

  if(myidState==0&&iperd>0){
    ewald3d_selfbgr_cp(clatoms_info,ewald,ptens,vol,
                      &vself,&vbgr,iperd);
    printf("vrecip %.8lg vself %.16lg vbgr %.16lg\n",stat_avg->vrecip,vself,vbgr);
    stat_avg->vrecip += vself+vbgr;
    stat_avg->vintert = stat_avg->vrecip;
    stat_avg->vcoul = stat_avg->vrecip;
  }//endif myid_state

/*======================================================================*/
/* VI) Calculate real space nuclei-nuclei interaction	                */
  /*
  energy_control_inter_real(class,bonded,general_data);

  if(numProcStates==1){
    memcpy(&fxNuclei[0],&fx[1],numAtomTot*sizeof(double));
    memcpy(&fyNuclei[0],&fy[1],numAtomTot*sizeof(double));
    memcpy(&fzNuclei[0],&fz[1],numAtomTot*sizeof(double));
  }
  else{
    Reduce(&fx[1],&fxNuclei[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&fy[1],&fyNuclei[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&fz[1],&fzNuclei[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
  }
  */

  if(myidState==0){
    //debug
    /*
    for(iAtom=1;iAtom<=numAtomTot;iAtom++){
      fx[iAtom] = 0.0;
      fy[iAtom] = 0.0;
      fz[iAtom] = 0.0;
    }
    printf("stat_avg->vintert %lg\n",stat_avg->vintert);
    */
    printf("stat_avg->vintert %lg\n",stat_avg->vintert);
    energy_control_inter_real(class,bonded,general_data);

    //debug
    /*
    for(iAtom=1;iAtom<=numAtomTot;iAtom++){
      printf("aaaall nuclei %.16lg %.16lg %.16lg\n",fx[iAtom],fy[iAtom],fz[iAtom]);
    }
    */
    memcpy(&fxNuclei[0],&fx[1],numAtomTot*sizeof(double));
    memcpy(&fyNuclei[0],&fy[1],numAtomTot*sizeof(double));
    memcpy(&fzNuclei[0],&fz[1],numAtomTot*sizeof(double));

    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      fx[iAtom+1] = fxNlTrue[iAtom]+fxLoc[iAtom]+fxNuclei[iAtom];
      fy[iAtom+1] = fyNlTrue[iAtom]+fyLoc[iAtom]+fyNuclei[iAtom];
      fz[iAtom+1] = fzNlTrue[iAtom]+fzLoc[iAtom]+fzNuclei[iAtom];
      if(calcFragFlag==1){
        fxUnCor[iAtom] += fxLoc[iAtom]+fxNuclei[iAtom];
        fyUnCor[iAtom] += fyLoc[iAtom]+fyNuclei[iAtom];
        fzUnCor[iAtom] += fzLoc[iAtom]+fzNuclei[iAtom];
      }//endif calcFragFlag
    }//endfor iAtom
  }//endif myidState
  
/*======================================================================*/ 
/* VII) Output the energy Term                                          */


  if(myidState==0){
    energyTotElec = energyKe+energyPnl+energyHartTemp
		    +energyExtTemp+energyExcTemp;
    vInter = stat_avg->vintert;
    energyTot = energyTotElec+vInter;
    printf("==============================================\n");
    printf("Total Energy\n");
    printf("==============================================\n");
    printf("Electron Kinetic Energy:      %.20lg\n",energyKe);
    printf("Electron NLPP:                %.20lg\n",energyPnl);
    printf("Electron Hartree Energy:      %.20lg\n",energyHartTemp);
    printf("Electron Ext Energy:          %.20lg\n",energyExtTemp);
    printf("Electron Ex-Cor Energy:       %.20lg\n",energyExcTemp);
    printf("Electron Total Elec Energy:   %.20lg\n",energyTotElec);
    printf("Atom Energy:		  %.20lg\n",vInter);
    printf("Total Energy:		  %.20lg\n",energyTot);
    printf("==============================================\n");


    FILE *fileForce = fopen("atom-force","w");
    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      fprintf(fileForce,"atom %i cor %.16lg %.16lg %.16lg Uncor %.8lg %.8lg %.8lg loc %.8lg %.8lg %.8lg\n",
	     iAtom,fx[iAtom+1],fy[iAtom+1],fz[iAtom+1],
	     fxUnCor[iAtom],fyUnCor[iAtom],fzUnCor[iAtom],
             fxLoc[iAtom],fyLoc[iAtom],fzLoc[iAtom]);
    }
    fclose(fileForce);
    /*
    FILE *fileForceClass = fopen("atom-force-classic","w");
    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      fprintf(fileForceClass,"%i %.16lg %.16lg %.16lg\n",
              fxCl[iAtom+1],fyCl[iAtom+1],fzCl[iAtom+1]);
    }
    */
  }

/*======================================================================*/
/* VII) Free all temp vectors                                           */

  free(&fxNuclei[0]);  
  free(&fyNuclei[0]);
  free(&fzNuclei[0]);
  free(&fxLoc[0]);
  free(&fyLoc[0]);
  free(&fzLoc[0]);
  free(&fxUnCor[0]);
  free(&fyUnCor[0]);
  free(&fzUnCor[0]);

  free(&fxNlTrue[0]);
  free(&fyNlTrue[0]);
  free(&fzNlTrue[0]);
  for(iChem=0;iChem<numChemPot;iChem++){
    free(&fxNl[iChem][0]);
    free(&fyNl[iChem][0]);
    free(&fzNl[iChem][0]);
  }
  free(&fxNl[0]);
  free(&fyNl[0]);
  free(&fzNl[0]);
  free(&lagFunValue[0]);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/





