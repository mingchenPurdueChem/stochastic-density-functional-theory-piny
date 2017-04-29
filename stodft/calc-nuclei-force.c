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
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  COMMUNICATE *communicate      = &(cp->communicate);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  CPEWALD *cpewald              = &(cp->cpewald);
  CELL *cell                    = &(general_data->cell);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  EWALD *ewald			= &(general_data->ewald);
  PTENS *ptens			= &(general_data->ptens);
  FRAGINFO *fragInfo            = stodftInfo->fragInfo;
  
  int cpLsda         = cpopts->cp_lsda;
  int numStateStoUp  = stodftInfo->numStateStoUp;
  int atomForceFlag  = stodftInfo->atomForceFlag;
  int chemPotOpt     = stodftInfo->chemPotOpt;
  int calcFragFlag   = stodftInfo->calcFragFlag;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int iperd	      = cell->iperd;
  int numChemPot = stodftInfo->numChemPot;
  int occNumber = stodftInfo->occNumber;
  int myidState         = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int numAtomTot = clatoms_info->natm_tot;
  int iState,iCoeff,iChem,iAtom;
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
  double energyTotElec,energyTot;
  double energyExtTemp,energyExcTemp,energyHartTemp;
  double vInter;
  double vrecip;
  double vself,vbgr;

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
  double *ak2_sm  =  cpewald->ak2_sm;
  double *chemPot = stodftCoefPos->chemPot;
  double *fx = clatoms_pos->fx;
  double *fy = clatoms_pos->fy;
  double *fz = clatoms_pos->fz;
  double *vnlFxCor = fragInfo->vnlFxCor;
  double *vnlFyCor = fragInfo->vnlFyCor;
  double *vnlFzCor = fragInfo->vnlFzCor;
  double *lagFunValue = (double*)cmalloc(numChemPot*sizeof(double));

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
  double *fxBackup,*fyBackup,*fzBackup;
  double *fxNlTrue,*fyNlTrue,*fzNlTrue;
  double *fxUnCor,*fyUnCor,*fzUnCor;

  MPI_Comm commStates = communicate->comm_states;

/*======================================================================*/
/* I) Calculate Local pp	                                        */

  fxBackup = (double *)cmalloc(numAtomTot*sizeof(double));
  fyBackup = (double *)cmalloc(numAtomTot*sizeof(double));
  fzBackup = (double *)cmalloc(numAtomTot*sizeof(double));
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
    //fscanf(fp_rhok,"%lg",&(cp->cpscr.cpscr_rho.rhocr_up[iCoeff]));
    //fscanf(fp_rhok,"%lg",&(cp->cpscr.cpscr_rho.rhoci_up[iCoeff]));
    //printf("rho k %lg %lg\n",cp->cpscr.cpscr_rho.rhocr_up[1],cp->cpscr.cpscr_rho.rhoci_up[1]);
  }
  fclose(fp_rhok);
  */

  calcLocExtPostScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

  if(numProcStates==1){
    memcpy(&fxBackup[0],&fx[1],numAtomTot*sizeof(double));
    memcpy(&fyBackup[0],&fy[1],numAtomTot*sizeof(double));
    memcpy(&fzBackup[0],&fz[1],numAtomTot*sizeof(double));
  }
  else{
    Reduce(&fx[1],&fxBackup[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&fy[1],&fyBackup[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&fz[1],&fzBackup[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&(stat_avg->vrecip),&vrecip,1,MPI_DOUBLE,MPI_SUM,0,commStates);
  }

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

    //pp 
    for(iAtom=1;iAtom<=numAtomTot;iAtom++){
      fx[iAtom] = 0.0;
      fy[iAtom] = 0.0;
      fz[iAtom] = 0.0;
    }

/*--------------------------------------------------------------------------*/
/* ii) Calculate nl pp force						    */

    //calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    calcNlPseudoPostScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    //calcCoefForceExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    for(iAtom=1;iAtom<=numAtomTot;iAtom++){
      fx[iAtom] *= occNumber;
      fy[iAtom] *= occNumber;
      fz[iAtom] *= occNumber;
    }

/*--------------------------------------------------------------------------*/
/* iii) Reduce forces to the master proc                                    */

    if(numProcStates>1){
      if(atomForceFlag==1){
        Reduce(fxNl[iChem],&fx[1],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
        Reduce(fyNl[iChem],&fy[1],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
        Reduce(fzNl[iChem],&fz[1],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
      }
    }
    else{
      memcpy(fxNl[iChem],&fx[1],numAtomTot*sizeof(double));
      memcpy(fyNl[iChem],&fy[1],numAtomTot*sizeof(double));
      memcpy(fzNl[iChem],&fz[1],numAtomTot*sizeof(double));
      printf("fx[1] %lg fy[1] %lg fz[1] %lg\n",fx[1],fy[1],fz[1]);
    }

/*--------------------------------------------------------------------------*/
/* iv) Calculate the average values                                         */

    if(myidState==0){
      for(iAtom=0;iAtom<numAtomTot;iAtom++){
	fxNl[iChem][iAtom] /= numStateStoUp;
	fyNl[iChem][iAtom] /= numStateStoUp;
	fzNl[iChem][iAtom] /= numStateStoUp;
      }
      //printf("iChem %i chemPot %lg K %lg NL %lg\n",iChem,chemPot[iChem],energyKineticTemp,energyNLTemp);
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
	printf("fxNlTrue %lg fyNlTrue %lg fzNlTrue %lg\n",fxNlTrue[iAtom],fyNlTrue[iAtom],fzNlTrue[iAtom]);
      }//endfor iAtom
    }//endif chemPotOpt
    //Correct the non local force by fragment       
  }//endif myidState
 
/*======================================================================*/
/* IV) Add fragmentation correction and local contribution              */

  if(calcFragFlag==1&&myidState==0){
    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      // Add contribution from local pp
      fxNlTrue[iAtom] += fxBackup[iAtom];
      fyNlTrue[iAtom] += fyBackup[iAtom];
      fzNlTrue[iAtom] += fzBackup[iAtom];
      // backup the uncorrected	version
      fxUnCor[iAtom] = fxNlTrue[iAtom];
      fyUnCor[iAtom] = fyNlTrue[iAtom];
      fzUnCor[iAtom] = fzNlTrue[iAtom];
      // fragmentation correction
      fxNlTrue[iAtom] += vnlFxCor[iAtom];
      fyNlTrue[iAtom] += vnlFyCor[iAtom];
      fzNlTrue[iAtom] += vnlFzCor[iAtom];
    }
  }

/*======================================================================*/
/* V)   Ewald self and background terms                                */

  vself       = 0.0;
  vbgr        = 0.0;

  for(iAtom=1;iAtom<=numAtomTot;iAtom++){
    fx[iAtom] = 0.0;
    fy[iAtom] = 0.0;
    fz[iAtom] = 0.0;
  }

  if(myidState==0&&iperd>0){
    ewald3d_selfbgr_cp(clatoms_info,ewald,ptens,vol,
                      &vself,&vbgr,iperd);
    stat_avg->vrecip += vself+vbgr;
    stat_avg->vintert = stat_avg->vrecip;
    stat_avg->vcoul = stat_avg->vrecip;
  }//endif myid_state

/*======================================================================*/
/* VI) Calculate real space nuclei-nuclei interaction	                */
  class->energy_ctrl.iget_full_inter = 1;
  class->energy_ctrl.iget_res_inter = 0;
  energy_control_inter_real(class,bonded,general_data);

  if(numProcStates==1){
    memcpy(&fxBackup[0],&fx[1],numAtomTot*sizeof(double));
    memcpy(&fyBackup[0],&fy[1],numAtomTot*sizeof(double));
    memcpy(&fzBackup[0],&fz[1],numAtomTot*sizeof(double));
  }
  else{
    Reduce(&fx[1],&fxBackup[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&fy[1],&fyBackup[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&fz[1],&fzBackup[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
  }

  if(myidState==0){
    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      fx[iAtom+1] = fxNlTrue[iAtom]+fxBackup[iAtom];
      fy[iAtom+1] = fyNlTrue[iAtom]+fyBackup[iAtom];
      fz[iAtom+1] = fzNlTrue[iAtom]+fzBackup[iAtom];
      if(calcFragFlag==1){
        fxUnCor[iAtom] += fxBackup[iAtom];
        fyUnCor[iAtom] += fyBackup[iAtom];
        fzUnCor[iAtom] += fzBackup[iAtom];
      }//endif calcFragFlag
    }//endfor iAtom
  }//endif myidState
  
/*======================================================================*/ 
/* VII) Output the energy Term                                          */


  if(myidState==0){
    energyHartTemp = stat_avg->cp_ehart;
    energyExtTemp = stat_avg->cp_eext;
    energyExcTemp = stat_avg->cp_exc;

    energyTotElec = energyKe+energyPnl+energyHartTemp
		    +energyExtTemp+energyExcTemp;
    vInter = stat_avg->vintert;
    energyTot = energyTotElec+vInter;
    printf("==============================================\n");
    printf("Total Energy\n");
    printf("==============================================\n");
    printf("Electron Kinetic Energy:      %.20lg\n",stat_avg->cp_eke);
    printf("Electron NLPP:                %.20lg\n",energyPnl);
    printf("Electron Hartree Energy:      %.20lg\n",energyHartTemp);
    printf("Electron Ext Energy:          %.20lg\n",energyExtTemp);
    printf("Electron Ex-Cor Energy:       %.20lg\n",energyExcTemp);
    printf("Electron Total Elec Energy:   %.20lg\n",energyTotElec);
    printf("Atom Energy:		  %.20lg\n",vInter);
    printf("Total Energy:		  %.20lg\n",energyTot);
    printf("==============================================\n");

    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      printf("atom %i uncor %.8lg %.8lg %.8lg cor %.8lg %.8lg %.8lg\n",
	     iAtom,fx[iAtom+1],fy[iAtom+1],fz[iAtom+1],
	     fxUnCor[iAtom],fyUnCor[iAtom],fzUnCor[iAtom]);
    }
  }

/*======================================================================*/
/* VII) Free all temp vectors                                           */

  free(&fxBackup[0]);  
  free(&fyBackup[0]);
  free(&fzBackup[0]);
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





