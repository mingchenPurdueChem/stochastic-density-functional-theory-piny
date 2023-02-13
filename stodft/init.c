/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: init.c                                         */
/*                                                                          */
/* This routine initialize the stochastic dft calculation.                  */
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
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_stodft_local.h"
#include "../proto_defs/proto_frag_entry.h"

#define TIME_CP_OFF
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void commStodft(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate expanCoeff for different chemical potentials                */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  COMMUNICATE   *communicate      = &(cp->communicate);
  STODFTINFO    *stodftInfo;
  STODFTCOEFPOS *stodftCoefPos;

  int numProcStates             = communicate->np_states;
  int myidState                 = communicate->myid_state;
  int myid			= communicate->myid;
  MPI_Comm world		= communicate->world;
  
  if(myid==0){
    stodftInfo = cp->stodftInfo;
    stodftCoefPos = cp->stodftCoefPos;
  }
  else{
    cp->stodftInfo = (STODFTINFO*)cmalloc(sizeof(STODFTINFO));
    cp->stodftCoefPos = (STODFTCOEFPOS*)cmalloc(sizeof(STODFTCOEFPOS));
    stodftInfo = cp->stodftInfo;
    stodftCoefPos = cp->stodftCoefPos;
    stodftInfo->metallic = (METALLIC*)cmalloc(sizeof(METALLIC));
    stodftInfo->fragInfo = (FRAGINFO*)cmalloc(sizeof(FRAGINFO));
  }
 
  if(numProcStates>1){ 
    Bcast(&(stodftInfo->stodftOn),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->missionType),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->numScf),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->expanType),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->filterFunType),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->numOrbital),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->numChemPot),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->energyWindowOn),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->printChebyMoment),1,MPI_INT,0,world);
    //Bcast(&(stodftInfo->readCoeffFlag),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->numStateStoUp),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->numStateStoDn),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->chemPotOpt),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->filterDiagFlag),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->densityMixFlag),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->numDiis),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->numStepMix),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->checkpointWriteFreq),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->checkpointParFlag),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->numStatePrintUp),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->numStatePrintDn),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->smearOpt),1,MPI_INT,0,world);

    //frag
    Bcast(&(stodftInfo->calcFragFlag),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->fragOpt),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->fragCellOpt),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->fragDFTMethod),1,MPI_INT,0,world);

    Bcast(&(stodftInfo->fitErrTol),1,MPI_DOUBLE,0,world);
    Bcast(&(stodftInfo->beta),1,MPI_DOUBLE,0,world);
    Bcast(&(stodftInfo->numElecTrue),1,MPI_DOUBLE,0,world);
    Bcast(&(stodftInfo->numElecTrueUp),1,MPI_DOUBLE,0,world);
    Bcast(&(stodftInfo->numElecTrueDn),1,MPI_DOUBLE,0,world);
    Bcast(&(stodftInfo->numElecTrue),1,MPI_DOUBLE,0,world);
    Bcast(&(stodftInfo->chemPotInit),1,MPI_DOUBLE,0,world);
    Bcast(&(stodftInfo->gapInit),1,MPI_DOUBLE,0,world);
    Bcast(&(stodftInfo->mixRatioSM),1,MPI_DOUBLE,0,world);
    Bcast(&(stodftInfo->mixRatioSM2),1,MPI_DOUBLE,0,world);
    Bcast(&(stodftInfo->energyTol),1,MPI_DOUBLE,0,world);
    Bcast(&(stodftInfo->smearTemperature),1,MPI_DOUBLE,0,world);

    Bcast(stodftInfo->densityFileName,MAXWORD,MPI_CHAR,0,world);
    Bcast(stodftInfo->densityFinalFileName,MAXWORD,MPI_CHAR,0,world);
    Bcast(stodftInfo->densityReadFileName,MAXWORD,MPI_CHAR,0,world);

    Bcast(&(stodftInfo->metallic->electronFricFlag),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->metallic->numAtomFric),1,MPI_INT,0,world);
    Bcast(&(stodftInfo->metallic->sigma),1,MPI_DOUBLE,0,world);
    Bcast(&(stodftInfo->calcLocalTraceOpt),1,MPI_INT,0,world);
  }

}/*end Routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initStodft(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp,
		int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate expanCoeff for different chemical potentials                */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS  *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  ATOMMAPS     *atommaps     = &(class->atommaps);
  CELL         *cell         = &(general_data->cell);
  FOR_SCR      *for_scr      = &(class->for_scr);
  EWD_SCR      *ewd_scr      = &(class->ewd_scr);
  PTENS        *ptens        = &(general_data->ptens);

  CPOPTS        *cpopts           = &(cp->cpopts);
  PSEUDO        *pseudo           = &(cp->pseudo);
  CPCOEFFS_INFO *cpcoeffs_info    = &(cp->cpcoeffs_info);
  COMMUNICATE   *communicate      = &(cp->communicate);
  CPCOEFFS_POS  *cpcoeffs_pos     = &(cp->cpcoeffs_pos[ip_now]);
  STODFTINFO    *stodftInfo       = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos    = cp->stodftCoefPos;
  CPSCR         *cpscr            = &(cp->cpscr);  
  NEWTONINFO    *newtonInfo;
  CHEBYSHEVINFO *chebyshevInfo;
  FRAGINFO	*fragInfo;
  //PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d; // use lg or sparse
  
  int iperd          = cell->iperd;
  int cpLsda         = cpopts->cp_lsda;
  int cpGga          = cpopts->cp_gga;
  int cpParaOpt      = cpopts->cp_para_opt;
  int realSparseOpt  = cpopts->realSparseOpt;
  int expanType      = stodftInfo->expanType;
  int numOrbital     = stodftInfo->numOrbital;
  int polynormLength = stodftInfo->polynormLength;
  int energyWindowOn = stodftInfo->energyWindowOn;
  int numChemPot     = stodftInfo->numChemPot;
  int readCoeffFlag  = cpopts->readCoeffFlag;
  int densityMixFlag = stodftInfo->densityMixFlag;
  int numDiis	     = stodftInfo->numDiis;
  int numStepMix     = stodftInfo->numStepMix;
  int numScf	     = stodftInfo->numScf;
  int numElecTrue    = stodftInfo->numElecTrue;
  int chemPotOpt     = stodftInfo->chemPotOpt;
  int filterDiagFlag = stodftInfo->filterDiagFlag;
  int calcFragFlag   = stodftInfo->calcFragFlag;
  int fragOpt	     = stodftInfo->fragOpt;
  int fragCellOpt    = stodftInfo->fragCellOpt;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numStateUpTot  = numStateUpProc*numCoeff;
  int numStateDnTot  = numStateDnProc*numCoeff;
  int totalPoly	     = polynormLength*numChemPot;
  int filterFunType   = stodftInfo->filterFunType;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int checkPerdSize             = cpopts->icheck_perd_size;
  int checkDualSize             = cpopts->icheck_dual_size;
  int coefFormUp                = cpcoeffs_pos->icoef_form_up;
  int coefOrthUp                = cpcoeffs_pos->icoef_orth_up;
  int forceCoefFormUp           = cpcoeffs_pos->ifcoef_form_up;
  int forceCoefOrthUp           = cpcoeffs_pos->ifcoef_orth_up;
  int coefFormDn                = cpcoeffs_pos->icoef_form_dn;
  int coefOrthDn                = cpcoeffs_pos->icoef_orth_dn;
  int forceCoefFormDn           = cpcoeffs_pos->ifcoef_form_dn;
  int forceCoefOrthDn           = cpcoeffs_pos->ifcoef_orth_dn;
  int numProcStates             = communicate->np_states;
  int myidState                 = communicate->myid_state;
  //int numFFTProc        = cp_para_fft_pkg3d_lg->nfft_proc;
  //int numFFT            = cp_para_fft_pkg3d_lg->nfft;
  //int numFFT2           = numFFT/2;
  //int numFFT2Proc       = numFFTProc/2;
  int numFFTProc,numFFT,numFFT2,numFFT2Proc;
  int iChem,iSamp,iCell,iProc,iCoeff,iMol,iDiis;
  int div,res;
  int count,numChemProc,rhoRealGridNum,rhoRealGridTot;
  int numChemProcMalloc;
  int numSendNoise;
  // frag
  int numMolTot;
  int numMolType	= atommaps->nmol_typ;
  int numFragTot;
  int numFragProc;
  int iFrag;

  MPI_Comm comm_states   =    communicate->comm_states;

  int *pcoefFormUp                   = &(cpcoeffs_pos->icoef_form_up);
  int *pcoefOrthUp                   = &(cpcoeffs_pos->icoef_orth_up);
  int *pforceCoefFormUp              = &(cpcoeffs_pos->ifcoef_form_up);
  int *pforceCoefOrthUp              = &(cpcoeffs_pos->ifcoef_orth_up);
  int *pcoefFormDn                   = &(cpcoeffs_pos->icoef_form_dn);
  int *pcoefOrthDn                   = &(cpcoeffs_pos->icoef_orth_dn);
  int *pforceCoefFormDn              = &(cpcoeffs_pos->icoef_form_dn);
  int *pforceCoefOrthDn              = &(cpcoeffs_pos->icoef_orth_dn);
  int *rhoRealSendCounts;
  int *rhoRealDispls;
  int *noiseSendCounts;
  int *noiseDispls;
  //frag
  int *numMolJmolType		     = atommaps->nmol_jmol_typ;

  char *ggaxTyp     = pseudo->ggax_typ;
  char *ggacTyp     = pseudo->ggac_typ;
  
  double Smin = -2.0;
  double Smax = 2.0;
  double energyDiff;
  double tolEdgeDist            = cpopts->tol_edge_dist;

  double *coeffReUp        = cpcoeffs_pos->cre_up;
  double *coeffImUp        = cpcoeffs_pos->cim_up;
  double *forceCoeffReUp   = cpcoeffs_pos->fcre_up;
  double *forceCoeffImUp   = cpcoeffs_pos->fcre_up;
  double *coeffReDn        = cpcoeffs_pos->cre_dn;
  double *coeffImDn        = cpcoeffs_pos->cim_dn;
  double *cpScrCoeffReUp   = cpscr->cpscr_wave.cre_up;
  double *cpScrCoeffImUp   = cpscr->cpscr_wave.cim_up;
  double *cpScrCoeffReDn   = cpscr->cpscr_wave.cre_dn;
  double *cpScrCoeffImDn   = cpscr->cpscr_wave.cim_dn;
  double *ptensPvtenTmp    = ptens->pvten_tmp;

  FILE *densityFile;
  
  if(numProcStates>1)Barrier(comm_states);

/*==========================================================================*/
/* I) General parameters and malloc					    */

  if(realSparseOpt==0){
    cp_para_fft_pkg3d = &(cp->cp_para_fft_pkg3d_lg);
  }
  else{
    cp_para_fft_pkg3d = &(cp->cp_para_fft_pkg3d_sparse);
  }
  numFFTProc = cp_para_fft_pkg3d->nfft_proc;
  numFFT = cp_para_fft_pkg3d->nfft;
  numFFT2 = numFFT/2;
  numFFT2Proc = numFFTProc/2;
  
  stodftInfo->vpsAtomListFlag = 0;
  stodftInfo->filterFlag = 0;
  stodftInfo->numThreads = communicate->numThreads;
  //stodftInfo->fragInfo->iFrag = 0; //TODO sagar
  // Chebyshev way to calculate chem pot (if we do not use energy window)
  if(chemPotOpt==2&&energyWindowOn==0)stodftInfo->numChemPot = 1;
  
  // The following defines the fragment+energy window
  // In this combination, we have to generate energy windows and 
  // another window covering the unoccupied space w.r.t. fragment density. 
  // After that, we will use the regular energy widows. Some functions will 
  // be involved in both cases therefore we need to use this flag to control 
  // it. 
  if(calcFragFlag==1&&energyWindowOn==1){
    stodftInfo->fragWindowFlag = 1;
    // numChemPot+1 at reading paramter file step
    //stodftInfo->numChemPot += 1;
    //numChemPot += 1;
  } 

  stodftInfo->orbRealPrintFlag = 0;

  stodftInfo->numElecSys = stodftInfo->numElecTrue;

  stodftCoefPos->chemPot = (double*)cmalloc(numChemPot*sizeof(double));
  stodftCoefPos->chemPotBackUp = (double*)cmalloc(numChemPot*sizeof(double));
  stodftInfo->energyKe = (double*)cmalloc(numChemPot*sizeof(double));
  stodftInfo->energyPNL = (double*)cmalloc(numChemPot*sizeof(double));
  stodftCoefPos->testWfMaxRe = (double*)cmalloc(numCoeff*sizeof(double));
  stodftCoefPos->testWfMaxIm = (double*)cmalloc(numCoeff*sizeof(double));
  stodftCoefPos->testWfMinRe = (double*)cmalloc(numCoeff*sizeof(double));
  stodftCoefPos->testWfMinIm = (double*)cmalloc(numCoeff*sizeof(double));

  /*
  if(expanType==2&&filterFunType==1){
    if(energyWindowOn==0){
      stodftInfo->fermiFunctionReal = &fermiExpReal;
    }
    else{
      stodftInfo->fermiFunctionReal = &fermiExpReal;
      stodftInfo->fermiFunctionLongDouble = &fermiExpLongDouble;
    }
  }
  if(expanType==2&&filterFunType==2)stodftInfo->fermiFunctionReal = &fermiErfcReal;
  if(expanType==2&&filterFunType==3)stodftInfo->fermiFunctionReal = &gaussianReal;
  if(expanType==3&&filterFunType==1)stodftInfo->fermiFunctionComplex = &fermiExpComplex;
  */
  switch(expanType){
    case 1:
      switch(filterFunType){
        case 1:
          if(energyWindowOn==0){
            stodftInfo->fermiFunctionReal = &fermiExpReal;
          }
          else{
            stodftInfo->fermiFunctionReal = &fermiExpReal;
            stodftInfo->fermiFunctionLongDouble = &fermiExpLongDouble;
          }
          break;
        case 2:
          stodftInfo->fermiFunctionReal = &fermiErfcReal;
          break;
        case 3:
          stodftInfo->fermiFunctionReal = &gaussianReal;
          break;
        case 4:
          stodftInfo->fermiFunctionReal = &entropyReal;
          break;
        default:
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          printf("Internal Error! Bad filter type!\n");
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          exit(1);
      }  
      break;
    case 2:
      switch(filterFunType){
        case 1:
          if(energyWindowOn==0){
            stodftInfo->fermiFunctionReal = &fermiExpReal;
          }
          else{
            stodftInfo->fermiFunctionReal = &fermiExpReal;
            stodftInfo->fermiFunctionLongDouble = &fermiExpLongDouble;
          }
          break;
        case 2:
          stodftInfo->fermiFunctionReal = &fermiErfcReal;
          break;
        case 3:
          stodftInfo->fermiFunctionReal = &gaussianReal;
          break;
        case 4:
          stodftInfo->fermiFunctionReal = &entropyReal;
          break;
        default:
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          printf("Internal Error! Bad filter type!\n");
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          exit(1);
      }
      break;
    case 3:
      switch(filterFunType){
        case 1:
          stodftInfo->fermiFunctionComplex = &fermiExpComplex;
          break;
        default:
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n"); 
          printf("Only support Fermi function for Non-Hermitian Hamiltonian!\n");
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          exit(1);
      }
      break;
    case 4:
      switch(filterFunType){
        case 1:
          stodftInfo->fermiFunctionReal = &fermiExpReal;
          break;
        case 2:
          stodftInfo->fermiFunctionReal = &fermiErfcReal;
          break;
        case 3:
          stodftInfo->fermiFunctionReal = &gaussianReal;
          break;
        case 4:
          stodftInfo->fermiFunctionReal = &entropyReal;
          break;
        default:
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          printf("Internal Error! Bad filter type!\n");
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          exit(1);
      }  
      break;
    default:
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Bad expansion type!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }

  if(myidState==0){
    if(expanType==3&&filterFunType==2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("We haven't implement erfc type of Fermi function \n");
      printf("for non-Hermitian KS Hamiltonian. Please use the \n");
      printf("exponential type in this case!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }
    if(expanType==3&&filterFunType==3){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("We haven't implement gaussian type of Fermi function \n");
      printf("for non-Hermitian KS Hamiltonian. Please use the \n");
      printf("exponential type in this case!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }
    if(energyWindowOn==1&&filterFunType==3){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("sDFT with energy window does not support gaussian\n");
      printf("type window. Please use the exponential type in this case\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);

    }
  }

/*==========================================================================*/
/* II) Malloc by expension type						    */


  stodftInfo->newtonInfo = (NEWTONINFO *)cmalloc(sizeof(NEWTONINFO));
  stodftInfo->chebyshevInfo = (CHEBYSHEVINFO *)cmalloc(sizeof(CHEBYSHEVINFO));

  if(expanType==2){
    //stodftCoefPos->expanCoeff = (double *)cmalloc(totalPoly*sizeof(double));
    newtonInfo = stodftInfo->newtonInfo;
    //newtonInfo->sampPoint = (double *)cmalloc(polynormLength*sizeof(double));
    //newtonInfo->sampPointUnscale = (double *)cmalloc(polynormLength*sizeof(double));
    newtonInfo->Smin = Smin;
    newtonInfo->Smax = Smax;
    stodftCoefPos->expanCoeff = NULL;
    newtonInfo->sampPoint = NULL;
    newtonInfo->sampPointUnscale = NULL;
    //newtonInfo->scale = (Smax-Smin)/energyDiff;      
  }
  if(expanType==1||chemPotOpt==2){
    chebyshevInfo = stodftInfo->chebyshevInfo;
    chebyshevInfo->Smin = -1.0;
    chebyshevInfo->Smax = 1.0;
  }

  //For debug only
  //stodftCoefPos->chemPot[0] = 0.4990160113690864;
  //stodftCoefPos->chemPot[0] = -0.17435045;
  //stodftCoefPos->chemPot[0] = 0.5045818049407941;
  //stodftCoefPos->chemPot[1] = 0.5045818049407941;
 
/*==========================================================================*/
/* III) Initialize utility data						    */
  /*
  FILE *fileSampPoint = fopen("samp-point","r");
  double *sampLocal = (double*)(newtonInfo->sampPoint);

  switch(expanType){
    case 2:
      
      for(iSamp=0;iSamp<polynormLength;iSamp++){
	fscanf(fileSampPoint,"%lg",&(sampLocal[iSamp]));
	//printf("samp %lg\n",sampLocal[iSamp]);
      }
      
      //genSampNewtonHermit(stodftInfo,stodftCoefPos);
      break;
    // I'll do chebyshev and non-Hermitain Newtonian later    
  }
  fclose(fileSampPoint);
  */
/*==========================================================================*/
/* IV) Initialize Flags							    */


  general_data->stat_avg.count_diag_srot      = 0.0;
  general_data->stat_avg.fatm_mag = 10000.0;
  general_data->stat_avg.fatm_max = 10000.0;

  cpopts->cp_becke=0;
  cpopts->cp_pw91x=0;
  cpopts->cp_fila_1x=0;
  cpopts->cp_fila_2x=0;
  cpopts->cp_pbe_x=0;
  cpopts->cp_revpbe_x=0;
  cpopts->cp_rpbe_x=0;
  cpopts->cp_xpbe_x=0;
  cpopts->cp_brx89=0;
  cpopts->cp_brx2k=0;
  cpopts->cp_lyp=0;
  cpopts->cp_lypm1=0;
  cpopts->cp_pw91c=0;
  cpopts->cp_pbe_c=0;
  cpopts->cp_xpbe_c=0;
  cpopts->cp_tau1_c=0;
  cpopts->cp_debug_xc=0;
  if(cpGga==1){
    if(strcasecmp(ggaxTyp,"becke"   )==0){cpopts->cp_becke=1;}
    if(strcasecmp(ggaxTyp,"pw91x"   )==0){cpopts->cp_pw91x=1;}
    if(strcasecmp(ggaxTyp,"fila_1x" )==0){cpopts->cp_fila_1x=1;}
    if(strcasecmp(ggaxTyp,"fila_2x" )==0){cpopts->cp_fila_2x=1;}
    if(strcasecmp(ggaxTyp,"pbe_x"   )==0){cpopts->cp_pbe_x=1;}
    if(strcasecmp(ggaxTyp,"revpbe_x")==0){cpopts->cp_revpbe_x=1;}
    if(strcasecmp(ggaxTyp,"rpbe_x"  )==0){cpopts->cp_rpbe_x=1;}
    if(strcasecmp(ggaxTyp,"xpbe_x"  )==0){cpopts->cp_xpbe_x=1;}
    if(strcasecmp(ggaxTyp,"brx89"   )==0){cpopts->cp_brx89=1;}
    if(strcasecmp(ggaxTyp,"brx2k"   )==0){cpopts->cp_brx2k=1;}
    if(strcasecmp(ggacTyp,"lyp"     )==0){cpopts->cp_lyp=1;  }
    if(strcasecmp(ggacTyp,"lypm1"   )==0){cpopts->cp_lypm1=1;  }
    if(strcasecmp(ggacTyp,"pw91c"   )==0){cpopts->cp_pw91c=1;}
    if(strcasecmp(ggacTyp,"pbe_c"   )==0){cpopts->cp_pbe_c=1;}
    if(strcasecmp(ggacTyp,"xpbe_c"  )==0){cpopts->cp_xpbe_c=1;}
    if(strcasecmp(ggacTyp,"tau1_c"  )==0){cpopts->cp_tau1_c=1;}
    if(strcasecmp(ggacTyp,"debug97x")==0){cpopts->cp_debug_xc=1;}
  }/*endif*/

  stodftInfo->readCoeffFlag = readCoeffFlag;
  if(readCoeffFlag==1||readCoeffFlag<0)stodftInfo->reInitFlag = 0;
  else stodftInfo->reInitFlag = 1;

  stodftInfo->energyElecTot = 0.0;
  stodftInfo->energyElecTotOld = 0.0;

  //printf("readCoeffFlag %i reInitFlag %i\n",readCoeffFlag,stodftInfo->reInitFlag);

/*==========================================================================*/
/* V) Calculate the non-local pseudopotential list                          */
  /*
  if(stodftInfo->vpsAtomListFlag==0||cpDualGridOptOn>= 1){
    control_vps_atm_list(pseudo,cell,clatoms_pos,clatoms_info,
                         atommaps,ewd_scr,for_scr,cpDualGridOptOn,
                         stodftInfo->vpsAtomListFlag);
    stodftInfo->vpsAtomListFlag = 1;
  }
  */

/*==========================================================================*/
/* VI) Initialize noise orbital scattering	                            */



  stodftInfo->randSeedTot = (double*)cmalloc(numProcStates*sizeof(double));

  /*
  stodftInfo->noiseSendCounts = (int*)cmalloc(numProcStates*sizeof(int));
  stodftInfo->noiseDispls     = (int*)cmalloc(numProcStates*sizeof(int));
  noiseSendCounts = stodftInfo->noiseSendCounts;
  noiseDispls = stodftInfo->noiseDispls;

  if(numProcStates>1){
    if(myidState==0){
      if(cpLsda==0)noiseSendCounts[0] = numStateUpTot*2;
      else noiseSendCounts[0] = (numStateUpTot+numStateDnTot)*2;
      for(iProc=1;iProc<numProcStates;iProc++){
	Recv(&noiseSendCounts[iProc],1,MPI_INT,iProc,iProc,comm_states);
      }//endfor iProc
    }//endif
    else{
      if(cpLsda==0){
	numSendNoise = numStateUpTot*2;
	Send(&numStateUpTot,1,MPI_INT,0,myidState,comm_states);
      }
      else{
	numSendNoise = (numStateUpTot+numStateDnTot)*2;
	Send(&numSendNoise,1,MPI_INT,0,myidState,comm_states);
      } 
    }
    Barrier(comm_states);
    Bcast(noiseSendCounts,numProcStates,MPI_INT,0,comm_states);
    Barrier(comm_states);
    printf("0 %i 1 %i 2 %i 3 %i\n",noiseSendCounts[0],noiseSendCounts[1],noiseSendCounts[2],noiseSendCounts[3]);
    noiseDispls[0] = 0;
    for(iProc=1;iProc<numProcStates;iProc++){
      noiseDispls[iProc] = noiseDispls[iProc-1]+noiseSendCounts[iProc-1];
    }
    stodftInfo->numRandTot = 0;
    for(iProc=0;iProc<numProcStates;iProc++){
      stodftInfo->numRandTot += noiseSendCounts[iProc];
    }
  }//endif 
  else{
    stodftInfo->numRandTot = numStateUpTot*2;
    if(cpLsda==1)stodftInfo->numRandTot += numStateDnTot*2;
  }
  */
  


/*==========================================================================*/
/* V) Initialize for the density calculation				    */

  // I need to do this for both deterministic/stochastic density calculation since
  // I use similiar functions. I don't need this if I read in density or generate 
  // density from fragmentation.


  //printf("coefFormUp %i forceCoefFormUp %i\n",coefFormUp,forceCoefFormUp);
  if(readCoeffFlag>=0){
    if(numProcStates>1){
      if((coefFormUp+forceCoefFormUp)!=2){
	printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	printf("Up CP vectors are not in transposed form \n");
	printf("on state processor %d in min_STD_cp \n",myidState);
	printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	fflush(stdout);
	exit(1);
      }/*endif*/
      if(cpLsda==1){
	if((coefFormDn+forceCoefFormDn)!=2){
	  printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	  printf("Up CP vectors are not in transposed form \n");
	  printf("on state processor %d in min_STD_cp \n",myidState);
	  printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	  fflush(stdout);
	  exit(1);
	}/*endif*/
      }/*endif*/
    }/*endif*/
    if(numProcStates>1)Barrier(comm_states);

    if((iperd<3||iperd==4)&&checkPerdSize==1){
      cp_boundary_check(cell,clatoms_info,clatoms_pos,tolEdgeDist);
    }/*endif*/
    if(cpDualGridOptOn>=1&&checkDualSize==1){
      cp_dual_check(cell,clatoms_info,clatoms_pos,
		    atommaps->cp_atm_lst,tolEdgeDist);
    }/*endif*/

    if(numProcStates>1){
      cp_transpose_bck(coeffReUp,coeffImUp,pcoefFormUp,
		      cpScrCoeffReUp,cpScrCoeffImUp,&(cp->cp_comm_state_pkg_up));
      if(cpLsda==1&&numStateDnProc>0){
	cp_transpose_bck(coeffReDn,coeffImDn,pcoefFormDn,
		       cpScrCoeffReDn,cpScrCoeffImDn,&(cp->cp_comm_state_pkg_dn));
      }/*endif*/
    }/*endif*/
  }

  cpcoeffs_pos->ifcoef_form_up = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1;
  if(cpLsda==1&&numStateDnProc>0){
    cpcoeffs_pos->ifcoef_form_dn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }/*endif*/

  for(iCell=1;iCell<=9;iCell++){ptensPvtenTmp[iCell] = 0.0;}
  gethinv(cell->hmat_cp,cell->hmati_cp,&(cell->vol_cp),iperd);
  gethinv(cell->hmat,cell->hmati,&(cell->vol),iperd);

/*==========================================================================*/
/* VI) Initialize density interpolation                                     */


  //Set occupitation number
  //This is only used in reading wave functions. occNumber will be 
  //changed to normal after reading wave functions. 
  stodftInfo->occNumber = 1;
  if(readCoeffFlag==1&&cpLsda==0)stodftInfo->occNumber = 2;
  
  //set real space density grid number and chemical potential map
  //Now we just use hybrid before interpolation and fullg afterward
  //This will make life easier.
  //For very big system, we may have few stochastic orbitals and 
  //may need to use fullg all the time. That's sounds like a future 
  //implementation

  //if(cpParaOpt==0)stodftInfo->rhoRealGridNum = numFFT2;
  //else stodftInfo->rhoRealGridNum = numFFT2Proc;
  if(chemPotOpt==1){
    stodftInfo->rhoRealGridNum = numFFT2Proc;
    rhoRealGridNum = stodftInfo->rhoRealGridNum;
    stodftInfo->rhoRealGridTot = numFFT2;
    rhoRealGridTot = numFFT2;
    if(cpParaOpt==0){//hybrid case
      div = numChemPot/numProcStates;
      res = numChemPot%numProcStates;
      if(myidState<res)stodftInfo->numChemProc = div+1;
      else stodftInfo->numChemProc = div;
      numChemProc = stodftInfo->numChemProc;
      stodftInfo->densityMap = (int*)cmalloc(numChemPot*sizeof(int));
      stodftInfo->indexChemProc = (int*)cmalloc(numChemPot*sizeof(int));
      stodftInfo->chemProcIndexInv = (int*)cmalloc(numChemProc*sizeof(int));
      for(iChem=0;iChem<numChemPot;iChem++){
	if(iChem<(div+1)*res){
	  stodftInfo->densityMap[iChem] = iChem/(div+1);
	  stodftInfo->indexChemProc[iChem] = iChem%(div+1);
	}
	else{
	  stodftInfo->densityMap[iChem] = (iChem-(div+1)*res)/div+res;
	  stodftInfo->indexChemProc[iChem] = (iChem-(div+1)*res)%div;
	}//endif
      }//endfor iChem
      if(myidState<res)count = myidState*(div+1);
      else count = (div+1)*res+(myidState-res)*div;
      for(iChem=0;iChem<numChemProc;iChem++){
	stodftInfo->chemProcIndexInv[iChem] = count+iChem;
      }//endfor iChem
    }//endfor cpParaOpt
    else{//full g case
      stodftInfo->numChemProc = numChemPot;
      numChemProc = stodftInfo->numChemProc;
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Please use hybrid option!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }
    
    stodftInfo->rhoRealSendCounts = (int*)cmalloc(numProcStates*sizeof(int));
    stodftInfo->rhoRealDispls     = (int*)cmalloc(numProcStates*sizeof(int));
    rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
    rhoRealDispls = stodftInfo->rhoRealDispls;

    if(numProcStates>1){
      if(myidState==0){
	rhoRealSendCounts[0] = numFFT2Proc;
	for(iProc=1;iProc<numProcStates;iProc++){
	  Recv(&rhoRealSendCounts[iProc],1,MPI_INT,iProc,iProc,comm_states);
	}//endfor iProc
      }//endif
      else{
	Send(&numFFT2Proc,1,MPI_INT,0,myidState,comm_states);
      }
      Barrier(comm_states);
      Bcast(rhoRealSendCounts,numProcStates,MPI_INT,0,comm_states);
      Barrier(comm_states);
      rhoRealDispls[0] = 0;
      for(iProc=1;iProc<numProcStates;iProc++){
	rhoRealDispls[iProc] = rhoRealDispls[iProc-1]+rhoRealSendCounts[iProc-1];
      }
    }//endif
    
    // We do this because MPI_Reduce will segfault if the malloc space is not large enough.
    numChemProcMalloc = (res==0?div:div+1);
    //for(iChem=0;iChem<numChemPot;iChem++){
    //printf("myidState %i densityMap[0] %i indexchemproc[0] %i numChemProc %i rhoRealGridTot %i numChemProcMalloc %i\n",
    //    myidState,stodftInfo->densityMap[0],stodftInfo->indexChemProc[0],numChemProc,rhoRealGridTot,numChemProcMalloc);
    //}
    stodftCoefPos->rhoUpChemPot = (double**)cmalloc(numChemProcMalloc*sizeof(double*));
    for(iChem=0;iChem<numChemProcMalloc;iChem++){
      stodftCoefPos->rhoUpChemPot[iChem] = (double*)cmalloc(rhoRealGridTot*sizeof(double));
    }
    stodftCoefPos->rhoUpCorrect = (double*)cmalloc(rhoRealGridNum*sizeof(double));
    if(cpLsda==1&&numStateDnProc>0){
      stodftCoefPos->rhoDnChemPot = (double**)cmalloc(numChemProcMalloc*sizeof(double*));
      for(iChem=0;iChem<numChemProcMalloc;iChem++){
	stodftCoefPos->rhoDnChemPot[iChem] = (double*)cmalloc(rhoRealGridTot*sizeof(double));
      }
      stodftCoefPos->rhoDnCorrect = (double*)cmalloc(rhoRealGridNum*sizeof(double));    
    }
    stodftCoefPos->numElectron = (double*)cmalloc(numChemPot*sizeof(double));
    
    genChemPotInterpPoints(stodftInfo,stodftCoefPos);
    //stodftCoefPos->chemPot[0] = 0.5045818049407941;
    /*
    if(myidState==0){
      for(iChem=0;iChem<numChemPot;iChem++)printf("myid %i numChemProc %i densityMap %i indexChemProc %i\n",myidState,numChemProc,stodftInfo->densityMap[iChem],stodftInfo->indexChemProc[iChem]);
      fflush(stdout);
    }
    exit(0);
    if(myidState==0){
      printf("myid %i pointer rhoUpChemPot[0] %p\n",myidState,stodftCoefPos->rhoUpChemPot[0]);
    }
    */
  }else{//don't do interpolation
    stodftInfo->rhoRealGridNum = numFFT2Proc;
    rhoRealGridNum = stodftInfo->rhoRealGridNum;
    stodftInfo->rhoRealGridTot = numFFT2;
    rhoRealGridTot = numFFT2;
    //if(myidState==0){
    stodftCoefPos->rhoUpChemPot = (double**)cmalloc(sizeof(double*));
    if(myidState==0){
      stodftCoefPos->rhoUpChemPot[0] = (double*)cmalloc(rhoRealGridTot*sizeof(double));
    }
    stodftCoefPos->rhoUpCorrect = (double*)cmalloc(rhoRealGridNum*sizeof(double));
    if(cpLsda==1&&numStateDnProc>0){
      if(myidState==0){
	stodftCoefPos->rhoDnChemPot = (double**)cmalloc(sizeof(double*));
	stodftCoefPos->rhoDnChemPot[0] = (double*)cmalloc(rhoRealGridTot*sizeof(double));
      }
      stodftCoefPos->rhoDnCorrect = (double*)cmalloc(rhoRealGridNum*sizeof(double));
    }
    stodftInfo->rhoRealSendCounts = (int*)cmalloc(numProcStates*sizeof(int));
    stodftInfo->rhoRealDispls     = (int*)cmalloc(numProcStates*sizeof(int));
    rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
    rhoRealDispls = stodftInfo->rhoRealDispls;
    if(numProcStates>1){
      if(myidState==0){
	rhoRealSendCounts[0] = numFFT2Proc;
	for(iProc=1;iProc<numProcStates;iProc++){
	  Recv(&rhoRealSendCounts[iProc],1,MPI_INT,iProc,iProc,comm_states);
	}//endfor iProc
      }//endif
      else{
	Send(&numFFT2Proc,1,MPI_INT,0,myidState,comm_states);
      }
      Barrier(comm_states);
      Bcast(rhoRealSendCounts,numProcStates,MPI_INT,0,comm_states);
      Barrier(comm_states);
      rhoRealDispls[0] = 0;
      for(iProc=1;iProc<numProcStates;iProc++){
	rhoRealDispls[iProc] = rhoRealDispls[iProc-1]+rhoRealSendCounts[iProc-1];
      }
    }//endif numProcStates
  }

/*==========================================================================*/
/* VII) Initialize density mixing                                           */
 
  
  if(densityMixFlag==1){//mixing only, overwirte some parameters
    stodftInfo->numStepMix = numScf+1;
    stodftInfo->numDiis = 1;
  }
  if(densityMixFlag==2){//diis only
    stodftInfo->numStepMix = -1;
  }
  
  if(densityMixFlag>0){//do mix
    stodftCoefPos->rhoUpBank = (double**)cmalloc((numDiis+1)*sizeof(double*));
    stodftCoefPos->rhoDnBank = (double**)cmalloc((numDiis+1)*sizeof(double*));
    stodftCoefPos->rhoUpErr = (double**)cmalloc(numDiis*sizeof(double*));
    stodftCoefPos->rhoDnErr = (double**)cmalloc(numDiis*sizeof(double*));
    stodftCoefPos->rhoUpOld = (double*)cmalloc(rhoRealGridNum*sizeof(double));
    stodftCoefPos->rhoDnOld = (double*)cmalloc(rhoRealGridNum*sizeof(double));

    for(iDiis=0;iDiis<numDiis;iDiis++){
      stodftCoefPos->rhoUpBank[iDiis] = NULL;
      stodftCoefPos->rhoDnBank[iDiis] = NULL;
      stodftCoefPos->rhoUpErr[iDiis] = NULL;
      stodftCoefPos->rhoDnErr[iDiis] = NULL;
    }

    stodftInfo->mixRatioBig = 0.9; // Can tune this parameter
    if(densityMixFlag>1){
      stodftCoefPos->diisMatrix = (double*)cmalloc((numDiis+1)*(numDiis+1)*sizeof(double));
      stodftCoefPos->diisCoeff  = (double*)cmalloc((numDiis+1)*sizeof(double));
      stodftInfo->diisMatrixCalcFullFlag = 1;
    }
  }

/*==========================================================================*/
/* VII) Initialize dynamic density                                          */


  stodftInfo->chemPotHistory = (double*)cmalloc(numScf*sizeof(double));


/*==========================================================================*/
/* VIII) Initialize density output                                          */


  // master proc only check existance
  if(myidState==0){
    printf("%s\n",stodftInfo->densityFileName);
    densityFile = NULL;
    densityFile = cfopen(stodftInfo->densityFileName,"w");
    if(densityFile!=NULL)fclose(densityFile);
  }

/*==========================================================================*/
/* IX) Initialize backup determ wf                                          */


  if(filterDiagFlag==1){
    stodftInfo->numStatesDet = numElecTrue/2;
    int numStatesDet = stodftInfo->numStatesDet;
    if(myidState==0){
      stodftCoefPos->wfUpReDet = (double*)cmalloc(numStatesDet*numCoeff*sizeof(double));
      stodftCoefPos->wfUpImDet = (double*)cmalloc(numStatesDet*numCoeff*sizeof(double));
    }
    if(numProcStates>1){
      stodftInfo->numStatesAllDet = (int*)cmalloc(numProcStates*sizeof(int));
      int *numStatesAllDet = stodftInfo->numStatesAllDet;
      Barrier(comm_states);
      Allgather(&numStateUpTot,1,MPI_INT,numStatesAllDet,1,MPI_INT,0,comm_states);
      Barrier(comm_states);
    
      stodftInfo->dsplStatesAllDet = (int*)cmalloc(numProcStates*sizeof(int));
      int *dsplStatesAllDet = stodftInfo->dsplStatesAllDet;
      dsplStatesAllDet[0] = 0;
      for(iProc=1;iProc<numProcStates;iProc++){
	dsplStatesAllDet[iProc] = dsplStatesAllDet[iProc-1]+numStatesAllDet[iProc-1];
      }
      Barrier(comm_states);
    }
  }

/*==========================================================================*/
/* X) Other allocations                                                     */

  /*
  if(calcFragFlag==1){// We don't initialize frag scf here
    initFrag(class,bonded,general_data,cp,ip_now);
  }
  */

/*==========================================================================*/
/* XI) Initialize Electron Friction                                         */

 
  METALLIC *metallic = stodftInfo->metallic;
  int electronFricFlag = metallic->electronFricFlag;
  int numAtomFric = metallic->numAtomFric;
  int numAtomFricProc;
  FILE *fileAtomFric = NULL;
  int *atomFricInd;
  int *numAtomFricAllProc;
  int *numAtomFricDspl;
  int iAtom;

  if(electronFricFlag==1){
    if(myidState==0){
      fileAtomFric = fopen("atom_index_friction","r");

      if(fileAtomFric==NULL){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Can not find a file named as atom_index_friction\n");
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(0);
      }
      atomFricInd = (int*)cmalloc(numAtomFric*sizeof(int));
      printf("The total number of Atoms to calculate Frictions is %i\n",numAtomFric);
      printf("The list of atoms:\n");
      for(iAtom=0;iAtom<numAtomFric;iAtom++){
        fscanf(fileAtomFric,"%i",&atomFricInd[iAtom]);
	printf("%i ",atomFricInd[iAtom]);
      }
      printf("\n");
      printf("Finish printing the list of atoms\n");
      numAtomFricAllProc = (int*)cmalloc(numProcStates*sizeof(int));
      div = numAtomFric/numProcStates;
      res = numAtomFric%numProcStates;
      for(iProc=0;iProc<numProcStates;iProc++){
        if(iProc<res)numAtomFricAllProc[iProc] = div+1;
        else numAtomFricAllProc[iProc] = div;
      }
      numAtomFricDspl = (int*)cmalloc(numProcStates*sizeof(int));
      numAtomFricDspl[0] = 0;
      for(iProc=1;iProc<numProcStates;iProc++){
        numAtomFricDspl[iProc] = numAtomFricDspl[iProc-1]+numAtomFricAllProc[iProc-1];
      }
    }
    if(numProcStates>1){
      Scatter(numAtomFricAllProc,1,MPI_INT,&(metallic->numAtomFricProc),
              1,MPI_INT,0,comm_states);
    }
    else{
      metallic->numAtomFricProc = numAtomFric;
    }
    
    numAtomFricProc = metallic->numAtomFricProc;
    metallic->atomFricIndProc = (int*)cmalloc(numAtomFricProc*sizeof(int));
    if(numProcStates>1){
      Scatterv(atomFricInd,numAtomFricAllProc,numAtomFricDspl,MPI_INT,metallic->atomFricIndProc,
               numAtomFricProc,MPI_INT,0,comm_states);
    }
    else{
      memcpy(metallic->atomFricIndProc,atomFricInd,numAtomFric*sizeof(int));
    }
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void reInitWaveFunMin(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                   CP *cp,int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the routine to reset number of wave functions after           */
/* calculating the initial density. This is the sp/opt version.		 */
/* check files:#coords_cp/mall_properties.c				 */
/*	       #cp_ewald/control_set_cp_ewald.c				 */
/*	       #parse/parse.c						 */
/*	       #parse/zero_cp.c						 */
/*	       scratch/mall_scratch.c					 */
/*	       				 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  #include "../typ_defs/typ_mask.h"

  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  COMMUNICATE  *communicate     = &(cp->communicate);
  CPOPTS       *cpopts          = &(cp->cpopts);
  CPSCR        *cpscr           = &(cp->cpscr);
  //PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_sclr_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d;

  int numStateStoUp = stodftInfo->numStateStoUp;
  int numStateStoDn = stodftInfo->numStateStoDn;
  int realSparseOpt = cpopts->realSparseOpt;
  int numCoeff      = cpcoeffs_info->ncoef;
  int numStateUpProc,numStateDnProc;
  int numStateUpTot,numStateDnTot;
  int numProcStates = communicate->np_states;
  int myidState	    = communicate->myid;
  int cpLsda        = cpopts->cp_lsda;
  int piBeadsProc   = cpcoeffs_info->pi_beads_proc;
  int hessCalc      = class->clatoms_info.hess_calc;
  int numChemPot    = stodftInfo->numChemPot;
  int numSendNoise;
  int iState,iChem,iProc;
  int nfft,nfft2;
  int reInitFlag    = stodftInfo->reInitFlag;
  int smearOpt      = stodftInfo->smearOpt;

  //int nfft             = cp_para_fft_pkg3d_lg->nfft;
  //int nfft2	       = nfft/2; 
  MPI_Comm comm_states = communicate->comm_states;

  int *noiseSendCounts;
  int *noiseDispls;
  int *randSeedSendCounts;
  int *randSeedDispls;

/*==========================================================================*/
/* I) Initialize Check                                                      */

  if(realSparseOpt==0){
    cp_para_fft_pkg3d = &(cp->cp_sclr_fft_pkg3d_lg);
  }
  else{
    cp_para_fft_pkg3d = &(cp->cp_sclr_fft_pkg3d_sparse);
  }
  nfft = cp_para_fft_pkg3d->nfft;
  nfft2 = nfft/2;
     
  cpcoeffs_info->nstate_up = numStateStoUp;
  cpcoeffs_info->nstate_dn = numStateStoDn;

  if(cpLsda==0){
    if(numStateStoUp<numProcStates){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("Number of states less than number of processors\n");
      printf("If possible, reduce number of processors to be\n");
      printf("less than the number of states or run a bigger system.\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }//endif
  }else{
    if(numStateStoUp+numStateStoDn<numProcStates){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("Number of states less than number of processors\n");
      printf("If possible, reduce number of processors to be\n");
      printf("less than the number of states or run a bigger system.\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }//endif
  }//endif

/*==========================================================================*/
/* II) Reinit communication group for new number of wave function           */
/*     (//N means don't remalloc these arrays)				    */

  
  if(reInitFlag==1)reInitComm(cp,cpcoeffs_pos);

  numStateUpProc = cpcoeffs_info->nstate_up_proc;
  numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  numStateUpTot  = numStateDnProc*numCoeff;
  numStateDnTot  = numStateDnProc*numCoeff;

/*==========================================================================*/
/* III) Reinit arrays, except scratch					    */
  
  if(reInitFlag==1)stoRealloc(cp,cpcoeffs_pos);

/*==========================================================================*/
/* IV) Reinit scratch		                                            */

  if(reInitFlag==1)reallocScratch(cp,hessCalc);


/*==========================================================================*/
/* V) Malloc stochastic wave function                                       */

  stodftCoefPos->stoWfUpRe = (double**)cmalloc(numChemPot*sizeof(double*));
  stodftCoefPos->stoWfUpIm = (double**)cmalloc(numChemPot*sizeof(double*));

  for(iChem=0;iChem<numChemPot;iChem++){
    stodftCoefPos->stoWfUpRe[iChem] = (double*)cmalloc(numStateUpTot*sizeof(double))-1;
    stodftCoefPos->stoWfUpIm[iChem] = (double*)cmalloc(numStateUpTot*sizeof(double))-1;
  }
  if(cpLsda==1&&numStateDnProc!=0){
    //printf("11111111111111111\n");
    stodftCoefPos->stoWfDnRe = (double**)cmalloc(numChemPot*sizeof(double*));
    stodftCoefPos->stoWfDnIm = (double**)cmalloc(numChemPot*sizeof(double*));
    for(iChem=0;iChem<numChemPot;iChem++){
      stodftCoefPos->stoWfDnRe[iChem] = (double*)cmalloc(numStateDnTot*sizeof(double))-1;
      stodftCoefPos->stoWfDnIm[iChem] = (double*)cmalloc(numStateDnTot*sizeof(double))-1;
    }//endfor iChem
  }//endif

  if(smearOpt>0){
    stodftCoefPos->entropyUpRe = (double*)cmalloc(numStateUpTot*sizeof(double))-1;
    stodftCoefPos->entropyUpIm = (double*)cmalloc(numStateUpTot*sizeof(double))-1;
    if(cpLsda==1&&numStateDnProc!=0){
      stodftCoefPos->entropyDnRe = (double*)cmalloc(numStateDnTot*sizeof(double))-1;
      stodftCoefPos->entropyDnIm = (double*)cmalloc(numStateDnTot*sizeof(double))-1;
    }
  }

/*==========================================================================*/
/* VI) Initialize noise orbital scattering                                  */

  stodftInfo->noiseSendCounts = (int*)cmalloc(numProcStates*sizeof(int));
  stodftInfo->noiseDispls     = (int*)cmalloc(numProcStates*sizeof(int));
  noiseSendCounts = stodftInfo->noiseSendCounts;
  noiseDispls = stodftInfo->noiseDispls;

  if(numProcStates>1){
    if(myidState==0){
      if(cpLsda==0)noiseSendCounts[0] = numStateUpTot*2;
      else noiseSendCounts[0] = (numStateUpTot+numStateDnTot)*2;
      for(iProc=1;iProc<numProcStates;iProc++){
        Recv(&noiseSendCounts[iProc],1,MPI_INT,iProc,iProc,comm_states);
      }//endfor iProc
    }//endif
    else{
      if(cpLsda==0){ 
        numSendNoise = numStateUpTot*2;
        Send(&numSendNoise,1,MPI_INT,0,myidState,comm_states);
      }
      else{
        numSendNoise = (numStateUpTot+numStateDnTot)*2;
        Send(&numSendNoise,1,MPI_INT,0,myidState,comm_states);
      }
    }
    Barrier(comm_states);
    Bcast(noiseSendCounts,numProcStates,MPI_INT,0,comm_states);
    Barrier(comm_states);
    //printf("0 %i 1 %i 2 %i 3 %i\n",noiseSendCounts[0],noiseSendCounts[1],noiseSendCounts[2],noiseSendCounts[3]); 
    noiseDispls[0] = 0;
    for(iProc=1;iProc<numProcStates;iProc++){
      noiseDispls[iProc] = noiseDispls[iProc-1]+noiseSendCounts[iProc-1];
    }
    stodftInfo->numRandTot = 0;
    for(iProc=0;iProc<numProcStates;iProc++){
      stodftInfo->numRandTot += noiseSendCounts[iProc];
    }
  }//endif 
  else{
    stodftInfo->numRandTot = numStateUpTot*2;
    if(cpLsda==1)stodftInfo->numRandTot += numStateDnTot*2;
  }

/*==========================================================================*/
/* VI) Reset some flags so that the program will not crash                  */

  cpcoeffs_pos->icoef_orth_up = 1;
  if(cpLsda==0)cpcoeffs_pos->icoef_orth_dn = 1;

  cpcoeffs_pos->icoef_form_up = 0;
  cpcoeffs_pos->ifcoef_form_up = 0;
  if(cpLsda==0){
    cpcoeffs_pos->icoef_form_dn = 0;
    cpcoeffs_pos->ifcoef_form_dn = 0;    
  }
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void reInitComm(CP *cp,CPCOEFFS_POS *cpCoeffsPos)
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*    Local Variables   */
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPCOEFFS_INFO *cpCoeffsInfo  = &(cp->cpcoeffs_info);
  COMMUNICATE  *communicate     = &(cp->communicate);
  CPOPTS       *cpopts          = &(cp->cpopts);
  CPSCR        *cpscr           = &(cp->cpscr);
  CP_COMM_STATE_PKG *cpCommStatePkgUp  = &(cp->cp_comm_state_pkg_up);
  CP_COMM_STATE_PKG *cpCommStatePkgDn  = &(cp->cp_comm_state_pkg_up);

  int irem,idiv,iii;
  int numStateUp = cpCoeffsInfo->nstate_up;
  int numStateDn = cpCoeffsInfo->nstate_dn;
  int numProcStates = communicate->np_states;
  int myidState = communicate->myid_state;
  int numStateUpProc,iStateUpSt;
  int numStateDnProc,iStateDnSt;

/*==========================================================================*/
/* I) Up states, state per process                                          */

  idiv = numStateUp/numProcStates;
  irem = numStateUp%numProcStates;
  numStateUpProc = idiv;

  if(myidState<irem)numStateUpProc = idiv+1;
  if(myidState<=irem)iStateUpSt = myidState*(idiv+1)+1;
  else iStateUpSt = irem*(idiv+1)+(myidState-irem)*idiv+1;

  cpCoeffsInfo->nstate_up_proc = numStateUpProc;
  cpCoeffsInfo->istate_up_st = iStateUpSt;
  cpCoeffsInfo->istate_up_end = iStateUpSt+numStateUpProc-1;

  cp->cp_comm_state_pkg_up.nstate      = numStateUp;
  cp->cp_comm_state_pkg_up.nstate_proc = numStateUpProc;

  // I have no idea what is this used for, since it is initialized in the 
  // begining, I shall reinit this to be safe.

  irem = numStateUp%numProcStates;
  if(irem>0){
    cp->cp_comm_state_pkg_up.nstate_proc_max = idiv+1;
    cp->cp_comm_state_pkg_up.nstate_max = (idiv+1)*numProcStates;
  }
  else{
    cp->cp_comm_state_pkg_up.nstate_proc_max = idiv;
    cp->cp_comm_state_pkg_up.nstate_max = idiv*numProcStates;
  }

/*==========================================================================*/
/* I) Down states, state per process                                        */

  idiv = numStateDn/numProcStates;
  irem = numStateDn%numProcStates;
  numStateDnProc = idiv;

  if(myidState<irem)numStateDnProc = idiv+1;
  if(myidState<=irem)iStateDnSt = myidState*(idiv+1)+1;
  else iStateDnSt = irem*(idiv+1)+(myidState-irem)*idiv+1;

  cpCoeffsInfo->nstate_dn_proc = numStateDnProc;
  cpCoeffsInfo->istate_dn_st = iStateDnSt;
  cpCoeffsInfo->istate_dn_end = iStateDnSt+numStateDnProc-1;

  cp->cp_comm_state_pkg_dn.nstate      = numStateDn;
  cp->cp_comm_state_pkg_dn.nstate_proc = numStateDnProc;

  irem = numStateDn%numProcStates;
  if(irem>0){
    cp->cp_comm_state_pkg_dn.nstate_proc_max = idiv+1;
    cp->cp_comm_state_pkg_dn.nstate_max = (idiv+1)*numProcStates;
  }
  else{
    cp->cp_comm_state_pkg_dn.nstate_proc_max = idiv;
    cp->cp_comm_state_pkg_dn.nstate_max = idiv*numProcStates;
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void stoRealloc(CP *cp,CPCOEFFS_POS *cpcoeffs_pos)
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*    Local Variables   */
#include "../typ_defs/typ_mask.h"
  CPOPTS *cpopts = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);

  int is,i;
  int par_size_up,par_size_dn,ncoef_up_tot;
  int nstate,nstate2,ncoef_dn_tot;
  double *cre,*cim,*vcre,*vcim,*fcre,*fcim;
  int *ioff_up,*ioff_upt,*ioff_dn,*ioff_dnt;

/*  Local Pointers */
  int nstate_up_proc = cp->cpcoeffs_info.nstate_up_proc;
  int nstate_dn_proc = cp->cpcoeffs_info.nstate_dn_proc;
  int nstate_up      = cp->cpcoeffs_info.nstate_up;
  int nstate_dn      = cp->cpcoeffs_info.nstate_dn;
  int cp_lda         = cp->cpopts.cp_lda;
  int cp_lsda        = cp->cpopts.cp_lsda;
  int ncoef_up       = cp->cpcoeffs_info.ncoef;
  int ncoef_dn       = cp->cpcoeffs_info.ncoef;
  int cp_norb        = cp->cpopts.cp_norb;

  int np_states               = cp->communicate.np_states;
  int ncoef_up_proc           = cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
  int ncoef_dn_proc           = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;
  int nstate_max_up           =cp->cp_comm_state_pkg_up.nstate_max;
  int nstate_ncoef_proc_max_up=cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
  int nstate_max_dn           =cp->cp_comm_state_pkg_dn.nstate_max;
  int nstate_ncoef_proc_max_dn=cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;

  if(cp_lda == 1){ncoef_dn = 0;}
/*==========================================================================*/
/* I) Calculate the sizes */

  par_size_up = nstate_max_up*(nstate_ncoef_proc_max_up);
  ncoef_up_tot = nstate_up_proc*ncoef_up;
  ncoef_up_tot = MAX(ncoef_up_tot,par_size_up);

  par_size_dn = nstate_max_dn*(nstate_ncoef_proc_max_dn);
  if(cp_lda ==1){par_size_dn = 0;}
  ncoef_dn_tot = MAX(nstate_dn_proc*ncoef_dn,1);
  ncoef_dn_tot = MAX(ncoef_dn_tot,par_size_dn);

  nstate  = MAX(nstate_up,nstate_dn);
  nstate2 = nstate*nstate;

/*==========================================================================*/
/* II) Free everything */

  free(&(cpcoeffs_pos->cre_up[1]));
  free(&(cpcoeffs_pos->cim_up[1]));
  free(&(cpcoeffs_pos->cre_dn[1]));
  free(&(cpcoeffs_pos->cim_dn[1]));
  free(&(cpcoeffs_pos->vcre_up[1]));//N
  free(&(cpcoeffs_pos->vcim_up[1]));//N
  free(&(cpcoeffs_pos->vcre_dn[1]));//N
  free(&(cpcoeffs_pos->vcim_dn[1]));//N
  free(&(cpcoeffs_pos->fcre_up[1]));
  free(&(cpcoeffs_pos->fcim_up[1]));
  free(&(cpcoeffs_pos->fcre_dn[1]));
  free(&(cpcoeffs_pos->fcim_dn[1]));
  free(&(cpcoeffs_pos->kfcre_up[1]));//N
  free(&(cpcoeffs_pos->kfcim_up[1]));//N
  free(&(cpcoeffs_pos->ksmat_up[1]));//N
  free(&(cpcoeffs_pos->ksmat_dn[1]));//N
  free(&(cpcoeffs_pos->ksmat_eig_up[1]));//N
  free(&(cpcoeffs_pos->ksmat_eig_dn[1]));//N
  free(&(cpcoeffs_pos->norbmat_up[1]));//N
  free(&(cpcoeffs_pos->norbmat_dn[1]));//N
  free(&(cpcoeffs_pos->norbmati_up[1]));//N
  free(&(cpcoeffs_pos->norbmati_dn[1]));//N
  free(&(cpcoeffs_pos->ovmat_eigv_up[1]));//N
  free(&(cpcoeffs_pos->ovmat_eigv_dn[1]));//N

  // cp_min_on>0 I need fake cp_min_on to cheat the code 
  // when I call coef_force_control
  // I still need to realloc cp_hess to avoid segfault
  free(&(cpcoeffs_pos->cp_hess_re_up[1]));//N
  free(&(cpcoeffs_pos->cp_hess_im_up[1]));//N
  free(&(cpcoeffs_pos->cp_hess_re_dn[1]));//N
  free(&(cpcoeffs_pos->cp_hess_im_dn[1]));//N

  free(&(cpcoeffs_info->ioff_up[1]));
  free(&(cpcoeffs_info->ioff_dn[1]));
  free(&(cpcoeffs_info->ioff_upt[1]));
  free(&(cpcoeffs_info->ioff_dnt[1]));

  free(&(cpopts->occ_up[1]));//N
  free(&(cpopts->occ_dn[1]));//N
  free(&(cpopts->rocc_sum_up[1]));//N
  free(&(cpopts->rocc_sum_dn[1]));//N



/*==========================================================================*/
/* III) Malloc the variables */

  cpcoeffs_pos->cre_up =(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cpcoeffs_pos->cim_up =(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cpcoeffs_pos->cre_dn =(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  cpcoeffs_pos->cim_dn =(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  cpcoeffs_pos->fcre_up=(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cpcoeffs_pos->fcim_up=(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cpcoeffs_pos->fcre_dn=(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  cpcoeffs_pos->fcim_dn=(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;

  cp->cpcoeffs_info.ioff_up  = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_dn  = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_upt = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_dnt = (int *)cmalloc(nstate*sizeof(int))-1;

/*==========================================================================*/
/* IV) Assign the offsets */

  ioff_up  = cp->cpcoeffs_info.ioff_up;
  ioff_upt = cp->cpcoeffs_info.ioff_upt;
  ioff_dn  = cp->cpcoeffs_info.ioff_dn;
  ioff_dnt = cp->cpcoeffs_info.ioff_dnt;

  for(is=1;is<=nstate;is++){ioff_up[is]=(is-1)*ncoef_up;}
  for(is=1;is<=nstate;is++){ioff_dn[is]=(is-1)*ncoef_dn;}

  if(np_states==1){
    for(is=1;is<=nstate;is++){ioff_upt[is]=(is-1)*ncoef_up;}
    for(is=1;is<=nstate;is++){ioff_dnt[is]=(is-1)*ncoef_dn;}
  }else{
    for(is=1;is<=nstate;is++){ioff_upt[is]=(is-1)*ncoef_up_proc;}
    for(is=1;is<=nstate;is++){ioff_dnt[is]=(is-1)*ncoef_dn_proc;}
  }/*endif*/

/*========================================================================*/
/* V) Initialize coeficient arrays and form flags                       */


  cpcoeffs_pos->icoef_form_up  = 0;
  cpcoeffs_pos->ivcoef_form_up = 0;
  cpcoeffs_pos->ifcoef_form_up = 1;
  cpcoeffs_pos->icoef_orth_up  = 1;
  if(cp_norb>0){cpcoeffs_pos->icoef_orth_up = 0;}
  cpcoeffs_pos->ivcoef_orth_up  = 1;
  if(cp_norb>0){cpcoeffs_pos->ivcoef_orth_up = 0;}
  cre = cpcoeffs_pos->cre_up;
  cim = cpcoeffs_pos->cim_up;
  fcre = cpcoeffs_pos->fcre_up;
  fcim = cpcoeffs_pos->fcim_up;
  for(i=1;i<=ncoef_up_tot;i++){
    cre[i] = 0.0;
    cim[i] = 0.0;
    fcre[i] = 0.0;
    fcim[i] = 0.0;
  }/*endfor*/
  if(cp_lsda==1){
    cpcoeffs_pos->icoef_form_dn  = 0;
    cpcoeffs_pos->ivcoef_form_dn = 0;
    cpcoeffs_pos->ifcoef_form_dn = 1;
    cpcoeffs_pos->icoef_orth_dn  = 1;
    if(cp_norb>0)cpcoeffs_pos->icoef_orth_dn = 0;
    cpcoeffs_pos->ivcoef_orth_dn  = 1;
    if(cp_norb>0)cpcoeffs_pos->ivcoef_orth_dn = 0;

    cre = cpcoeffs_pos->cre_dn;
    cim = cpcoeffs_pos->cim_dn;
    fcre = cpcoeffs_pos->fcre_dn;
    fcim = cpcoeffs_pos->fcim_dn;
    for(i=1;i<=ncoef_dn_tot;i++){
      cre[i] = 0.0;
      cim[i] = 0.0;
      fcre[i] = 0.0;
      fcim[i] = 0.0;
    }/*endfor*/
  }/*endif:lsda*/


/*-----------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void reallocScratch(CP *cp,int hess_calc)
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*    Local Variables   */
  CPOPTS *cpopts = &(cp->cpopts);
  CPSCR *cpscr = &(cp->cpscr);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPTHERM_INFO *cptherm_info = &(cp->cptherm_info);
  PSEUDO *pseudo = &(cp->pseudo);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_up = &(cp->cp_comm_state_pkg_up);
  //PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d;
 
  CPSCR_LOC    *cpscr_loc = &(cpscr->cpscr_loc);
  CPSCR_NONLOC *cpscr_nonloc = &(cpscr->cpscr_nonloc);
  CPSCR_RHO    *cpscr_rho = &(cpscr->cpscr_rho);
  CPSCR_GRHO   *cpscr_grho = &(cpscr->cpscr_grho);
  CPSCR_OVMAT  *cpscr_ovmat = &(cpscr->cpscr_ovmat);
  CPSCR_WAVE   *cpscr_wave = &(cpscr->cpscr_wave);
  CPSCR_THERM  *cpscr_therm = &(cpscr->cpscr_therm);
  CPSCR_DUAL_PME *cpscr_dual_pme = &(cpscr->cpscr_dual_pme);
  CPSCR_ATOM_PME *cpscr_atom_pme = &(cpscr->cpscr_atom_pme);
  CPSCR_WANNIER *cpscr_wannier = &(cpscr->cpscr_wannier);

  int cp_ptens_calc     = cpopts->cp_ptens_calc;
  int cp_dual_grid_opt_on = cpopts->cp_dual_grid_opt;
  int nstate_up         = cpcoeffs_info->nstate_up_proc;
  int nstate_dn         = cpcoeffs_info->nstate_dn_proc;
  int nstate_up_tot     = cpcoeffs_info->nstate_up;
  int nstate_dn_tot     = cpcoeffs_info->nstate_dn;
  int ncoef             = cpcoeffs_info->ncoef;
  int ncoef_l           = cpcoeffs_info->ncoef_l;
  int num_c_nhc_proc    = cptherm_info->num_c_nhc_proc;
  int massiv_flag       = cptherm_info->massiv_flag;
  int n_ang_max         = pseudo->n_ang_max;
  int n_rad_max         = pseudo->n_rad_max;
  int natm_nls_max      = cpscr_nonloc->natm_nls_max;
  int cp_lsda           = cpopts->cp_lsda;
  int cp_norb           = cpopts->cp_norb;
  int realSparseOpt	= cpopts->realSparseOpt;
  int np_states         = cp_comm_state_pkg_up->num_proc;
  int nstate_max_up     = cp_comm_state_pkg_up->nstate_max;
  int nstate_ncoef_proc_max_up = cp_comm_state_pkg_up->nstate_ncoef_proc_max;
  int nstate_max_dn     = cp_comm_state_pkg_up->nstate_max;
  int nstate_ncoef_proc_max_dn = cp_comm_state_pkg_up->nstate_ncoef_proc_max;

  int ncoef_l_pme_dual,ncoef_l_pme_dual_proc;
  int ncoef_l_proc_max_mall;

/*--------------------------------------------------------------------------*/
/*         Local variable declarations                                      */

  int i,iii,irem;
  int nlscr_up,nlscr_dn,nlscr_up_pv,nlscr_dn_pv,ncoef_l_pv,ncoef_l_proc_max;
  int ncoef_up,ncoef_dn;
  int par_size_up,par_size_dn;
  int ncoef2_up_c,ncoef2_dn_c;
  int ncoef2_up,ncoef2_dn,ncoef2_up_spec,ncoef2_dn_spec;
  int ncoef2_up_par,nstate,nstate2,nstate_tot,nstate2_tot;

  int cp_wan_opt        = cpopts->cp_wan_opt;
  int cp_wan_min_opt    = cpopts->cp_wan_min_opt;
  int cp_wan_init_opt   = cpopts->cp_wan_init_opt;
  int ndim_wannier;
  int mm=5;

/*==========================================================================*/
/* Free everything*/

  if(realSparseOpt==0){
    cp_para_fft_pkg3d = &(cp->cp_para_fft_pkg3d_lg);
  }
  else{
    cp_para_fft_pkg3d = &(cp->cp_para_fft_pkg3d_sparse);
  }

  free(&(cpscr_nonloc->vnlre_up[1]));
  free(&(cpscr_nonloc->vnlim_up[1]));
  free(&(cpscr_nonloc->vnlre_dn[1]));
  free(&(cpscr_nonloc->vnlim_dn[1]));

  free(&(cpscr_nonloc->dvnlre_x_up[1]));
  free(&(cpscr_nonloc->dvnlre_y_up[1]));
  free(&(cpscr_nonloc->dvnlre_z_up[1]));
  free(&(cpscr_nonloc->dvnlim_x_up[1]));
  free(&(cpscr_nonloc->dvnlim_y_up[1]));
  free(&(cpscr_nonloc->dvnlim_z_up[1]));

  free(&(cpscr_nonloc->dvnlre_x_dn[1]));
  free(&(cpscr_nonloc->dvnlre_y_dn[1]));
  free(&(cpscr_nonloc->dvnlre_z_dn[1]));
  free(&(cpscr_nonloc->dvnlim_x_dn[1]));
  free(&(cpscr_nonloc->dvnlim_y_dn[1]));
  free(&(cpscr_nonloc->dvnlim_z_dn[1]));

  free(&(cpscr_nonloc->dvnlre_gxgx_up[1]));
  free(&(cpscr_nonloc->dvnlre_gygy_up[1]));
  free(&(cpscr_nonloc->dvnlre_gzgz_up[1]));
  free(&(cpscr_nonloc->dvnlre_gxgy_up[1]));
  free(&(cpscr_nonloc->dvnlre_gygz_up[1]));
  free(&(cpscr_nonloc->dvnlre_gxgz_up[1]));
  free(&(cpscr_nonloc->dvnlim_gxgx_up[1]));
  free(&(cpscr_nonloc->dvnlim_gygy_up[1]));
  free(&(cpscr_nonloc->dvnlim_gzgz_up[1]));
  free(&(cpscr_nonloc->dvnlim_gxgy_up[1]));
  free(&(cpscr_nonloc->dvnlim_gygz_up[1]));
  free(&(cpscr_nonloc->dvnlim_gxgz_up[1]));

  free(&(cpscr_nonloc->dvnlre_gxgx_dn[1]));
  free(&(cpscr_nonloc->dvnlre_gygy_dn[1]));
  free(&(cpscr_nonloc->dvnlre_gzgz_dn[1]));
  free(&(cpscr_nonloc->dvnlre_gxgy_dn[1]));
  free(&(cpscr_nonloc->dvnlre_gygz_dn[1]));
  free(&(cpscr_nonloc->dvnlre_gxgz_dn[1]));
  free(&(cpscr_nonloc->dvnlim_gxgx_dn[1]));
  free(&(cpscr_nonloc->dvnlim_gygy_dn[1]));
  free(&(cpscr_nonloc->dvnlim_gzgz_dn[1]));
  free(&(cpscr_nonloc->dvnlim_gxgy_dn[1]));
  free(&(cpscr_nonloc->dvnlim_gygz_dn[1]));
  free(&(cpscr_nonloc->dvnlim_gxgz_dn[1]));

  free(&(cpscr_wave->cre_up[1]));
  free(&(cpscr_wave->cim_up[1]));
  free(&(cpscr_wave->cre_dn[1]));
  free(&(cpscr_wave->cim_dn[1]));

  free(&(cpscr_ovmat->ovlap1[1]));
  free(&(cpscr_ovmat->ovlap2[1]));
  free(&(cpscr_ovmat->ovlap3[1]));
  free(&(cpscr_ovmat->ovlap4[1]));
  free(&(cpscr_ovmat->ovlap5[1]));
  free(&(cpscr_ovmat->ovlap6[1]));
  free(&(cpscr_ovmat->ovlap7[1]));
  free(&(cpscr_ovmat->ovlap8[1]));

  free(&(cpscr_ovmat->state_vec1[1]));
  free(&(cpscr_ovmat->state_vec2[1]));
  free(&(cpscr_ovmat->state_vec3[1]));
  free(&(cpscr_ovmat->state_vec4[1]));
  free(&(cpscr_ovmat->state_vec5[1]));
  free(&(cpscr_ovmat->state_vec6[1]));
  free(&(cpscr_ovmat->state_vec7[1]));

  if(cpscr_wannier->cp_wannier_on==1){
    free(cpscr_wannier->A[1]);
    free(cpscr_wannier->R_real[1]);
    free(cpscr_wannier->R_imag[1]);
    free(cpscr_wannier->U_real[1]);
    free(cpscr_wannier->U_tmp1[1]);
    free(cpscr_wannier->U_tmp2[1]);
    free(cpscr_wannier->U_tmp3[1]);
    free(cpscr_wannier->U_tmp4[1]);
    free(cpscr_wannier->Bt_real[1]);
    free(cpscr_wannier->Bt_imag[1]);
    free(cpscr_wannier->M_real[1]);
    free(cpscr_wannier->phi[1]);
    free(&(cpscr_wannier->A[1]));
    free(&(cpscr_wannier->R_real[1]));
    free(&(cpscr_wannier->R_imag[1]));
    free(&(cpscr_wannier->U_real[1]));
    free(&(cpscr_wannier->U_tmp1[1]));
    free(&(cpscr_wannier->U_tmp2[1]));
    free(&(cpscr_wannier->U_tmp3[1]));
    free(&(cpscr_wannier->U_tmp4[1]));
    free(&(cpscr_wannier->Bt_real[1]));
    free(&(cpscr_wannier->Bt_imag[1]));
    free(&(cpscr_wannier->M_real[1]));
    free(&(cpscr_wannier->phi[1]));

    free(&(cpscr_wannier->real));
    free(&(cpscr_wannier->imag));
    free(&(cpscr_wannier->D));
    free(&(cpscr_wannier->norm));
  }

/*=========================================================================*/
/* I) Malloc size calculation */

 /*-------------------------------------------------------------------------*/
 /* i) Dual grid CP : Define the small dense grid sizes */

 if(cp_dual_grid_opt_on == 2){
   ncoef_l_pme_dual = cp_para_fft_pkg3d->ncoef;
 }/*endif cp_dual_grid_opt_on*/

 /*-------------------------------------------------------------------------*/
 /* ii) Normal CP : Define the grid size              */
 /*     Dual   CP : Define the large sparse grid size */


 /*-------------------------------------------------------------------------*/
 /* iii) Choose the correct size for your application                       */
 /*      This is always the small dense grid                                */

 /*-------------------------------------------------------------------------*/
 /* iv) Wave function size (spherically cutoff small g-space for dense box) */

  ncoef_up  = ncoef;
  ncoef_dn  = (cp_lsda == 1 ? ncoef : 0);

 /*-------------------------------------------------------------------------*/
 /* v) Normal CP: The sphere cut large g-space for the dense box            */
 /*    Dual CP  : The sphere cut large g-space for the large sparse box     */

  ncoef_l_proc_max = ncoef_l/np_states;
  irem = (ncoef_l % np_states);
  if(irem>0){ncoef_l_proc_max++;}

  if(cp_dual_grid_opt_on == 2){
    ncoef_l_pme_dual_proc = ncoef_l_pme_dual/np_states;
    irem = (ncoef_l_pme_dual_proc % np_states);
    if(irem>0){ncoef_l_pme_dual_proc++;}
  }/*endif cp_dual_grid_opt_on */

  ncoef_l_proc_max_mall = (cp_dual_grid_opt_on==2?ncoef_l_pme_dual_proc
                                                    :ncoef_l_proc_max);

 /*-------------------------------------------------------------------------*/
 /* vi) Choose the large g-space malloc size based on the dual or normal opt*/
 /*     The malloc size is the always the dense grid.                       */
 /*     Set special dual malloc sizes to zero to avoid mallocing extra      */
 /*     memory during normal CP.                                            */

 /*-------------------------------------------------------------------*/
 /* vii) Wavefunction scratch sizes : always on small dense grid     */

  par_size_up = nstate_max_up*nstate_ncoef_proc_max_up;
  par_size_dn = nstate_max_dn*nstate_ncoef_proc_max_dn;

  if(massiv_flag==0){
    ncoef2_up_c = MAX(ncoef_up*nstate_up,num_c_nhc_proc);
  }else{
    ncoef2_up_c = MAX(ncoef_up*nstate_up,2*ncoef_up);
  }/*endif*/

  if(np_states>1){ncoef2_up_c = MAX(ncoef2_up_c,par_size_up);}
  ncoef2_up      = ncoef_up*nstate_up;
  ncoef2_dn      = ncoef_dn*nstate_dn;
  ncoef2_up      = MAX(ncoef2_up,par_size_up);
  ncoef2_dn      = MAX(ncoef2_dn,par_size_dn);
  ncoef2_up_spec = MAX(ncoef2_up,ncoef2_dn);
  ncoef2_dn_spec = ncoef2_up_spec;
  if(np_states == 1 &&  cp_norb==0){ncoef2_up_spec = 1;}
  ncoef2_up_par = ncoef2_up_spec;
  if(np_states == 1){ncoef2_up_par = 1;}

  nstate        = MAX(nstate_up, nstate_dn);
  nstate2       = nstate*nstate;
  nstate_tot    = MAX(nstate_up_tot,nstate_dn_tot);
  nstate2_tot   = nstate_tot*nstate_tot;

 /*-------------------------------------------------------------------*/
 /* viii) Nonlocal sizes : always on small dense grid */

  nlscr_up  = nstate_up*(n_ang_max+1)*(n_ang_max+1)
                        *(n_rad_max)*natm_nls_max;
  nlscr_dn  = 0;
  if(cp_lsda==1){
     nlscr_dn = nstate_dn*(n_ang_max+1)*(n_ang_max+1)
                         *(n_rad_max)*(natm_nls_max);
  }/*endif*/
  nlscr_up_pv = 0;
  nlscr_dn_pv = 0;
  ncoef_l_pv  = 0;
  if(cp_ptens_calc == 1){
    nlscr_up_pv = nlscr_up;
    nlscr_dn_pv = nlscr_dn;
    ncoef_l_pv  = ncoef_l_proc_max_mall;
  }/* endif */
  if(hess_calc == 3){
    nlscr_up_pv = nlscr_up;
    nlscr_dn_pv = nlscr_dn;
  }/* endif */

 /*-------------------------------------------------------------------*/
 /* ix) GGA sizes : Always the small dense grid                       */


  /*------------------------------------------------------------------------*/
  /* x) Wannier scratch size   */

  if(cp_wan_opt ==1 || cp_wan_min_opt==1 || cp_wan_init_opt==1){
    cpscr_wannier->cp_wannier_on=1;
    ndim_wannier=nstate_up_tot*nstate_up_tot;
  }else{
    cpscr_wannier->cp_wannier_on=0;
    ndim_wannier=0;
  }


/*===========================================================================*/
/* II) Malloc the vectors  */

/*------------------------------------------------------------------*/
/* Non_local  : Always small dense grid                             */

  cpscr_nonloc->vnlre_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr_nonloc->vnlim_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;

  cpscr_nonloc->vnlre_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr_nonloc->vnlim_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;

  cpscr_nonloc->dvnlre_x_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr_nonloc->dvnlre_y_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr_nonloc->dvnlre_z_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr_nonloc->dvnlim_x_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr_nonloc->dvnlim_y_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr_nonloc->dvnlim_z_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;

  cpscr_nonloc->dvnlre_x_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr_nonloc->dvnlre_y_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr_nonloc->dvnlre_z_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr_nonloc->dvnlim_x_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr_nonloc->dvnlim_y_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr_nonloc->dvnlim_z_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;

  cpscr_nonloc->dvnlre_gxgx_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlre_gygy_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlre_gzgz_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlre_gxgy_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlre_gygz_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlre_gxgz_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlim_gxgx_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlim_gygy_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlim_gzgz_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlim_gxgy_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlim_gygz_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlim_gxgz_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;


  cpscr_nonloc->dvnlre_gxgx_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlre_gygy_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlre_gzgz_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlre_gxgy_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlre_gygz_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlre_gxgz_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;

  cpscr_nonloc->dvnlim_gxgx_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlim_gygy_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlim_gzgz_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlim_gxgy_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlim_gygz_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr_nonloc->dvnlim_gxgz_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;

/*------------------------------------------------------------------*/
/* wave */
  cpscr_wave->cre_up
                    = (double *)cmalloc(ncoef2_up_c*sizeof(double))-1;
  cpscr_wave->cim_up
                    = (double *)cmalloc(ncoef2_up*sizeof(double))-1;
  cpscr_wave->cre_dn
                    = (double *)cmalloc(ncoef2_dn*sizeof(double))-1;
  cpscr_wave->cim_dn
                    = (double *)cmalloc(ncoef2_dn*sizeof(double))-1;


/*------------------------------------------------------------------*/
/* ovmat */
  cpscr_ovmat->ovlap1
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr_ovmat->ovlap2
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr_ovmat->ovlap3
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr_ovmat->ovlap4
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr_ovmat->ovlap5
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr_ovmat->ovlap6
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr_ovmat->ovlap7
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr_ovmat->ovlap8
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;

  cpscr_ovmat->state_vec1
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr_ovmat->state_vec2
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr_ovmat->state_vec3
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr_ovmat->state_vec4
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr_ovmat->state_vec5
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr_ovmat->state_vec6
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr_ovmat->state_vec7
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;

/*-----------------------------------------------------------------------*/
/* Wannier scratch                                                       */

  if(cpscr_wannier->cp_wannier_on==1){
    cpscr_wannier->U_final = (double *) cmalloc(ndim_wannier*sizeof(double))-1;
    cpscr_wannier->g       = (double *) cmalloc(ndim_wannier*sizeof(double))-1;
    cpscr_wannier->gg      = (double *) cmalloc(ndim_wannier*sizeof(double))-1;
    cpscr_wannier->diag    = (double *) cmalloc(ndim_wannier*sizeof(double))-1;
    cpscr_wannier->scr     =
                   (double *) cmalloc((ndim_wannier*(2*mm+1)+2*mm)*sizeof(double))-1;
    cpscr_wannier->iprint  = (int *)cmalloc(2*sizeof(int))-1;

    cpscr_wannier->A       = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr_wannier->R_real  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr_wannier->R_imag  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr_wannier->U_real  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr_wannier->U_imag  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr_wannier->U_tmp1  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr_wannier->U_tmp2  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr_wannier->U_tmp3  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr_wannier->U_tmp4  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr_wannier->Bt_real = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr_wannier->Bt_imag = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr_wannier->M_real  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);

    cpscr_wannier->real    = (double *)cmalloc(nstate_up_tot*sizeof(double))-1;
    cpscr_wannier->imag    = (double *)cmalloc(nstate_up_tot*sizeof(double))-1;
    cpscr_wannier->D       = (double *)cmalloc(nstate_up_tot*sizeof(double))-1;
    cpscr_wannier->norm    = (double *)cmalloc(nstate_up_tot*sizeof(double))-1;
    cpscr_wannier->phi     = cmall_mat(1,nstate_up_tot,1,3);
    cpscr_wannier->HMatrix = cmall_mat(1,3,1,3);
    printf("CP_WANNIER memory allocation \n");
  } 

  
/*-----------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initFilterDiag(CP *cp)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  CPCOEFFS_INFO *cpcoeffs_info    = &(cp->cpcoeffs_info);
  STODFTINFO    *stodftInfo       = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos    = cp->stodftCoefPos;
  COMMUNICATE *communicate      = &(cp->communicate);

  int iChem,iCoeff,iState,jState,iProc;
  int index1,index2;
  int numChemPot     = stodftInfo->numChemPot;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numStateUpTotal = stodftInfo->numStateStoUp;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int numCoeffUpAllChemPot = numCoeffUpTotal*numChemPot;
  int numCoeffUpAllProc = numChemPot*numStateUpTotal*numCoeff;
  int numStateUpAllProc = numChemPot*numStateUpTotal;
  int numStateUpIdp = stodftInfo->numStateUpIdp;
  int numElecTrue = stodftInfo->numElecTrue;
  int numOccState = numElecTrue/2;
  int numStateProcTot = numStateUpProc*numChemPot;
  int stateIndex;
  int myidState = communicate->myid_state;
  int numProcStates = communicate->np_states;
  MPI_Comm comm_states = communicate->comm_states;
  
  int *stowfRecvCounts,*stowfDispls;
  int *numStates;
  int *dsplStates;

  double pre = sqrt(2.0);

  double *numOccDetAll,*numOccDetProc;

/*===========================================================================*/
/* I) Malloc the vectors  */
  //printf("numStateUpProc %i\n",numStateUpProc);

  if(myidState==0){
    printf("Start Initializing Filter Diagonalization!\n");
    stodftCoefPos->wfBfOrthUp = (double *)cmalloc(2*numCoeffUpAllProc*
				sizeof(double));
    //stodftCoefPos->wfOrthUpRe = (double *)cmalloc(numCoeffUpAllProc*
    //	 			sizeof(double))-1;
    //stodftCoefPos->wfOrthUpIm = (double *)cmalloc(numCoeffUpAllProc*
    //				sizeof(double))-1;
    //stodftCoefPos->KSMatrix = (double *)cmalloc(numStateUpAllProc*numStateUpAllProc*
    //				sizeof(double));
  }
  stodftCoefPos->KSMatrix = NULL;
  //stodftCoefPos->energyLevel = (double*)cmalloc(numStateUpAllProc*sizeof(double));
  stodftCoefPos->energyLevel = NULL;
  //stodftCoefPos->moUpRe = (double*)cmalloc(numCoeffUpTotal*numChemPot*sizeof(double))-1;
  //stodftCoefPos->moUpIm = (double*)cmalloc(numCoeffUpTotal*numChemPot*sizeof(double))-1;
  stodftCoefPos->moUpRe = NULL;
  stodftCoefPos->moUpIm = NULL;

/*===========================================================================*/
/* II) MPI things  */

  stodftInfo->stowfRecvCounts = (int*)cmalloc(numProcStates*sizeof(int));
  stodftInfo->stowfDispls = (int*)cmalloc(numProcStates*sizeof(int));
  stodftInfo->stowfRecvCountsComplex = (int*)cmalloc(numProcStates*sizeof(int));
  stodftInfo->stowfDisplsComplex = (int*)cmalloc(numProcStates*sizeof(int));
  stodftInfo->stowfRecvCountsComplex2 = (int*)cmalloc(numProcStates*sizeof(int));
  stodftInfo->stowfDisplsComplex2 = (int*)cmalloc(numProcStates*sizeof(int));

  stodftInfo->numStates = (int*)cmalloc(numProcStates*sizeof(int));
  stodftInfo->dsplStates = (int*)cmalloc(numProcStates*sizeof(int));
  stodftInfo->numStates2 = (int*)cmalloc(numProcStates*sizeof(int));
  stodftInfo->dsplStates2 = (int*)cmalloc(numProcStates*sizeof(int));


  stowfRecvCounts = stodftInfo->stowfRecvCounts;
  stowfDispls = stodftInfo->stowfDispls;
  numStates = stodftInfo->numStates;
  dsplStates = stodftInfo->dsplStates;
  if(numProcStates>1)Barrier(comm_states);

  if(numProcStates>1){
    Allgather(&numCoeffUpTotal,1,MPI_INT,stowfRecvCounts,1,MPI_INT,0,comm_states); 
    //Allgather(&numCoeffUpAllChemPot,1,MPI_INT,stowfRecvCounts,1,MPI_INT,0,comm_states);
  }
  //Gather(&numCoeffUpAllProc,1,MPI_INT,stowfRecvCounts,numProcStates,MPI_INT,0,comm_states); 
  //Bcast(stowfRecvCounts,numProcStates,MPI_INT,0,comm_states);
  if(numProcStates>1)Barrier(comm_states);

  stowfDispls[0] = 0;
  for(iProc=1;iProc<numProcStates;iProc++){
    stowfDispls[iProc] = stowfDispls[iProc-1]+stowfRecvCounts[iProc-1];
  }
  for(iProc=0;iProc<numProcStates;iProc++){
    stodftInfo->stowfRecvCountsComplex[iProc] = stowfRecvCounts[iProc]*2;
    stodftInfo->stowfDisplsComplex[iProc] = stowfDispls[iProc]*2;
  }

  if(myidState==0){
    stodftInfo->numOccDetAll = (double*)cmalloc(numStateUpAllProc*sizeof(double));
    for(iState=0;iState<numStateUpAllProc;iState++){
      if(iState<numOccState)stodftInfo->numOccDetAll[iState] = pre;
      else stodftInfo->numOccDetAll[iState] = 0;
    }
  }
  if(numProcStates>1){
    Allgather(&numStateProcTot,1,MPI_INT,numStates,1,MPI_INT,0,comm_states); 
  }
  dsplStates[0] = 0;
  for(iProc=1;iProc<numProcStates;iProc++){
    dsplStates[iProc] = dsplStates[iProc-1]+numStates[iProc-1];
  }

  if(myidState==0){
    stodftInfo->stateStartIndex = 0;
    stateIndex = 0;
    for(iProc=1;iProc<numProcStates;iProc++){
      stateIndex += numStates[iProc];
      Send(&stateIndex,1,MPI_INT,iProc,iProc,comm_states);
    }
  }
  if(myidState!=0){ 
    Recv(&(stodftInfo->stateStartIndex),1,MPI_INT,0,myidState,comm_states);
  }
  if(numProcStates>1)Barrier(comm_states);
  stateIndex = stodftInfo->stateStartIndex;
  stodftInfo->numOccDetProc = (double*)cmalloc(numStateProcTot*sizeof(double));
  for(iState=0;iState<numStateProcTot;iState++){
    if(iState+stateIndex<numOccState)stodftInfo->numOccDetProc[iState] = pre;
    else stodftInfo->numOccDetProc[iState] = 0;
  }
  if(myidState==0)printf("Finish Initializing Filter Diagonalization!\n");

  //exit(0);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/




