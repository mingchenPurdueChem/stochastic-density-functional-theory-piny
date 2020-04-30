/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: min-CP-stodft.c                              */
/*                                                                          */
/* Subroutine for SCF calculation			                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_vps_params_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_stodft_local.h"

//#define EW_FAKE
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
#ifdef FAST_FILTER
void scfStodftInterp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                    CP *cp,
                    CLASS *class2,BONDED *bonded2,GENERAL_DATA *general_data2,CP *cp2,
                    int ip_now)
#else
void scfStodftInterp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                    CP *cp,int ip_now)
#endif
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  CELL *cell			 = &(general_data->cell);  
  PTENS *ptens			 = &(general_data->ptens);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  
  STODFTINFO *stodftInfo         = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos   = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CPEWALD *cpewald              = &(cp->cpewald);
  CPSCR *cpscr                  = &(cp->cpscr);
  PSEUDO *pseudo                = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal	= &(pseudo->pseudoReal);
  COMMUNICATE *communicate      = &(cp->communicate);

  FOR_SCR      *for_scr         = &(class->for_scr);
  EWD_SCR      *ewd_scr         = &(class->ewd_scr);
  CLATOMS_INFO *clatoms_info	= &(class->clatoms_info);
  CLATOMS_POS *clatoms_pos	= &(class->clatoms_pos[ip_now]);
  ATOMMAPS *atommaps		= &(class->atommaps);

  int iperd            		= cell->iperd;
  int iScf,iCell,iCoeff,iState,iChem;
  int numScf			= stodftInfo->numScf; //Need claim this in cp
  int numChemPot		= stodftInfo->numChemPot;
  int cpLsda 			= cpopts->cp_lsda;
  int cpParaOpt			= cpopts->cp_para_opt;

  int checkPerdSize 		= cpopts->icheck_perd_size;
  int checkDualSize 		= cpopts->icheck_dual_size;
  int cpDualGridOptOn 	 	= cpopts->cp_dual_grid_opt;
  int numProcStates 		= communicate->np_states; 
  int myidState 		= communicate->myid_state;
  int coefFormUp 		= cpcoeffs_pos->icoef_form_up;
  int coefOrthUp                = cpcoeffs_pos->icoef_orth_up;
  int forceCoefFormUp 		= cpcoeffs_pos->ifcoef_form_up;
  int forceCoefOrthUp           = cpcoeffs_pos->ifcoef_orth_up;
  int coefFormDn                = cpcoeffs_pos->icoef_form_dn;
  int coefOrthDn                = cpcoeffs_pos->icoef_orth_dn;
  int forceCoefFormDn           = cpcoeffs_pos->ifcoef_form_dn;
  int forceCoefOrthDn           = cpcoeffs_pos->ifcoef_orth_dn;
  int numStateUp		= cpcoeffs_info->nstate_up_proc;
  int numStateDn                = cpcoeffs_info->nstate_dn_proc;
  int numCoeff			= cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUp*numCoeff;
  int numCoeffDnTotal = numStateDn*numCoeff;
  int scfStopFlag = 0; // 0=do scf 1=stop scf
  int pseudoRealFlag = pseudoReal->pseudoRealFlag;
  MPI_Comm commStates = communicate->comm_states;
  int calcFragFlag              = stodftInfo->calcFragFlag;


  int *pcoefFormUp		     = &(cpcoeffs_pos->icoef_form_up);
  int *pcoefOrthUp		     = &(cpcoeffs_pos->icoef_orth_up);
  int *pforceCoefFormUp              = &(cpcoeffs_pos->ifcoef_form_up);
  int *pforceCoefOrthUp              = &(cpcoeffs_pos->ifcoef_orth_up);
  int *pcoefFormDn		     = &(cpcoeffs_pos->icoef_form_dn);
  int *pcoefOrthDn		     = &(cpcoeffs_pos->icoef_orth_dn);
  int *pforceCoefFormDn              = &(cpcoeffs_pos->icoef_form_dn);
  int *pforceCoefOrthDn              = &(cpcoeffs_pos->icoef_orth_dn);


  double tolEdgeDist 		= cpopts->tol_edge_dist;
  double energyDiff		= -1.0;
  double energyTol		= stodftInfo->energyTol;
  double timeStart,timeEnd;

  double *coeffReUp        = cpcoeffs_pos->cre_up;
  double *coeffImUp        = cpcoeffs_pos->cim_up;
  double *forceCoeffReUp   = cpcoeffs_pos->fcre_up;
  double *forceCoeffImUp   = cpcoeffs_pos->fcre_up;
  double *coeffReDn        = cpcoeffs_pos->cre_dn;
  double *coeffImDn        = cpcoeffs_pos->cim_dn;
  double *forceCoeffReDn   = cpcoeffs_pos->fcre_dn;
  double *forceCoeffImDn   = cpcoeffs_pos->fcre_dn;
  //double *cpScrCoeffReUp   = cpscr->cpscr_wave.cre_up;
  //double *cpScrCoeffImUp   = cpscr->cpscr_wave.cim_up;
  //double *cpScrCoeffReDn   = cpscr->cpscr_wave.cre_dn;
  //double *cpScrCoeffImDn   = cpscr->cpscr_wave.cim_dn;
  double *ptensPvtenTmp    = ptens->pvten_tmp;
  double *chemPot          = stodftCoefPos->chemPot;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;


/*======================================================================*/
/* 0.05) Check the approximations in the methods                        */


  if((iperd<3||iperd==4)&&checkPerdSize==1){
    cp_boundary_check(cell,clatoms_info,clatoms_pos,tolEdgeDist);
  }/*endif*/
  if(cpDualGridOptOn>=1&&checkDualSize==1){
    cp_dual_check(cell,clatoms_info,clatoms_pos,
                  atommaps->cp_atm_lst,tolEdgeDist);
  }/*endif*/

/*======================================================================*/
/* I) Orthogonalize the coefs if norbing                                */

/*======================================================================*/
/* II) In parallel, transpose coefs back to normal form                 */

  /*
  // I do this once in init.c, in case we read in deterministic wf. For 
  // stochastic wf, we don't need this
  if(numProcStates>1){
    cp_transpose_bck(coeffReUp,coeffImUp,pcoefFormUp,
                    cpScrCoeffReUp,cpScrCoeffImUp,&(cp->cp_comm_state_pkg_up));
    if(cpLsda==1&&numStateDn>0){
      cp_transpose_bck(coeffReDn,coeffImDn,pcoefFormDn,
                     cpScrCoeffReDn,cpScrCoeffImDn,&(cp->cp_comm_state_pkg_dn));
    }//endif
  }//endif
  */

/*======================================================================*/
/* III) Initialize forces, pressure tensor, inverse hmat                */

  cpcoeffs_pos->ifcoef_form_up = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1;
  if(cpLsda==1&&numStateDn>0){
    cpcoeffs_pos->ifcoef_form_dn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }/*endif*/

  for(iCell=1;iCell<=9;iCell++){ptensPvtenTmp[iCell] = 0.0;}
  gethinv(cell->hmat_cp,cell->hmati_cp,&(cell->vol_cp),iperd);
  gethinv(cell->hmat,cell->hmati,&(cell->vol),iperd);

/*======================================================================*/
/* IV) Check the forms                                                   */

/*======================================================================*/
/* V) Initial KS potential calculation			                */

  stat_avg->cp_ehart = 0.0;
  stat_avg->cp_eext = 0.0;
  stat_avg->cp_exc = 0.0;
  stodftInfo->energyElecTot = 0.0;
  stodftInfo->energyElecTotOld = 0.0;
  if(myidState==0)printf("**Calculating Initial Kohn-Sham Potential...\n");
  calcLocalPseudoScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  calcKSPot(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  if(pseudoRealFlag==1){
    pseudoReal->forceCalcFlag = 1;
    initRealNlppWf(cp,class,general_data);
    allocRealNl(cp,class);
    pseudoReal->forceCalcFlag = 0;
  }

  if(myidState==0)printf("**Finish Calculating Initial Kohn-Sham Potential\n");

#ifdef FAST_FILTER
  cp2->cpcoeffs_pos[1].ifcoef_form_up = 0;
  cp2->cpcoeffs_pos[1].ifcoef_orth_up = 1;

  gethinv(general_data2->cell.hmat_cp,general_data2->cell.hmati_cp,
          &(general_data2->cell.vol_cp),iperd);
  gethinv(general_data2->cell.hmat,general_data2->cell.hmati,
          &(general_data2->cell.vol),iperd);
  if(calcFragFlag==0){
    calcLocalPseudoScf(class2,general_data2,cp2,
                       &(cp2->cpcoeffs_pos[1]),&(class2->clatoms_pos[1]));
  }

  if(cp2->pseudo.pseudoReal.pseudoRealFlag==1){
    cp2->pseudo.pseudoReal.forceCalcFlag = 1;
    initRealNlppWf(cp2,class2,general_data2);
    allocRealNl(cp2,class2);
    cp2->pseudo.pseudoReal.forceCalcFlag = 0;
  }
  cp->stodftInfo->numStatePrintUp = cp2->stodftInfo->numStatePrintUp;
  cp->stodftCoefPos->moUpRePrint = (double*)cmalloc(cp->stodftInfo->numStatePrintUp*numCoeff*sizeof(double))-1;
  cp->stodftCoefPos->moUpImPrint = (double*)cmalloc(cp->stodftInfo->numStatePrintUp*numCoeff*sizeof(double))-1;
  // attention: size of cp->stodftCoefPos.energyLevel is different from 
  // the size of cp2->stodftCoefPos.energyLevel
  cp->stodftCoefPos->energyLevel = (double*)cmalloc(cp->stodftInfo->numStatePrintUp*sizeof(double));

#endif

  
  //exit(0);
/*======================================================================*/
/* V) SCF loop						                */


  if(myidState==0){
    printf("===============================================================================\n");
    printf("Runing SCF Calculation\n");
    printf("-------------------------------------------------------------------------------\n");
  }
  //for(iScf=1;iScf<=numScf;iScf++){
  iScf = 0;
  while(scfStopFlag==0){ 
    timeStart = omp_get_wtime();
    iScf += 1;
    stodftInfo->iScf = iScf;
    if(myidState==0){
      printf("********************************************************\n");
      printf("SCF Step %i\n",iScf);
      printf("--------------------------------------------------------\n");
    }

/*----------------------------------------------------------------------*/
/* i) Generate KS potential                                             */    
   
    /*
    stat_avg->cp_ehart = 0.0;
    stat_avg->cp_eext = 0.0;
    stat_avg->cp_exc = 0.0;
    if(myidState==0)printf("**Calculating Kohn-Sham Potential...\n");
    calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    calcKSForceControlWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

    if(myidState==0)printf("**Finish Calculating Kohn-Sham Potential\n");
    */

/*----------------------------------------------------------------------*/
/* i) Generate stochastic WF for different chemical potentials          */

    if(myidState==0)printf("**Generating Stochastic Orbitals...\n");
#ifdef FAST_FILTER
    genStoOrbitalInterpFake(class,general_data,cp,class2,general_data2,cp2,ip_now);
#else
    genStoOrbitalInterp(class,general_data,cp,ip_now);
#endif
 
    if(stodftInfo->filterDiagFlag==1){
      orthDiagDriver(cp,class,general_data,ip_now);
    }
    //exit(0);
    //debug
    /*
    for(iChem=0;iChem<numChemPot;iChem++){
      if(checkNanArray(&stoWfUpRe[iChem][1],numCoeffUpTotal)==1){
	printf("iChem %i myid %i Bad number in upper sate real part!\n",iChem,myidState);
      }
      if(checkNanArray(&stoWfUpIm[iChem][1],numCoeffUpTotal)==1){
        printf("iChem %i myid %i Bad number in upper sate imag part!\n",iChem,myidState);
      }
    }
    */
    if(myidState==0)printf("**Finish Generating Stochastic Orbitals\n");

    //exit(0);   
    
    /*   
    char wfname[100];
    sprintf(wfname,"sto-wf-save-%i",myidState);

    FILE *filePrintWF = fopen(wfname,"r");
    for(iChem=0;iChem<numChemPot;iChem++){
      for(iState=0;iState<numStateUp;iState++){
	for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
	  fprintf(filePrintWF,"%.16lg %.16lg",&stoWfUpRe[iChem][iState*numCoeff+iCoeff],
		  &stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
	  //fscanf(filePrintWF,"%lg",&stoWfUpRe[iChem][iState*numCoeff+iCoeff]);
	  //fscanf(filePrintWF,"%lg",&stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
	}//endfor iCoeff
      }//endfor iState
    }//endfor iChem
    fclose(filePrintWF);
    printf("myid %i finish reading in WF.\n",myidState);
    if(numProcStates>1)Barrier(commStates);
    exit(0);
    */
    
    /*
    printf("Start Readin WF\n");

    char wfname[100];
    
    FILE *filePrintWF;
    for(iState=0;iState<numStateUp;iState++){
      printf("iState %i\n",iState);
      sprintf(wfname,"./all-data/sto-wf-save-%i",iState);
      filePrintWF = fopen(wfname,"r");
      for(iChem=0;iChem<numChemPot;iChem++){
	for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
	  fscanf(filePrintWF,"%lg",&stoWfUpRe[iChem][iState*numCoeff+iCoeff]);
	  fscanf(filePrintWF,"%lg",&stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
	}//endif iCoeff
      }//endfor iChem
      fclose(filePrintWF);
    }//endfor iState
    printf("Finish Readin WF\n");
    */
    
    if(myidState==0)printf("**Calculating KE and NLPPE...\n");
    calcEnergyChemPot(cp,class,general_data,cpcoeffs_pos,clatoms_pos);   
    if(myidState==0)printf("**Finish Calculating KE and NLP E\n");

    
//debug 
    /*
    double norm;    
    double repart,impart;
    for(iState=0;iState<numStateUp;iState++){
      norm = 0.0;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	repart = stoWfUpRe[0][iState*numCoeff+iCoeff];
	impart = stoWfUpIm[0][iState*numCoeff+iCoeff];
	norm += repart*repart+impart*impart;
      }
      norm *= 2.0;
      norm += stoWfUpRe[0][iState*numCoeff+numCoeff]*stoWfUpRe[0][iState*numCoeff+numCoeff];
      printf("iState %i norm %lg\n",iState,norm);
    }
    */
    

/*----------------------------------------------------------------------*/
/* ii)  Get the total density, for each chemical potential and get      */
/*     total electron number for each chemical potential	        */
/*     Interpolate the chemical potential w.r.t                         */
/*     Correct electron number and use the interpolation coefficients   */
/*     to generate the density w.r.t. correct number of electrons.      */

    if(myidState==0)printf("**Calculating Density...\n");
    if(cpParaOpt==0) calcRhoStoHybridInterp(class,bonded,general_data,cp,ip_now); 
    if(myidState==0)printf("**Finish Calculating Density\n");

/*----------------------------------------------------------------------*/
/* iii) Rerun if necessary	                                        */

    while(stodftInfo->reRunFlag==1){
      if(myidState==0)printf("**Generating Stochastic Orbitals...\n");
      adjChemPot(stodftInfo,stodftCoefPos);
      genStoOrbitalInterp(class,general_data,cp,ip_now);
      if(myidState==0)printf("**Finish Generating Stochastic Orbitals\n");
      if(myidState==0)printf("**Calculating KE and NLPPE...\n");
      calcEnergyChemPot(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
      if(myidState==0)printf("**Finish Calculating KE and NLP E\n");
      if(myidState==0)printf("**Calculating Density...\n");
      if(cpParaOpt==0) calcRhoStoHybridInterp(class,bonded,general_data,cp,ip_now);
      if(myidState==0)printf("**Finish Calculating Density\n");
    }//endwhile reRunFlag

/*----------------------------------------------------------------------*/
/* iv) Generate KS potential                                            */

    stat_avg->cp_ehart = 0.0;
    stat_avg->cp_eext = 0.0;
    stat_avg->cp_exc = 0.0;
    if(myidState==0)printf("**Calculating Kohn-Sham Potential...\n");
    //calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    calcKSPot(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

    if(myidState==0)printf("**Finish Calculating Kohn-Sham Potential\n");

/*----------------------------------------------------------------------*/
/* v) Calculate the total energy		                        */


    if(myidState==0)printf("**Calculating Total Energy...\n");
    stodftInfo->energyElecTotOld = stodftInfo->energyElecTot;
    calcTotEnergy(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
    energyDiff = fabs(stodftInfo->energyElecTot-stodftInfo->energyElecTotOld);
    if(myidState==0)printf("**Finish Calculating Total Energy\n");

    //exit(0);   

/*----------------------------------------------------------------------*/
/* iv) Generate chemical potentials for the next step                   */

    if(myidState==0)printf("**Prepare Chemical Potentials for next SCF...\n");
    updateChemPot(stodftInfo,stodftCoefPos);
    if(myidState==0){
      printf("The min chem pot for next step is %.6lg\n",chemPot[0]);
      printf("The max chem pot for next step is %.6lg\n",chemPot[numChemPot-1]);
    }
    if(myidState==0)printf("**Finish Preparing Chemical Potentials for next SCF\n");
    
    
/*----------------------------------------------------------------------*/
/* v) Finish this SCF step			                        */

    timeEnd = omp_get_wtime();

    if(myidState==0){
      printf("The master process spend %lgs in this SCF step.\n",timeEnd-timeStart);
      printf("--------------------------------------------------------\n");
      printf("Finish SCF Step %i\n",iScf);
      printf("********************************************************\n");
      printf("\n");
    }
    if(energyDiff<energyTol||iScf>=numScf)scfStopFlag = 1;
    //exit(0);
  }//endwhile scfStopFlag

/*======================================================================*/
/* VI) Calcualte energy wih nuclei forces                 		*/

  

  

/*======================================================================*/
/* VI) In parallel, transpose coefs and coef forces fwd                 */

  

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
#ifdef FAST_FILTER
void scfStodftCheby(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                    CP *cp,
                    CLASS *class2,BONDED *bonded2,GENERAL_DATA *general_data2,
                    CP *cp2,
                    int ip_now)   
#else
void scfStodftCheby(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                    CP *cp,int ip_now)
#endif
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  CELL *cell			 = &(general_data->cell);  
  PTENS *ptens			 = &(general_data->ptens);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  
  STODFTINFO *stodftInfo         = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos   = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CPEWALD *cpewald              = &(cp->cpewald);
  CPSCR *cpscr                  = &(cp->cpscr);
  PSEUDO *pseudo                = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal	= &(pseudo->pseudoReal);
  COMMUNICATE *communicate      = &(cp->communicate);

  FOR_SCR      *for_scr         = &(class->for_scr);
  EWD_SCR      *ewd_scr         = &(class->ewd_scr);
  CLATOMS_INFO *clatoms_info	= &(class->clatoms_info);
  CLATOMS_POS *clatoms_pos	= &(class->clatoms_pos[ip_now]);
  ATOMMAPS *atommaps		= &(class->atommaps);

  int iperd			= cell->iperd;
  int iScf,iCell,iCoeff,iState,iChem;
  int numScf			= stodftInfo->numScf; //Need claim this in cp
  int numChemPot		= stodftInfo->numChemPot;
  int cpLsda			= cpopts->cp_lsda;
  int cpParaOpt			= cpopts->cp_para_opt;
  int pseudoRealFlag		= pseudoReal->pseudoRealFlag;
  int checkPerdSize		= cpopts->icheck_perd_size;
  int checkDualSize		= cpopts->icheck_dual_size;
  int cpDualGridOptOn		= cpopts->cp_dual_grid_opt;
  int numProcStates		= communicate->np_states; 
  int myidState			= communicate->myid_state;
  int coefFormUp		= cpcoeffs_pos->icoef_form_up;
  int coefOrthUp                = cpcoeffs_pos->icoef_orth_up;
  int forceCoefFormUp		= cpcoeffs_pos->ifcoef_form_up;
  int forceCoefOrthUp           = cpcoeffs_pos->ifcoef_orth_up;
  int coefFormDn                = cpcoeffs_pos->icoef_form_dn;
  int coefOrthDn                = cpcoeffs_pos->icoef_orth_dn;
  int forceCoefFormDn           = cpcoeffs_pos->ifcoef_form_dn;
  int forceCoefOrthDn           = cpcoeffs_pos->ifcoef_orth_dn;
  int numStateUp		= cpcoeffs_info->nstate_up_proc;
  int numStateDn                = cpcoeffs_info->nstate_dn_proc;
  int numCoeff			= cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUp*numCoeff;
  int numCoeffDnTotal = numStateDn*numCoeff;
  int scfStopFlag     = 0;
  int checkpointWriteFreq = stodftInfo->checkpointWriteFreq;
  int checkpointParFlag = stodftInfo->checkpointParFlag;
  int readCoeffFlag = stodftInfo->readCoeffFlag;
  MPI_Comm commStates   =    communicate->comm_states;

  int *pcoefFormUp		     = &(cpcoeffs_pos->icoef_form_up);
  int *pcoefOrthUp		     = &(cpcoeffs_pos->icoef_orth_up);
  int *pforceCoefFormUp              = &(cpcoeffs_pos->ifcoef_form_up);
  int *pforceCoefOrthUp              = &(cpcoeffs_pos->ifcoef_orth_up);
  int *pcoefFormDn		     = &(cpcoeffs_pos->icoef_form_dn);
  int *pcoefOrthDn		     = &(cpcoeffs_pos->icoef_orth_dn);
  int *pforceCoefFormDn              = &(cpcoeffs_pos->icoef_form_dn);
  int *pforceCoefOrthDn              = &(cpcoeffs_pos->icoef_orth_dn);


  double numElecTrue = stodftInfo->numElecTrue;
  double tolEdgeDist		= cpopts->tol_edge_dist;
  double energyDiff		= -1.0;
  double energyTol		= stodftInfo->energyTol*numElecTrue;
  double timeStart,timeEnd;

  double *coeffReUp        = cpcoeffs_pos->cre_up;
  double *coeffImUp        = cpcoeffs_pos->cim_up;
  double *forceCoeffReUp   = cpcoeffs_pos->fcre_up;
  double *forceCoeffImUp   = cpcoeffs_pos->fcre_up;
  double *coeffReDn        = cpcoeffs_pos->cre_dn;
  double *coeffImDn        = cpcoeffs_pos->cim_dn;
  double *forceCoeffReDn   = cpcoeffs_pos->fcre_dn;
  double *forceCoeffImDn   = cpcoeffs_pos->fcre_dn;
  //double *cpScrCoeffReUp   = cpscr->cpscr_wave.cre_up;
  //double *cpScrCoeffImUp   = cpscr->cpscr_wave.cim_up;
  //double *cpScrCoeffReDn   = cpscr->cpscr_wave.cre_dn;
  //double *cpScrCoeffImDn   = cpscr->cpscr_wave.cim_dn;
  double *ptensPvtenTmp    = ptens->pvten_tmp;
  double *chemPot          = stodftCoefPos->chemPot;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  //debug
  FILE *fileRhoRecip;
  int numCoeffLargeProc = cp->cp_para_fft_pkg3d_lg.ncoef_proc;
  double *rhoCoeffReUp = cpscr->cpscr_rho.rhocr_up;
  double *rhoCoeffImUp = cpscr->cpscr_rho.rhoci_up;

  //timeing
  double timeStart1,timeEnd1;
  double timeStart2,timeEnd2;
  double timeStart3,timeEnd3;
  double timeStart4,timeEnd4;
  double timeStart5,timeEnd5;
  double timeStart6,timeEnd6;
  double diffTime1 = 0.0;
  double diffTime2 = 0.0;
  double diffTime3 = 0.0;
  double diffTime4 = 0.0;
  double diffTime5 = 0.0;
  double diffTime6 = 0.0;

/*======================================================================*/
/* I) Check the approximations in the methods				*/

  if((iperd<3||iperd==4)&&checkPerdSize==1){
    cp_boundary_check(cell,clatoms_info,clatoms_pos,tolEdgeDist);
  }/*endif*/
  if(cpDualGridOptOn>=1&&checkDualSize==1){
    cp_dual_check(cell,clatoms_info,clatoms_pos,
                  atommaps->cp_atm_lst,tolEdgeDist);
  }/*endif*/

/*======================================================================*/
/* II) Initialize forces, pressure tensor, inverse hmat                 */

  cpcoeffs_pos->ifcoef_form_up = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1;
  if(cpLsda==1&&numStateDn>0){
    cpcoeffs_pos->ifcoef_form_dn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }/*endif*/

  for(iCell=1;iCell<=9;iCell++){ptensPvtenTmp[iCell] = 0.0;}
  gethinv(cell->hmat_cp,cell->hmati_cp,&(cell->vol_cp),iperd);
  gethinv(cell->hmat,cell->hmati,&(cell->vol),iperd);

/*======================================================================*/
/* III) Initial KS potential calculation			        */
/*	Initialize real sapce nlpp*/

  stat_avg->cp_ehart = 0.0;
  stat_avg->cp_eext = 0.0;
  stat_avg->cp_exc = 0.0;
  if(myidState==0)printf("**Calculating Initial Kohn-Sham Potential...\n");
  calcLocalPseudoScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  calcKSPot(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

  if(pseudoRealFlag==1){
    if(myidState==0)printf("**Calculating Real Space Non-local Pseudopotential...\n");
    pseudoReal->forceCalcFlag = 1;
    initRealNlppWf(cp,class,general_data);
    allocRealNl(cp,class);
    pseudoReal->forceCalcFlag = 0;
  }

  if(myidState==0)printf("**Finish Calculating Initial Kohn-Sham Potential\n");

  //exit(0);
/*======================================================================*/
/* IV) SCF loop						                */
#ifdef FAST_FILTER  
  cp2->cpcoeffs_pos[1].ifcoef_form_up = 0;
  cp2->cpcoeffs_pos[1].ifcoef_orth_up = 1;

  gethinv(general_data2->cell.hmat_cp,general_data2->cell.hmati_cp,
          &(general_data2->cell.vol_cp),iperd);
  gethinv(general_data2->cell.hmat,general_data2->cell.hmati,
          &(general_data2->cell.vol),iperd);
  calcLocalPseudoScf(class2,general_data2,cp2,
                     &(cp2->cpcoeffs_pos[1]),&(class2->clatoms_pos[1]));

  if(cp2->pseudo.pseudoReal.pseudoRealFlag==1){
    cp2->pseudo.pseudoReal.forceCalcFlag = 1;
    initRealNlppWf(cp2,class2,general_data2);
    allocRealNl(cp2,class2);
    cp2->pseudo.pseudoReal.forceCalcFlag = 0;
  }

  cp->stodftInfo->numStatePrintUp = cp2->stodftInfo->numStatePrintUp;
  cp->stodftCoefPos->moUpRePrint = (double*)cmalloc(cp->stodftInfo->numStatePrintUp*numCoeff*sizeof(double))-1;
  cp->stodftCoefPos->moUpImPrint = (double*)cmalloc(cp->stodftInfo->numStatePrintUp*numCoeff*sizeof(double))-1;
  // attention: size of cp->stodftCoefPos.energyLevel is different from 
  // the size of cp2->stodftCoefPos.energyLevel
  cp->stodftCoefPos->energyLevel = (double*)cmalloc(cp->stodftInfo->numStatePrintUp*sizeof(double));
#endif

  if(myidState==0){
    printf("===============================================================================\n");
    printf("Runing SCF Calculation\n");
    printf("-------------------------------------------------------------------------------\n");
  }
  //printf("numStateUp %i\n",numStateUp);
  
  //for(iScf=1;iScf<=numScf;iScf++){

  if(readCoeffFlag==-3){
    iScf = stodftInfo->iScf;
  }
  else{
    iScf = 0;
  }
  while(scfStopFlag==0){
    timeStart = omp_get_wtime();
    iScf += 1;
    stodftInfo->iScf = iScf;

    if(myidState==0){
      printf("********************************************************\n");
      printf("SCF Step %i\n",iScf);
      printf("--------------------------------------------------------\n");
    }
/*----------------------------------------------------------------------*/
/* i) Generate stochastic WF for different chemical potentials          */

    timeStart1 = omp_get_wtime();
    if(myidState==0)printf("**Generating Stochastic Orbitals...\n");
#ifdef FAST_FILTER  
    genStoOrbitalChebyFakeNew(class,general_data,cp,class2,general_data2,cp2,ip_now);
#else
    //genStoOrbitalCheby(class,general_data,cp,ip_now);
#endif
    genStoOrbitalFake(class,general_data,cp,ip_now);
 
    //exit(0);
    //debug
    /*
    for(iChem=0;iChem<numChemPot;iChem++){
      if(checkNanArray(&stoWfUpRe[iChem][1],numCoeffUpTotal)==1){
	printf("iChem %i myid %i Bad number in upper sate real part!\n",iChem,myidState);
      }
      if(checkNanArray(&stoWfUpIm[iChem][1],numCoeffUpTotal)==1){
        printf("iChem %i myid %i Bad number in upper sate imag part!\n",iChem,myidState);
      }
    }
    */
    if(myidState==0)printf("**Finish Generating Stochastic Orbitals\n");
    timeEnd1 = omp_get_wtime();
    diffTime1 += timeEnd1-timeStart1;

    //exit(0);
    
    /*     
    char wfname[100];
    //sprintf(wfname,"/scratch/mingchen/tmp/sto-wf-save-%i",myidState);
    printf("Read in stochastic orbitals...\n");
    sprintf(wfname,"sto-wf-save-new-%i",myidState);
    
    
    FILE *filePrintWF = fopen(wfname,"r");
    //FILE *filePrintWF = fopen("sto-wf-save-all-2","r");
    for(iChem=0;iChem<numChemPot;iChem++){
      for(iState=0;iState<numStateUp;iState++){
	for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
	  fscanf(filePrintWF,"%lg",&stoWfUpRe[iChem][iState*numCoeff+iCoeff]);
	  fscanf(filePrintWF,"%lg",&stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
	  //fprintf(filePrintWF,"%.16lg %.16lg\n",
	  //	    stoWfUpRe[iChem][iState*numCoeff+iCoeff],
	  //	    stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
	}//endfor iCoeff
      }//endfor iState
    }//endfor iChem
    fclose(filePrintWF);
    printf("myid %i finish reading in WF.\n",myidState);
    printf("%lg %lg\n",stoWfUpRe[0][1],stoWfUpIm[0][1]);
    fflush(stdout);
    Barrier(commStates);
    //exit(0);
    */
    
    /*
    printf("Start Readin WF\n");

    char wfname[100];
    
    FILE *filePrintWF;
    numChemPot = stodftInfo->numChemPot;
    for(iState=0;iState<numStateUp;iState++){
      //printf("iState %i\n",iState);
      //sprintf(wfname,"./all-data/sto-wf-save-%i",iState);
      //filePrintWF = fopen(wfname,"r");
      filePrintWF = fopen("sto-wf-save","r");
      for(iChem=0;iChem<numChemPot;iChem++){
	for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
	  fscanf(filePrintWF,"%lg",&stoWfUpRe[iChem][iState*numCoeff+iCoeff]);
	  fscanf(filePrintWF,"%lg",&stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
	}//endif iCoeff
      }//endfor iChem
      fclose(filePrintWF);
    }//endfor iState
    printf("Finish Readin WF\n");
    */
    
    timeStart2 = omp_get_wtime();
    if(myidState==0)printf("**Calculating KE and NLPPE...\n");
    calcEnergyChemPot(cp,class,general_data,cpcoeffs_pos,clatoms_pos);   
    if(myidState==0)printf("**Finish Calculating KE and NLP E\n");
    timeEnd2 = omp_get_wtime();
    diffTime2 += timeEnd2-timeStart2;
    

    
//debug 
    /*
    double norm;    
    double repart,impart;
    for(iState=0;iState<numStateUp;iState++){
      norm = 0.0;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	repart = stoWfUpRe[0][iState*numCoeff+iCoeff];
	impart = stoWfUpIm[0][iState*numCoeff+iCoeff];
	norm += repart*repart+impart*impart;
      }
      norm *= 2.0;
      norm += stoWfUpRe[0][iState*numCoeff+numCoeff]*stoWfUpRe[0][iState*numCoeff+numCoeff];
      printf("iState %i norm %lg\n",iState,norm);
    }
    */
    

/*----------------------------------------------------------------------*/
/* ii)  Get the total density.						*/

    timeStart3 = omp_get_wtime();
    if(myidState==0)printf("**Calculating Density...\n");
    if(cpParaOpt==0) calcRhoStoHybridCheby(class,bonded,general_data,cp,ip_now); 
    if(myidState==0)printf("**Finish Calculating Density\n");
    timeEnd3 = omp_get_wtime();
    diffTime3 += timeEnd3-timeStart3;

/*----------------------------------------------------------------------*/
/* iv) Generate KS potential                                            */
    
    if(numProcStates>1)Barrier(commStates);

    timeStart4 = omp_get_wtime();
    stat_avg->cp_ehart = 0.0;
    stat_avg->cp_eext = 0.0;
    stat_avg->cp_exc = 0.0;
    if(myidState==0)printf("**Calculating Kohn-Sham Potential...\n");
    //calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    calcKSPot(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    //debug
    /*
    if(myidState==0){
      int nfft = cp->cp_sclr_fft_pkg3d_lg.nfft;
      int nfft2 = nfft/2;
      int iGrid;
      for(iGrid=1;iGrid<=nfft2;iGrid++){
	printf("vkssssss %lg\n",cp->cpscr.cpscr_rho.v_ks_up[iGrid]);
      }
    }
    */
    if(myidState==0)printf("**Finish Calculating Kohn-Sham Potential\n");
    timeEnd4 = omp_get_wtime();
    diffTime4 += timeEnd4-timeStart4;

/*----------------------------------------------------------------------*/
/* v) Calculate the total energy		                        */


    timeStart5 = omp_get_wtime();
    if(myidState==0)printf("**Calculating Total Energy...\n");
    stodftInfo->energyElecTotOld = stodftInfo->energyElecTot;
    calcTotEnergy(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
    energyDiff = fabs(stodftInfo->energyElecTot-stodftInfo->energyElecTotOld);
    if(myidState==0)printf("**Finish Calculating Total Energy\n");
    timeEnd5 = omp_get_wtime();
    diffTime5 += timeEnd5-timeStart5;


    //exit(0);   
/*----------------------------------------------------------------------*/
/* v) Write checkpoint file if necessary                                */

    if(iScf%checkpointWriteFreq==0){
      if(checkpointParFlag==0)checkpointOutput(cp,general_data);
      else checkpointOutputDist(cp,general_data);
    }

/*----------------------------------------------------------------------*/
/* v) Finish this SCF step			                        */

    timeEnd = omp_get_wtime();
    if(myidState==0){
      printf("**The master process spend %lgs in this SCF step.\n",timeEnd-timeStart);
      printf("--------------------------------------------------------\n");
      printf("Energy difference is %.16lg Hartree\n",energyDiff);
      printf("Finish SCF Step %i\n",iScf);
      printf("********************************************************\n");
      printf("\n");
    }
    if(energyDiff<energyTol||iScf>=numScf)scfStopFlag = 1;
    //exit(0);
  }//endfor iScf

  if(numProcStates>1)Barrier(commStates);

  if(myidState==0){
    printf("SCF time myid %i gen-stowf %.8lg energy %.8lg density %.8lg KS potential %.8lg total-energy %.8lg\n",myidState,diffTime1,diffTime2,diffTime3,diffTime4,diffTime5);
  }

  /*  
  char wfname[100];
  //sprintf(wfname,"/scratch/mingchen/tmp/sto-wf-save-%i",myidState);
  printf("Read in stochastic orbitals...\n");
  sprintf(wfname,"sto-wf-save-%i",myidState);

  FILE *filePrintWF = fopen(wfname,"w");
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iState=0;iState<numStateUp;iState++){
      for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
	//fscanf(filePrintWF,"%lg",&stoWfUpRe[iChem][iState*numCoeff+iCoeff]);
	//fscanf(filePrintWF,"%lg",&stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
	fprintf(filePrintWF,"%.16lg %.16lg\n",stoWfUpRe[iChem][iState*numCoeff+iCoeff],
	        stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
      }//endfor iCoeff
    }//endfor iState
  }//endfor iChem
  fclose(filePrintWF);
 
  FILE *fp_rhok = fopen("rho_conv","w");
  int ncoef_l = cp->cp_para_fft_pkg3d_lg.ncoef;
  for(iCoeff=1;iCoeff<=ncoef_l;iCoeff++){
    //fscanf(fp_rhok,"%lg",&(cp->cpscr.cpscr_rho.rhocr_up[iCoeff]));
    //fscanf(fp_rhok,"%lg",&(cp->cpscr.cpscr_rho.rhoci_up[iCoeff]));
    fprintf(fp_rhok,"%.16lg %.16lg\n",cp->cpscr.cpscr_rho.rhocr_up[iCoeff],cp->cpscr.cpscr_rho.rhoci_up[iCoeff]);
  }
  fclose(fp_rhok);
  */
/*======================================================================*/
/* VI) Calculate nuclei forces after SCF loop	                        */

  calcEnergyForce(class,general_data,cp,bonded,cpcoeffs_pos,clatoms_pos);

/*======================================================================*/
/* VI) In parallel, transpose coefs and coef forces fwd                 */


/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void scfStodftFilterDiag(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                    CP *cp,int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  CELL *cell		 = &(general_data->cell);  
  PTENS *ptens		 = &(general_data->ptens);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  
  STODFTINFO *stodftInfo         = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos   = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CPEWALD *cpewald              = &(cp->cpewald);
  CPSCR *cpscr                  = &(cp->cpscr);
  PSEUDO *pseudo                = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal       = &(pseudo->pseudoReal);
  COMMUNICATE *communicate      = &(cp->communicate);

  FOR_SCR      *for_scr         = &(class->for_scr);
  EWD_SCR      *ewd_scr         = &(class->ewd_scr);
  CLATOMS_INFO *clatoms_info	= &(class->clatoms_info);
  CLATOMS_POS *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  ATOMMAPS *atommaps	    = &(class->atommaps);

  int iperd		    = cell->iperd;
  int iScf,iCell,iCoeff,iState,iChem;
  int numScf		= stodftInfo->numScf; //Need claim this in cp
  int numChemPot	= stodftInfo->numChemPot;
  int cpLsda		= cpopts->cp_lsda;
  int cpParaOpt		= cpopts->cp_para_opt;
  int pseudoRealFlag	= pseudoReal->pseudoRealFlag;

  int checkPerdSize	    = cpopts->icheck_perd_size;
  int checkDualSize	    = cpopts->icheck_dual_size;
  int cpDualGridOptOn	    = cpopts->cp_dual_grid_opt;
  int numProcStates	    = communicate->np_states; 
  int myidState		= communicate->myid_state;
  int coefFormUp	= cpcoeffs_pos->icoef_form_up;
  int coefOrthUp                = cpcoeffs_pos->icoef_orth_up;
  int forceCoefFormUp	    = cpcoeffs_pos->ifcoef_form_up;
  int forceCoefOrthUp           = cpcoeffs_pos->ifcoef_orth_up;
  int coefFormDn                = cpcoeffs_pos->icoef_form_dn;
  int coefOrthDn                = cpcoeffs_pos->icoef_orth_dn;
  int forceCoefFormDn           = cpcoeffs_pos->ifcoef_form_dn;
  int forceCoefOrthDn           = cpcoeffs_pos->ifcoef_orth_dn;
  int numStateUp	= cpcoeffs_info->nstate_up_proc;
  int numStateDn                = cpcoeffs_info->nstate_dn_proc;
  int numCoeff		= cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUp*numCoeff;
  int numCoeffDnTotal = numStateDn*numCoeff;
  int scfStopFlag     = 0;
  int checkpointWriteFreq = stodftInfo->checkpointWriteFreq;
  int checkpointParFlag = stodftInfo->checkpointParFlag;
  int readCoeffFlag = stodftInfo->readCoeffFlag;
  MPI_Comm commStates   =    communicate->comm_states;

  int *pcoefFormUp	     = &(cpcoeffs_pos->icoef_form_up);
  int *pcoefOrthUp	     = &(cpcoeffs_pos->icoef_orth_up);
  int *pforceCoefFormUp              = &(cpcoeffs_pos->ifcoef_form_up);
  int *pforceCoefOrthUp              = &(cpcoeffs_pos->ifcoef_orth_up);
  int *pcoefFormDn	     = &(cpcoeffs_pos->icoef_form_dn);
  int *pcoefOrthDn	     = &(cpcoeffs_pos->icoef_orth_dn);
  int *pforceCoefFormDn              = &(cpcoeffs_pos->icoef_form_dn);
  int *pforceCoefOrthDn              = &(cpcoeffs_pos->icoef_orth_dn);

  double numElecTrue = stodftInfo->numElecTrue;
  double tolEdgeDist	= cpopts->tol_edge_dist;
  double energyDiff	= -1.0;
  double energyTol	= stodftInfo->energyTol*numElecTrue;
  double timeStart,timeEnd;

  double *coeffReUp        = cpcoeffs_pos->cre_up;
  double *coeffImUp        = cpcoeffs_pos->cim_up;
  double *forceCoeffReUp   = cpcoeffs_pos->fcre_up;
  double *forceCoeffImUp   = cpcoeffs_pos->fcre_up;
  double *coeffReDn        = cpcoeffs_pos->cre_dn;
  double *coeffImDn        = cpcoeffs_pos->cim_dn;
  double *forceCoeffReDn   = cpcoeffs_pos->fcre_dn;
  double *forceCoeffImDn   = cpcoeffs_pos->fcre_dn;
  //double *cpScrCoeffReUp   = cpscr->cpscr_wave.cre_up;
  //double *cpScrCoeffImUp   = cpscr->cpscr_wave.cim_up;
  //double *cpScrCoeffReDn   = cpscr->cpscr_wave.cre_dn;
  //double *cpScrCoeffImDn   = cpscr->cpscr_wave.cim_dn;
  double *ptensPvtenTmp    = ptens->pvten_tmp;
  double *chemPot          = stodftCoefPos->chemPot;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  //timeing
  double timeStart1,timeEnd1;
  double timeStart2,timeEnd2;
  double timeStart3,timeEnd3;
  double timeStart4,timeEnd4;
  double timeStart5,timeEnd5;
  double timeStart6,timeEnd6;
  double diffTime1 = 0.0;
  double diffTime2 = 0.0;
  double diffTime3 = 0.0;
  double diffTime4 = 0.0;
  double diffTime5 = 0.0;
  double diffTime6 = 0.0;

/*======================================================================*/
/* 0.05) Check the approximations in the methods                        */


  if((iperd<3||iperd==4)&&checkPerdSize==1){
    cp_boundary_check(cell,clatoms_info,clatoms_pos,tolEdgeDist);
  }/*endif*/
  if(cpDualGridOptOn>=1&&checkDualSize==1){
    cp_dual_check(cell,clatoms_info,clatoms_pos,
                  atommaps->cp_atm_lst,tolEdgeDist);
  }/*endif*/

/*======================================================================*/
/* I) Orthogonalize the coefs if norbing                                */

/*======================================================================*/
/* II) In parallel, transpose coefs back to normal form                 */

/*======================================================================*/
/* III) Initialize forces, pressure tensor, inverse hmat                */

  cpcoeffs_pos->ifcoef_form_up = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1;
  if(cpLsda==1&&numStateDn>0){
    cpcoeffs_pos->ifcoef_form_dn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }/*endif*/

  for(iCell=1;iCell<=9;iCell++){ptensPvtenTmp[iCell] = 0.0;}
  gethinv(cell->hmat_cp,cell->hmati_cp,&(cell->vol_cp),iperd);
  gethinv(cell->hmat,cell->hmati,&(cell->vol),iperd);

/*======================================================================*/
/* IV) Check the forms                                                   */


/*======================================================================*/
/* V) Initial KS potential calculation                                  */

  stat_avg->cp_ehart = 0.0;
  stat_avg->cp_eext = 0.0;
  stat_avg->cp_exc = 0.0;
  if(myidState==0)printf("**Calculating Initial Kohn-Sham Potential...\n");
  calcLocalPseudoScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  calcKSPot(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  if(pseudoRealFlag==1){
    if(myidState==0)printf("**Calculating Real Space Non-local Pseudopotential...\n");
    pseudoReal->forceCalcFlag = 1;
    initRealNlppWf(cp,class,general_data);
    allocRealNl(cp,class);
    pseudoReal->forceCalcFlag = 0;
  }

  //calcKSPotExtRecipWrapPreScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  //calcKSForceControlWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

  if(myidState==0)printf("**Finish Calculating Initial Kohn-Sham Potential\n");

  
  //exit(0);
/*======================================================================*/
/* V) SCF loop			                    */

  if(myidState==0){
    printf("===============================================================================\n");
    printf("Runing SCF Calculation\n");
    printf("-------------------------------------------------------------------------------\n");
  }
  if(readCoeffFlag==-3){
    iScf = stodftInfo->iScf;
  }
  else{
    iScf = 0;
  }

  while(scfStopFlag==0){
    timeStart = omp_get_wtime();
    iScf += 1;
    stodftInfo->iScf = iScf;

    if(myidState==0){
      printf("********************************************************\n");
      printf("SCF Step %i\n",iScf);
      printf("--------------------------------------------------------\n");
    }

/*----------------------------------------------------------------------*/
/* ii) Generate stochastic WF for different chemical potentials         */

    timeStart1 = omp_get_wtime();
    if(myidState==0)printf("**Generating Stochastic Orbitals...\n");
    genStoOrbitalInterp(class,general_data,cp,ip_now);
 
    //exit(0);
    //debug
    /*
    for(iChem=0;iChem<numChemPot;iChem++){
      if(checkNanArray(&stoWfUpRe[iChem][1],numCoeffUpTotal)==1){
    printf("iChem %i myid %i Bad number in upper sate real part!\n",iChem,myidState);
      }
      if(checkNanArray(&stoWfUpIm[iChem][1],numCoeffUpTotal)==1){
        printf("iChem %i myid %i Bad number in upper sate imag part!\n",iChem,myidState);
      }
    }
    */

  /*  
  char wfname[100];
  //sprintf(wfname,"/scratch/mingchen/tmp/sto-wf-save-%i",myidState);
  printf("Read in stochastic orbitals...\n");
  sprintf(wfname,"sto-wf-save-%i",myidState);

  FILE *filePrintWF = fopen(wfname,"r");
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iState=0;iState<numStateUp;iState++){
      for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
        fscanf(filePrintWF,"%lg",&stoWfUpRe[iChem][iState*numCoeff+iCoeff]);
        fscanf(filePrintWF,"%lg",&stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
        //fprintf(filePrintWF,"%.16lg %.16lg\n",stoWfUpRe[iChem][iState*numCoeff+iCoeff],
        //        stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
      }//endfor iCoeff
    }//endfor iState
  }//endfor iChem
  fclose(filePrintWF);
  */


    if(myidState==0){
      printf("**Finish Generating Stochastic Orbitals\n");
      fflush(stdout);
    }
    timeEnd1 = omp_get_wtime();
    diffTime1 += timeEnd1-timeStart1;
    
    timeStart2 = omp_get_wtime();
    if(myidState==0)printf("**Filter Diagonalization...\n");
    orthDiagDriver(cp,class,general_data,ip_now);
    if(myidState==0)printf("**Finish Filter Diagonalization\n");
    timeEnd2 = omp_get_wtime();
    diffTime2 += timeEnd2-timeStart2;

    //exit(0);   
    
    timeStart3 = omp_get_wtime();  
    if(myidState==0)printf("**Calculating KE and NLPPE...\n");
    //calcEnergyChemPot(cp,class,general_data,cpcoeffs_pos,clatoms_pos);   
    calcKNEEnergyFilterDiag(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
    if(myidState==0)printf("**Finish Calculating KE and NLP E\n");
    timeEnd3 = omp_get_wtime();
    diffTime3 += timeEnd3-timeStart3;

/*----------------------------------------------------------------------*/
/* iii)  Get the total density, for each chemical potential and get     */
/*     total electron number for each chemical potential            */
/*     Interpolate the chemical potential w.r.t                         */
/*     Correct electron number and use the interpolation coefficients   */
/*     to generate the density w.r.t. correct number of electrons.      */

    timeStart4 = omp_get_wtime();
    if(myidState==0)printf("**Calculating Density...\n");
    //if(cpParaOpt==0) calcRhoStoHybrid(class,bonded,general_data,cp,ip_now); 
    if(cpParaOpt==0) calcRhoFilterDiagHybrid(class,bonded,general_data,cp,ip_now);
    if(myidState==0)printf("**Finish Calculating Density\n");
    timeEnd4 = omp_get_wtime();
    diffTime4 += timeEnd3-timeStart3;

/*----------------------------------------------------------------------*/
/* iv) Rerun if necessary                                           */

    /*
    while(stodftInfo->reRunFlag==1){
      if(myidState==0)printf("**Generating Stochastic Orbitals...\n");
      adjChemPot(stodftInfo,stodftCoefPos);
      genStoOrbital(class,bonded,general_data,cp,ip_now);
      if(myidState==0)printf("**Finish Generating Stochastic Orbitals\n");
      if(myidState==0)printf("**Calculating KE and NLPPE...\n");
      calcEnergyChemPot(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
      if(myidState==0)printf("**Finish Calculating KE and NLP E\n");
      if(myidState==0)printf("**Calculating Density...\n");
      if(cpParaOpt==0) calcRhoStoHybrid(class,bonded,general_data,cp,ip_now);
      if(myidState==0)printf("**Finish Calculating Density\n");
    }//endwhile reRunFlag
    */

/*----------------------------------------------------------------------*/
/* iv) Generate KS potential                                            */

    if(numProcStates>1)Barrier(commStates);

    timeStart5 = omp_get_wtime();
    stat_avg->cp_ehart = 0.0;
    stat_avg->cp_eext = 0.0;
    stat_avg->cp_exc = 0.0;
    if(myidState==0)printf("**Calculating Kohn-Sham Potential...\n");
    //calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    calcKSPot(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

    if(myidState==0)printf("**Finish Calculating Kohn-Sham Potential\n");
    timeEnd5 = omp_get_wtime();
    diffTime5 += timeEnd5-timeStart5;

/*----------------------------------------------------------------------*/
/* v) Calculate the total energy	                        */

    timeStart6 = omp_get_wtime();
    if(myidState==0)printf("**Calculating Total Energy...\n");
    stodftInfo->energyElecTotOld = stodftInfo->energyElecTot;
    calcTotEnergyFilterDiag(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
    energyDiff = fabs(stodftInfo->energyElecTot-stodftInfo->energyElecTotOld);
    if(myidState==0)printf("**Finish Calculating Total Energy\n");
    timeEnd6 = omp_get_wtime();
    diffTime6 += timeEnd6-timeStart6;

    //exit(0);   

/*----------------------------------------------------------------------*/
/* vi)  Write checkpoint file if necessary                              */

    if(iScf%checkpointWriteFreq==0){
      if(checkpointParFlag==0)checkpointOutput(cp,general_data);
      else checkpointOutputDist(cp,general_data);
    }

    /*
    if(myidState==0)printf("**Prepare Chemical Potentials for next SCF...\n");
    updateChemPot(stodftInfo,stodftCoefPos);
    if(myidState==0){
      printf("The min chem pot for next step is %.6lg\n",chemPot[0]);
      printf("The max chem pot for next step is %.6lg\n",chemPot[numChemPot-1]);
    }
    if(myidState==0)printf("**Finish Preparing Chemical Potentials for next SCF\n");
    */

    
/*----------------------------------------------------------------------*/
/* v) Finish this SCF step	                            */

    timeEnd = omp_get_wtime();
    if(myidState==0){
      printf("**The master process spend %lgs in this SCF step.\n",timeEnd-timeStart);
      printf("--------------------------------------------------------\n");
      printf("Finish SCF Step %i\n",iScf);
      printf("********************************************************\n");
      printf("\n");
    }
    if(energyDiff<energyTol||iScf>=numScf)scfStopFlag = 1;
    //exit(0);

  }//endfor iScf

  printf("SCF time myid %i gen-stowf %.8lg Diag %.8lg energy %.8lg density %.8lg KS potential %.8lg total-energy %.8lg\n",myidState,diffTime1,diffTime2,diffTime3,diffTime4,diffTime5,diffTime6);


/*======================================================================*/
/* VI) Calculate nuclei forces after SCF loop                           */

  // TODO: a filter diag version of force calculation
  calcEnergyForceFilterDiag(class,general_data,cp,bonded,cpcoeffs_pos,clatoms_pos);

/*======================================================================*/
/* VI) In parallel, transpose coefs and coef forces fwd                 */


/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
#ifdef FAST_FILTER
void scfStodftEnergyWindow(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                    CP *cp,
                    CLASS *class2,BONDED *bonded2,GENERAL_DATA *general_data2,CP *cp2,
                    int ip_now)
#else
void scfStodftEnergyWindow(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                    CP *cp, int ip_now)
#endif
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  CELL *cell		 = &(general_data->cell);  
  PTENS *ptens		 = &(general_data->ptens);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  
  STODFTINFO *stodftInfo         = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos   = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CPEWALD *cpewald              = &(cp->cpewald);
  CPSCR *cpscr                  = &(cp->cpscr);
  PSEUDO *pseudo                = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal       = &(pseudo->pseudoReal);  
  COMMUNICATE *communicate      = &(cp->communicate);

  FOR_SCR      *for_scr         = &(class->for_scr);
  EWD_SCR      *ewd_scr         = &(class->ewd_scr);
  CLATOMS_INFO *clatoms_info	= &(class->clatoms_info);
  CLATOMS_POS *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  ATOMMAPS *atommaps	    = &(class->atommaps);

  int iperd		    = cell->iperd;
  int iScf,iCell,iCoeff,iState,iChem;
  int numScf			= stodftInfo->numScf; //Need claim this in cp
  int numChemPot		= stodftInfo->numChemPot;
  int cpLsda			= cpopts->cp_lsda;
  int cpParaOpt			= cpopts->cp_para_opt;
  int pseudoRealFlag            = pseudoReal->pseudoRealFlag;
  int checkPerdSize		= cpopts->icheck_perd_size;
  int checkDualSize		= cpopts->icheck_dual_size;
  int cpDualGridOptOn		= cpopts->cp_dual_grid_opt;
  int numProcStates		= communicate->np_states; 
  int myidState			= communicate->myid_state;
  int coefFormUp		= cpcoeffs_pos->icoef_form_up;
  int coefOrthUp                = cpcoeffs_pos->icoef_orth_up;
  int forceCoefFormUp		= cpcoeffs_pos->ifcoef_form_up;
  int forceCoefOrthUp           = cpcoeffs_pos->ifcoef_orth_up;
  int coefFormDn                = cpcoeffs_pos->icoef_form_dn;
  int coefOrthDn                = cpcoeffs_pos->icoef_orth_dn;
  int forceCoefFormDn           = cpcoeffs_pos->ifcoef_form_dn;
  int forceCoefOrthDn           = cpcoeffs_pos->ifcoef_orth_dn;
  int numStateUp		= cpcoeffs_info->nstate_up_proc;
  int numStateDn                = cpcoeffs_info->nstate_dn_proc;
  int numCoeff			= cpcoeffs_info->ncoef;
  int numCoeffUpTotal		= numStateUp*numCoeff;
  int numCoeffDnTotal		= numStateDn*numCoeff;
  int scfStopFlag     = 0;
  int checkpointWriteFreq = stodftInfo->checkpointWriteFreq;
  int checkpointParFlag = stodftInfo->checkpointParFlag;
  int readCoeffFlag = stodftInfo->readCoeffFlag;
  int calcFragFlag              = stodftInfo->calcFragFlag;

  MPI_Comm commStates		= communicate->comm_states;

  int *pcoefFormUp	     = &(cpcoeffs_pos->icoef_form_up);
  int *pcoefOrthUp	     = &(cpcoeffs_pos->icoef_orth_up);
  int *pforceCoefFormUp              = &(cpcoeffs_pos->ifcoef_form_up);
  int *pforceCoefOrthUp              = &(cpcoeffs_pos->ifcoef_orth_up);
  int *pcoefFormDn	     = &(cpcoeffs_pos->icoef_form_dn);
  int *pcoefOrthDn	     = &(cpcoeffs_pos->icoef_orth_dn);
  int *pforceCoefFormDn              = &(cpcoeffs_pos->icoef_form_dn);
  int *pforceCoefOrthDn              = &(cpcoeffs_pos->icoef_orth_dn);

  double numElecTrue = stodftInfo->numElecTrue;
  double tolEdgeDist	    = cpopts->tol_edge_dist;
  double energyDiff	= -1.0;
  double energyTol	= stodftInfo->energyTol*numElecTrue;

  double *coeffReUp        = cpcoeffs_pos->cre_up;
  double *coeffImUp        = cpcoeffs_pos->cim_up;
  double *forceCoeffReUp   = cpcoeffs_pos->fcre_up;
  double *forceCoeffImUp   = cpcoeffs_pos->fcre_up;
  double *coeffReDn        = cpcoeffs_pos->cre_dn;
  double *coeffImDn        = cpcoeffs_pos->cim_dn;
  double *forceCoeffReDn   = cpcoeffs_pos->fcre_dn;
  double *forceCoeffImDn   = cpcoeffs_pos->fcre_dn;
  //double *cpScrCoeffReUp   = cpscr->cpscr_wave.cre_up;
  //double *cpScrCoeffImUp   = cpscr->cpscr_wave.cim_up;
  //double *cpScrCoeffReDn   = cpscr->cpscr_wave.cre_dn;
  //double *cpScrCoeffImDn   = cpscr->cpscr_wave.cim_dn;
  double *ptensPvtenTmp    = ptens->pvten_tmp;
  double *chemPot          = stodftCoefPos->chemPot;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  //debug
  FILE *fileRhoRecip;
  int numCoeffLargeProc = cp->cp_para_fft_pkg3d_lg.ncoef_proc;
  double *rhoCoeffReUp = cpscr->cpscr_rho.rhocr_up;
  double *rhoCoeffImUp = cpscr->cpscr_rho.rhoci_up;

  //timeing
  double timeStart,timeEnd;
  double timeStart1,timeEnd1;
  double timeStart2,timeEnd2;
  double timeStart3,timeEnd3;
  double timeStart4,timeEnd4;
  double timeStart5,timeEnd5;
  double timeStart6,timeEnd6;
  double diffTime1 = 0.0;
  double diffTime2 = 0.0;
  double diffTime3 = 0.0;
  double diffTime4 = 0.0;
  double diffTime5 = 0.0;
  double diffTime6 = 0.0;

/*======================================================================*/
/* 0.05) Check the approximations in the methods                        */


  if((iperd<3||iperd==4)&&checkPerdSize==1){
    cp_boundary_check(cell,clatoms_info,clatoms_pos,tolEdgeDist);
  }/*endif*/
  if(cpDualGridOptOn>=1&&checkDualSize==1){
    cp_dual_check(cell,clatoms_info,clatoms_pos,
                  atommaps->cp_atm_lst,tolEdgeDist);
  }/*endif*/

/*======================================================================*/
/* I) Orthogonalize the coefs if norbing                                */

/*======================================================================*/
/* II) In parallel, transpose coefs back to normal form                 */

/*======================================================================*/
/* III) Initialize forces, pressure tensor, inverse hmat                */

  cpcoeffs_pos->ifcoef_form_up = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1;
  if(cpLsda==1&&numStateDn>0){
    cpcoeffs_pos->ifcoef_form_dn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }/*endif*/

  for(iCell=1;iCell<=9;iCell++){ptensPvtenTmp[iCell] = 0.0;}
  gethinv(cell->hmat_cp,cell->hmati_cp,&(cell->vol_cp),iperd);
  gethinv(cell->hmat,cell->hmati,&(cell->vol),iperd);

/*======================================================================*/
/* IV) Check the forms                                                   */


/*======================================================================*/
/* V) Initial KS potential calculation                                  */

  stat_avg->cp_ehart = 0.0;
  stat_avg->cp_eext = 0.0;
  stat_avg->cp_exc = 0.0;
  if(myidState==0)printf("**Calculating Initial Kohn-Sham Potential...\n");
  if(calcFragFlag==0){
    calcLocalPseudoScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  }

  calcKSPot(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

  if(calcFragFlag==0){
    if(pseudoRealFlag==1){
      if(myidState==0)printf("**Calculating Real Space Non-local Pseudopotential...\n");
      pseudoReal->forceCalcFlag = 1;
      initRealNlppWf(cp,class,general_data);
      allocRealNl(cp,class);
      pseudoReal->forceCalcFlag = 0;
    }
  }

#ifdef FAST_FILTER
  cp2->cpcoeffs_pos[1].ifcoef_form_up = 0;
  cp2->cpcoeffs_pos[1].ifcoef_orth_up = 1;

  gethinv(general_data2->cell.hmat_cp,general_data2->cell.hmati_cp,
          &(general_data2->cell.vol_cp),iperd);
  gethinv(general_data2->cell.hmat,general_data2->cell.hmati,
          &(general_data2->cell.vol),iperd);
  if(calcFragFlag==0){
    calcLocalPseudoScf(class2,general_data2,cp2,
                       &(cp2->cpcoeffs_pos[1]),&(class2->clatoms_pos[1]));  
  }

  if(cp2->pseudo.pseudoReal.pseudoRealFlag==1){
    cp2->pseudo.pseudoReal.forceCalcFlag = 1;
    initRealNlppWf(cp2,class2,general_data2);
    allocRealNl(cp2,class2);
    cp2->pseudo.pseudoReal.forceCalcFlag = 0;
  }
  cp->stodftInfo->numStatePrintUp = cp2->stodftInfo->numStatePrintUp;
  cp->stodftCoefPos->moUpRePrint = (double*)cmalloc(cp->stodftInfo->numStatePrintUp*numCoeff*sizeof(double))-1;
  cp->stodftCoefPos->moUpImPrint = (double*)cmalloc(cp->stodftInfo->numStatePrintUp*numCoeff*sizeof(double))-1;
  // attention: size of cp->stodftCoefPos.energyLevel is different from 
  // the size of cp2->stodftCoefPos.energyLevel
  cp->stodftCoefPos->energyLevel = (double*)cmalloc(cp->stodftInfo->numStatePrintUp*sizeof(double));

#endif
  //calcKSPotExtRecipWrapPreScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  //calcKSForceControlWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

  if(myidState==0)printf("**Finish Calculating Initial Kohn-Sham Potential\n");

  //exit(0);
/*======================================================================*/
/* V) SCF loop			                    */

  if(myidState==0){
    printf("===============================================================================\n");
    printf("Runing SCF Calculation\n");
    printf("-------------------------------------------------------------------------------\n");
  }

  if(readCoeffFlag==-3){
    iScf = stodftInfo->iScf;
  }
  else{
    iScf = 0;
  }


  while(scfStopFlag==0){
    timeStart = omp_get_wtime();
    iScf += 1;
    stodftInfo->iScf = iScf;


    if(myidState==0){
      printf("********************************************************\n");
      printf("SCF Step %i\n",iScf);
      printf("--------------------------------------------------------\n");
    }

/*----------------------------------------------------------------------*/
/* ii) Generate stochastic WF for different chemical potentials         */

    timeStart1 = omp_get_wtime();
    if(myidState==0)printf("**Generating Stochastic Orbitals...\n");
/*
#ifndef EW_FAKE
    genStoOrbitalEnergyWindow(class,general_data,cp,ip_now);
#else
    genStoOrbitalEnergyWindowFake(class,general_data,cp,ip_now);
#endif
*/
#ifdef FAST_FILTER
      genStoOrbitalEnergyWindowFakeNew(class,general_data,cp,
                             class2,general_data2,cp2,ip_now);
#else
      genStoOrbitalEnergyWindow(class,general_data,cp,ip_now);
#endif
 
    //exit(0);
    //debug
    /*
    for(iChem=0;iChem<numChemPot;iChem++){
      if(checkNanArray(&stoWfUpRe[iChem][1],numCoeffUpTotal)==1){
    printf("iChem %i myid %i Bad number in upper sate real part!\n",iChem,myidState);
      }
      if(checkNanArray(&stoWfUpIm[iChem][1],numCoeffUpTotal)==1){
        printf("iChem %i myid %i Bad number in upper sate imag part!\n",iChem,myidState);
      }
    }
    */

  /*  
  char wfname[100];
  //sprintf(wfname,"/scratch/mingchen/tmp/sto-wf-save-%i",myidState);
  printf("Read in stochastic orbitals...\n");
  sprintf(wfname,"sto-wf-save-%i",myidState);

  FILE *filePrintWF = fopen(wfname,"r");
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iState=0;iState<numStateUp;iState++){
      for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
        fscanf(filePrintWF,"%lg",&stoWfUpRe[iChem][iState*numCoeff+iCoeff]);
        fscanf(filePrintWF,"%lg",&stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
        //fprintf(filePrintWF,"%.16lg %.16lg\n",stoWfUpRe[iChem][iState*numCoeff+iCoeff],
        //        stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
      }//endfor iCoeff
    }//endfor iState
  }//endfor iChem
  fclose(filePrintWF);
  */


    if(myidState==0)printf("**Finish Generating Stochastic Orbitals\n");
    fflush(stdout);
    timeEnd1 = omp_get_wtime();
    diffTime1 += timeEnd1-timeStart1;


/*----------------------------------------------------------------------*/
/* iii) Calculate the energy (KE,NL PP) for each energy window          */


    timeStart2 = omp_get_wtime();
    if(myidState==0)printf("**Calculating KE and NLPPE...\n");
    calcEnergyChemPot(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
    if(myidState==0)printf("**Finish Calculating KE and NLP E\n");
    timeEnd2 = omp_get_wtime();
    diffTime2 += timeEnd2-timeStart2;

    //exit(0);   

/*----------------------------------------------------------------------*/
/* iv) Calcualte the density for each energy window. The total density  */
/*     the summation of each energy window density.                     */

    timeStart3 = omp_get_wtime();
    if(myidState==0)printf("**Calculating Density...\n");
    if(cpParaOpt==0) calcRhoStoHybridEnergyWindow(class,bonded,general_data,cp,ip_now);
    if(myidState==0)printf("**Finish Calculating Density\n");
    timeEnd3 = omp_get_wtime();
    diffTime3 += timeEnd3-timeStart3;

/*----------------------------------------------------------------------*/
/* v) Generate KS potential                                             */

    if(numProcStates>1)Barrier(commStates);

    stat_avg->cp_ehart = 0.0;
    stat_avg->cp_eext = 0.0;
    stat_avg->cp_exc = 0.0;
    if(myidState==0)printf("**Calculating Kohn-Sham Potential...\n");
    //calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    calcKSPot(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

    if(myidState==0)printf("**Finish Calculating Kohn-Sham Potential\n");

/*----------------------------------------------------------------------*/
/* vi) Calculate the total energy	                                */

    if(myidState==0)printf("**Calculating Total Energy...\n");
    stodftInfo->energyElecTotOld = stodftInfo->energyElecTot;
    calcTotEnergy(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
    energyDiff = fabs(stodftInfo->energyElecTot-stodftInfo->energyElecTotOld);
    if(myidState==0)printf("**Finish Calculating Total Energy\n");

    //exit(0);   

/*----------------------------------------------------------------------*/
/* vii) Write checkpoint file if necessary                              */

    if(iScf%checkpointWriteFreq==0){
      if(checkpointParFlag==0)checkpointOutput(cp,general_data);
      else checkpointOutputDist(cp,general_data);
    }
    
/*----------------------------------------------------------------------*/
/* viii) Finish this SCF step	                                        */

    timeEnd = omp_get_wtime();

    if(myidState==0){
      printf("--------------------------------------------------------\n");
      printf("Finish SCF Step %i\n",iScf);
      printf("********************************************************\n");
      printf("\n");
    }
    if(energyDiff<energyTol||iScf>=numScf)scfStopFlag = 1;

  }//endwhile iScf
  printf("SCF time myid %i gen-stowf %.8lg energy %.8lg density %.8lg KS potential %.8lg total-energy %.8lg\n",
         myidState,diffTime1,diffTime2,diffTime3,diffTime4,diffTime5);


/*======================================================================*/
/* ix) Calculate nuclei forces after SCF loop                           */

  calcEnergyForce(class,general_data,cp,bonded,cpcoeffs_pos,clatoms_pos);

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
#ifdef FAST_FILTER
void scfStodftEnergyWindowFrag(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                    CP *cp,
                    CLASS *class2,BONDED *bonded2,GENERAL_DATA *general_data2,CP *cp2,
                    CP *cpMini, GENERAL_DATA *generalDataMini, CLASS *classMini,
                    int ip_now)
#else
void scfStodftEnergyWindowFrag(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                    CP *cp, CP *cpMini, GENERAL_DATA *generalDataMini, CLASS *classMini,
                    int ip_now)
#endif
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  CELL *cell		 = &(general_data->cell);  
  PTENS *ptens		 = &(general_data->ptens);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  
  STODFTINFO *stodftInfo         = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos   = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CPEWALD *cpewald              = &(cp->cpewald);
  CPSCR *cpscr                  = &(cp->cpscr);
  PSEUDO *pseudo                = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal       = &(pseudo->pseudoReal);  
  COMMUNICATE *communicate      = &(cp->communicate);

  FOR_SCR      *for_scr         = &(class->for_scr);
  EWD_SCR      *ewd_scr         = &(class->ewd_scr);
  CLATOMS_INFO *clatoms_info	= &(class->clatoms_info);
  CLATOMS_POS *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  ATOMMAPS *atommaps	    = &(class->atommaps);

  int iperd		    = cell->iperd;
  int iScf,iCell,iCoeff,iState,iChem;
  int numScf			= stodftInfo->numScf; //Need claim this in cp
  int numChemPot		= stodftInfo->numChemPot;
  int cpLsda			= cpopts->cp_lsda;
  int cpParaOpt			= cpopts->cp_para_opt;
  int pseudoRealFlag            = pseudoReal->pseudoRealFlag;
  int checkPerdSize		= cpopts->icheck_perd_size;
  int checkDualSize		= cpopts->icheck_dual_size;
  int cpDualGridOptOn		= cpopts->cp_dual_grid_opt;
  int numProcStates		= communicate->np_states; 
  int myidState			= communicate->myid_state;
  int coefFormUp		= cpcoeffs_pos->icoef_form_up;
  int coefOrthUp                = cpcoeffs_pos->icoef_orth_up;
  int forceCoefFormUp		= cpcoeffs_pos->ifcoef_form_up;
  int forceCoefOrthUp           = cpcoeffs_pos->ifcoef_orth_up;
  int coefFormDn                = cpcoeffs_pos->icoef_form_dn;
  int coefOrthDn                = cpcoeffs_pos->icoef_orth_dn;
  int forceCoefFormDn           = cpcoeffs_pos->ifcoef_form_dn;
  int forceCoefOrthDn           = cpcoeffs_pos->ifcoef_orth_dn;
  int numStateUp		= cpcoeffs_info->nstate_up_proc;
  int numStateDn                = cpcoeffs_info->nstate_dn_proc;
  int numCoeff			= cpcoeffs_info->ncoef;
  int numCoeffUpTotal		= numStateUp*numCoeff;
  int numCoeffDnTotal		= numStateDn*numCoeff;
  int scfStopFlag     = 0;
  int checkpointWriteFreq = stodftInfo->checkpointWriteFreq;
  int checkpointParFlag = stodftInfo->checkpointParFlag;
  int readCoeffFlag = stodftInfo->readCoeffFlag;
  int calcFragFlag              = stodftInfo->calcFragFlag;
  int fragWindowFlag = stodftInfo->fragWindowFlag;

  MPI_Comm commStates		= communicate->comm_states;

  int *pcoefFormUp	     = &(cpcoeffs_pos->icoef_form_up);
  int *pcoefOrthUp	     = &(cpcoeffs_pos->icoef_orth_up);
  int *pforceCoefFormUp              = &(cpcoeffs_pos->ifcoef_form_up);
  int *pforceCoefOrthUp              = &(cpcoeffs_pos->ifcoef_orth_up);
  int *pcoefFormDn	     = &(cpcoeffs_pos->icoef_form_dn);
  int *pcoefOrthDn	     = &(cpcoeffs_pos->icoef_orth_dn);
  int *pforceCoefFormDn              = &(cpcoeffs_pos->icoef_form_dn);
  int *pforceCoefOrthDn              = &(cpcoeffs_pos->icoef_orth_dn);

  double numElecTrue = stodftInfo->numElecTrue;
  double tolEdgeDist	    = cpopts->tol_edge_dist;
  double energyDiff	= -1.0;
  double energyTol	= stodftInfo->energyTol*numElecTrue;

  double *coeffReUp        = cpcoeffs_pos->cre_up;
  double *coeffImUp        = cpcoeffs_pos->cim_up;
  double *forceCoeffReUp   = cpcoeffs_pos->fcre_up;
  double *forceCoeffImUp   = cpcoeffs_pos->fcre_up;
  double *coeffReDn        = cpcoeffs_pos->cre_dn;
  double *coeffImDn        = cpcoeffs_pos->cim_dn;
  double *forceCoeffReDn   = cpcoeffs_pos->fcre_dn;
  double *forceCoeffImDn   = cpcoeffs_pos->fcre_dn;
  //double *cpScrCoeffReUp   = cpscr->cpscr_wave.cre_up;
  //double *cpScrCoeffImUp   = cpscr->cpscr_wave.cim_up;
  //double *cpScrCoeffReDn   = cpscr->cpscr_wave.cre_dn;
  //double *cpScrCoeffImDn   = cpscr->cpscr_wave.cim_dn;
  double *ptensPvtenTmp    = ptens->pvten_tmp;
  double *chemPot          = stodftCoefPos->chemPot;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  //debug
  FILE *fileRhoRecip;
  int numCoeffLargeProc = cp->cp_para_fft_pkg3d_lg.ncoef_proc;
  double *rhoCoeffReUp = cpscr->cpscr_rho.rhocr_up;
  double *rhoCoeffImUp = cpscr->cpscr_rho.rhoci_up;

  //timeing
  double timeStart,timeEnd;
  double timeStart1,timeEnd1;
  double timeStart2,timeEnd2;
  double timeStart3,timeEnd3;
  double timeStart4,timeEnd4;
  double timeStart5,timeEnd5;
  double timeStart6,timeEnd6;
  double diffTime1 = 0.0;
  double diffTime2 = 0.0;
  double diffTime3 = 0.0;
  double diffTime4 = 0.0;
  double diffTime5 = 0.0;
  double diffTime6 = 0.0;

/*======================================================================*/
/* 0.05) Check the approximations in the methods                        */


  if((iperd<3||iperd==4)&&checkPerdSize==1){
    cp_boundary_check(cell,clatoms_info,clatoms_pos,tolEdgeDist);
  }/*endif*/
  if(cpDualGridOptOn>=1&&checkDualSize==1){
    cp_dual_check(cell,clatoms_info,clatoms_pos,
                  atommaps->cp_atm_lst,tolEdgeDist);
  }/*endif*/

/*======================================================================*/
/* I) Orthogonalize the coefs if norbing                                */

/*======================================================================*/
/* II) In parallel, transpose coefs back to normal form                 */

/*======================================================================*/
/* III) Initialize forces, pressure tensor, inverse hmat                */

  cpcoeffs_pos->ifcoef_form_up = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1;
  if(cpLsda==1&&numStateDn>0){
    cpcoeffs_pos->ifcoef_form_dn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }/*endif*/

  for(iCell=1;iCell<=9;iCell++){ptensPvtenTmp[iCell] = 0.0;}
  gethinv(cell->hmat_cp,cell->hmati_cp,&(cell->vol_cp),iperd);
  gethinv(cell->hmat,cell->hmati,&(cell->vol),iperd);

/*======================================================================*/
/* IV) Check the forms                                                   */


/*======================================================================*/
/* V) Initial KS potential calculation                                  */

  stat_avg->cp_ehart = 0.0;
  stat_avg->cp_eext = 0.0;
  stat_avg->cp_exc = 0.0;
  if(myidState==0)printf("**Calculating Initial Kohn-Sham Potential...\n");
  calcLocalPseudoScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

  calcKSPot(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

  printf("1111111111111 pseudoRealFlag %i\n",pseudoRealFlag);
  if(pseudoRealFlag==1){
    if(myidState==0)printf("**Calculating Real Space Non-local Pseudopotential...\n");
    pseudoReal->forceCalcFlag = 1;
    initRealNlppWf(cp,class,general_data);
    allocRealNl(cp,class);
    pseudoReal->forceCalcFlag = 0;
  }

#ifdef FAST_FILTER
  cp2->cpcoeffs_pos[1].ifcoef_form_up = 0;
  cp2->cpcoeffs_pos[1].ifcoef_orth_up = 1;

  gethinv(general_data2->cell.hmat_cp,general_data2->cell.hmati_cp,
          &(general_data2->cell.vol_cp),iperd);
  gethinv(general_data2->cell.hmat,general_data2->cell.hmati,
          &(general_data2->cell.vol),iperd);
  calcLocalPseudoScf(class2,general_data2,cp2,
                     &(cp2->cpcoeffs_pos[1]),&(class2->clatoms_pos[1]));  

  if(cp2->pseudo.pseudoReal.pseudoRealFlag==1){
    cp2->pseudo.pseudoReal.forceCalcFlag = 1;
    initRealNlppWf(cp2,class2,general_data2);
    allocRealNl(cp2,class2);
    cp2->pseudo.pseudoReal.forceCalcFlag = 0;
  }
  cp->stodftInfo->numStatePrintUp = cp2->stodftInfo->numStatePrintUp;
  cp->stodftCoefPos->moUpRePrint = (double*)cmalloc(cp->stodftInfo->numStatePrintUp*numCoeff*sizeof(double))-1;
  cp->stodftCoefPos->moUpImPrint = (double*)cmalloc(cp->stodftInfo->numStatePrintUp*numCoeff*sizeof(double))-1;
  // attention: size of cp->stodftCoefPos.energyLevel is different from 
  // the size of cp2->stodftCoefPos.energyLevel
  cp->stodftCoefPos->energyLevel = (double*)cmalloc(cp->stodftInfo->numStatePrintUp*sizeof(double));

#endif
  //calcKSPotExtRecipWrapPreScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  //calcKSForceControlWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

  if(myidState==0)printf("**Finish Calculating Initial Kohn-Sham Potential\n");

  //exit(0);
/*======================================================================*/
/* V) SCF loop			                    */

  if(myidState==0){
    printf("===============================================================================\n");
    printf("Runing SCF Calculation\n");
    printf("-------------------------------------------------------------------------------\n");
  }

  if(readCoeffFlag==-3){
    iScf = stodftInfo->iScf;
  }
  else{
    iScf = 0;
  }


  while(scfStopFlag==0){
    timeStart = omp_get_wtime();
    iScf += 1;
    stodftInfo->iScf = iScf;


    if(myidState==0){
      printf("********************************************************\n");
      printf("SCF Step %i\n",iScf);
      printf("--------------------------------------------------------\n");
    }

/*----------------------------------------------------------------------*/
/* ii) Generate stochastic WF for different chemical potentials         */

    timeStart1 = omp_get_wtime();
    if(myidState==0)printf("**Generating Stochastic Orbitals...\n");
/*
#ifndef EW_FAKE
    genStoOrbitalEnergyWindow(class,general_data,cp,ip_now);
#else
    genStoOrbitalEnergyWindowFake(class,general_data,cp,ip_now);
#endif
*/
    if(fragWindowFlag==0){
#ifdef FAST_FILTER
      genStoOrbitalEnergyWindowFakeNew(class,general_data,cp,
                             class2,general_data2,cp2,ip_now);
#else
      genStoOrbitalEnergyWindow(class,general_data,cp,ip_now);
#endif
    }
    else{
#ifdef FAST_FILTER
      printf("I haven't prepared FAST_FILTER for energy_window fragment\n");
      fflush(stdout);
      exit(0);
#else
      genStoOrbitalEnergyWindowFragFake(class,general_data,
                       cp,generalDataMini,cpMini,classMini,ip_now);

#endif
    }
 
    //exit(0);
    //debug
    /*
    for(iChem=0;iChem<numChemPot;iChem++){
      if(checkNanArray(&stoWfUpRe[iChem][1],numCoeffUpTotal)==1){
    printf("iChem %i myid %i Bad number in upper sate real part!\n",iChem,myidState);
      }
      if(checkNanArray(&stoWfUpIm[iChem][1],numCoeffUpTotal)==1){
        printf("iChem %i myid %i Bad number in upper sate imag part!\n",iChem,myidState);
      }
    }
    */

  /*  
  char wfname[100];
  //sprintf(wfname,"/scratch/mingchen/tmp/sto-wf-save-%i",myidState);
  printf("Read in stochastic orbitals...\n");
  sprintf(wfname,"sto-wf-save-%i",myidState);

  FILE *filePrintWF = fopen(wfname,"r");
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iState=0;iState<numStateUp;iState++){
      for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
        fscanf(filePrintWF,"%lg",&stoWfUpRe[iChem][iState*numCoeff+iCoeff]);
        fscanf(filePrintWF,"%lg",&stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
        //fprintf(filePrintWF,"%.16lg %.16lg\n",stoWfUpRe[iChem][iState*numCoeff+iCoeff],
        //        stoWfUpIm[iChem][iState*numCoeff+iCoeff]);
      }//endfor iCoeff
    }//endfor iState
  }//endfor iChem
  fclose(filePrintWF);
  */


    if(myidState==0)printf("**Finish Generating Stochastic Orbitals\n");
    fflush(stdout);
    timeEnd1 = omp_get_wtime();
    diffTime1 += timeEnd1-timeStart1;


/*----------------------------------------------------------------------*/
/* iii) Calculate the energy (KE,NL PP) for each energy window          */


    timeStart2 = omp_get_wtime();
    if(myidState==0)printf("**Calculating KE and NLPPE...\n");
    calcEnergyChemPot(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
    if(myidState==0)printf("**Finish Calculating KE and NLP E\n");
    timeEnd2 = omp_get_wtime();
    diffTime2 += timeEnd2-timeStart2;

    //exit(0);   

/*----------------------------------------------------------------------*/
/* iv) Calcualte the density for each energy window. The total density  */
/*     the summation of each energy window density.                     */

    timeStart3 = omp_get_wtime();
    if(myidState==0)printf("**Calculating Density...\n");
    if(cpParaOpt==0) calcRhoStoHybridEnergyWindow(class,bonded,general_data,cp,ip_now);
    if(myidState==0)printf("**Finish Calculating Density\n");
    timeEnd3 = omp_get_wtime();
    diffTime3 += timeEnd3-timeStart3;

/*----------------------------------------------------------------------*/
/* v) Generate KS potential                                             */

    if(numProcStates>1)Barrier(commStates);

    stat_avg->cp_ehart = 0.0;
    stat_avg->cp_eext = 0.0;
    stat_avg->cp_exc = 0.0;
    if(myidState==0)printf("**Calculating Kohn-Sham Potential...\n");
    //calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    calcKSPot(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

    if(myidState==0)printf("**Finish Calculating Kohn-Sham Potential\n");

/*----------------------------------------------------------------------*/
/* vi) Calculate the total energy	                                */

    if(myidState==0)printf("**Calculating Total Energy...\n");
    stodftInfo->energyElecTotOld = stodftInfo->energyElecTot;
    calcTotEnergy(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
    energyDiff = fabs(stodftInfo->energyElecTot-stodftInfo->energyElecTotOld);
    if(myidState==0)printf("**Finish Calculating Total Energy\n");

    //exit(0);   

/*----------------------------------------------------------------------*/
/* vii) Write checkpoint file if necessary                              */

    if(iScf%checkpointWriteFreq==0){
      if(checkpointParFlag==0)checkpointOutput(cp,general_data);
      else checkpointOutputDist(cp,general_data);
    }
    
/*----------------------------------------------------------------------*/
/* viii) Finish this SCF step	                                        */

    timeEnd = omp_get_wtime();

    if(myidState==0){
      printf("--------------------------------------------------------\n");
      printf("Finish SCF Step %i\n",iScf);
      printf("********************************************************\n");
      printf("\n");
    }
    if(energyDiff<energyTol||iScf>=numScf)scfStopFlag = 1;

  }//endwhile iScf
  printf("SCF time myid %i gen-stowf %.8lg energy %.8lg density %.8lg KS potential %.8lg total-energy %.8lg\n",
         myidState,diffTime1,diffTime2,diffTime3,diffTime4,diffTime5);


/*======================================================================*/
/* ix) Calculate nuclei forces after SCF loop                           */

  calcEnergyForce(class,general_data,cp,bonded,cpcoeffs_pos,clatoms_pos);

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/
