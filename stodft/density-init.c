/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: density-init.c                                 */
/*                                                                          */
/* This routine calculate initial density from either deterministic or      */
/* random wave functions.                                                   */
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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcRhoInit(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                   CP *cp,int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the routine to calculate the initial density                  */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  CELL *cell                     = &(general_data->cell);
  CLATOMS_POS *clatoms_pos      = &(class->clatoms_pos[ip_now]);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  ATOMMAPS *atommaps            = &(class->atommaps);
  EWD_SCR      *ewd_scr         = &(class->ewd_scr);
  FOR_SCR      *for_scr         = &(class->for_scr);
  COMMUNICATE  *communicate	= &(class->communicate);

  STODFTINFO *stodftInfo       = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos = cp->stodftCoefPos;
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CPOPTS *cpopts                = &(cp->cpopts);
  PSEUDO *pseudo                = &(cp->pseudo);


  int readCoeffFlag = stodftInfo->readCoeffFlag;
  int vpsAtomListFlag = stodftInfo->vpsAtomListFlag;
  int cpDualGridOptOn           = cpopts->cp_dual_grid_opt;
  int cpLsda			= cpopts->cp_lsda;
  int myid			= communicate->myid;


/*==========================================================================*/
/* I) Generate initial density, stochastic Case                             */

  if(myid==0){
    PRINT_LINE_STAR;
    printf("Start Calculating Initial Density\n");
    PRINT_LINE_DASH;
  }

  if(readCoeffFlag==1) calcRhoStoInit(class,bonded,general_data,cp,cpcoeffs_pos);
  if(readCoeffFlag==2) calcRhoDetInit(class,bonded,general_data,cp,cpcoeffs_pos);
  if(readCoeffFlag==3) readRho(class,bonded,general_data,cp,cpcoeffs_pos);
  

  if(myid==0){
    PRINT_LINE_DASH;
    printf("Finish Calculating Initial Density\n");
    PRINT_LINE_STAR;
  }

/*==========================================================================*/
/* II) Calculate the non-local pseudopotential list                          */

  if(myid==0){
    PRINT_LINE_STAR;
    printf("Start Generating Pseudopotential List\n");
    PRINT_LINE_DASH;
  }

  if(stodftInfo->vpsAtomListFlag==0||cpDualGridOptOn>= 1){
    control_vps_atm_list(pseudo,cell,clatoms_pos,clatoms_info,
                         atommaps,ewd_scr,for_scr,cpDualGridOptOn,
                         stodftInfo->vpsAtomListFlag);
    stodftInfo->vpsAtomListFlag = 1;
  }

  if(myid==0){
    PRINT_LINE_DASH;
    printf("Finish Generating Pseudopotential List\n");
    PRINT_LINE_STAR;
  }

/*==========================================================================*/
/* III) Change the occupation number after reading initial deisity          */
/*      density calculation                                                 */

  if(cpLsda==1)stodftInfo->occNumber = 1;
  else stodftInfo->occNumber = 2;


/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcRhoDetInit(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
	   CP *cp,CPCOEFFS_POS  *cpcoeffs_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the routine to calculate the density from deterministic orbitals */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  #include "../typ_defs/typ_mask.h"
  
  EWALD        *ewald        = &(general_data->ewald);
  CELL         *cell         = &(general_data->cell);
  CPOPTS       *cpopts       = &(cp->cpopts);  
  CPSCR        *cpscr	     = &(cp->cpscr);
  CPEWALD      *cpewald      = &(cp->cpewald);
  PSEUDO       *pseudo	     = &(cp->pseudo);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  COMMUNICATE   *commCP         = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);

  int densityMixFlag = stodftInfo->densityMixFlag;
  int cpParaOpt = cpopts->cp_para_opt;
  int cpLsda = cpopts->cp_lsda;
  int cpGGA  = cpopts->cp_gga;
  int numCoeffLarge       = cpcoeffs_info->ncoef_l;
  int numCoeffLargeProc   = cp->cp_para_fft_pkg3d_lg.ncoef_proc;
  int numCoeffLargeProcDensCpBox = cp->cp_para_fft_pkg3d_dens_cp_box.ncoef_proc;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numInterpPmeDual = pseudo->n_interp_pme_dual;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int rhoRealGridNum = stodftInfo->rhoRealGridNum;
  int filterDiagFlag = stodftInfo->filterDiagFlag;
  int iCoeff,iGrid;
  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);
  MPI_Comm comm_states   =    commCP->comm_states;

  double volCP;

  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;
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
  double *hmatCP	 = cell->hmat_cp;
  double *rhoUpCorrect   = stodftCoefPos->rhoUpCorrect;
  double *rhoDnCorrect   = stodftCoefPos->rhoDnCorrect;

/*======================================================================*/
/* I) Calculate the density			                        */


  //Debug Flag: we temp do this for debug
  if(cpParaOpt==0){
    cp_rho_calc_hybrid(cpewald,cpscr,cpcoeffs_info,
	       ewald,cell,coeffReUp,coeffImUp,*coefFormUp,*coefOrthUp,
	       rhoCoeffReUp,rhoCoeffImUp,rhoUp,rhoCoeffReUpDensCpBox,
	       rhoCoeffImUpDensCpBox,divRhoxUp,divRhoyUp,
	       divRhozUp,d2RhoUp,numStateUpProc,numCoeff,
	       cpGGA,cpDualGridOptOn,numInterpPmeDual,commCP,
	       &(cp->cp_para_fft_pkg3d_lg),&(cp->cp_sclr_fft_pkg3d_lg),
	       &(cp->cp_para_fft_pkg3d_dens_cp_box),
	       &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
	       &(cp->cp_sclr_fft_pkg3d_sm));
  }
  if(cpParaOpt==1){
    cp_rho_calc_full_g(cpewald,cpscr,cpcoeffs_info,
                       ewald,cell,coeffReUp,coeffImUp,*coefFormUp,*coefOrthUp,
                       rhoCoeffReUp,rhoCoeffImUp,rhoUp,rhoCoeffReUpDensCpBox,
                       rhoCoeffImUpDensCpBox,divRhoxUp,divRhoyUp,
                       divRhozUp,d2RhoUp,numStateUpProc,numCoeff,
                       cpGGA,cpDualGridOptOn,numInterpPmeDual,commCP,
                       &(cp->cp_para_fft_pkg3d_lg),
                       &(cp->cp_para_fft_pkg3d_dens_cp_box),
                       &(cp->cp_para_fft_pkg3d_sm));
  }
  if(cpLsda==1&&numStateDnProc>0){
    if(cpParaOpt==0){
      cp_rho_calc_hybrid(cpewald,cpscr,cpcoeffs_info,
                     ewald,cell,coeffReDn,coeffImUp,*coefFormDn,*coefOrthDn,
                     rhoCoeffReDn,rhoCoeffImDn,rhoDn,rhoCoeffReUpDensCpBox,
                     rhoCoeffImDnDensCpBox,divRhoxDn,divRhoyDn,
                     divRhozDn,d2RhoDn,numStateDnProc,numCoeff,
                     cpGGA,cpDualGridOptOn,numInterpPmeDual,commCP,
                     &(cp->cp_para_fft_pkg3d_lg),
                     &(cp->cp_sclr_fft_pkg3d_lg),
                     &(cp->cp_para_fft_pkg3d_dens_cp_box),
                     &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
                     &(cp->cp_sclr_fft_pkg3d_sm));
    }
    if(cpParaOpt==1){
      cp_rho_calc_full_g(cpewald,cpscr,cpcoeffs_info,
	       ewald,cell,coeffReDn,coeffImUp,*coefFormDn,*coefOrthDn,
	       rhoCoeffReDn,rhoCoeffImDn,rhoDn,rhoCoeffReUpDensCpBox,
	       rhoCoeffImDnDensCpBox,divRhoxDn,divRhoyDn,
	       divRhozDn,d2RhoDn,numStateDnProc,numCoeff,
	       cpGGA,cpDualGridOptOn,numInterpPmeDual,commCP,
                       &(cp->cp_para_fft_pkg3d_lg),
                       &(cp->cp_para_fft_pkg3d_dens_cp_box),
                       &(cp->cp_para_fft_pkg3d_sm));
    }
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

  // Backup wavefunctions if necessary
  if(filterDiagFlag==1){
    double *wfUpReDet = stodftCoefPos->wfUpReDet;
    double *wfUpImDet = stodftCoefPos->wfUpImDet;
    int numStatesDet = stodftInfo->numStatesDet;
    int *numStatesAllDet = stodftInfo->numStatesAllDet;
    int *dsplStatesAllDet = stodftInfo->dsplStatesAllDet;

    Gatherv(&coeffReUp[1],numStateUpProc*numCoeff,MPI_DOUBLE,
	    wfUpReDet,numStatesAllDet,dsplStatesAllDet,MPI_DOUBLE,
	    0,comm_states);

    Gatherv(&coeffImUp[1],numStateUpProc*numCoeff,MPI_DOUBLE,
	    wfUpImDet,numStatesAllDet,dsplStatesAllDet,MPI_DOUBLE,
	    0,comm_states);
  }

  //debug
  /*
  stodftCoefPos->creTest = (double*)cmalloc(numCoeff*sizeof(double))-1;
  stodftCoefPos->cimTest = (double*)cmalloc(numCoeff*sizeof(double))-1;
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    stodftCoefPos->creTest[iCoeff] = coeffReUp[numCoeff+iCoeff];
    stodftCoefPos->cimTest[iCoeff] = coeffImUp[numCoeff+iCoeff];
  }
  */

/*======================================================================*/
/* II) Store the density for mixing			                */
  
  //printf("densityMixFlag %i\n",densityMixFlag);
  if(densityMixFlag>0){
    volCP  = getdeth(hmatCP);
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
      rhoUpCorrect[iGrid] = rhoUp[iGrid+1]*volCP;
    }
    if(cpLsda==1&&numStateDnProc>0){
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
	rhoDnCorrect[iGrid] = rhoDn[iGrid+1]*volCP;
      }
    }
    //debug
    genDensityMix(cp,0);
    /*
    int itest;
    for(itest=0;itest<100;itest++){
      genDensityMix(cp,itest);
      rhoUpCorrect[0] += 0.1/(itest+1);
    }
    */
    
  }
  //exit(0);
  /*
  //debug
  FILE *fileRhoInit = fopen("density-init","w");
  FILE *fileRhoInitRecip = fopen("density-recip","w");
  int iGrid;
  for(iGrid=1;iGrid<=rhoRealGridNum;iGrid++)fprintf(fileRhoInit,"%.10lg\n",rhoUp[iGrid]);
  for(iCoeff=1;iCoeff<=numCoeffLargeProc;iCoeff++)fprintf(fileRhoInitRecip,"%.10lg %.10lg\n",rhoCoeffReUp[iCoeff],rhoCoeffImUp[iCoeff]);
  fclose(fileRhoInit);
  fclose(fileRhoInitRecip);
  */

/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcRhoStoInit(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                   CP *cp,CPCOEFFS_POS  *cpcoeffs_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the routine to calculate the density from stochastic orbitals */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  EWALD        *ewald        = &(general_data->ewald);
  CELL         *cell         = &(general_data->cell);
  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR        *cpscr        = &(cp->cpscr);
  CPEWALD      *cpewald      = &(cp->cpewald);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  COMMUNICATE   *commCP         = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  PSEUDO       *pseudo       = &(cp->pseudo);

  int cpParaOpt = cpopts->cp_para_opt;
  int cpLsda = cpopts->cp_lsda;
  int cpGGA  = cpopts->cp_gga;
  int numProcStates = commCP->np_states;
  int numCoeffLarge       = cpcoeffs_info->ncoef_l;
  int numCoeffLargeProc   = cp->cp_para_fft_pkg3d_lg.ncoef_proc;
  int numCoeffLargeProcDensCpBox = cp->cp_para_fft_pkg3d_dens_cp_box.ncoef_proc;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numInterpPmeDual = pseudo->n_interp_pme_dual;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numStateStoUp  = stodftInfo->numStateStoUp;
  int numStateStoDn  = stodftInfo->numStateStoDn;
  int numCoeff       = cpcoeffs_info->ncoef;
  int rhoRealGridNum    = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot    = stodftInfo->rhoRealGridTot;
  int occNumber         = stodftInfo->occNumber;

  int iCoeff;
  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);

  double numGridTotInv = 1.0/rhoRealGridTot;
  double aveFactUp = occNumber/numStateStoUp;

  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;
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
  double *rhoScr;
  double *rhoTemp = (double*)cmalloc(rhoRealGridNum*sizeof(double));


/*======================================================================*/
/* III) Calculate the density                                           */


  //Debug Flag: we temp do this for debug
  /*
  if(cpParaOpt==0){
    if(numProcStates>1)rhoScr = cpscr->cpscr_rho.v_ks_up;
    else rhoScr = rhoUp;

    rhoCalcRealStoHybrid(cpscr,cpcoeffs_info,
                   cell,stodftInfo,stoWfUpRe[iChem],
                   stoWfUpIm[iChem],*coefFormUp,*coefOrthUp,rhoScr,
                   numStateUpProc,numCoeff,cpDualGridOptOn,commCP,
                   &(cp->cp_sclr_fft_pkg3d_lg),
                   &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
                   &(cp->cp_sclr_fft_pkg3d_sm));
    //Reduce(&rhoScr[1],rhoRealGridNum,rhoTemp,);
    //calcRhoStoRecipHybrid()
    

  }
  */
/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void readRho(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
           CP *cp,CPCOEFFS_POS  *cpcoeffs_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the routine to calculate the density from deterministic orbitals */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  #include "../typ_defs/typ_mask.h"
  EWALD        *ewald        = &(general_data->ewald);
  CELL         *cell         = &(general_data->cell);
  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR        *cpscr        = &(cp->cpscr);
  CPEWALD      *cpewald      = &(cp->cpewald);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  COMMUNICATE   *commCP         = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  PSEUDO        *pseudo         = &(cp->pseudo);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

  
  int cpParaOpt = cpopts->cp_para_opt;
  int cpLsda = cpopts->cp_lsda;
  int cpGGA  = cpopts->cp_gga;  
  int numCoeffLarge       = cpcoeffs_info->ncoef_l;
  int numCoeffLargeProc   = cp->cp_para_fft_pkg3d_lg.ncoef_proc;
  int numCoeffLargeProcDensCpBox = cp->cp_para_fft_pkg3d_dens_cp_box.ncoef_proc;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numInterpPmeDual = pseudo->n_interp_pme_dual;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numChemPot     = stodftInfo->numChemPot;
  int numFFTProc        = cp_para_fft_pkg3d_lg->nfft_proc;
  int numFFT            = cp_para_fft_pkg3d_lg->nfft;
  int numFFT2           = numFFT/2;
  int numFFT2Proc       = numFFTProc/2;

  int rhoRealGridNum    = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot    = stodftInfo->rhoRealGridTot;
  int numChemProc       = stodftInfo->numChemProc;
  int numStateStoUp     = stodftInfo->numStateStoUp;
  int numStateStoDn     = stodftInfo->numStateStoDn;
  int occNumber         = stodftInfo->occNumber;
  int densityMixFlag    = stodftInfo->densityMixFlag;
  int iScf              = stodftInfo->iScf;
  int myidState         = commCP->myid_state;
  int numProcStates     = commCP->np_states;

  int iCoeff,iChem,iGrid;
  int index;
  int i,j,k;
  int reRunFlag;
  MPI_Comm comm_states = commCP->comm_states;

  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);
  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;

  char *densityReadFileName = stodftInfo->densityReadFileName; 
  FILE *densityReadFile;

  double volCP,rvolCP;

  double *rhoUpRead,*rhoDnRead;
  double *hmatCP    = cell->hmat_cp;
  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;
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
  double *rhoUpCorrect = stodftCoefPos->rhoUpCorrect;
  double *rhoDnCorrect = stodftCoefPos->rhoDnCorrect;
  double *divRhoxDn       = cpscr->cpscr_grho.d_rhox_dn;
  double *divRhoyDn       = cpscr->cpscr_grho.d_rhoy_dn;
  double *divRhozDn       = cpscr->cpscr_grho.d_rhoz_dn;
  double *d2RhoDn        = cpscr->cpscr_grho.d2_rho_dn;

/*==========================================================================*/
/* I) Read in Real Space Density                */

  if(myidState==0){
    volCP = getdeth(hmatCP);

    rhoUpRead = (double*)cmalloc(rhoRealGridTot*sizeof(double));
    rhoDnRead = (double*)cmalloc(rhoRealGridTot*sizeof(double));
    densityReadFile = fopen(densityReadFileName,"r");
    //printf("%s\n",densityReadFileName);
    if(densityReadFile==NULL){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("Density input file %s does not exist!\n",densityReadFileName);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }
    else{
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	fscanf(densityReadFile,"%lg",&rhoUpRead[iGrid]);
        //printf("111 %lg\n",rhoUpRead[iGrid]);
	rhoUpRead[iGrid] *= volCP;
      }
      //printf("volCP %lg rhoUp %lg\n",volCP,rhoUpRead[0]);
      if(cpLsda==1&&numStateDnProc!=0){
	for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	  fscanf(densityReadFile,"%lg",&rhoDnRead[iGrid]);
	  rhoDnRead[iGrid] *= volCP;
	}
      }
      fclose(densityReadFile);
    }    
  }

/*==========================================================================*/
/* II) Scatter the density                */

  if(numProcStates>1){
    Barrier(comm_states);
    Scatterv(rhoUpRead,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
             rhoUpCorrect,rhoRealGridNum,MPI_DOUBLE,0,comm_states);
  }
  else{
    memcpy(rhoUpCorrect,rhoUpRead,rhoRealGridNum*sizeof(double));
  }
  if(cpLsda==1&&numStateDnProc!=0){
    rhoDnRead = (double*)cmalloc(rhoRealGridTot*sizeof(double));
    if(numProcStates>1){
      Scatterv(rhoDnRead,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
               rhoDnCorrect,rhoRealGridNum,MPI_DOUBLE,0,comm_states);
    }
    else{
      memcpy(rhoDnCorrect,rhoDnRead,rhoRealGridNum*sizeof(double));
    }
  }
  
  if(myidState==0){
    free(rhoUpRead);
    free(rhoDnRead);
  }

/*==========================================================================*/
/* II) Initial step of density mixing                */

  genDensityMix(cp,0);
  memcpy(&rhoUp[1],rhoUpCorrect,rhoRealGridNum*sizeof(double));
  if(cpLsda==1&&numStateDnProc!=0){
    memcpy(&rhoDn[1],rhoDnCorrect,rhoRealGridNum*sizeof(double));
  }

/*==========================================================================*/
/* III) Calculate Reciprocal Density                */

  if(numProcStates>1)Barrier(comm_states);


  calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
                     rhoCoeffReUp,rhoCoeffImUp,rhoUp,rhoCoeffReUpDensCpBox,rhoCoeffImUpDensCpBox,
                     divRhoxUp,divRhoyUp,divRhozUp,d2RhoUp,cpGGA,cpDualGridOptOn,numInterpPmeDual,
                     commCP,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));
  if(cpLsda==1&&numStateDnProc!=0){
    calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
                       rhoCoeffReDn,rhoCoeffImDn,rhoDn,rhoCoeffReDnDensCpBox,rhoCoeffImDnDensCpBox,
                       divRhoxDn,divRhoyDn,divRhozDn,d2RhoDn,cpGGA,cpDualGridOptOn,numInterpPmeDual,
                       commCP,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));
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


/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/

