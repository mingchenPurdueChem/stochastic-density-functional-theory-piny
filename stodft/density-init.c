/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: density-init.c                                 */
/*                                                                          */
/* This routine costruct P_N(H)|phi> where P_N is some                      */
/* polynomial.                                                              */
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
#define TIME_CP_OFF

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

  STODFTINFO *stodftInfo       = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos = cp->stodftCoefPos;
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CPOPTS *cpopts                = &(cp->cpopts);
  PSEUDO *pseudo                = &(cp->pseudo);


  int reInitFlag = stodftInfo->reInitFlag;
  int vpsAtomListFlag = stodftInfo->vpsAtomListFlag;
  int cpDualGridOptOn           = cpopts->cp_dual_grid_opt;


/*==========================================================================*/
/* I) Generate initial density, stochastic Case                             */

  if(reInitFlag==0) calcRhoStoInit(class,bonded,general_data,cp,cpcoeffs_pos);
  if(reInitFlag==1) calcRhoDetInit(class,bonded,general_data,cp,cpcoeffs_pos);

/*==========================================================================*/
/* II) Calculate the non-local pseudopotential list                          */

  if(stodftInfo->vpsAtomListFlag==0||cpDualGridOptOn>= 1){
    control_vps_atm_list(pseudo,cell,clatoms_pos,clatoms_info,
                         atommaps,ewd_scr,for_scr,cpDualGridOptOn,
                         stodftInfo->vpsAtomListFlag);
    stodftInfo->vpsAtomListFlag = 1;
  }

/*==========================================================================*/
/* III) Change the occupation number after reading initial deisity          */
/*      density calculation                                                 */

  if(cpLsda==1)stodftInfo->occNumber = 1;
  else stodftInfo->occNumber = 2;

  printf("Finish generating Pseudopotential list.\n");


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
  int iCoeff;
  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);


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


/*======================================================================*/
/* III) Calculate the density		                        */


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

  for(iCoeff=1;iCoeff<=10;iCoeff++)printf("i %i rhocr_up %lg rhoci_up %lg\n",iCoeff,rhoCoeffReUp[iCoeff],rhoCoeffImUp[iCoeff]);
  

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
  int numCoeffLarge       = cpcoeffs_info->ncoef_l;
  int numCoeffLargeProc   = cp->cp_para_fft_pkg3d_lg.ncoef_proc;
  int numCoeffLargeProcDensCpBox = cp->cp_para_fft_pkg3d_dens_cp_box.ncoef_proc;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numInterpPmeDual = pseudo->n_interp_pme_dual;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
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
  double *rhoTemp = (double*)cmalloc();


/*======================================================================*/
/* III) Calculate the density                                           */


  //Debug Flag: we temp do this for debug
  if(cpParaOpt==0){
    if(np_states>1)rhoScr = cpscr->cpscr_rho.v_ks_up;
    else rhoScr = rho;

    rhoCalcRealStoHybrid(cpscr,cpcoeffs_info,
                   cell,stodftInfo,stoWfUpRe[iChem],
                   stoWfUpIm[iChem],*coefFormUp,*coefOrthUp,rhoScr,
                   numStateUpProc,numCoeff,cpDualGridOptOn,commCP,
                   &(cp->cp_sclr_fft_pkg3d_lg),
                   &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
                   &(cp->cp_sclr_fft_pkg3d_sm));
    Reduce(&rhoScr[1],rhoRealGridNum,rhoTemp,);
    calcRhoStoRecipHybrid()

  }
  if(cpLsda==1&&numStateDnProc>0){
    if(cpParaOpt==0){
      cp_rho_calc_sto_hybrid(cpewald,cpscr,cpcoeffs_info,
	       ewald,cell,stodftInfo,coeffReDn,
	       coeffImUp,*coefFormDn,*coefOrthDn,
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
      cp_rho_calc_sto_full_g(cpewald,cpscr,cpcoeffs_info,
                           ewald,cell,stodftInfo,coeffReDn,
                           coeffImUp,*coefFormDn,*coefOrthDn,
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

/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/

