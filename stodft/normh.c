/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: normh.c                                        */
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
void normHNewtonHerm(CP *cp,CLASS *class,GENERAL_DATA *general_data,
		 CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos,double zn)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Let's first build a simple version. We will first apply this on       */
/* Hermitian operator. So all the sample points will locate on the real  */
/* axis. zn is the interpolation point. For details please read		 */
/* Ashkenazi et.al. JChemPhys 103, 10005(1995).				 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  COMMUNICATE *communicate      = &(cp->communicate);

  NEWTONINFO *newtonInfo = stodftInfo->newtonInfo;
  double Smin		= newtonInfo->Smin;
  double Smax		= newtonInfo->Smax;
  double scale		= newtonInfo->scale;
  double energyMean	= stodftInfo->energyMean;
  double energyDiff	= stodftInfo->energyDiff;
  double prefact	= -scale*energyMean-zn;
  double scale1		= -scale*0.5;
  double scale2		= -scale;

  int cpLsda         = cpopts->cp_lsda;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int cpDualGridOptOn  = cpopts->cp_dual_grid_opt;
  int numCoeffM1     = numCoeff-1;
  int incx = 1;
  int incy = 1;
  int iState,iCoeff,iCoeffStart,index1,index2;
  
  double *expanCoeff = (double*)stodftCoefPos->expanCoeff;
  
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  double dot;

  
/*==========================================================================*/
/* 1) Calculate the H/sigma|phi> */
  //control_vps_atm_list will be done somewhere else (perhaps in density calculation?)

  calcCoefForceWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

  for(iState=0;iState<numStateUpProc;iState++){
    iCoeffStart = iState*numCoeff;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      index1 = iCoeffStart+iCoeff;
      fcre_up[index1] *= scale1;
      fcim_up[index1] *= scale1;	
    }//endfor iCoeff
    index1 = iCoeffStart+numCoeff;
    fcre_up[index1] *= scale2;
  }//endfor iState
  if(cpLsda==1&&numStateDnProc!=0){
    for(iState=0;iState<numStateDnProc;iState++){
      iCoeffStart = iState*numCoeff;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	index1 = iCoeffStart+iCoeff;
	fcre_dn[index1] *= scale1;
	fcim_dn[index1] *= scale1;
      }//endfor iCoeff
      index1 = iCoeffStart+numCoeff;
      fcim_dn[index1] *= scale2;
    }//endfor iState
  }

/*==========================================================================*/
/* 2) Calculate P_(n+1)(H)|phi>. Everything store in fcre(im)_up(dn) */

  DAXPY(&numCoeffUpTotal,&prefact,&cre_up[1],&incx,&fcre_up[1],&incy);
  DAXPY(&numCoeffUpTotal,&prefact,&cim_up[1],&incx,&fcim_up[1],&incy);
  if(cpLsda==1&&numStateDnProc!=0){
    DAXPY(&numCoeffDnTotal,&prefact,&cre_dn[1],&incx,&fcre_dn[1],&incy);
    DAXPY(&numCoeffDnTotal,&prefact,&cim_dn[1],&incx,&fcim_dn[1],&incy);
  }

/*==========================================================================*/
/* 3) Copy the force back to the coefficients */
  //debug
  /*
  for(iState=0;iState<numStateUpProc;iState++){
    dot = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      dot += cre_up[iState*numCoeff+iCoeff]*fcre_up[iState*numCoeff+iCoeff]+
	    cim_up[iState*numCoeff+iCoeff]*fcim_up[iState*numCoeff+iCoeff];
    }
    dot *= 2.0;
    dot += cre_up[iState*numCoeff+numCoeff]*fcre_up[iState*numCoeff+numCoeff];
    dot *= 0.5;
    printf("zn %lg dot %lg\n",zn,dot);
  } 
  */
  for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
    cre_up[iCoeff] = fcre_up[iCoeff];
    cim_up[iCoeff] = fcim_up[iCoeff];
  }
  if(cpLsda==1&&numStateDnProc!=0){
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      cre_dn[iCoeff] = fcre_dn[iCoeff];
      cim_dn[iCoeff] = fcim_dn[iCoeff];
    }
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcRhoDet(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
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
/* III) Calculate the density			                        */


  //Debug Flag: we temp do this for debug
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

  if(cpLsda==1&&numStateDnProc>0){
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
void calcRhoSto(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
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
/* III) Calculate the density                                           */


  //Debug Flag: we temp do this for debug
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

  if(cpLsda==1&&numStateDnProc>0){
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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcCoefForceWrap(CLASS *class,GENERAL_DATA *general_data,
                   CP *cp,CPCOEFFS_POS  *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the wrapper to calculate H|phi> given |phi>			 */
/* |phi> are stored in cre(im)_up(dn) and H|phi> are stored in		 */
/* fcre(im)_up(dn), for correct H|phi>, you need to scale all coeff	 */
/* with k!=0 by -0.5, and k=0 term by -1. Since this step can be	 */
/* combined with some other scalings, I just output the raw force.	 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  CELL *cell                    = &(general_data->cell);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  EWALD *ewald                  = &(general_data->ewald);
  EWD_SCR *ewd_scr              = &(class->ewd_scr);
  ATOMMAPS *atommaps            = &(class->atommaps);
  FOR_SCR *for_scr              = &(class->for_scr);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  PTENS *ptens                  = &(general_data->ptens);
  SIMOPTS *simopts              = &(general_data->simopts);


  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPEWALD *cpewald              = &(cp->cpewald);
  CPSCR *cpscr                  = &(cp->cpscr);
  PSEUDO *pseudo                = &(cp->pseudo);
  COMMUNICATE *communicate      = &(cp->communicate);

  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm             = &(cp->cp_sclr_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_sm             = &(cp->cp_para_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_dens_cp_box    = &(cp->cp_sclr_fft_pkg3d_dens_cp_box);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box    = &(cp->cp_para_fft_pkg3d_dens_cp_box);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_lg             = &(cp->cp_sclr_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg             = &(cp->cp_para_fft_pkg3d_lg);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_up          = &(cp->cp_comm_state_pkg_up);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_dn          = &(cp->cp_comm_state_pkg_dn);

  int cpLsda         = cpopts->cp_lsda;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int cpDualGridOptOn  = cpopts->cp_dual_grid_opt;
  int iState,iCoeff,iCoeffStart,index1,index2;
  int cpMinOn = 0; //I don't want to calculate cp_hess

  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;

/*==========================================================================*/
/* 0) Copy the input wave function to CP coeff and zero the force */

  for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
    fcre_up[iCoeff] = 0.0;
    fcim_up[iCoeff] = 0.0;
  }
  if(cpLsda==1&&numStateDnProc!=0){
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      fcre_dn[iCoeff] = 0.0;
      fcim_dn[iCoeff] = 0.0;
    }
  }

/*==========================================================================*/
/* 1) Calculate the H/sigma|phi> */
  //control_vps_atm_list will be done somewhere else (perhaps in density calculation?)

  control_cp_eext_recip(clatoms_info,clatoms_pos,cpcoeffs_info,
                       cpcoeffs_pos,cpewald,cpscr,cpopts,pseudo,
                       ewd_scr,atommaps,cell,ewald,ptens,&(stat_avg->vrecip),
                       &(stat_avg->cp_enl),communicate,for_scr,cpDualGridOptOn,
                       cp_para_fft_pkg3d_lg);

  coef_force_control(cpopts,cpcoeffs_info,cpcoeffs_pos,cpscr,ewald,cpewald,
                    cell,stat_avg,pseudo->vxc_typ,ptens->pvten_tmp,pseudo->gga_cut,
                    pseudo->alpha_conv_dual,pseudo->n_interp_pme_dual,cpMinOn,
                    communicate,cp_comm_state_pkg_up,
                     cp_comm_state_pkg_dn,cp_para_fft_pkg3d_lg,cp_sclr_fft_pkg3d_lg,
                     cp_para_fft_pkg3d_dens_cp_box,cp_sclr_fft_pkg3d_dens_cp_box,
                     cp_para_fft_pkg3d_sm,cp_sclr_fft_pkg3d_sm,cpDualGridOptOn);

/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/

