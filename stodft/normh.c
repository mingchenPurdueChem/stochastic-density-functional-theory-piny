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

#define TIME_CP_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void normHNewton(STODFTCOEFPOS *stodftCoeffPos,STODFTINFO *stodftInfo,
		 CPCOEFFS_POS *cpcoeffs_pos,CPCOEFFS_INFO *cpcoeffs_info,
		 CELL *cell,CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
		 CPOPTS *cpopts,CPEWALD *cpewald,CPSCR *cpscr,PSEUDO *pseudo,
		 EWALD *ewald,EWD_SCR *ewd_scr,ATOMMAPS *atommaps)


/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Let's first build a simple version. We will first apply this on       */
/* Hermitian operator. So all the sample points will locate on the real  */
/* axis.								 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int cpLsda         = cpopts->cp_lsda;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int iCoeff;
  
  double *expanCoeff = stodftCoefPos->expanCoeff;
  double *wfInReUp   = stodftCoefPos->wfInReUp;
  double *wfInImUp   = stodftCoefPos->wfInImUp;
  double *wfInReDn   = stodftCoefPos->wfInReDn;
  double *wfInImDn   = stodftCoefPos->wfInImDn;
  double *wfOutReUp  = stodftCoefPos->wfOutReUp;
  double *wfOutImUp  = stodftCoefPos->wfOutImUp;
  double *wfOutReDn  = stodftCoefPos->wfOutReDn;
  double *wfOutImDn  = stodftCoefPos->wfOutImDn;
  
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;

/*==========================================================================*/
/* 0) Copy the input wave function to CP coeff and zero the force */
 
  memcpy(&(cre_up[1]),wfInReUp,numCoeffUpTotal);
  memcpy(&(cim_up[1]),wfInImUp,numCoeffUpTotal);
  for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
    fcre_up[iCoeff] = 0.0;
    fcim_up[iCoeff] = 0.0;
  }
  if(cpLsda==1&&numStateDnProc!=0){
    memcpy(&(cre_dn[1]),wfInReDn,numCoeffDnTotal);
    memcpy(&(cim_dn[1]),wfInImDn,numCoeffDnTotal);
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      fcre_up[iCoeff] = 0.0;
      fcim_up[iCoeff] = 0.0;
    }
  }
  
/*==========================================================================*/
/* 1) Calculate the H|phi> */
  //control_vps_atm_list will be done somewhere else

  control_cp_eext_recip(clatoms_info,clatoms_pos,cpcoeffs_info,
                       cpcoeffs_pos,cpewald,cpscr,
                       cpopts,pseudo,ewd_scr,atommaps,cell,
                       ewald,ptens,&(stat_avg->vrecip),
                       &(stat_avg->cp_enl),&(cp->communicate),for_scr,
                       cp_dual_grid_opt_on,
                       &(cp->cp_para_fft_pkg3d_lg));

  coef_force_control(&(cp->cpopts),&(cp->cpcoeffs_info),
                            &(cp->cpcoeffs_pos[ip_now]),
                            &(cp->cpscr),ewald,&(cp->cpewald),cell,stat_avg,
                            cp->pseudo.vxc_typ,ptens->pvten_tmp,
                            cp->pseudo.gga_cut,cp->pseudo.alpha_conv_dual,
                            cp->pseudo.n_interp_pme_dual,cp_min_on,
                            &(cp->communicate),
                            &(cp->cp_comm_state_pkg_up),
                            &(cp->cp_comm_state_pkg_dn),
                            &(cp->cp_para_fft_pkg3d_lg),
                            &(cp->cp_sclr_fft_pkg3d_lg),
                            &(cp->cp_para_fft_pkg3d_dens_cp_box),
                            &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
                            &(cp->cp_para_fft_pkg3d_sm),
                            &(cp->cp_sclr_fft_pkg3d_sm),
                            cp_dual_grid_opt_on);
  
  
  
 
  
/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/

