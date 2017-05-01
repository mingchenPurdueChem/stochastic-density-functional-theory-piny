/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: energy-wrapper-post-scf.c                      */
/*                                                                          */
/* This routine wrapps all functions used within SCF. Nuclei forces are not */
/* calculated.								    */
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


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcLocExtPostScf(CLASS *class,GENERAL_DATA *general_data,
		 CP *cp,CPCOEFFS_POS  *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the wrapper to calculate Kohn-Sham potential	after the SCF    */
/* loop. Neuclei forces are also evaluated here.			 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"
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

  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg             = &(cp->cp_para_fft_pkg3d_lg);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_up          = &(cp->cp_comm_state_pkg_up);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_dn          = &(cp->cp_comm_state_pkg_dn);

  int idual_switch;
  int pme_on = cpscr->cpscr_atom_pme.pme_on;
  int cp_dual_grid_opt  = cpopts->cp_dual_grid_opt;

  // keep these for debug
  int myid_state	= communicate->myid_state;
  int np_states		= communicate->np_states;
  int np_forc		= communicate->np_forc;

  MPI_Comm comm_states = communicate->comm_states;

  double vrecip; //keep
  double *vrecip_ret = &(stat_avg->vrecip); //keep
/*======================================================================*/
/* VI) Perform the ewald sum/ cp local pseudopotential calculation      */

  idual_switch = 0; // cp_dual_grid_opt<=1 : get vrecip vext on dense grid
                    // cp_dual_grid_opt==2 : get vrecip vext on sparse grid
  vrecip = 0.0;
  if(pme_on==1&&cp_dual_grid_opt==2){
    //Need to change this in future
    control_ewd_loc_pme(clatoms_info, clatoms_pos, cell, ptens, ewald, cpewald,
                    cpscr, pseudo, ewd_scr, cpopts, atommaps,
                    &vrecip, &(cpcoeffs_info->pseud_hess_loc),communicate,
                    for_scr,cp_dual_grid_opt,idual_switch,
                    cp_para_fft_pkg3d_lg);
  }
  else{
    control_ewd_loc(clatoms_info,clatoms_pos,cell,ptens,ewald,cpewald,
			cpscr,pseudo,ewd_scr,cpopts,atommaps,
			&vrecip,&(cpcoeffs_info->pseud_hess_loc),communicate,
			for_scr,cp_dual_grid_opt,idual_switch);
  }//endif
  if(cp_dual_grid_opt== 2){
    idual_switch = 1; /*get vext on small dense grid */
    control_ewd_loc(clatoms_info,clatoms_pos,cell,ptens,ewald,cpewald,
		        cpscr,pseudo,ewd_scr,cpopts,atommaps,
                        &vrecip,&(cpcoeffs_info->pseud_hess_loc),communicate,
		        for_scr,cp_dual_grid_opt,idual_switch);
  }//endif cp_dual_grid_opt
  *vrecip_ret += vrecip;
/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcNlPseudoPostScf(CLASS *class,GENERAL_DATA *general_data,
                 CP *cp,CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the wrapper to calculate the H|phi> without calculating K-S   */
/* potential. This part comes from control_cp_eext_recip	         */
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

  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg             = &(cp->cp_para_fft_pkg3d_lg);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_up          = &(cp->cp_comm_state_pkg_up);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_dn          = &(cp->cp_comm_state_pkg_dn);


  int i,j,iii,igh;
  int nlmtot,ntot_up,ntot_dn;
  int nl_max_kb,np_nlmax_kb,nl_max_gh,np_nlmax_gh,nl_max_all,np_nlmax_all;
  int iperd     = cell->iperd;
  int cp_ptens  = cpopts->cp_ptens_calc;
  int n_ang_max         = pseudo->n_ang_max;
  int n_ang_max_kb      = pseudo->n_ang_max_kb;
  int n_ang_max_gh      = pseudo->n_ang_max_gh;
  int n_rad_max         = pseudo->n_rad_max;
  int natm_tot          = clatoms_info->natm_tot;
  int hess_size;
  int hess_calc = clatoms_info->hess_calc;
  int cp_lsda   = cpopts->cp_lsda;
  int icoef_orth_up  = cpcoeffs_pos->icoef_orth_up;
  int icoef_form_up  = cpcoeffs_pos->icoef_form_up;
  int ifcoef_form_up = cpcoeffs_pos->ifcoef_form_up;
  int icoef_orth_dn  = cpcoeffs_pos->icoef_orth_dn;
  int icoef_form_dn  = cpcoeffs_pos->icoef_form_dn;
  int ifcoef_form_dn = cpcoeffs_pos->ifcoef_form_dn;
  int nstate_up = cpcoeffs_info->nstate_up_proc;
  int nstate_dn = cpcoeffs_info->nstate_dn_proc;
  int myid_state  = communicate->myid_state;
  int np_states   = communicate->np_states;
  int np_forc     = communicate->np_forc;
  int ncoef	  = cpcoeffs_info->ncoef;

  MPI_Comm comm_states = communicate->comm_states;

  int *np_nl            = pseudo->np_nl;
  int *np_nl_gh         = pseudo->np_nl_gh;

  double vrecip,cp_enl,cp_enl_gh;

  double *pvten = ptens->pvten_tmp;
  double *vnlreal_up    = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up    = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn    = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn    = cpscr->cpscr_nonloc.vnlim_dn;
  double *dvnlreal_x_up = cpscr->cpscr_nonloc.dvnlre_x_up;
  double *dvnlreal_y_up = cpscr->cpscr_nonloc.dvnlre_y_up;
  double *dvnlreal_z_up = cpscr->cpscr_nonloc.dvnlre_z_up;
  double *dvnlimag_x_up = cpscr->cpscr_nonloc.dvnlim_x_up;
  double *dvnlimag_y_up = cpscr->cpscr_nonloc.dvnlim_y_up;
  double *dvnlimag_z_up = cpscr->cpscr_nonloc.dvnlim_z_up;
  double *dvnlreal_x_dn = cpscr->cpscr_nonloc.dvnlre_x_dn;
  double *dvnlreal_y_dn = cpscr->cpscr_nonloc.dvnlre_y_dn;
  double *dvnlreal_z_dn = cpscr->cpscr_nonloc.dvnlre_z_dn;
  double *dvnlimag_x_dn = cpscr->cpscr_nonloc.dvnlim_x_dn;
  double *dvnlimag_y_dn = cpscr->cpscr_nonloc.dvnlim_y_dn;
  double *dvnlimag_z_dn = cpscr->cpscr_nonloc.dvnlim_z_dn;
  double *dvnlreal_gxgx_up = cpscr->cpscr_nonloc.dvnlre_gxgx_up;
  double *dvnlimag_gxgx_up = cpscr->cpscr_nonloc.dvnlim_gxgx_up;
  double *dvnlreal_gygy_up = cpscr->cpscr_nonloc.dvnlre_gygy_up;
  double *dvnlimag_gygy_up = cpscr->cpscr_nonloc.dvnlim_gygy_up;
  double *dvnlreal_gzgz_up = cpscr->cpscr_nonloc.dvnlre_gzgz_up;
  double *dvnlimag_gzgz_up = cpscr->cpscr_nonloc.dvnlim_gzgz_up;
  double *dvnlreal_gxgy_up = cpscr->cpscr_nonloc.dvnlre_gxgy_up;
  double *dvnlimag_gxgy_up = cpscr->cpscr_nonloc.dvnlim_gxgy_up;
  double *dvnlreal_gygz_up = cpscr->cpscr_nonloc.dvnlre_gygz_up;
  double *dvnlimag_gygz_up = cpscr->cpscr_nonloc.dvnlim_gygz_up;
  double *dvnlreal_gxgz_up = cpscr->cpscr_nonloc.dvnlre_gxgz_up;
  double *dvnlimag_gxgz_up = cpscr->cpscr_nonloc.dvnlim_gxgz_up;
  double *dvnlreal_gxgx_dn = cpscr->cpscr_nonloc.dvnlre_gxgx_dn;
  double *dvnlimag_gxgx_dn = cpscr->cpscr_nonloc.dvnlim_gxgx_dn;
  double *dvnlreal_gygy_dn = cpscr->cpscr_nonloc.dvnlre_gygy_dn;
  double *dvnlimag_gygy_dn = cpscr->cpscr_nonloc.dvnlim_gygy_dn;
  double *dvnlreal_gzgz_dn = cpscr->cpscr_nonloc.dvnlre_gzgz_dn;
  double *dvnlimag_gzgz_dn = cpscr->cpscr_nonloc.dvnlim_gzgz_dn;
  double *dvnlreal_gxgy_dn = cpscr->cpscr_nonloc.dvnlre_gxgy_dn;
  double *dvnlimag_gxgy_dn = cpscr->cpscr_nonloc.dvnlim_gxgy_dn;
  double *dvnlreal_gygz_dn = cpscr->cpscr_nonloc.dvnlre_gygz_dn;
  double *dvnlimag_gygz_dn = cpscr->cpscr_nonloc.dvnlim_gygz_dn;
  double *dvnlreal_gxgz_dn = cpscr->cpscr_nonloc.dvnlre_gxgz_dn;
  double *dvnlimag_gxgz_dn = cpscr->cpscr_nonloc.dvnlim_gxgz_dn;
  double *fx_tmp           = ewd_scr->fx;
  double *fy_tmp           = ewd_scr->fy;
  double *fz_tmp           = ewd_scr->fz;
  double *fx               = clatoms_pos->fx;
  double *fy               = clatoms_pos->fy;
  double *fz               = clatoms_pos->fz;
  double *x                = clatoms_pos->x;
  double *y                = clatoms_pos->y;
  double *z                = clatoms_pos->z;
  double *q                = clatoms_info->q;
  double *hess_xx         = clatoms_pos->hess_xx;
  double *hess_xy         = clatoms_pos->hess_xy;
  double *hess_xz         = clatoms_pos->hess_xz;
  double *hess_yy         = clatoms_pos->hess_yy;
  double *hess_yz         = clatoms_pos->hess_yz;
  double *hess_zz         = clatoms_pos->hess_zz;

  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;

  double *cp_enl_ret = &(stat_avg->cp_enl);

/*======================================================================*/
/* II) Malloc a hessian scratch vector if necessary                     */

  if(hess_calc == 3 && np_states > 1){
    hess_size = natm_tot*natm_tot;
  }

/*======================================================================*/
/* III) Determine the maximum open non-local angular momentum channel   */
/*      for Kleinman-Bylander and Goedecker pseudo potentials           */

  nl_max_kb = -1;
  for(i=1;i<=(n_ang_max_kb+1);i++){
   if(np_nl[i]>0){nl_max_kb=i-1;}
  }/*endfor*/

  nl_max_gh = -1;
  for(i=1;i<=(n_ang_max_gh+1);i++){
   if(np_nl_gh[i]>0){nl_max_gh=i-1;}
  }/*endfor*/

/*======================================================================*/
/* IV) Determine the maximum number of atoms in any                     */
/*       open angular momentum channel                                  */

  np_nlmax_kb = 1;
  for(i = 1;i<=(nl_max_kb+1);i++){
   np_nlmax_kb = MAX(np_nlmax_kb,np_nl[i]);
  }/*endfor*/

  np_nlmax_gh = 1;
  for(i = 1;i<=(nl_max_gh+1);i++){
   np_nlmax_gh = MAX(np_nlmax_gh,np_nl_gh[i]);
  }/*endfor*/

  np_nlmax_all = (np_nlmax_gh > np_nlmax_kb ? np_nlmax_gh : np_nlmax_kb);

/*======================================================================*/
/* V) Zero the non-local tensors                                        */

  nl_max_all = (nl_max_gh > nl_max_kb ? nl_max_gh : nl_max_kb);
  nlmtot = (nl_max_all+1)*(nl_max_all+1);

  ntot_up = nstate_up*np_nlmax_all*nlmtot*n_rad_max;
  ntot_dn = 0;

  if(nl_max_all>= 0){
    for(i=1;i<=ntot_up;i++){
      vnlreal_up[i] = 0.0;
      vnlimag_up[i] = 0.0;
      dvnlreal_x_up[i] = 0.0;
      dvnlreal_y_up[i] = 0.0;
      dvnlreal_z_up[i] = 0.0;
      dvnlimag_x_up[i] = 0.0;
      dvnlimag_y_up[i] = 0.0;
      dvnlimag_z_up[i] = 0.0;
    }//endfor
    if(cp_ptens==1||hess_calc==3){
      for(i=1;i<=ntot_up;i++){
        dvnlreal_gxgx_up[i] = 0.0;
        dvnlreal_gzgz_up[i] = 0.0;
        dvnlreal_gygy_up[i] = 0.0;
        dvnlreal_gxgy_up[i] = 0.0;
        dvnlreal_gxgz_up[i] = 0.0;
        dvnlreal_gygz_up[i] = 0.0;

        dvnlimag_gxgx_up[i] = 0.0;
        dvnlimag_gxgy_up[i] = 0.0;
        dvnlimag_gygy_up[i] = 0.0;
        dvnlimag_gxgz_up[i] = 0.0;
        dvnlimag_gygz_up[i] = 0.0;
        dvnlimag_gzgz_up[i] = 0.0;
      }//endfor
    }//endif:ptens
    if(cp_lsda==1){
      ntot_dn = nstate_dn*np_nlmax_all*nlmtot*n_rad_max;
      for(i=1;i<=ntot_dn;i++){
        vnlreal_dn[i]    = 0.0;
        vnlimag_dn[i]    = 0.0;
        dvnlreal_x_dn[i] = 0.0;
        dvnlreal_y_dn[i] = 0.0;
        dvnlreal_z_dn[i] = 0.0;
        dvnlimag_x_dn[i] = 0.0;
        dvnlimag_y_dn[i] = 0.0;
        dvnlimag_z_dn[i] = 0.0;
      }//endfor
      if(cp_ptens==1 || hess_calc == 3){
        for(i=1;i<=ntot_dn;i++){
	  dvnlreal_gxgx_dn[i] = 0.0;
	  dvnlreal_gxgy_dn[i] = 0.0;
	  dvnlreal_gxgz_dn[i] = 0.0;
	  dvnlreal_gygy_dn[i] = 0.0;
	  dvnlreal_gygz_dn[i] = 0.0;
	  dvnlreal_gzgz_dn[i] = 0.0;

	  dvnlimag_gxgx_dn[i] = 0.0;
	  dvnlimag_gxgy_dn[i] = 0.0;
	  dvnlimag_gxgz_dn[i] = 0.0;
	  dvnlimag_gygy_dn[i] = 0.0;
	  dvnlimag_gygz_dn[i] = 0.0;
	  dvnlimag_gzgz_dn[i] = 0.0;
        }//endfor
      }//endif:ptens
    }//endif:lsda
  }//endif : non-local potential on*/
/*======================================================================*/
/* VII) Get the nl pe, pvten and particle forces then the coef forces   */

   cp_enl    = 0.0;
   cp_enl_gh = 0.0;

 if((nl_max_all >= 0) && (n_rad_max>1)){
      non_loc_chng_ord(clatoms_pos,clatoms_info,
                       atommaps, pseudo,ewd_scr,for_scr,1);
 }/*endif*/

/*-------------------------------------------------------------------------*/
/* A) KB/Goedecker NLs */

  if( (nl_max_kb >= 0) && ((ntot_up+ntot_dn)>0) ){
    control_ewd_non_loc(clatoms_info,clatoms_pos,cpcoeffs_info,cpcoeffs_pos,
                        cell,ptens,cpewald,cpscr,pseudo,ewd_scr,
                        cpopts,atommaps,communicate,for_scr);
  }else{
    get_ak2_sm(cpewald,cell);
  }/*endif*/
  
  if((nl_max_kb >= 0)&&((ntot_up+ntot_dn)>0)&&(pseudo->np_nonloc_cp_box_kb>0) ){
    getnl_pot_pv_fatm(clatoms_info,clatoms_pos,cell,cpcoeffs_info,
                      cpscr,ewd_scr,cpopts,pseudo,atommaps,&cp_enl,
                      np_nlmax_kb,pvten);
    getnl_fcoef(clatoms_info,clatoms_pos,cpcoeffs_info,cpcoeffs_pos,
                  cpscr,ewd_scr,cpopts,pseudo,cpewald,atommaps,
                  cell,np_nlmax_kb,pvten,for_scr);

  }/*endif*/

/*-------------------------------------------------------------------------*/
/* B) Gauss-Hermite NLs                                                    */

  if((nl_max_gh >= 0)&&((ntot_up+ntot_dn)>0)&&(pseudo->np_nonloc_cp_box_gh>0)){
      if(cp_ptens==1){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf(" CP-PTENS is not implemented for Gauss-Hermite nonlocality \n");
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
  }/* endif */

 for(igh=1;igh<=pseudo->ngh;igh++){

  cp_enl_gh = 0.0;

  if(nl_max_all >= 0){
   for(i=1;i<=ntot_up;i++){
    vnlreal_up[i] = 0.0;
    vnlimag_up[i] = 0.0;
    dvnlreal_x_up[i] = 0.0;
    dvnlreal_y_up[i] = 0.0;
    dvnlreal_z_up[i] = 0.0;
    dvnlimag_x_up[i] = 0.0;
    dvnlimag_y_up[i] = 0.0;
    dvnlimag_z_up[i] = 0.0;
   }/*endfor*/

   if(cp_ptens==1 || hess_calc == 3){
    for(i=1;i<=ntot_up;i++){
     dvnlreal_gxgx_up[i] = 0.0;
     dvnlreal_gzgz_up[i] = 0.0;
     dvnlreal_gygy_up[i] = 0.0;
     dvnlreal_gxgy_up[i] = 0.0;
     dvnlreal_gxgz_up[i] = 0.0;
     dvnlreal_gygz_up[i] = 0.0;

     dvnlimag_gxgx_up[i] = 0.0;
     dvnlimag_gxgy_up[i] = 0.0;
     dvnlimag_gygy_up[i] = 0.0;
     dvnlimag_gxgz_up[i] = 0.0;
     dvnlimag_gygz_up[i] = 0.0;
     dvnlimag_gzgz_up[i] = 0.0;
    }/*endfor*/
   }/*endif:ptens*/

   if(cp_lsda==1){

    ntot_dn = nstate_dn*np_nlmax_all*nlmtot*n_rad_max;
    for(i=1;i<=ntot_dn;i++){
     vnlreal_dn[i]    = 0.0;
     vnlimag_dn[i]    = 0.0;
     dvnlreal_x_dn[i] = 0.0;
     dvnlreal_y_dn[i] = 0.0;
     dvnlreal_z_dn[i] = 0.0;
     dvnlimag_x_dn[i] = 0.0;
     dvnlimag_y_dn[i] = 0.0;
     dvnlimag_z_dn[i] = 0.0;
    }/*endfor*/

    if(cp_ptens==1 || hess_calc == 3){
     for(i=1;i<=ntot_dn;i++){
      dvnlreal_gxgx_dn[i] = 0.0;
      dvnlreal_gxgy_dn[i] = 0.0;
      dvnlreal_gxgz_dn[i] = 0.0;
      dvnlreal_gygy_dn[i] = 0.0;
      dvnlreal_gygz_dn[i] = 0.0;
      dvnlreal_gzgz_dn[i] = 0.0;

      dvnlimag_gxgx_dn[i] = 0.0;
      dvnlimag_gxgy_dn[i] = 0.0;
      dvnlimag_gxgz_dn[i] = 0.0;
      dvnlimag_gygy_dn[i] = 0.0;
      dvnlimag_gygz_dn[i] = 0.0;
      dvnlimag_gzgz_dn[i] = 0.0;
     }/*endfor*/
    }/*endif:ptens*/

   }/*endif:lsda*/
  }/*endif : non-local potential on*/

   /* Create Zilm(alpha) and Xilm(alpha)= dZilm/dR  */

       control_nonloc_gh(clatoms_info,clatoms_pos,cpcoeffs_info,cpcoeffs_pos,
                         cell,ptens,cpewald,cpscr,pseudo,ewd_scr,
                         cpopts,atommaps,communicate,for_scr,
                         pseudo->rgh[igh]);

 /*i) Calculate the nonlocal Energy  */
 /*ii) and the contribution to the forces on the atoms due to the non-local term*/

       getnl_pot_pv_fatm_gh(clatoms_info,clatoms_pos,cell,cpcoeffs_info,
                            cpscr,ewd_scr,cpopts,pseudo,atommaps,
                            &cp_enl_gh,np_nlmax_gh,pvten,pseudo->wgh,igh);
   /*Calculate nonlocal contribution to the coefficient forces */

       getnl_fcoef_gh(clatoms_info,clatoms_pos,cpcoeffs_info,cpcoeffs_pos,
                      cpscr,ewd_scr,cpopts,pseudo,cpewald,atommaps,
                      cell,np_nlmax_gh,pvten,for_scr,
                      pseudo->rgh[igh],pseudo->wgh,igh);

        *cp_enl_ret += cp_enl_gh;
      }/*endfor igh gauss-hermite integration points */

    }/*endif gauss-hermit*/

    if((nl_max_all >= 0)&&((ntot_up+ntot_dn)>0)&&(n_rad_max>1)){
      non_loc_restore_ord(clatoms_pos,clatoms_info,
                          atommaps, pseudo,ewd_scr,for_scr);
    }/*endif*/

/*======================================================================*/
/* VIII) Assign the potential energy                                    */

  *cp_enl_ret += cp_enl;

/*======================================================================*/
/* IX) Reduce particle forces if necessary                              */

  if(np_states>1 && np_forc == 1){
    reduce_cp_atm_forc(natm_tot,fx,fy,fz,fx_tmp,fy_tmp,fz_tmp,
                       comm_states,myid_state);
  }/* endif:npstates */

/*======================================================================*/
/* X) Reduce particle hessian if necessary                              */

  if(hess_calc == 3 && np_states>1){
    reduce_cp_hess_stuff(hess_xx,hess_yy,hess_zz,
                         hess_xy,hess_xz,hess_yz,hess_size,
                         myid_state,comm_states);
  }/*endif*/
 
 //for(i=1;i<=ncoef;i++)printf("coef_force_eext %.10lg %.10lg\n",fcre_up[i],fcim_up[i]);


/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcCoefForcePosScf(CLASS *class,GENERAL_DATA *general_data,
                   CP *cp,CPCOEFFS_POS  *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the wrapper to calculate the H|phi> without calculating K-S   */
/* potential. The potential is calculated in calcKSPotWrap so that we    */
/* would have one less FFT for every step H|phi>.                        */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
/*==========================================================================*/
/* III) get the force on the states (up and down)                           */
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

  int i,iii;
  int cp_lda            = cpopts->cp_lda;
  int cp_lsda           = cpopts->cp_lsda;
  int cp_sic            = cpopts->cp_sic;
  int cp_nonint         = cpopts->cp_nonint;
  int cp_ptens_calc     = cpopts->cp_ptens_calc;
  int cp_gga            = cpopts->cp_gga;
  int cp_para_opt       = cpopts->cp_para_opt;
  int nstate_up         = cpcoeffs_info->nstate_up_proc;
  int nstate_dn         = cpcoeffs_info->nstate_dn_proc;
  int nstate_up_tot     = cpcoeffs_info->nstate_up;
  int nstate_dn_tot     = cpcoeffs_info->nstate_dn;
  int laplacian_on      = cpcoeffs_info->cp_laplacian_on;
  int cp_tau_functional = cpcoeffs_info->cp_tau_functional;
  int alpha_conv_dual   = pseudo->alpha_conv_dual;
  int n_interp_pme_dual = pseudo->n_interp_pme_dual;
  int cp_dual_grid_opt  = cpopts->cp_dual_grid_opt;
  int cp_min_on = 0; //I don't want to calculate cp_hess
  int myid_state           =  communicate->myid_state;
  int np_states            =  communicate->np_states;
  int nstate_ncoef_proc_max_up = cpcoeffs_info->nstate_ncoef_proc_max_up;
  int nstate_ncoef_proc_max_dn = cpcoeffs_info->nstate_ncoef_proc_max_dn;
  int nstate_ncoef_proc_up     = cpcoeffs_info->nstate_ncoef_proc_up;
  int nstate_ncoef_proc_dn     = cpcoeffs_info->nstate_ncoef_proc_dn;
  int icoef_orth_up          =  cpcoeffs_pos->icoef_orth_up;
  int icoef_form_up          =  cpcoeffs_pos->icoef_form_up;
  int ifcoef_form_up         =  cpcoeffs_pos->ifcoef_form_up;
  int icoef_orth_dn          =  cpcoeffs_pos->icoef_orth_dn;
  int icoef_form_dn          =  cpcoeffs_pos->icoef_form_dn;
  int ifcoef_form_dn         =  cpcoeffs_pos->ifcoef_form_dn;
  char *vxc_typ              =  pseudo->vxc_typ;

  double gc_cut         = pseudo->gga_cut;
  double cp_eke_dn,cp_eke;

  double *creal_up         =  cpcoeffs_pos->cre_up;
  double *cimag_up         =  cpcoeffs_pos->cim_up;
  double *cimag_dn         =  cpcoeffs_pos->cim_dn;
  double *creal_dn         =  cpcoeffs_pos->cre_dn;
  double *fcreal_up        =  cpcoeffs_pos->fcre_up;
  double *fcimag_up        =  cpcoeffs_pos->fcim_up;
  double *fcimag_dn        =  cpcoeffs_pos->fcim_dn;
  double *fcreal_dn        =  cpcoeffs_pos->fcre_dn;
  double *kfcre_up         =  cpcoeffs_pos->kfcre_up;
  double *kfcim_up         =  cpcoeffs_pos->kfcim_up;
  double *cp_hess_re_up    =  cpcoeffs_pos->cp_hess_re_up;
  double *cp_hess_im_up    =  cpcoeffs_pos->cp_hess_im_up;
  double *cp_hess_re_dn    =  cpcoeffs_pos->cp_hess_re_dn;
  double *cp_hess_im_dn    =  cpcoeffs_pos->cp_hess_im_dn;
  double *ak2_sm           =  cpewald->ak2_sm;
  double *v_ks_up          =  cpscr->cpscr_rho.v_ks_up;
  double *v_ks_dn          =  cpscr->cpscr_rho.v_ks_dn;
  double *v_ks_tau_up      =  cpscr->cpscr_rho.v_ks_tau_up;
  double *v_ks_tau_dn      =  cpscr->cpscr_rho.v_ks_tau_dn;
  double *zfft             =  cpscr->cpscr_wave.zfft;
  double *zfft_tmp         =  cpscr->cpscr_wave.zfft_tmp;
  double *cre_scr          =  cpscr->cpscr_wave.cre_up;
  double *cim_scr          =  cpscr->cpscr_wave.cim_up;
  double *hmati_cp         =  cell->hmati_cp;
  double *pvten_cp         =  ptens->pvten_tmp;
  double *cp_eke_ret       =  &(stat_avg->cp_eke);
  double *ks_offset        =  &(cpcoeffs_pos->ks_offset);

  MPI_Comm comm_states = communicate->comm_states;
  MPI_Comm world       = communicate->world;

  switch(cp_para_opt){

/*-------------------------------------------------------------------------*/

    case 0: /* hybrid */
 /*-----------------------------------------*/
 /* i)  Up states                           */
      coef_force_calc_hybrid(cpewald,nstate_up,creal_up,cimag_up,
                          fcreal_up,fcimag_up,cre_scr,cim_scr,cp_hess_re_up,cp_hess_im_up,
                          zfft,zfft_tmp,v_ks_up,v_ks_tau_up,ak2_sm,&cp_eke,pvten_cp,
                          cp_ptens_calc,hmati_cp,communicate,icoef_form_up,
                          icoef_orth_up,ifcoef_form_up,cp_tau_functional,cp_min_on,
                          cp_sclr_fft_pkg3d_sm);
      *cp_eke_ret += cp_eke;

 /*--------------------------------------------*/
 /* ii) down states (if necessary)             */

      if(cp_lsda == 1 && nstate_dn != 0){
        coef_force_calc_hybrid(cpewald,nstate_dn,creal_dn,cimag_dn,
                          fcreal_dn,fcimag_dn,cre_scr,cim_scr,cp_hess_re_dn,cp_hess_im_dn,
                          zfft,zfft_tmp,v_ks_dn,v_ks_tau_dn,ak2_sm,&cp_eke_dn,pvten_cp,
                          cp_ptens_calc,hmati_cp,communicate,icoef_form_dn,
                          icoef_orth_dn,ifcoef_form_dn,cp_tau_functional,cp_min_on,
                          cp_sclr_fft_pkg3d_sm);
        *cp_eke_ret += cp_eke_dn;
      }/*endif*/

    break;
/*-------------------------------------------------------------------------*/

    case 1: /* full g */
 /*-----------------------------------------*/
 /* i)  Up states                           */

      coef_force_calc_full_g(cpewald,nstate_up_tot,nstate_ncoef_proc_up,
                          nstate_ncoef_proc_max_up,
                          creal_up,cimag_up,fcreal_up,fcimag_up,cre_scr,cim_scr,
                          cp_hess_re_up,cp_hess_im_up,
                          zfft,zfft_tmp,v_ks_up,v_ks_tau_up,ak2_sm,&cp_eke,pvten_cp,
                          cp_ptens_calc,hmati_cp,communicate,icoef_form_up,
                          icoef_orth_up,ifcoef_form_up,cp_tau_functional,cp_min_on,
                          cp_para_fft_pkg3d_sm);
      *cp_eke_ret += cp_eke;

 /*--------------------------------------------*/
 /* ii) down states (if necessary)             */

      if(cp_lsda == 1 && nstate_dn != 0){
        coef_force_calc_full_g(cpewald,nstate_dn_tot,nstate_ncoef_proc_dn,
                          nstate_ncoef_proc_max_dn,
                          creal_dn,cimag_dn,fcreal_dn,fcimag_dn,cre_scr,cim_scr,
                          cp_hess_re_dn,cp_hess_im_dn,
                          zfft,zfft_tmp,v_ks_dn,v_ks_tau_dn,ak2_sm,&cp_eke_dn,pvten_cp,
                          cp_ptens_calc,hmati_cp,communicate,icoef_form_dn,
                          icoef_orth_dn,ifcoef_form_dn,cp_tau_functional,cp_min_on,
                          cp_para_fft_pkg3d_sm);
        *cp_eke_ret += cp_eke_dn;
      }/*endif*/
    break;
  }/* end switch cp_para_opt */

/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcCoefForceWrapReduce(CLASS *class,GENERAL_DATA *general_data,
                   CP *cp,CPCOEFFS_POS  *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the wrapper to calculate H|phi> given |phi>                   */
/* |phi> are stored in cre(im)_up(dn) and H|phi> are stored in           */
/* fcre(im)_up(dn), for correct H|phi>, you need to scale all coeff      */
/* with k!=0 by -0.5, and k=0 term by -1. Since this step can be         */
/* combined with some other scalings, I just output the raw force.       */
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
    for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
      fcre_dn[iCoeff] = 0.0;
      fcim_dn[iCoeff] = 0.0;
    }
  }

/*==========================================================================*/
/* 1) Calculate the H/sigma|phi> */
  //control_vps_atm_list will be done somewhere else (perhaps in density calculation?)

  //calcCoefForceExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  //calcCoefForceForceControlWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/

