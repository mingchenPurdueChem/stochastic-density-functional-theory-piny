/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: cp_norb.c                                    */
/*                                                                          */
/* File contains functions necessary for computing norb force               */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*  Control which orthogonalization scheme is used                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void orthog_control_cp_dvr(CP *cp,int ip)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

 int iii;
 int icoef_orth;
 double max_off_diag=0.0,max_diag=0.0;

/* Local pointers */
  int nstate_up            = cp->cpcoeffs_info.nstate_up;
  int nstate_dn            = cp->cpcoeffs_info.nstate_dn;
  int cp_gs                = cp->cpopts.cp_gs;
  int cp_low               = cp->cpopts.cp_low;
  int cp_lsda              = cp->cpopts.cp_lsda;
  int cp_norm              = cp->cpopts.cp_normalize;
  int np_states            = cp->communicate.np_states;

  int icoef_form_up        = cp->cpcoeffs_pos_dvr[ip].icoef_form_up;
  double *dvrc_up          = cp->cpcoeffs_pos_dvr[ip].dvrc_up;
  int icoef_form_dn        = cp->cpcoeffs_pos_dvr[ip].icoef_form_dn;
  double *dvrc_dn          = cp->cpcoeffs_pos_dvr[ip].dvrc_dn;

  double *cpscr_cre_up     = cp->cpscr.cpscr_wave.cre_up;
  double *cpscr_cre_dn     = cp->cpscr.cpscr_wave.cre_dn;

  double *occ_up           = cp->cpopts.occ_up;
  double *occ_dn           = cp->cpopts.occ_dn;

  double *omat             = cp->cpscr.cpscr_ovmat.ovlap7;
  double *omat_tmp         = cp->cpscr.cpscr_ovmat.ovlap8;
  double *anorm            = cp->cpscr.cpscr_ovmat.state_vec1;
  double *anorm_tmp        = cp->cpscr.cpscr_ovmat.state_vec2;
  double *norbmat          = cp->cpscr.cpscr_ovmat.ovlap6;
  /*WARNING CPSCR_OVMAT ALSO PASSED TO ROUTINE*/
  double *norbmati         = cp->cpscr.cpscr_ovmat.ovlap7;
  double *ovmat_eigv       = cp->cpscr.cpscr_ovmat.ovlap8;

  int *ioff_upt            = cp->cpcoeffs_info.ioff_upt;
  int *ioff_dnt            = cp->cpcoeffs_info.ioff_dnt;

  double *rocc_sum_up      = cp->cpopts.rocc_sum_up;
  double *rocc_sum_dn      = cp->cpopts.rocc_sum_dn;

  double scale_fact        = cp->cpcoeffs_info.scale_fact;

/*========================================================================*/
/*  Up states */

  if(cp_gs== 1) {
    control_cp_gram_schmidt_dvr(dvrc_up,icoef_form_up,occ_up,omat,
                                omat_tmp,ioff_upt,scale_fact,
                                &(cp->cp_comm_state_pkg_dvr_up));
  }/*endif:gs*/

  if(cp_low== 1) {

    icoef_orth = 0;

    cp_rotate_coef_ortho_dvr(dvrc_up,icoef_form_up,&icoef_orth,
                             norbmat,norbmati,ovmat_eigv,
                             cpscr_cre_up,occ_up,ioff_upt,
                             &max_off_diag,&max_diag,scale_fact,
                             &(cp->cpscr.cpscr_ovmat),
                             &(cp->cp_comm_state_pkg_dvr_up));

    rotate_occ_shuffle_dvr(dvrc_up,icoef_form_up,omat,omat_tmp,
                           occ_up,rocc_sum_up,ioff_upt, scale_fact,
                           &(cp->cp_comm_state_pkg_dvr_up));

  }/*endif:lowdin*/
  if(cp_norm== 1) {
    cp_normalize_dvr(dvrc_up,icoef_form_up,occ_up,ioff_upt, 
                     anorm,anorm_tmp,&(cp->cp_comm_state_pkg_dvr_up));
  }/*endif:norm*/

/*========================================================================*/
/*  Down states */

  if( (cp_lsda == 1) && (nstate_dn != 0) ){

    if(cp_gs== 1) {
      control_cp_gram_schmidt_dvr(dvrc_dn,icoef_form_dn,occ_dn,omat,
                                  omat_tmp,ioff_dnt,scale_fact,
                                  &(cp->cp_comm_state_pkg_dvr_dn));
    }/*endif:gs*/

    if(cp_low== 1) {
      icoef_orth = 0;
      cp_rotate_coef_ortho_dvr(dvrc_dn,icoef_form_dn,&icoef_orth,
                               norbmat,norbmati,ovmat_eigv,
                               cpscr_cre_dn,occ_dn,ioff_dnt,
                               &max_off_diag,&max_diag, scale_fact,
                               &(cp->cpscr.cpscr_ovmat),
                               &(cp->cp_comm_state_pkg_dvr_dn));

      rotate_occ_shuffle_dvr(dvrc_dn,icoef_form_dn,omat,omat_tmp,
                             occ_dn,rocc_sum_dn,ioff_dnt,scale_fact,
                             &(cp->cp_comm_state_pkg_dvr_dn));
    }/*endif:low*/

    if(cp_norm== 1) {
      cp_normalize_dvr(dvrc_dn,icoef_form_dn,occ_dn,ioff_dnt,
                       anorm,anorm_tmp,&(cp->cp_comm_state_pkg_dvr_dn));
    }/*endif:norm*/

  }/* endif lsda */

/*========================================================================*/
}/* end routine */
/*========================================================================*/

