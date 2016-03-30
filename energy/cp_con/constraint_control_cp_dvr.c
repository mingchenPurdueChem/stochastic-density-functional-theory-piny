/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: constraint_control_cp                        */
/*                                                                          */
/* This subprogram controls shake and rattle calls for CP                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_cpcon_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void shake_control_cp_dvr(CP *cp,int *iter_shake, double dt,int ip,double dfact)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */


  /* Local pointers */

  int cp_norb              = cp->cpopts.cp_norb;
  int cp_lsda              = cp->cpopts.cp_lsda;
  double c_tolshake        = cp->cpconstrnt.c_tolshake;

  int    *cpcoeffs_ioff_upt = cp->cpcoeffs_info.ioff_upt;
  int    *cpcoeffs_ioff_dnt = cp->cpcoeffs_info.ioff_dnt;

  double *cpcoeffs_dvrc_up  = cp->cpcoeffs_pos_dvr[ip].dvrc_up;
  double *cpcoeffs_dvrc_dn  = cp->cpcoeffs_pos_dvr[ip].dvrc_dn;
  double *cpcoeffs_dvrvc_up = cp->cpcoeffs_pos_dvr[ip].dvrvc_up;
  double *cpcoeffs_dvrvc_dn = cp->cpcoeffs_pos_dvr[ip].dvrvc_dn;

  int icoef_form_up        = cp->cpcoeffs_pos_dvr[ip].icoef_form_up;
  int icoef_form_old_up    = cp->cpcoeffs_pos_dvr[ip].icoef_form_up;
  int icoef_orth_up        = cp->cpcoeffs_pos_dvr[ip].icoef_orth_up;
  int ivcoef_form_up       = cp->cpcoeffs_pos_dvr[ip].ivcoef_form_up;
  int icoef_form_dn        = cp->cpcoeffs_pos_dvr[ip].icoef_form_dn;
  int icoef_form_old_dn    = cp->cpcoeffs_pos_dvr[ip].icoef_form_dn;
  int icoef_orth_dn        = cp->cpcoeffs_pos_dvr[ip].icoef_orth_dn;
  int ivcoef_form_dn       = cp->cpcoeffs_pos_dvr[ip].ivcoef_form_dn;
  int nstate_dn            = cp->cpcoeffs_info.nstate_dn;

  double *cpcoeffs_cmass   = cp->cpcoeffs_info.cmass;
  int icmass_unif          = cp->cpcoeffs_info.icmass_unif;

  double *occ_up           = cp->cpopts.occ_up;
  double *occ_dn           = cp->cpopts.occ_dn;
  double *rocc_sum_up      = cp->cpopts.rocc_sum_up;
  double *rocc_sum_dn      = cp->cpopts.rocc_sum_dn;

  double *cpscr_dvrc_up    = cp->cpscr.cpscr_wave.cre_up;
  double *cpscr_dvrc_dn    = cp->cpscr.cpscr_wave.cre_dn;


/*========================================================================*/
/* I) Norb off or norb with full ortho */

  if(icmass_unif !=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("   CP MASS has to be uniform for DVR basis\n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(cp_norb <= 1) {

    cp_shake_dvr(cpcoeffs_dvrc_up,icoef_form_up,icoef_orth_up,
                 cpcoeffs_dvrvc_up,ivcoef_form_up,
                 cpscr_dvrc_up,icoef_form_old_up,
                 c_tolshake,dt,&(cp->cpscr.cpscr_ovmat),
                 cpcoeffs_ioff_upt,occ_up,rocc_sum_up,cp_norb,
                 iter_shake,&(cp->cp_comm_state_pkg_dvr_up),dfact);

    if(cp_lsda == 1 && nstate_dn != 0) {
      cp_shake_dvr(cpcoeffs_dvrc_dn,icoef_form_dn,icoef_orth_dn,
                   cpcoeffs_dvrvc_dn,ivcoef_form_dn,
                   cpscr_dvrc_dn,icoef_form_old_dn,
                   c_tolshake,dt,&(cp->cpscr.cpscr_ovmat),
                   cpcoeffs_ioff_dnt,occ_dn,rocc_sum_dn,cp_norb,
                   iter_shake,&(cp->cp_comm_state_pkg_dvr_dn),dfact);
    }/* endif: lsda */

  }/* endif cp_norb */

/*========================================================================*/
/* II) Norb with norm only */

  if(cp_norb== 2){

    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("         DEBUG cp_shake_norb_dvr!!! \n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);

    /*
    cp_shake_norb_dvr(cpcoeffs_dvrc_up,icoef_form_up,icoef_orth_up,
                      cpcoeffs_dvrvc_up,ivcoef_form_up,
                      cpscr_dvrc_up,icoef_form_old_up,
                      dt,cpcoeffs_ioff_upt,occ_up,cpcoeffs_cmass,icmass_unif,
                      &(cp->cpscr.cpscr_ovmat),
                      &(cp->cp_comm_state_pkg_dvr_dn),dfact);

    if(cp_lsda == 1 && nstate_dn != 0) {

      cp_shake_norb_dvr(cpcoeffs_dvrc_dn,icoef_form_dn,icoef_orth_dn,
                        cpcoeffs_dvrvc_dn,ivcoef_form_dn,
                        cpscr_dvrc_dn,icoef_form_old_dn,
                        dt,cpcoeffs_ioff_dnt,occ_dn,cpcoeffs_cmass,icmass_unif,
                        &(cp->cpscr.cpscr_ovmat),
                        &(cp->cp_comm_state_pkg_dvr_dn),dfact);

    }
    */

  }/* endif cp_norb norm only*/

/*========================================================================*/
}/* end routine */
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rattle_control_cp_dvr(CP *cp,int *iter_rattle,double dt,int ip,double dfact)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  /*  Local pointers */

  int cp_norb              = cp->cpopts.cp_norb;
  int cp_lsda              = cp->cpopts.cp_lsda;
  double c_tolratl         = cp->cpconstrnt.c_tolratl;

  int    *cpcoeffs_ioff_upt = cp->cpcoeffs_info.ioff_upt;
  int    *cpcoeffs_ioff_dnt = cp->cpcoeffs_info.ioff_dnt;
  double *cpcoeffs_dvrc_up  = cp->cpcoeffs_pos_dvr[ip].dvrc_up;
  double *cpcoeffs_dvrc_dn  = cp->cpcoeffs_pos_dvr[ip].dvrc_dn;
  double *cpcoeffs_dvrvc_up = cp->cpcoeffs_pos_dvr[ip].dvrvc_up;
  double *cpcoeffs_dvrvc_dn = cp->cpcoeffs_pos_dvr[ip].dvrvc_dn;

  int icoef_form_up         = cp->cpcoeffs_pos_dvr[ip].icoef_form_up;
  int icoef_orth_up         = cp->cpcoeffs_pos_dvr[ip].icoef_orth_up;
  int ivcoef_form_up        = cp->cpcoeffs_pos_dvr[ip].ivcoef_form_up;
  int icoef_form_dn         = cp->cpcoeffs_pos_dvr[ip].icoef_form_dn;
  int icoef_orth_dn         = cp->cpcoeffs_pos_dvr[ip].icoef_orth_dn;
  int ivcoef_form_dn        = cp->cpcoeffs_pos_dvr[ip].ivcoef_form_dn;
  int nstate_dn             = cp->cpcoeffs_info.nstate_dn;


  double *cpcoeffs_cmass   = cp->cpcoeffs_info.cmass;
  int icmass_unif          = cp->cpcoeffs_info.icmass_unif;

  double *occ_up           = cp->cpopts.occ_up;
  double *occ_dn           = cp->cpopts.occ_dn;
  double *rocc_sum_up      = cp->cpopts.rocc_sum_up;
  double *rocc_sum_dn      = cp->cpopts.rocc_sum_dn;

  double *cpscr_dvrc_up     = cp->cpscr.cpscr_wave.cre_up;
  double *cpscr_dvrc_dn     = cp->cpscr.cpscr_wave.cre_dn;

/*========================================================================*/
/* I) Norb off or full ortho */

  if(cp_norb <= 1) {

    if(icmass_unif !=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("   CP MASS has to be uniform for DVR basis\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/

    cp_rattle_dvr(cpcoeffs_dvrc_up,icoef_form_up,icoef_orth_up,
                  cpcoeffs_dvrvc_up,ivcoef_form_up,&(cp->cpscr.cpscr_ovmat),
                  cpcoeffs_ioff_upt,rocc_sum_up,&(cp->cp_comm_state_pkg_dvr_up),
                  cp_norb,dfact);

    if(cp_lsda== 1 && nstate_dn != 0){

      cp_rattle_dvr(cpcoeffs_dvrc_dn,icoef_form_dn,icoef_orth_dn,
                    cpcoeffs_dvrvc_dn,ivcoef_form_dn,&(cp->cpscr.cpscr_ovmat),
                    cpcoeffs_ioff_dnt,rocc_sum_dn,&(cp->cp_comm_state_pkg_dvr_dn),
                    cp_norb,dfact);

    }/* endif lsda */

  }/* endif cp_norb */

/*========================================================================*/
/* I) Norb on */

  if(cp_norb== 2){

    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("         DEBUG cp_shake_norb_dvr!!! \n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);

    /*
    cp_rattle_norb_dvr(cpcoeffs_dvrc_up,icoef_form_up,icoef_orth_up,
                       cpcoeffs_dvrvc_up,ivcoef_form_up,cpscr_dvrc_up,
                       cpcoeffs_ioff_upt,occ_up,cpcoeffs_cmass,icmass_unif,
                       &(cp->cpscr.cpscr_ovmat),&(cp->cp_comm_state_pkg_dvr_up));

    if(cp_lsda== 1 && nstate_dn != 0) {

      cp_rattle_norb_dvr(cpcoeffs_dvrc_dn,icoef_form_dn,icoef_orth_dn,
                         cpcoeffs_dvrvc_dn,ivcoef_form_dn,cpscr_dvrc_dn,
                         cpcoeffs_ioff_dnt,occ_dn,cpcoeffs_cmass,icmass_unif,
                         &(cp->cpscr.cpscr_ovmat),&(cp->cp_comm_state_pkg_dvr_dn));

    } */
  }/*endif: cp_norb */

/*========================================================================*/
}/* end routine */
/*========================================================================*/

