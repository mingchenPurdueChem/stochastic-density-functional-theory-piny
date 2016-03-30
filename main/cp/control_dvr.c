/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_cp                                   */
/*                                                                          */
/* This subprogram performs Minization on a classical+ab-initio potential   */
/* energy surface (GGLSDA/GGLDA-PES)                                        */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_main_local.h"
#include "../proto_defs/proto_main_cp_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_integrate_cpmin_entry.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_output_entry.h"
#include "../proto_defs/proto_output_cp_entry.h"
#include "../proto_defs/proto_output_cp_local.h"
#include "../proto_defs/proto_analysis_md_entry.h"
#include "../proto_defs/proto_analysis_cp_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void check_coef_grad_mag_dvr(CP *cp,SIMOPTS *simopts,
                             double *fc_mag_up_ret,double *fc_mag_dn_ret,
                             int *ireset_ret,int *idone_ret, double tol_coef,
                             int ip_start,int ip_end,STAT_AVG *stat_avg)
/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

 int iii;
 double fc_mag_up_old  = *fc_mag_up_ret;
 double fc_mag_dn_old  = *fc_mag_dn_ret;
 double fc_mag_up,fc_mag_dn;
 double fc_max_up,fc_max_dn;
 double fc_mag_up_tmp,fc_mag_dn_tmp;
 double fc_max_up_tmp,fc_max_dn_tmp;
 int pi_beads       = cp->cpcoeffs_info.pi_beads;
 int pi_beads_proc  = cp->cpcoeffs_info.pi_beads_proc;
 int ip;
 int i,idone,ireset;
 int ncoef_tot;
 double *dvrfc_up,*dvrfc_dn;
 double *dvrc_up, *dvrc_dn;
 double *ksmat_up,*ksmat_dn;

 double *occ_up     = cp->cpopts.occ_up;
 double *occ_dn     = cp->cpopts.occ_dn;
 double *ksmat_scr  = cp->cpscr.cpscr_ovmat.ovlap1;
 int *ioff_upt      = cp->cpcoeffs_info.ioff_upt;
 int *ioff_dnt      = cp->cpcoeffs_info.ioff_dnt;
 int cp_norb        = cp->cpopts.cp_norb;
 int cp_min;
 int ncoef_up       = cp->cp_comm_state_pkg_dvr_up.nstate_ncoef_proc;
 int ncoef_dn       = cp->cp_comm_state_pkg_dvr_dn.nstate_ncoef_proc;
 int nstate_up      = cp->cpcoeffs_info.nstate_up;
 int nstate_dn      = cp->cpcoeffs_info.nstate_dn;
 int cp_lsda        = cp->cpopts.cp_lsda;
 int myid_state     = cp->communicate.myid_state;
 int myid_bead      = cp->communicate.myid_bead;
 int np_states      = cp->communicate.np_states;
 int np_beads       = cp->communicate.np_beads;

 int icoef_orth_up,icoef_orth_dn;
 int icoef_form_up,icoef_form_dn;
 int ifcoef_orth_up,ifcoef_orth_dn;
 int ifcoef_form_up,ifcoef_form_dn;

 MPI_Comm comm_states = cp->communicate.comm_states;
 MPI_Comm comm_beads  = cp->communicate.comm_beads;

/*=======================================================================*/
/* 0) Parallel checks  */

  if(np_states>1){
    for(ip=ip_start;ip<=ip_end;ip++){
      ifcoef_form_up = cp->cpcoeffs_pos_dvr[ip].ifcoef_form_up;
      ifcoef_form_dn = cp->cpcoeffs_pos_dvr[ip].ifcoef_form_dn;
      if(ifcoef_form_up!=1){
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("Up Coef forces are not in transposed form \n");
         printf("on state processor %d in check_coef_grad \n",myid_state);
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
         exit(1);
       }/*endif*/
       if(cp_lsda==1){
         if(ifcoef_form_dn!=1){
           printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
           printf("Dn Coef forces are not in transposed form \n");
           printf("on state processor %d in check_coef_grad \n",myid_state);
           printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
           fflush(stdout);
           exit(1);
         }/*endif*/
       }/*endif*/
     }/* endfor */
   }/*endif*/

/*=======================================================================*/
/* 0.1) Initialize and calculate some useful quantities                  */

  ncoef_tot = ncoef_up*nstate_up;

  cp_min = simopts->cp_min + simopts->cp_wave_min + simopts->cp_wave_min_pimd;

  fc_mag_up = 0.0;
  fc_max_up = 0.0;
  fc_mag_dn = 0.0;
  fc_max_dn = 0.0;

/*=======================================================================*/
/* Begin for loop over beads */

  for(ip=ip_start;ip<=ip_end;ip++){
/*-----------------------------------------------------------------------*/
/* Assign and project forces if necessary */

    dvrc_up        = cp->cpcoeffs_pos_dvr[ip].dvrc_up;
    dvrc_dn        = cp->cpcoeffs_pos_dvr[ip].dvrc_dn;
    icoef_orth_up  = cp->cpcoeffs_pos_dvr[ip].icoef_orth_up;
    icoef_orth_dn  = cp->cpcoeffs_pos_dvr[ip].icoef_orth_dn;
    icoef_form_up  = cp->cpcoeffs_pos_dvr[ip].icoef_form_up;
    icoef_form_dn  = cp->cpcoeffs_pos_dvr[ip].icoef_form_dn;
    ifcoef_orth_up = cp->cpcoeffs_pos_dvr[ip].ifcoef_orth_up;
    ifcoef_orth_dn = cp->cpcoeffs_pos_dvr[ip].ifcoef_orth_dn;
    ifcoef_form_up = cp->cpcoeffs_pos_dvr[ip].ifcoef_form_up;
    ifcoef_form_dn = cp->cpcoeffs_pos_dvr[ip].ifcoef_form_dn;
    ksmat_up       = cp->cpcoeffs_pos_dvr[ip].ksmat_up;
    ksmat_dn       = cp->cpcoeffs_pos_dvr[ip].ksmat_dn;

    if(cp_min == 1 || cp_norb >= 1) {
      dvrfc_up = cp->cpcoeffs_pos_dvr[ip].dvrfc_up;
      dvrfc_dn = cp->cpcoeffs_pos_dvr[ip].dvrfc_dn;
    }else { /* CHECK PARALLEL */
      dvrfc_up = cp->cpscr.cpscr_wave.cre_up;
      dvrfc_dn = cp->cpscr.cpscr_wave.cre_dn;
      for(i=1;i<=ncoef_tot; i++){
        dvrfc_up[i] = cp->cpcoeffs_pos_dvr[ip].dvrfc_up[i];
      }
      if(cp->cpopts.cp_lsda == 1 && nstate_dn != 0){
        for(i=1;i<=ncoef_tot; i++){
          dvrfc_dn[i] = cp->cpcoeffs_pos_dvr[ip].dvrfc_dn[i];
        }
      }/* endif */

      cp_add_ksmat_force_dvr(dvrc_up,icoef_form_up,icoef_orth_up,
                            dvrfc_up,ifcoef_form_up,ifcoef_orth_up,
                            ksmat_up,ksmat_scr,ioff_upt,cp_lsda,cp_min,occ_up,
                            cp->cpcoeffs_info.scale_fact,
                            &(cp->cp_comm_state_pkg_dvr_up));

      if( (cp_lsda==1) && (nstate_dn!=0) ){
        cp_add_ksmat_force_dvr(dvrc_dn,icoef_form_dn,icoef_orth_dn,
                               dvrfc_dn,ifcoef_form_dn,ifcoef_orth_dn,
                               ksmat_dn,ksmat_scr,ioff_dnt,cp_lsda,cp_min,occ_dn,
                               cp->cpcoeffs_info.scale_fact,
                               &(cp->cp_comm_state_pkg_dvr_dn));
      }/*endif*/
    }/* endif */

/*-----------------------------------------------------------------------*/
/* I) Up tolerence */

    for(i=1;i <= ncoef_tot; i++) {
      fc_mag_up += (dvrfc_up[i]*dvrfc_up[i]);
      fc_max_up = MAX(fabs(dvrfc_up[i]),fc_max_up);
    }/*endfor*/

/*-----------------------------------------------------------------------*/
/* II) Dn tolerence */

    if(cp_lsda == 1){
      ncoef_tot = ncoef_dn*nstate_dn;
      for(i=1;i <= ncoef_tot; i++) {
        fc_mag_dn += dvrfc_dn[i] * dvrfc_dn[i];
        fc_max_dn = MAX(fabs(dvrfc_dn[i]),fc_max_dn);
      }/*endfor*/
    }/*endif*/

  }/* endfor ip */

/*=======================================================================*/
/* II.V) Parallel reductions */

/*------------------------------------------------------------------------*/
/* i) First state level  */

  if(np_states > 1){
    fc_mag_up_tmp = fc_mag_up;
    fc_mag_dn_tmp = fc_mag_dn;
    Allreduce(&(fc_mag_up_tmp),&(fc_mag_up),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
    Allreduce(&(fc_mag_dn_tmp),&(fc_mag_dn),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
    fc_max_up_tmp = fc_max_up;
    fc_max_dn_tmp = fc_max_dn;
    Allreduce(&(fc_max_up_tmp),&(fc_max_up),1,MPI_DOUBLE,MPI_MAX,0,comm_states);
    Allreduce(&(fc_max_dn_tmp),&(fc_max_dn),1,MPI_DOUBLE,MPI_MAX,0,comm_states);
  }/*endif*/

/*------------------------------------------------------------------------*/
/* ii) Next, bead level  */

  if(np_beads > 1 && cp_min == 0){
    fc_mag_up_tmp = fc_mag_up;
    fc_mag_dn_tmp = fc_mag_dn;
    Allreduce(&(fc_mag_up_tmp),&(fc_mag_up),1,MPI_DOUBLE,MPI_SUM,0,comm_beads);
    Allreduce(&(fc_mag_dn_tmp),&(fc_mag_dn),1,MPI_DOUBLE,MPI_SUM,0,comm_beads);
    fc_max_up_tmp = fc_max_up;
    fc_max_dn_tmp = fc_max_dn;
    Allreduce(&(fc_max_up_tmp),&(fc_max_up),1,MPI_DOUBLE,MPI_MAX,0,comm_beads);
    Allreduce(&(fc_max_dn_tmp),&(fc_max_dn),1,MPI_DOUBLE,MPI_MAX,0,comm_beads);
  }/*endif*/

/*=======================================================================*/
/* II.VI) Calculate the magnitude */

  fc_mag_up = sqrt(fc_mag_up/((double)pi_beads));
  fc_mag_dn = sqrt(fc_mag_dn/((double)pi_beads));

/*=======================================================================*/
/* III) Set the flags (used for minimization)                            */

  if(myid_state==0){
    idone = 1;
    if(fc_mag_up  > tol_coef){idone=0;}
    if((fc_mag_dn > tol_coef)&&(cp_lsda==1)){idone=0;}
    ireset = 0;
    if(fc_mag_up  > fc_mag_up_old){ireset=1;}
    if((fc_mag_dn > fc_mag_dn_old)&&(cp_lsda==1)){ireset=1;}
  }/*endif*/
  if(np_states>1){
    Bcast(&idone,1,MPI_INT,0,comm_states);
    Bcast(&ireset,1,MPI_INT,0,comm_states);
  }/*endif*/


/*=======================================================================*/
/* IV) Set return values and put things in structures                    */

  *fc_mag_up_ret = fc_mag_up;
  *fc_mag_dn_ret = fc_mag_dn;
  *idone_ret     = idone;
  *ireset_ret    = ireset;
  (stat_avg->fc_mag_up) = fc_mag_up;
  (stat_avg->fc_mag_dn) = fc_mag_dn;
  (stat_avg->fc_max_up) = fc_max_up;
  (stat_avg->fc_max_dn) = fc_max_dn;

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


