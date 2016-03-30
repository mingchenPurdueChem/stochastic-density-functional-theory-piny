/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NVE                                      */
/*                                                                          */
/* This subprogram integrates the system using Vel Verlet                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_integrate_cp_entry.h"
#include "../proto_defs/proto_integrate_cp_local.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_coords_local.h"

#define DEBUG_OFF
#define DEBUG_NHC_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_0_to_dt2_cp_dvr(CLASS *class,BONDED *bonded,
                         GENERAL_DATA *general_data,CP *cp)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

     /* local pointers*/
    double dt         = general_data->timeinfo.dt;
    double dt2        = general_data->timeinfo.dt/2.0;

    int natm_tot = class->clatoms_info.natm_tot;

    int *cpcoeffs_ioff_up    = cp->cpcoeffs_info.ioff_upt;
    int *cpcoeffs_ioff_dn    = cp->cpcoeffs_info.ioff_dnt;

    int nstate_up        = cp->cpcoeffs_info.nstate_up;
    int nstate_dn        = cp->cpcoeffs_info.nstate_dn;
    int ncoef_tot        = cp->cpcoeffs_info.ncoef;

    int massiv_flag          = cp->cptherm_info.massiv_flag;
    int myid_state           = class->communicate.myid_state;
    int np_states            = class->communicate.np_states;
    MPI_Comm comm_states     = class->communicate.comm_states;

    int icmoff_up            = cp->cpcoeffs_info.icoef_start_up-1;
    int icmoff_dn            = cp->cpcoeffs_info.icoef_start_dn-1;

    int icoef_form_up        = cp->cpcoeffs_pos_dvr[1].icoef_form_up;
    int icoef_orth_up        = cp->cpcoeffs_pos_dvr[1].icoef_orth_up;
    int ivcoef_form_up       = cp->cpcoeffs_pos_dvr[1].ivcoef_form_up;
    int ivcoef_orth_up       = cp->cpcoeffs_pos_dvr[1].ivcoef_orth_up;
    int ifcoef_form_up       = cp->cpcoeffs_pos_dvr[1].ifcoef_form_up;
    int ifcoef_orth_up       = cp->cpcoeffs_pos_dvr[1].ifcoef_orth_up;
    int icoef_form_dn        = cp->cpcoeffs_pos_dvr[1].icoef_form_dn;
    int icoef_orth_dn        = cp->cpcoeffs_pos_dvr[1].icoef_orth_dn;
    int ivcoef_form_dn       = cp->cpcoeffs_pos_dvr[1].ivcoef_form_dn;
    int ivcoef_orth_dn       = cp->cpcoeffs_pos_dvr[1].ivcoef_orth_dn;
    int ifcoef_form_dn       = cp->cpcoeffs_pos_dvr[1].ifcoef_form_dn;
    int ifcoef_orth_dn       = cp->cpcoeffs_pos_dvr[1].ifcoef_orth_dn;

    int cp_norb              = cp->cpopts.cp_norb;
    int cp_lsda              = cp->cpopts.cp_lsda;
    int cp_isok_opt          = cp->cpopts.cp_isok_opt;
    int pi_beads_proc        = class->clatoms_info.pi_beads_proc;

    double *cpscr_cre_up  = cp->cpscr.cpscr_wave.cre_up;
    double *cpscr_cre_dn  = cp->cpscr.cpscr_wave.cre_dn;

    double *cpcoeffs_cmass    = cp->cpcoeffs_info.cmass;

    /* local variables */

    int i,ip,icoef,is,js;

    int ncoef,ncoef_up_tot,ncoef_dn_tot;
    int ncoef_up,ncoef_dn,ncoef_up_max,ncoef_dn_max;
    int nscale_up, nscale_dn;

    double *cpcoeffs_dvrc_up;
    double *cpcoeffs_dvrc_dn;

    double *cpcoeffs_dvrvc_up;
    double *cpcoeffs_dvrvc_dn;

    double *cpcoeffs_dvrfc_up;
    double *cpcoeffs_dvrfc_dn;

    /* Inverse of volume factor */
    double dfact = (double)(ncoef_tot)/(general_data->cell.vol_cp);

    /* Isokinetic option */    
    double temp_ext = cp->cpopts.te_ext;
    double k0;
    double ac_up = 0.0;
    double bc_up = 0.0;
    double ac_dn = 0.0;
    double bc_dn = 0.0;
    double kinet_cp = 0.0;

    nscale_up = ncoef_tot*nstate_up;
    nscale_dn = ncoef_tot*nstate_dn;

    k0 = temp_ext*((double)nscale_up)/BOLTZ;

/*==========================================================================*/
/* 0) Checks                                                                */

  if(cp_norb>0){
    if((icoef_orth_up+ivcoef_orth_up+ifcoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are in orthonormal form \n");
      printf("on state processor %d in int_NVE_cp \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
      if((icoef_orth_dn+ivcoef_orth_dn+ifcoef_orth_dn)!=0){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Dn CP vectors are in orthonormal form \n");
        printf("on state processor %d in int_NVE_cp \n",myid_state);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
    }/*endif*/
  }/*endif*/

  if(np_states>1){
    if((icoef_form_up+ivcoef_form_up+ifcoef_form_up)!=3){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP vectors are not in transposed form \n");
     printf("on state processor %d in int_NVE_cp \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cp_lsda==1){
      if((icoef_form_dn+ivcoef_form_dn+ifcoef_form_dn)!=3){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Up CP vectors are not in transposed form \n");
        printf("on state processor %d in int_NVE_cp \n",myid_state);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
    }/*endif*/
  }/*endif*/

/*==========================================================================*/
/* 0) Useful constants                                                      */

    general_data->timeinfo.int_res_tra    = 0;
    general_data->timeinfo.int_res_ter    = 0;

    zero_constrt_iters(&(general_data->stat_avg));

    if(np_states==1){
      ncoef_up     = cp->cp_para_fft_pkg3d_sm.nfft/2;
      ncoef_dn     = cp->cp_para_fft_pkg3d_sm.nfft/2;
      ncoef_up_max = cp->cp_para_fft_pkg3d_sm.nfft/2;
      ncoef_dn_max = cp->cp_para_fft_pkg3d_sm.nfft/2;
    }else{
      ncoef_up     = cp->cp_comm_state_pkg_dvr_up.nstate_ncoef_proc;
      ncoef_dn     = cp->cp_comm_state_pkg_dvr_dn.nstate_ncoef_proc;
      ncoef_up_max = cp->cp_comm_state_pkg_dvr_up.nstate_ncoef_proc_max;
      ncoef_dn_max = cp->cp_comm_state_pkg_dvr_dn.nstate_ncoef_proc_max;
    }/*endif*/

    ncoef_up_tot = ncoef_up_max*nstate_up;
    ncoef_dn_tot = ncoef_dn_max*nstate_dn;

    if((cp->cptherm_info.num_c_nhc > 0) && (cp->cptherm_info.len_c_nhc > 0)) {
      cp->cptherm_info.dt_nhc  = dt;
      cp->cptherm_info.dti_nhc = dt/( (double)(cp->cptherm_info.nres_c_nhc) );
      set_yosh(cp->cptherm_info.nyosh_c_nhc,cp->cptherm_info.dti_nhc ,
               cp->cptherm_info.wdti,cp->cptherm_info.wdti2,
               cp->cptherm_info.wdti4,cp->cptherm_info.wdti8,
               cp->cptherm_info.wdti16);
    }/*endif:therms on*/

/*==========================================================================*/
/* I) If gaussian dynamics apply the velocity dependent part              */
/*      of the Liouville operator                                           */

  if(cp->cpopts.cp_gauss == 1){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Gaussian dynamics has not been implemented for DVR \n");
     printf("basis sets \n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
  }

/*==========================================================================*/
/* II) First coefficient thermostat application                               */

  if((cp->cptherm_info.num_c_nhc > 0) && (cp->cptherm_info.len_c_nhc > 0)) {

    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Coefficient thermostat has not been implemented \n");
    printf("for DVR basis sets \n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);

    /*
    for(ip=1;ip<=pi_beads_proc;ip++) {
      if(massiv_flag==0){
        apply_c_nhc_dvr(&(cp->cptherm_info),&(cp->cptherm_pos[ip]),
                        &(cp->cpcoeffs_info),&(cp->cpcoeffs_pos_dvr[ip]),
                        &(cp->cpscr),cp->cpopts.cp_lsda,&(class->communicate));
      }else{
        apply_c_nhc_massiv_dvr(&(cp->cptherm_info),&(cp->cptherm_pos[ip]),
                               &(cp->cpcoeffs_info),&(cp->cpcoeffs_pos_dvr[ip]),
                               &(cp->cpscr),cp->cpopts.cp_lsda,
                               class->communicate.myid_state,
                               class->communicate.np_states);
      }
    }
    */
  }/*endif*/

/*==========================================================================*/
/* III) Propagate the coefficients                                        */

  for(ip=1;ip<=pi_beads_proc;ip++) {

    /**************************************/
    /* a) Save positions and coefficients */
    /**************************************/

    cpcoeffs_dvrc_up   = cp->cpcoeffs_pos_dvr[ip].dvrc_up;
    cpcoeffs_dvrvc_up  = cp->cpcoeffs_pos_dvr[ip].dvrvc_up;
    cpcoeffs_dvrfc_up  = cp->cpcoeffs_pos_dvr[ip].dvrfc_up;

    for(i=1;i<=ncoef_up_tot;i++){
      cpscr_cre_up[i] = cpcoeffs_dvrc_up[i];
    }/*endfor*/

    if( (cp->cpopts.cp_lsda==1) && (nstate_dn != 0) ){
      cpcoeffs_dvrc_dn   = cp->cpcoeffs_pos_dvr[ip].dvrc_dn;
      cpcoeffs_dvrvc_dn  = cp->cpcoeffs_pos_dvr[ip].dvrvc_dn;
      cpcoeffs_dvrfc_dn  = cp->cpcoeffs_pos_dvr[ip].dvrfc_dn;
      for(i=1;i<=ncoef_dn_tot;i++){
        cpscr_cre_dn[i] = cpcoeffs_dvrc_dn[i];
      }/*endfor*/
    }/*endif*/

    /******************************/
    /* b) CP isokinetic goes here */
    /******************************/

    if(cp_isok_opt == 1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("isokinetic method has not been implemented for \n");
      printf("DVR basis sets \n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/* endif isok opt */
 
    /*****************************/
    /* c) Evolve Coef velocities */
    /*****************************/

    /* Up states, excluding isokinetic option */

    for(is=1;is<=nstate_up;is++) {
      for(i=1;i<=ncoef_up;i++) {
        icoef = i + cpcoeffs_ioff_up[is];
        cpcoeffs_dvrvc_up[icoef] +=
            dt2*cpcoeffs_dvrfc_up[icoef]/cpcoeffs_cmass[(i+icmoff_up)];
      }/*endfor*/
    }/*endfor*/

    /* Down states, excluding isokinetic option */

    if( (cp->cpopts.cp_lsda == 1) && (nstate_dn != 0) ){
      for(is=1;is<=nstate_dn;is++) {
        for(i=1;i<=ncoef_dn;i++) {
          icoef = i+cpcoeffs_ioff_dn[is];
          cpcoeffs_dvrvc_dn[icoef] +=
                dt2*cpcoeffs_dvrfc_dn[icoef]/cpcoeffs_cmass[(i+icmoff_dn)];
        }/*endfor*/
      }/*endfor*/
    }/* endif lsda */

    /**************************/
    /* d) Evolve coefficients */
    /**************************/

    for(i=1;i<=ncoef_up_tot;i++){
      cpcoeffs_dvrc_up[i] += dt*cpcoeffs_dvrvc_up[i];
    }/*endfor*/

    if( (cp->cpopts.cp_lsda == 1) && (nstate_dn != 0) ){
      for(i=1;i<=ncoef_dn_tot;i++){
        cpcoeffs_dvrc_dn[i] += dt*cpcoeffs_dvrvc_dn[i];
      }/*endfor*/
    }/* endif */

   /*****************/
   /* e) CP-Shake   */ 
   /*****************/

    if( (cp->cpopts.cp_gauss == 0) && (cp->cpopts.cp_norb < 3) ){
       shake_control_cp_dvr(cp,&(general_data->stat_avg.iter_shake_cp),
                            general_data->timeinfo.dt,ip,dfact);
    }/*endif*/

/*-------------------------------------------------------------------------*/
  }/*endfor : ip*/
/*-------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

