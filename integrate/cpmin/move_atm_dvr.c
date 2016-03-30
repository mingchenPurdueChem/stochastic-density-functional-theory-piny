/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: move_atm                                     */
/*                                                                          */
/* This subprogram minimizes atomic positions using either                  */
/* steepest descent, conjugate gradient or DIIS                             */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_integrate_cpmin_entry.h"
#include "../proto_defs/proto_integrate_cpmin_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void move_atm_cg_dvr(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                 CP *cp,int ifirst)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    /* local variables */
 
    int i,ip=1,icoef,is,natm_tot;
    int nstate_up,nstate_dn;
    int ncoef_up,ncoef_dn;
    static double fovlap;
    double fovlap_old,gamma;
    double dt,dts;
    double fc_mag_up, fc_mag_dn;
    double *zeta_up,*zeta_dn;

    /* local pointers */
 
    double mass_sc_fact  = class->clatoms_info.mass_sc_fact;
    double *clatoms_x    = class->clatoms_pos[1].x;
    double *clatoms_y    = class->clatoms_pos[1].y;
    double *clatoms_z    = class->clatoms_pos[1].z;
    double *clatoms_fx   = class->clatoms_pos[1].fx;
    double *clatoms_fy   = class->clatoms_pos[1].fy;
    double *clatoms_fz   = class->clatoms_pos[1].fz;
    double *grad_x       = class->ewd_scr.fx;
    double *grad_y       = class->ewd_scr.fy;
    double *grad_z       = class->ewd_scr.fz;
    double *hess_xx      = class->clatoms_pos[1].hess_xx;
    double *hess_xy      = class->clatoms_pos[1].hess_xy;
    double *hess_xz      = class->clatoms_pos[1].hess_xz;
    double *hess_yy      = class->clatoms_pos[1].hess_yy;
    double *hess_yz      = class->clatoms_pos[1].hess_yz;
    double *hess_zz      = class->clatoms_pos[1].hess_zz;
    double *chx          = class->clatoms_pos[1].vx;
    double *chy          = class->clatoms_pos[1].vy;
    double *chz          = class->clatoms_pos[1].vz;
    double *scr_x        = class->clatoms_info.xold;
    double *scr_y        = class->clatoms_info.yold;
    double *scr_z        = class->clatoms_info.zold;
    double *clatoms_mass = class->clatoms_info.mass;

    double *dvrc_up       = cp->cpcoeffs_pos_dvr[1].dvrc_up;
    double *dvrc_dn       = cp->cpcoeffs_pos_dvr[1].dvrc_dn;

    double *dvrfc_up      = cp->cpcoeffs_pos_dvr[1].dvrfc_up;
    double *dvrfc_dn      = cp->cpcoeffs_pos_dvr[1].dvrfc_up;

    double *cmass        = cp->cpcoeffs_info.cmass;
    int *ioff_up         = cp->cpcoeffs_info.ioff_upt;
    int *ioff_dn         = cp->cpcoeffs_info.ioff_dnt;

    int np_states        = class->communicate.np_states;
    int myid             = class->communicate.myid;
    int myid_state       = class->communicate.myid_state;
    MPI_Comm comm_states = class->communicate.comm_states;

    int hess_calc        = class->clatoms_info.hess_calc;
    int cp_norb          = cp->cpopts.cp_norb;
    int cp_lsda          = cp->cpopts.cp_lsda;

    int icmoff_up        = cp->cpcoeffs_info.icoef_start_up-1;
    int icmoff_dn        = cp->cpcoeffs_info.icoef_start_dn-1;

    int icoef_form_up    = cp->cpcoeffs_pos_dvr[1].icoef_form_up;
    int icoef_orth_up    = cp->cpcoeffs_pos_dvr[1].icoef_orth_up;
    int ifcoef_form_up   = cp->cpcoeffs_pos_dvr[1].ifcoef_form_up;
    int ifcoef_orth_up   = cp->cpcoeffs_pos_dvr[1].ifcoef_orth_up;

    int icoef_form_dn    = cp->cpcoeffs_pos_dvr[1].icoef_form_dn;
    int icoef_orth_dn    = cp->cpcoeffs_pos_dvr[1].icoef_orth_dn;
    int ifcoef_form_dn   = cp->cpcoeffs_pos_dvr[1].ifcoef_form_dn;
    int ifcoef_orth_dn   = cp->cpcoeffs_pos_dvr[1].ifcoef_orth_dn;

    int min_atm_com_fix_opt = general_data->minopts.min_atm_com_fix_opt;


/*==========================================================================*/
/* 0) Checks                                                                */

  if(cp_norb>0){
    if((icoef_orth_up+ifcoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are in orthonormal form \n");
      printf("on state processor %d in move_atm_cg \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_orth_dn+ifcoef_orth_dn)!=0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dn CP vectors are in orthonormal form \n");
       printf("on state processor %d in move_atm_cg \n",myid_state);
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/

  if(np_states>1){
    if((icoef_form_up+ifcoef_form_up)!=2){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP vectors are not in transposed form \n");
     printf("on state processor %d in move_atm_cg \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_form_dn+ifcoef_form_dn)!=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("on state processor %d in move_atm_cg \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/

/*==========================================================================*/
/* 0) Useful constants                                                      */


    natm_tot = class->clatoms_info.natm_tot;
    general_data->timeinfo.int_res_tra    = 0;
    general_data->timeinfo.int_res_ter    = 0;
    dt  = general_data->timeinfo.dt;
    dts = dt*mass_sc_fact;

    zero_constrt_iters(&(general_data->stat_avg));

    nstate_up    = cp->cpcoeffs_info.nstate_up;
    nstate_dn    = cp->cpcoeffs_info.nstate_dn;
    if(np_states==1){
     ncoef_up      = cp->cp_para_fft_pkg3d_lg.nfft/2;
     ncoef_dn      = cp->cp_para_fft_pkg3d_lg.nfft/2;
    }else{
     ncoef_up     =  cp->cp_para_fft_pkg3d_lg.nfft_proc/2;
     ncoef_dn     =  cp->cp_para_fft_pkg3d_lg.nfft_proc/2;
    }/*endif*/

/*==========================================================================*/
/* 0.1) Zero conjugate gradients                                            */

   if(ifirst == 1){
     for(i=1;i<=natm_tot; i++){
      chx[i] = 0.0;
      chy[i] = 0.0;
      chz[i] = 0.0;
     }
     gamma = 0.0;
     fovlap = 1.0;
   }/* endif */

/*==========================================================================*/
/* I) Get forces                                                            */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;

    cp_energy_control(class,bonded,general_data,cp);

/*==========================================================================*/
/* II) Calculate the gamma's                                                */

  if(myid_state==0){
    if(min_atm_com_fix_opt==1){
      proj_com_out(natm_tot,clatoms_fx,clatoms_fy,clatoms_fz);
    }/*endif*/

    fovlap_old = fovlap;
    fovlap = 0.0;
    for(i=1;i<=natm_tot;i++){
      fovlap += clatoms_fx[i]*clatoms_fx[i]
             + clatoms_fy[i]*clatoms_fy[i]
             + clatoms_fz[i]*clatoms_fz[i];
    }
    if(ifirst != 1) {gamma = fovlap/fovlap_old;}
/*==========================================================================*/
/* II.V) Evolve gradients                                                   */

    for(i=1;i<=natm_tot; i++){
      chx[i] = clatoms_fx[i] + gamma*chx[i];
      chy[i] = clatoms_fy[i] + gamma*chy[i];
      chz[i] = clatoms_fz[i] + gamma*chz[i];
    }/* endfor */

/*==========================================================================*/
/* III) Calculate the step length                                           */

    switch(hess_calc){
      case 1:
      act_hess_inv_on_grad(hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
                                grad_x,grad_y,grad_z,chx,chy,chz,
                                natm_tot);
      break;

      case 2:
      act_hess_inv_on_grad_diag(hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
                           grad_x,grad_y,grad_z,chx,chy,chz,
                           natm_tot);
      break;

      default:

      for(i=1;i<=natm_tot;i++){
        grad_x[i] = chx[i];
        grad_y[i] = chy[i];
        grad_z[i] = chz[i];
      }/* endfor */
      break;

    }/* end switch */


/*==========================================================================*/
/* IV) Evolve positions                                                     */

    for(i=1;i<=natm_tot;i++){
      clatoms_x[i] += dts*grad_x[i]/clatoms_mass[i];
      clatoms_y[i] += dts*grad_y[i]/clatoms_mass[i];
      clatoms_z[i] += dts*grad_z[i]/clatoms_mass[i];
    }/*endfor*/

    get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[1]),
                  &(class->ghost_atoms));

  }/* endif myid_state */

  if(np_states > 1){
    Bcast(&(clatoms_x[1]),natm_tot,MPI_DOUBLE,0,comm_states);
    Bcast(&(clatoms_y[1]),natm_tot,MPI_DOUBLE,0,comm_states);
    Bcast(&(clatoms_z[1]),natm_tot,MPI_DOUBLE,0,comm_states);
  }/* endif */

/*==========================================================================*/
/* IV.V) Do a steepest descent step for coefficients                        */

  for(is=1;is<=nstate_up;is++) {
    for(i=1;i<=ncoef_up;i++) {
      icoef = i+ioff_up[is];
      dvrc_up[icoef] += dt*dvrfc_up[icoef]/cmass[(i+icmoff_up)];
    }/*endfor*/
  }/*endfor*/

  if( (cp->cpopts.cp_lsda == 1) && (nstate_dn != 0) ){
    for(is=1;is<=nstate_dn;is++) {
      for(i=1;i<=ncoef_dn;i++) {
        icoef = i+ioff_dn[is];
        dvrc_dn[icoef] += dt*dvrfc_dn[icoef]/cmass[(i+icmoff_dn)];
      }/*endfor*/
    }/*endfor*/
  }/* endif */

/*==========================================================================*/
/* V) Orthogonalize wave functions                                          */

  orthog_control_cp_dvr(cp,ip);

/*==========================================================================*/
/* ii) Free memory and shuffle states                                       */

  cp_shuffle_states_dvr(cp,ip);

/*========================================================================*/
}/*end routine*/
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void move_atm_std_dvr(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                      CP *cp,int ifirst)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    /* local variables */

    int i,ip=1,icoef,is,natm_tot;
    double dt,dts;
    int nstate_up,nstate_dn,ncoef_tot;
    int ncoef_up,ncoef_dn;
    double fc_mag_up, fc_mag_dn;
    double Tmass,xcom,ycom,zcom;
    double *zeta_up,*zeta_dn;

    /* local pointers */

    int hess_calc        = class->clatoms_info.hess_calc;
    double mass_sc_fact  = class->clatoms_info.mass_sc_fact;
    double *clatoms_x    = class->clatoms_pos[1].x;
    double *clatoms_y    = class->clatoms_pos[1].y;
    double *clatoms_z    = class->clatoms_pos[1].z;
    double *clatoms_fx   = class->clatoms_pos[1].fx;
    double *clatoms_fy   = class->clatoms_pos[1].fy;
    double *clatoms_fz   = class->clatoms_pos[1].fz;
    double *grad_x       = class->ewd_scr.fx;
    double *grad_y       = class->ewd_scr.fy;
    double *grad_z       = class->ewd_scr.fz;
    double *hess_xx      = class->clatoms_pos[1].hess_xx;
    double *hess_xy      = class->clatoms_pos[1].hess_xy;
    double *hess_xz      = class->clatoms_pos[1].hess_xz;
    double *hess_yy      = class->clatoms_pos[1].hess_yy;
    double *hess_yz      = class->clatoms_pos[1].hess_yz;
    double *hess_zz      = class->clatoms_pos[1].hess_zz;
    double *chx          = class->clatoms_pos[1].vx;
    double *chy          = class->clatoms_pos[1].vy;
    double *chz          = class->clatoms_pos[1].vz;
    double *scr_x        = class->clatoms_info.xold;
    double *scr_y        = class->clatoms_info.yold;
    double *scr_z        = class->clatoms_info.zold;
    double *clatoms_mass = class->clatoms_info.mass;

    double *dvrc_up      = cp->cpcoeffs_pos_dvr[1].dvrc_up;
    double *dvrfc_up     = cp->cpcoeffs_pos_dvr[1].dvrfc_up;

    double *dvrc_dn      = cp->cpcoeffs_pos_dvr[1].dvrc_dn;
    double *dvrfc_dn     = cp->cpcoeffs_pos_dvr[1].dvrfc_dn;

    double *cmass        = cp->cpcoeffs_info.cmass;
    int *ioff_up         = cp->cpcoeffs_info.ioff_upt;
    int *ioff_dn         = cp->cpcoeffs_info.ioff_dnt;

    int np_states        = class->communicate.np_states;
    int myid_state       = class->communicate.myid_state;
    MPI_Comm comm_states = class->communicate.comm_states;

    int cp_norb          = cp->cpopts.cp_norb;
    int cp_lsda          = cp->cpopts.cp_lsda;

    int icmoff_up        = cp->cpcoeffs_info.icoef_start_up-1;
    int icmoff_dn        = cp->cpcoeffs_info.icoef_start_dn-1;

    int icoef_form_up    = cp->cpcoeffs_pos_dvr[1].icoef_form_up;
    int icoef_orth_up    = cp->cpcoeffs_pos_dvr[1].icoef_orth_up;
    int ifcoef_form_up   = cp->cpcoeffs_pos_dvr[1].ifcoef_form_up;
    int ifcoef_orth_up   = cp->cpcoeffs_pos_dvr[1].ifcoef_orth_up;

    int icoef_form_dn    = cp->cpcoeffs_pos_dvr[1].icoef_form_dn;
    int icoef_orth_dn    = cp->cpcoeffs_pos_dvr[1].icoef_orth_dn;
    int ifcoef_form_dn   = cp->cpcoeffs_pos_dvr[1].ifcoef_form_dn;
    int ifcoef_orth_dn   = cp->cpcoeffs_pos_dvr[1].ifcoef_orth_dn;

    int min_atm_com_fix_opt = general_data->minopts.min_atm_com_fix_opt;

/*==========================================================================*/
/* 0) Checks                                                                */

  if(cp_norb>0){
    if((icoef_orth_up+ifcoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are in orthonormal form \n");
      printf("on state processor %d in move_atm_std \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_orth_dn+ifcoef_orth_dn)!=0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dn CP vectors are in orthonormal form \n");
       printf("on state processor %d in move_atm_std \n",myid_state);
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/

  if(np_states>1){
    if((icoef_form_up+ifcoef_form_up)!=2){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP vectors are not in transposed form \n");
     printf("on state processor %d in move_atm_std \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_form_dn+ifcoef_form_dn)!=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("on state processor %d in move_atm_std \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/

/*==========================================================================*/
/* 0) Useful constants                                                      */

    natm_tot = class->clatoms_info.natm_tot;
    general_data->timeinfo.int_res_tra    = 0;
    general_data->timeinfo.int_res_ter    = 0;
    dt  = general_data->timeinfo.dt;
    dts = dt*mass_sc_fact;
    zero_constrt_iters(&(general_data->stat_avg));

    nstate_up    = cp->cpcoeffs_info.nstate_up;
    nstate_dn    = cp->cpcoeffs_info.nstate_dn;
    if(np_states==1){
     ncoef_up      = cp->cp_para_fft_pkg3d_lg.nfft/2;
     ncoef_dn      = cp->cp_para_fft_pkg3d_lg.nfft/2;
    }else{
     ncoef_up     =  cp->cp_para_fft_pkg3d_lg.nfft_proc/2;
     ncoef_dn     =  cp->cp_para_fft_pkg3d_lg.nfft_proc/2;
    }/*endif*/

/*==========================================================================*/
/* I) Get forces                                                            */

    if(myid_state==0){
      Tmass=0.0;
      for(i=1;i<=natm_tot;i++){
        Tmass += clatoms_mass[i];
      }
      xcom=0.0;
      ycom=0.0;
      zcom=0.0;
      for(i=1;i<=natm_tot;i++){
        xcom += clatoms_x[i]*clatoms_mass[i];
        ycom += clatoms_y[i]*clatoms_mass[i];
        zcom += clatoms_z[i]*clatoms_mass[i];
      }/* endfor */
      xcom /= Tmass;
      ycom /= Tmass;
      zcom /= Tmass;
    }

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;

    cp_energy_control(class,bonded,general_data,cp);

/*==========================================================================*/
/* II) Get forces                                                            */


  if(myid_state == 0){

    if(min_atm_com_fix_opt==1){
      proj_com_out(natm_tot,clatoms_fx,clatoms_fy,clatoms_fz);
    }/*endif*/

    switch(hess_calc){
      case 1:
      act_hess_inv_on_grad(hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
                                grad_x,grad_y,grad_z,clatoms_fx,clatoms_fy,clatoms_fz,
                                natm_tot);
      break;

      case 2:
      act_hess_inv_on_grad_diag(hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
                           grad_x,grad_y,grad_z,clatoms_fx,clatoms_fy,clatoms_fz,
                           natm_tot);
      break;

      default:

      for(i=1;i<=natm_tot;i++){
        grad_x[i] = clatoms_fx[i];
        grad_y[i] = clatoms_fy[i];
        grad_z[i] = clatoms_fz[i];
      }/* endfor */

      break;
    }/* end switch */

/*==========================================================================*/
/* III) Evolve positions                                                     */

    for(i=1;i<=natm_tot;i++){
      clatoms_x[i] += dts*grad_x[i]/clatoms_mass[i];
      clatoms_y[i] += dts*grad_y[i]/clatoms_mass[i];
      clatoms_z[i] += dts*grad_z[i]/clatoms_mass[i];
    }/*endfor*/

    get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[1]),
                  &(class->ghost_atoms));

  }/* endif myid_state */

  if(np_states > 1){
    Bcast(&(clatoms_x[1]),natm_tot,MPI_DOUBLE,0,comm_states);
    Bcast(&(clatoms_y[1]),natm_tot,MPI_DOUBLE,0,comm_states);
    Bcast(&(clatoms_z[1]),natm_tot,MPI_DOUBLE,0,comm_states);
  }/* endif */

/*==========================================================================*/
/* IV.V) Do a steepest descent step for coefficients                        */

  for(is=1;is<=nstate_up;is++) {
    for(i=1;i<=ncoef_up;i++) {
      icoef = i+ioff_up[is];
      dvrc_up[icoef] += dt*dvrfc_up[icoef]/cmass[(i+icmoff_up)];
    }/*endfor*/
  }/*endfor*/

  if( (cp->cpopts.cp_lsda == 1) && (nstate_dn != 0) ){
    for(is=1;is<=nstate_dn;is++) {
      for(i=1;i<=ncoef_dn;i++) {
        icoef = i+ioff_dn[is];
        dvrc_dn[icoef] += dt*dvrfc_dn[icoef]/cmass[(i+icmoff_dn)];
      }/*endfor*/
    }/*endfor*/
  }/* endif */

/*==========================================================================*/
/* V) Orthogonalize wave functions                                          */

  orthog_control_cp_dvr(cp,ip);

/*==========================================================================*/
/* ii) Free memory and shuffle states                                       */

  cp_shuffle_states_dvr(cp,ip);

/*========================================================================*/
}/*end routine*/
/*========================================================================*/

