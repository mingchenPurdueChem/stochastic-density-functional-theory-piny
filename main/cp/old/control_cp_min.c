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
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_cp_min(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                                 CP *cp,ANALYSIS *analysis)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  double fc_mag_up,fc_mag_dn,f_atm_mag,elec_e_old,elec_e;
  double elec_e_old_tmp,elec_e_tmp;
  int ncoef_tot,i,itime,iii,ip_now=1;         
  int pi_beads = class->clatoms_info.pi_beads; 
  int atm_step,ifirst_atm,atm_min;
  int iexit = 0;
  int ireset,idone;
  int idone_atm=0;
  int iatm_count;
  int idum=0;
  int iwrite_confp=general_data->filenames.iwrite_confp;
  int iwrite_confc=general_data->filenames.iwrite_confc;
  int myid = class->communicate.myid;
  int np_state = cp->communicate.np_states;
  int num_proc = cp->communicate.np;
  MPI_Comm comm_states = class->communicate.comm_states;
  

/*======================================================================*/
/* 0) Write to Screen                                                   */

 if(myid==0){
  PRINT_LINE_STAR;
  printf("Running CP-MINIMIZATION \n");
  PRINT_LINE_DASH;
 }/* endif */

/*======================================================================*/
/* I) Set the flags/counters                                           */


  general_data->stat_avg.write_cp_atm_flag = 0;

  ifirst_atm=1;
  iatm_count=0;
  atm_min=0;
  if(general_data->simopts.cp_min == 1){atm_min=1;}
  atm_step=0;
  class->clatoms_info.cg_reset_flag = 1;
  general_data->stat_avg.updates           = 0.0;
  f_atm_mag = 10000.0;

  general_data->simopts.cp_min           = 0;
  general_data->simopts.cp_wave_min      = 1;
  cp->cpcoeffs_info.cg_reset_flag   = 1;
  cp->cpcoeffs_info.diis_reset_flag = 1;
  fc_mag_up = 10000.0;
  fc_mag_dn = 10000.0;
  general_data->stat_avg.count_diag_srot      = 0.0;
  general_data->stat_avg.fatm_mag = 10000.0;
  general_data->stat_avg.fatm_max = 10000.0;

/*======================================================================*/
/* Initial call to output_cp_min: need to open confp file               */

 if( atm_min == 1){
    general_data->filenames.ifile_open = 1;
    general_data->timeinfo.itime = 0;
    output_cp_min(class,general_data,bonded,cp,idum);
  }/*endif*/

/*======================================================================*/
/* II) Loop over the specified number of time steps */

  for(itime = 1;itime<=(general_data->timeinfo.ntime);itime++){
    general_data->timeinfo.itime = itime;
    class->energy_ctrl.itime     = itime;
    cputime(&(general_data->stat_avg.cpu1)); 

  /*---------------------------------------------------------------------*/
  /* 1) atm minimization                                                 */
    if(atm_step==1){
      general_data->stat_avg.write_cp_atm_flag = 1;
      if(general_data->minopts.min_std==1){
         move_atm_std(class,bonded,general_data,cp,ifirst_atm);
         iatm_count++;
      }/*endif*/
      if(general_data->minopts.min_cg==1){
         move_atm_cg(class,bonded,general_data,cp,ifirst_atm);
         iatm_count++;
         ifirst_atm = 0;
         class->clatoms_info.cg_reset_flag = 0;
      }/*endif*/
      if(general_data->minopts.min_diis==1 && atm_step==1){
         move_atm_diis(class,bonded,general_data,cp,ifirst_atm);
         iatm_count++;
      }/*endif*/
      check_atm_grad_mag(class,general_data,&f_atm_mag,&ireset,&idone,
                          general_data->minopts.tol_atom);

      write_dump_file_cp(class,bonded,general_data,cp);

      if(idone==1){
        idone_atm=1;
        atm_step = 0;
        general_data->simopts.cp_min        = 0;
        general_data->simopts.cp_wave_min   = 1;
        if(myid==0){
         printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         printf("Tolerance on atomic forces reached\n");
         printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        }/*endif*/
        write_dump_file_cp(class,bonded,general_data,cp);
        iexit = 1;
      }else{
        atm_step = 0;
        general_data->simopts.cp_min        = 0;
        general_data->simopts.cp_wave_min   = 1;
        if((ireset==1)&&(general_data->minopts.min_cg==1)){
         class->clatoms_info.cg_reset_flag = 1;
         ifirst_atm = 1;
         if(myid==0){
          printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
          printf("Resetting Atomic Conjugate Gradient\n");
          printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         }/*endif*/
         cp_shuffle_states(cp,1);
        }/*endif:reset*/
      }/*endif:not done*/
    }/*endif:atm move*/

  /*---------------------------------------------------------------------*/
  /* 2) CP_wave minimization                                             */
    if(atm_step==0){
      elec_e_old     = general_data->stat_avg.cp_eke
                     + general_data->stat_avg.cp_enl
                     + general_data->stat_avg.cp_ehart
                     + general_data->stat_avg.cp_exc
                     + general_data->stat_avg.cp_eext;
      if(num_proc>1){
       elec_e_old_tmp = elec_e_old;
       Allreduce(&(elec_e_old_tmp),&(elec_e_old),1,MPI_DOUBLE,MPI_SUM,0,
                 comm_states);
      }/*endif*/
      if(general_data->minopts.cp_min_std==1){
         min_STD_cp(class,bonded,general_data,cp,ip_now);
         cp_shuffle_states(cp,1);
      }/*endif*/
      if((general_data->minopts.cp_min_cg==1)){
         min_CG_cp(class,bonded,general_data,cp,ip_now);
         cp->cpcoeffs_info.cg_reset_flag = 0;
      }/*endif*/
      if((general_data->minopts.cp_min_diis==1)){
         min_DIIS_cp(class,bonded,general_data,cp,iatm_count,ip_now);
         cp->cpcoeffs_info.diis_reset_flag = 0;
      }/*endif*/
      elec_e         = general_data->stat_avg.cp_eke
                     + general_data->stat_avg.cp_enl
                     + general_data->stat_avg.cp_ehart
                     + general_data->stat_avg.cp_exc
                     + general_data->stat_avg.cp_eext;
      if(num_proc>1){
       elec_e_tmp = elec_e;
       Allreduce(&(elec_e_tmp),&(elec_e),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
      }/*endif*/
      check_coef_grad_mag(cp,&(general_data->simopts),
                          &fc_mag_up,&fc_mag_dn,&ireset,&idone,
                          general_data->minopts.tol_coef,1,1,
                          &(general_data->stat_avg));
      if(idone==1){
       if(myid==0){
        printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        printf("Tolerance on coefficient forces reached\n");
        if(atm_min==1){printf("Moving atoms on next step\n");}
        printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       }/* endif */
        if(atm_min==1){
         atm_step = 1;
         general_data->simopts.cp_min         = 1;
         general_data->simopts.cp_wave_min    = 0;
         cp->cpcoeffs_info.cg_reset_flag = 1;
         cp->cpcoeffs_info.diis_reset_flag = 1;
         control_coef_transpose_bck(cp,1);
/*         write_dump_file_cp(class,bonded,general_data,cp); */
         control_coef_transpose_fwd(cp,1);
        }else{
         iexit = 1;
        }/*endif*/
      }else{
        if((ireset==1)&&(general_data->minopts.cp_min_cg==1)&&
            (elec_e-elec_e_old)>0){
         cp->cpcoeffs_info.cg_reset_flag = 1;
         if(myid==0){
          printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
          printf("Resetting CP Conjugate Gradient\n");
          printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         }/* endif */
         cp_shuffle_states(cp,1);
        }/*endif:reset cg*/
      }/*endif:not done*/
    }/*endif:cp_wave move*/

  /*----------------------------------------------------------------------*/
  /*   6)Calculate some simple averages                                   */
    cputime(&(general_data->stat_avg.cpu2)); 
    (general_data->stat_avg.cpu_now)=(general_data->stat_avg.cpu2)-
                                     (general_data->stat_avg.cpu1);
    (general_data->stat_avg.acpu) += (general_data->stat_avg.cpu2)-
                                     (general_data->stat_avg.cpu1);

    simpavg_cp(&(general_data->timeinfo),&(general_data->stat_avg),
              &(general_data->cell),&(bonded->constrnt),
              &(general_data->ensopts),&(general_data->simopts),
              &(general_data->ptens),cp,&(class->communicate),
              &(class->nbr_list.verlist),&(class->energy_ctrl));

  /*-----------------------------------------------------------------------*/
  /*   7)Produce the output specified by the user                          */
    if(  (itime % (general_data->filenames.iwrite_screen))==0 ||
         (itime % (general_data->filenames.iwrite_dump  ))==0 ||
         (itime % (general_data->filenames.iwrite_confp ))==0 ||
         (itime % (general_data->filenames.iwrite_confc ))==0 ||
         (itime % (general_data->filenames.iwrite_confv)) ==0 ||
         (itime % (general_data->filenames.iwrite_inst))  ==0 ||
          idone_atm == 1                                        ){
         general_data->filenames.ifile_open = 0;

         output_cp_min(class,general_data,bonded,cp,idum);
         general_data->stat_avg.write_cp_atm_flag = 0;
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*  8) Analysis Routine                                               */
    analysis_cp(class,bonded,general_data,cp,analysis); 

  /*---------------------------------------------------------------------*/
    if(iexit==1){break;}

  }/*endfor:itime */

  /*======================================================================*/
  /*  II) Final dump  : get all energyies and write EVERYTHING            */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;
    (general_data->simopts.cp_min)      = 1;
    (general_data->simopts.cp_wave_min) = 0;

    cp_energy_control(class,bonded,general_data,cp);

    general_data->filenames.ifile_open    = 0;
    general_data->filenames.iwrite_screen = itime;
    general_data->filenames.iwrite_dump   = itime;
    general_data->stat_avg.write_cp_atm_flag = 1;

    output_cp_min(class,general_data,bonded,cp,idum);

  /*======================================================================*/
  /*  III)Write to Screen                                                 */

 if(myid==0){
  PRINT_LINE_DASH;
  printf("Completed CP-MINIMIZATION run \n");
  PRINT_LINE_STAR;
 }/* endif */

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void check_coef_grad_mag(CP *cp,SIMOPTS *simopts,
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
 double fc_mag_up;
 double fc_mag_dn;
 double fc_max_up;
 double fc_max_dn;
 double fc_mag_up_tmp;
 double fc_mag_dn_tmp;
 double fc_max_up_tmp;
 double fc_max_dn_tmp;
 int pi_beads      = cp->cpcoeffs_info.pi_beads;
 int pi_beads_proc = cp->cpcoeffs_info.pi_beads_proc;
 int ip;
 int i,idone,ireset; 
 int ncoef_tot;
 double *fcre_up;
 double *fcim_up;
 double *fcre_dn;
 double *fcim_dn;
 double *cre_up;
 double *cim_up;
 double *cre_dn;
 double *cim_dn;
 double *occ_up     = cp->cpopts.occ_up;
 double *occ_dn     = cp->cpopts.occ_dn;
 double *ksmat_up;
 double *ksmat_dn;
 double *ksmat_scr    = cp->cpscr.cpscr_ovmat.ovlap1;
 int *ioff_upt      = cp->cpcoeffs_info.ioff_upt;
 int *ioff_dnt      = cp->cpcoeffs_info.ioff_dnt;
 int cp_norb     = cp->cpopts.cp_norb;
 int cp_min;
 int ncoef_up    = cp->cpcoeffs_info.ncoef;
 int ncoef_dn    = cp->cpcoeffs_info.ncoef;
 int nstate_up = cp->cpcoeffs_info.nstate_up;
 int nstate_dn = cp->cpcoeffs_info.nstate_dn;
 int cp_lsda = cp->cpopts.cp_lsda;
 int myid_state = cp->communicate.myid_state;
 int myid_bead = cp->communicate.myid_bead;
 int np_states  = cp->communicate.np_states;
 int np_beads  = cp->communicate.np_beads;
 int icoef_orth_up;
 int icoef_orth_dn;
 int icoef_form_up;
 int icoef_form_dn;
 int ifcoef_orth_up;
 int ifcoef_orth_dn;
 int ifcoef_form_up;
 int ifcoef_form_dn;
 MPI_Comm comm_states = cp->communicate.comm_states;
 MPI_Comm comm_beads = cp->communicate.comm_beads;

/*=======================================================================*/
/* 0) Parallel checks                                         */

  if(np_states>1){
    for(ip=ip_start;ip<=ip_end;ip++){
     ifcoef_form_up = cp->cpcoeffs_pos[ip].ifcoef_form_up;
     ifcoef_form_dn = cp->cpcoeffs_pos[ip].ifcoef_form_dn;
      ncoef_up = cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
      if(ifcoef_form_up!=1){
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("Up Coef forces are not in transposed form \n");
         printf("on state processor %d in check_coef_grad \n",myid_state);
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
         exit(1);
       }/*endif*/
       if(cp_lsda==1){
        ncoef_dn = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;
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

   cre_up         = cp->cpcoeffs_pos[ip].cre_up;
   cim_up         = cp->cpcoeffs_pos[ip].cim_up;
   cre_dn         = cp->cpcoeffs_pos[ip].cre_dn;
   cim_dn         = cp->cpcoeffs_pos[ip].cim_dn;
   icoef_orth_up  = cp->cpcoeffs_pos[ip].icoef_orth_up;
   icoef_orth_dn  = cp->cpcoeffs_pos[ip].icoef_orth_dn;
   icoef_form_up  = cp->cpcoeffs_pos[ip].icoef_form_up;
   icoef_form_dn  = cp->cpcoeffs_pos[ip].icoef_form_dn;
   ifcoef_orth_up = cp->cpcoeffs_pos[ip].ifcoef_orth_up;
   ifcoef_orth_dn = cp->cpcoeffs_pos[ip].ifcoef_orth_dn;
   ifcoef_form_up = cp->cpcoeffs_pos[ip].ifcoef_form_up;
   ifcoef_form_dn = cp->cpcoeffs_pos[ip].ifcoef_form_dn;
   ksmat_up       = cp->cpcoeffs_pos[ip].ksmat_up;
   ksmat_dn       = cp->cpcoeffs_pos[ip].ksmat_dn;
  if(cp_min == 1 || cp_norb >= 1) {
     fcre_up = cp->cpcoeffs_pos[ip].fcre_up;
     fcim_up = cp->cpcoeffs_pos[ip].fcim_up;
     fcre_dn = cp->cpcoeffs_pos[ip].fcre_dn;
     fcim_dn = cp->cpcoeffs_pos[ip].fcim_dn;
  } else {
     fcre_up = cp->cpscr.cpscr_wave.cre_up;
     fcim_up = cp->cpscr.cpscr_wave.cim_up;
     fcre_dn = cp->cpscr.cpscr_wave.cre_dn;
     fcim_dn = cp->cpscr.cpscr_wave.cim_dn;
     for(i=1;i<=ncoef_tot; i++){
       fcre_up[i] = cp->cpcoeffs_pos[ip].fcre_up[i];
       fcim_up[i] = cp->cpcoeffs_pos[ip].fcim_up[i];
     }
     if(cp->cpopts.cp_lsda == 1 && nstate_dn != 0){
      for(i=1;i<=ncoef_tot; i++){
        fcre_dn[i] = cp->cpcoeffs_pos[ip].fcre_dn[i];
        fcim_dn[i] = cp->cpcoeffs_pos[ip].fcim_dn[i];
      }
     }/* endif */
     cp_add_ksmat_force(cre_up,cim_up,icoef_form_up,icoef_orth_up,
                        fcre_up,fcim_up,ifcoef_form_up,ifcoef_orth_up,
                        ksmat_up,ksmat_scr,ioff_upt,cp_lsda,cp_min,occ_up,
                        &(cp->cp_comm_state_pkg_up));
     if( (cp_lsda==1) && (nstate_dn!=0) ){
       cp_add_ksmat_force(cre_dn,cim_dn,icoef_form_dn,icoef_orth_dn,
                          fcre_dn,fcim_dn,ifcoef_form_dn,ifcoef_orth_dn,
                          ksmat_dn,ksmat_scr,ioff_dnt,cp_lsda,cp_min,occ_dn,
                          &(cp->cp_comm_state_pkg_dn));
     }/*endif*/
   }/* endif */

 

/*-----------------------------------------------------------------------*/
/* I) Up tolerence */

  for(i=1;i <= ncoef_tot; i++) {
    fc_mag_up += fcre_up[i] * fcre_up[i] 
                + fcim_up[i] * fcim_up[i];
    fc_max_up = MAX(fabs(fcre_up[i]),fc_max_up);
    fc_max_up = MAX(fabs(fcim_up[i]),fc_max_up);
  }/*endfor*/

/*-----------------------------------------------------------------------*/
/* II) Dn tolerence */

  if(cp_lsda == 1){
    ncoef_tot = ncoef_dn*nstate_dn;
    for(i=1;i <= ncoef_tot; i++) {
     fc_mag_dn += fcre_dn[i] * fcre_dn[i] 
                 + fcim_dn[i]  * fcim_dn[i];
     fc_max_dn = MAX(fabs(fcre_dn[i]),fc_max_dn);
     fc_max_dn = MAX(fabs(fcim_dn[i]),fc_max_dn);
    }/*endfor*/
  }/*endif*/

/*-----------------------------------------------------------------------*/
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



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void check_atm_grad_mag(CLASS *class,GENERAL_DATA *general_data,
                        double *f_atm_mag_ret,
                        int *ireset_ret, int *idone_ret, double tol_atom)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
 
 double f_atm_mag_old = *f_atm_mag_ret;
 double f_atm_mag;
 double f_atm_max;
 int i,idone,ireset; 
 int pi_beads = class->clatoms_info.pi_beads;
 int natm_tot = class->clatoms_info.natm_tot;
 double *fx = class->clatoms_pos[1].fx;
 double *fy = class->clatoms_pos[1].fy;
 double *fz = class->clatoms_pos[1].fz;
 int myid = class->communicate.myid;
 int np_states = class->communicate.np_states;
 MPI_Comm comm_states= class->communicate.comm_states;

/*==========================================================================*/
/*==========================================================================*/

 if(myid==0){

 /*==========================================================================*/
 /* I) Check atm tolerence */

  f_atm_mag = 0.0;
  f_atm_max = 0.0;
  for(i=1;i <= natm_tot;i++){
    f_atm_mag += (fx[i]*fx[i] + fy[i]*fy[i] + fz[i]*fz[i]);
    f_atm_max  = MAX(fabs(fx[i]),f_atm_max);
    f_atm_max  = MAX(fabs(fy[i]),f_atm_max);
    f_atm_max  = MAX(fabs(fz[i]),f_atm_max);
  }/*endfor*/
  f_atm_mag = sqrt(f_atm_mag/((double) 3*natm_tot));

 /*==========================================================================*/
 /* II) Set the Flags */

  idone = 1;
  if(f_atm_mag > tol_atom){idone=0;}
  ireset = 0;
  if(f_atm_mag > f_atm_mag_old){ireset=1;} 

 }/*endif*/

/*=======================================================================*/
/*=======================================================================*/
/* III) Set return values                                                */
 
 if(np_states>1){
   Bcast(&idone,1,MPI_INT,0,comm_states);
   Bcast(&ireset,1,MPI_INT,0,comm_states);
 }/*endif*/
 *f_atm_mag_ret = f_atm_mag;
 *idone_ret     = idone;
 *ireset_ret    = ireset;
 general_data->stat_avg.fatm_mag = f_atm_mag;
 general_data->stat_avg.fatm_max = f_atm_max;

/*-----------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_diag_cp_hess(CP *cp,CELL *cell,double scale)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

 int *kastore = cp->cpewald.kastr_sm;
 int *kbstore = cp->cpewald.kbstr_sm;
 int *kcstore = cp->cpewald.kcstr_sm;
 int ncoef = cp->cpcoeffs_info.ncoef;
 int i;
 int iii;
 int   ncoef_l = cp->cpcoeffs_info.ncoef_l;
 int cp_lsda = cp->cpopts.cp_lsda;
 int nstate_dn = cp->cpcoeffs_info.nstate_dn;
 int myid_state = cp->communicate.myid_state;
 int np_states = cp->communicate.np_states;
 int icmoff_up        = cp->cpcoeffs_info.icoef_start_up-1;
 int icmoff_dn        = cp->cpcoeffs_info.icoef_start_dn-1;
 MPI_Comm comm_states = cp->communicate.comm_states;
 double tpi;
 double aka,akb,akc;
 double xk,yk,zk;
 double enow,cmass_cut_def;
 double deth,*hmati;
 double hess0;
 double pseud_hess_loc = cp->cpcoeffs_info.pseud_hess_loc;
 double hess_tau = cp->cpcoeffs_info.cp_hess_tau;
 double hess_cut = cp->cpcoeffs_info.cp_hess_cut;
 double *cp_hess_up = cp->cpcoeffs_info.cp_hess_up;
 double *cp_hess_dn = cp->cpcoeffs_info.cp_hess_dn;

/*==========================================================================*/
/* 0) Get some useful quantities                                            */

   hmati = (double *) cmalloc((size_t)9*sizeof(double))-1;
   gethinv(cell->hmat,hmati,&deth,cell->iperd);    

/*==========================================================================*/
/* -I) Nonlocal contribution comes in from cp_energy_eext_nonloc.c         */
 
 for(i=1;i<=ncoef;i++){
   cp_hess_up[i] *= (double) np_states;
 }

 if(cp_lsda == 1 && nstate_dn != 0){
  for(i=1;i<=ncoef;i++){
    cp_hess_dn[i] *= (double) np_states;
  }
 }/* endif */


/*==========================================================================*/
/* I) Kinetic contribution                                               */

  tpi = 2.0*M_PI;
  for(i=1;i<=ncoef-1;i++) {
    aka = (double)(kastore[i]);
    akb = (double)(kbstore[i]);
    akc = (double)(kcstore[i]);
    xk = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
    yk = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
    zk = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
    cp_hess_up[i] += 0.5*(xk*xk+yk*yk+zk*zk);
  } /* endfor */
  if(cp_lsda == 1 && nstate_dn != 0){
   for(i=1;i<=ncoef-1;i++) {
     aka = (double)(kastore[i]);
     akb = (double)(kbstore[i]);
     akc = (double)(kcstore[i]);
     xk = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
     yk = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
     zk = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
     cp_hess_dn[i] += 0.5*(xk*xk+yk*yk+zk*zk);
   } /* endfor */
  }/* endif */

/*==========================================================================*/
/* II) Local external potential contribution                                */
   
  if(cp->cp_comm_state_pkg_up.num_proc > 1){
   Bcast(&pseud_hess_loc,1,MPI_DOUBLE,np_states-1,comm_states); 
  }
  for(i=1;i<=ncoef;i++){
   cp_hess_up[i] += pseud_hess_loc;
  }
  if(cp_lsda == 1 && nstate_dn != 0){
   for(i=1;i<=ncoef;i++){
    cp_hess_dn[i] += pseud_hess_loc;
   }
  }/* endif */
   
/*==========================================================================*/
/* Cut it off at low g-values and set it to appropriate constant            */

#define CP_EMAGIC 0.00050

  hess0 = 4.0*hess_tau*hess_tau*CP_EMAGIC;
  for(i=1;i<=ncoef-1;i++) {
    if(cp_hess_up[i] > 0.5*hess_cut){
      cp_hess_up[i] *= hess0;
    } else {
      cp_hess_up[i] = hess0*0.5*hess_cut;
    }/* endif */
  }/* endfor */

 if(cp_lsda == 1 && nstate_dn != 0){
  for(i=1;i<=ncoef-1;i++) {
    if(cp_hess_dn[i] > 0.5*hess_cut){
      cp_hess_dn[i] *= hess0;
    } else {
      cp_hess_dn[i] = hess0*0.5*hess_cut;
    }/* endif */
  }/* endfor */
 }/* endif */

  cp_hess_up[ncoef] = hess0*0.25*hess_cut;
  if(cp_lsda == 1 && nstate_dn != 0){
   cp_hess_dn[ncoef] = hess0*0.25*hess_cut;
  }/* endif */

/*-----------------------------------------------------------------------*/
  cfree(&(hmati[1]));
/*-----------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/










