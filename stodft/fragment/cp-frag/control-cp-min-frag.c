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

void controlCpMinFrag(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
			CP *cp,ANALYSIS *analysis)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  double fc_mag_up,fc_mag_dn,f_atm_mag,elec_e_old,elec_e;
  double elec_e_old_tmp,elec_e_tmp;
  double Delta_E = 0.0;
  double delta_finite_diff;
  double fp,fm;
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
  int calcul_freq_on;
  int do_analysis;
  int displace_index=0;
  int sign_index = 1;
  int natm_tot = class->clatoms_info.natm_tot;
  MPI_Comm comm_states = class->communicate.comm_states;
  MPI_Comm world       = class->communicate.world;
  
  int iopt_cp_pw       = cp->cpcoeffs_info.iopt_cp_pw;
  int iopt_cp_dvr      = cp->cpcoeffs_info.iopt_cp_dvr;

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

/*======================================================================*/
/* II) Loop over the specified number of time steps */

  for(itime = 1;itime<=(general_data->timeinfo.ntime);itime++){
    general_data->timeinfo.itime = itime;
    class->energy_ctrl.itime     = itime;
    cputime(&(general_data->stat_avg.cpu1)); 
  /*---------------------------------------------------------------------*/
  /* 2) CP_wave minimization                                             */
    elec_e_old = general_data->stat_avg.cp_eke+general_data->stat_avg.cp_enl
		+general_data->stat_avg.cp_ehart+general_data->stat_avg.cp_exc
		+general_data->stat_avg.cp_eext;
    /*
    if(general_data->minopts.cp_min_std==1){
      if(iopt_cp_pw){
	min_STD_cp(class,bonded,general_data,cp,ip_now);
	cp_shuffle_states(cp,1);
      }
    }//endif
    */
    if((general_data->minopts.cp_min_cg==1)){//Use CG only
      if(iopt_cp_pw) min_CG_cp(class,bonded,general_data,cp,ip_now);
      cp->cpcoeffs_info.cg_reset_flag = 0;
    }//endif
    /*
    if((general_data->minopts.cp_min_diis==1)){
      if(iopt_cp_pw) min_DIIS_cp(class,bonded,general_data,cp,iatm_count,ip_now);
      cp->cpcoeffs_info.diis_reset_flag = 0;
    }//endif
    */
    elec_e = general_data->stat_avg.cp_eke+general_data->stat_avg.cp_enl
	    +general_data->stat_avg.cp_ehart+general_data->stat_avg.cp_exc
	    +general_data->stat_avg.cp_eext;
    Delta_E = fabs(elec_e - elec_e_old);
    if(iopt_cp_pw){
      check_coef_grad_mag(cp,&(general_data->simopts),
			  &fc_mag_up,&fc_mag_dn,&ireset,&idone,
			  general_data->minopts.tol_coef,1,1,
			  &(general_data->stat_avg));
    }
    if(idone==1) iexit = 1;
    else{
      if((ireset==1)&&(general_data->minopts.cp_min_cg==1)&&(elec_e-elec_e_old)>0){
	cp->cpcoeffs_info.cg_reset_flag = 1;
	if(iopt_cp_pw) cp_shuffle_states(cp,1);
      }/*endif:reset cg*/
    }/*endif:not done*/

  /*----------------------------------------------------------------------*/
  /*   6)Calculate some simple averages                                   */
    /*
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
    */

  /*-----------------------------------------------------------------------*/
  /*   7)Produce the output specified by the user                          */

    /*
    if(myid == 0) check_auto_exit(&(general_data->timeinfo.exit_flag));
    if(num_proc > 1) Bcast(&(general_data->timeinfo.exit_flag),1,MPI_INT,0,world);

    if((itime % (general_data->filenames.iwrite_screen))==0 ||
       (itime % (general_data->filenames.iwrite_dump  ))==0 ||
       (itime % (general_data->filenames.iwrite_confp ))==0 ||
       (itime % (general_data->filenames.iwrite_confc ))==0 ||
       (itime % (general_data->filenames.iwrite_confv)) ==0 ||
       (itime % (general_data->filenames.iwrite_inst))  ==0 ||
       (general_data->timeinfo.exit_flag == 1) ||idone_atm == 1){
      general_data->filenames.ifile_open = 0;
      output_cp_min(class,general_data,bonded,cp,Delta_E,idum);
      general_data->stat_avg.write_cp_atm_flag = 0;
    }//endif
    */

  /*---------------------------------------------------------------------*/
  /*  8) Analysis Routine                                               */

    /*
    if(calcul_freq_on == 1 && do_analysis == 1) {
       general_data->simopts.cp_min = 1;
       assign_hessian(class,bonded,general_data,cp,delta_finite_diff,displace_index,sign_index);
       general_data->simopts.cp_min = 0;
    }

    if(do_analysis == 1) analysis_cp(class,bonded,general_data,cp,analysis); 
    */

  /*---------------------------------------------------------------------*/
    if(iexit==1){break;}

 /*---------------------------------------------------------------------*/
 /*   Check for exit condition                                          */

    if(general_data->timeinfo.exit_flag == 1) itime = general_data->timeinfo.ntime;

  }/*endfor:itime */

  /*======================================================================*/
  /*  II) Final dump  : get all energyies and write EVERYTHING            */
    
    /*
    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;
    (general_data->simopts.cp_min)      = 1;
    (general_data->simopts.cp_wave_min) = 0;

    cp_energy_control(class,bonded,general_data,cp);

    simpavg_cp(&(general_data->timeinfo),&(general_data->stat_avg),
               &(general_data->cell),&(bonded->constrnt),
               &(general_data->ensopts),&(general_data->simopts),
               &(general_data->ptens),cp,&(class->communicate),
               &(class->nbr_list.verlist),&(class->energy_ctrl));

    general_data->filenames.ifile_open    = 0;
    general_data->filenames.iwrite_screen = itime;
    general_data->filenames.iwrite_dump   = itime;
    general_data->stat_avg.write_cp_atm_flag = 1;

    output_cp_min(class,general_data,bonded,cp,Delta_E,idum);
    */

  /*======================================================================*/
  /*  III)Write to Screen                                                 */

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/












