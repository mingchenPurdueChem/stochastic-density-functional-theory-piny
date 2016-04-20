/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_cp                                   */
/*                                                                          */
/* This subprogram performs Stochastic DFT calculation and Geometric        */
/* Minimization.		                                            */
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
void controlStodftMin(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
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
/* 0) Write to Screen                                                   */

 if(myid==0){
  PRINT_LINE_STAR;
  printf("Running Stochastic DFT.\n");
  PRINT_LINE_DASH;
 }/* endif */

/*======================================================================*/
/* I) Set the flags/counters                                           */

/*======================================================================*/
/* Initial call to output_cp_min: need to open confp file               */

/*======================================================================*/
/* II) Loop over the specified number of time steps */

  /*---------------------------------------------------------------------*/
  /* 1) atm minimization                                                 */

//Let's temperaly test the filter by given the density and chemical potential, test how the filter looks like on KS eigenfunctions. 

  /*---------------------------------------------------------------------*/
  /* 2) CP_wave minimization                                             */

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

  //Minimize with stochastic dft
  minStodft(class,bonded,general_data,cp,ip_now);
  elec_e         = general_data->stat_avg.cp_eke
		 + general_data->stat_avg.cp_enl
		 + general_data->stat_avg.cp_ehart
		 + general_data->stat_avg.cp_exc
		 + general_data->stat_avg.cp_eext;
  if(num_proc>1){
   elec_e_tmp = elec_e;
   Allreduce(&(elec_e_tmp),&(elec_e),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
  }/*endif*/
  Delta_E = fabs(elec_e - elec_e_old);

  /*----------------------------------------------------------------------*/
  /*   6)Calculate some simple averages                                   */

  /*-----------------------------------------------------------------------*/
  /*   7)Produce the output specified by the user                          */


  /*---------------------------------------------------------------------*/
  /*  8) Analysis Routine                                               */

 /*---------------------------------------------------------------------*/
 /*   Check for exit condition                                      */

  /*======================================================================*/
  /*  II) Final dump  : get all energyies and write EVERYTHING            */

  /*======================================================================*/
  /*  III)Write to Screen                                                 */

  if(myid==0){
    PRINT_LINE_DASH;
    printf("Completed Stochastic DFT run \n");
    PRINT_LINE_STAR;
  }/* endif */

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/






