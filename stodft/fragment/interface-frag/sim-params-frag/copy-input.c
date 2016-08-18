/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: copy-input.c	                            */
/*                                                                          */
/*                                                                          */
/* This subprogram copies the simulation parameters to mini structures.	    */
/*				                                            */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_sim_params_entry.h"
#include "../proto_defs/proto_sim_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void copySimParam(GENERAL_DATA *general_data,BONDED *bonded,CLASS *class,
		  CP *cp,GENERAL_DATA *generalDataMini,BONDED *bondedMini,
		  CLASS *classMini,CP *cpMini)
/*========================================================================*/
/*             Begin Routine                                              */
{/*Begin subprogram: */
/*************************************************************************/
/* Copy the simulation parameters (in input file) down to mini           */
/* structures. The order is exactly the same as control_sim_params       */
/*************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */

/*=======================================================================*/
/*   I) set_sim_params_gen						 */
  
  /* 1) \simulation_typ */
  generalDataMini->simopts.minimize	    = general_data->simopts.minimize;
  generalDataMini->simopts.md		    = general_data->simopts.md;
  generalDataMini->simopts.pimd		    = general_data->simopts.pimd;
  generalDataMini->simopts.cp		    = general_data->simopts.cp;
  generalDataMini->simopts.cp_wave	    = general_data->simopts.cp_wave;
  generalDataMini->simopts.cp_wave_pimd	    = general_data->simopts.cp_wave_pimd;
  generalDataMini->simopts.cp_pimd	    = general_data->simopts.cp_pimd;
  generalDataMini->simopts.cp_min	    = general_data->simopts.cp_min;
  generalDataMini->simopts.cp_wave_min	    = general_data->simopts.cp_wave_min;
  generalDataMini->simopts.cp_wve_min_pimd  = general_data->simopts.cp_wave_min_pimd;
  generalDataMini->simopts.debug	    = general_data->simopts.debug;
  generalDataMini->simopts.debug_pimd	    = general_data->simopts.debug_pimd;
  generalDataMini->simopts.debug_cp	    = general_data->simopts.debug_cp;
  generalDataMini->simopts.debug_cp_pimd    = general_data->simopts.debug_cp_pimd;
  /* 2) \ensemble_typ */
  generalDataMini->ensopts.nve = general_data->ensopts.nve;
  generalDataMini->ensopts.nvt = general_data->ensopts.nvt;
  generalDataMini->ensopts.npt_i = general_data->ensopts.npt_i;
  generalDataMini->ensopts.npt_f = general_data->ensopts.npt_f;
  generalDataMini->ensopts.nst = general_data->ensopts.nst;
  /* 3) \num_time_step */
  generalDataMini->timeinfo.ntime = general_data->timeinfo.ntime;
  /* 4) \time_step */ 
  generalDataMini->timeinfo.dt = general_data->timeinfo.dt;
  /* 5) \temperature */
  /* I don't need this. Pass some fake number to stop warnning messages */
  generalDataMini->statepoint.t_ext = general_data->statepoint.t_ext;
  /* 6) \pressure{#} */
  /* I may need this if I want to do npt optimization */
  generalDataMini->statepoint.pext = general_data->statepoint.pext;  
  /* 8)\minimize_typ{min_std,min_cg,min_diis} */
  /* I don't do coord minimization pass some fake number. */
  generalDataMini->minopts.min_std = general_data->minopts.min_std;
  generalDataMini->minopts.min_cg = general_data->minopts.min_cg;
  generalDataMini->minopts.min_diis = general_data->minopts.min_diis;
  /* 9)\annealing_rate{#} */
  /* Not needed, just pass the value */
  generalDataMini->simopts.ann_rate = general_data->simopts.ann_rate;
  /* 10)\num_proc_beads{#} */
  /* 1 */
  classMini->communicate.np_beads = 1;
  /* 11)\num_proc_states{#} */
  /* 1 This is important*/
  classMini->communicate.np_states = 1;
  /* 12)\num_proc_class_forc{#} */
  /* 1 */
  classMini->communicate.np_forc = 1;
  /* 13)\num_proc_tot{#} */
  /* 1 */ 
  classMini->communicate.np = 1;
  classMini->communicate.myid = 0;
  /* 14)\rndm_seed{#} */
  /* Not really need it since we don't have random stuff. */
  classMini->vel_samp_class.qseed = class->vel_samp_class.qseed;
  /* 15)\rndm_seed2{#} */
  classMini->vel_samp_class.iseed2 = class->vel_samp_class.iseed2;
  cpMini->vel_samp_cp.iseed2 = cp->vel_samp_cp.iseed2;
  /* 16)\generic_fft_opt{on,off} */
  /* You are on, I am on. */
  cpMini->cp_sclr_fft_pkg3d_sm.igeneric_opt = cp->cp_sclr_fft_pkg3d_sm.igeneric_opt;
  cpMini->cp_para_fft_pkg3d_sm.igeneric_opt = cp->cp_para_fft_pkg3d_sm.igeneric_opt;
  cpMini->cp_sclr_fft_pkg3d_lg.igeneric_opt = cp->cp_sclr_fft_pkg3d_lg.igeneric_opt; 
  cpMini->cp_para_fft_pkg3d_lg.igeneric_opt = cp->cp_para_fft_pkg3d_lg.igeneric_opt;
  cpMini->cp_sclr_fft_pkg3d_dens_cp_box.igeneric_opt = cp->cp_sclr_fft_pkg3d_dens_cp_box.igeneric_opt;
  cpMini->cp_para_fft_pkg3d_dens_cp_box.igeneric_opt = cp->cp_para_fft_pkg3d_dens_cp_box.igeneric_opt;
  generalDataMini->pme_res_fft_pkg.igeneric_opt = general_data->pme_res_fft_pkg.igeneric_opt;
  generalDataMini->pme_fft_pkg.igeneric_opt = general_data->pme_fft_pkg.igeneric_opt;
  /* 17)\alpha_clus{#}   */
  generalDataMini->ewald.alp_clus = general_data->ewald.alp_clus;
  /* 18)\ecut_clus{#}   */
  generalDataMini->ewald.ecut_clus = general_data->ewald.ecut_clus;
  /* 19)\surf_tens{#} */
  generalDataMini->statepoint.stens_ext = general_data->statepoint.stens_ext;
  /* 20)\min_num_atoms_per_proc{#}   */
  classMini->clatoms_info.natm_proc = class->clatoms_info.natm_proc;
  /* 21)\num_proc_class_forc_src{#} */
  classMini->communicate.np_forc_src = classMini->communicate.np_forc;
  classMini->class_comm_forc_pkg.plimpton_ind.num_proc_source 
				= classMini->communicate.np_forc;
  /* 22)\num_proc_class_forc_trg{#} */
  classMini->communicate.np_forc_trg = class->communicate.np_forc_trg;
  classMini->class_comm_forc_pkg.plimpton_ind.num_proc_target 
	    = class->class_comm_forc_pkg.plimpton_ind.num_proc_target;
  /* 23)\annealing_opt{on,off} */
  generalDataMini->simopts.anneal_opt = general_data->simopts.anneal_opt;
  /* 24)\ann_start_temperature{#} */
  generalDataMini->simopts.ann_start_temp = general_data->simopts.ann_start_temp;

  

/*------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/




