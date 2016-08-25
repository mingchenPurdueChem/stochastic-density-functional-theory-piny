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
#include "../proto_defs/proto_interface_frag_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void copySimParam(GENERAL_DATA *general_data,BONDED *bonded,CLASS *class,
		  CP *cp,GENERAL_DATA *generalDataMini,BONDED *bondedMini,
		  CLASS *classMini,CP *cpMini,CLASS_PARSE *classParse,
		  CP_PARSE *cpParse,FILENAME_PARSE *fileNameParse)
/*========================================================================*/
/*             Begin Routine                                              */
{/*Begin subprogram: */
/*************************************************************************/
/* Copy the simulation parameters (in input file) down to mini           */
/* structures. The order is exactly the same as control_sim_params       */
/*************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */

  STODFTINFO    *stodftInfo	  = cp->stodftInfo;
  FRAGINFO      *fragInfo         = stodftInfo->fragInfo;

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
  /* 7)\restart_type{initial,restart_pos,restart_posvel,restart_all}*/
  classParse->istart = 1;
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

/*=======================================================================*/
/*  II) set_sim_params_list                                              */

  /* 1)\verlist_pad{#} */
  classMini->nbr_list.verlist.jver_pad = class->nbr_list.verlist.jver_pad;
  /* 2)\verlist_mem_safe{#} */
  classMini->nbr_list.verlist.mem_safe = class->nbr_list.verlist.mem_safe;
  /* 3)\verlist_mem_min{#} */
  classMini->nbr_list.verlist.nmem_min_lst = class->nbr_list.verlist.nmem_min_lst;
  /* 4)\verlist_skin{#} */
  classMini->interact.skin = class->interact.skin;
  /* 8)\lnkcell_cell_divs{#} */
  classMini->nbr_list.lnklist.ncell_div_avg = class->nbr_list.lnklist.ncell_div_avg;
  /* 9)\lnk_cell_excl_safe{#} */
  classMini->nbr_list.lnklist.lnk_excl_safe = class->nbr_list.lnklist.lnk_excl_safe;
  /* 10)\lnk_cell_vol_safe{#} */
  classMini->nbr_list.lnklist.lnk_vol_safe = class->nbr_list.lnklist.lnk_vol_safe;
  /* 11)\lnkcell_force_odd{on,off} */
  classMini->nbr_list.lnklist.lnk_for_odd = class->nbr_list.lnklist.lnk_for_odd; 
  /* 12)\shave_skin_opt{#} */
  classMini->interact.ishave_opt = class->interact.ishave_opt;
  /*  13)\neighbor_list{ver_list,lnk_list,no_list} */
  classMini->nbr_list.iver = class->nbr_list.iver;
  classMini->nbr_list.ilnk = class->nbr_list.ilnk;
  classMini->nbr_list.nolst = class->nbr_list.nolst;
  /*  14)\update_type{no_list,lnk_list} */
  classMini->nbr_list.verlist.nolst_ver_update = class->nbr_list.verlist.nolst_ver_update;
  classMini->nbr_list.verlist.lnk_ver_update = class->nbr_list.verlist.lnk_ver_update;
  /* 5)\brnch_root_list_opt{#} */
  classMini->nbr_list.brnch_root_list_opt = class->nbr_list.brnch_root_list_opt;
  /* 6)\brnch_root_list_skin{#} */
  classMini->interact.brnch_root_skin = class->interact.brnch_root_skin;
  /* 7)\brnch_root_cutoff{#} */
  classMini->interact.brnch_root_cut = class->interact.brnch_root_cut;

/*=======================================================================*/
/*  III) set_sim_params_cp                                               */

  /*  8)\cp_dft_typ{lsda,lda,gga_lda,gga_lsda} */
  cpMini->cpopts.cp_lda = cp->cpopts.cp_lda;
  cpMini->cpopts.cp_lsda = cp->cpopts.cp_lsda;
  cpMini->cpopts.cp_gga = cp->cpopts.cp_gga;
  /* 1)\cp_vxc_typ{pz_lda,pz_lsda,pw_lda,pw_lsda,pade_lda,pade_lsda}       */
  sscanf(cp->pseudo.vxc_typ,"%s",cpMini->pseudo.vxc_typ);
  /* 2)\cp_ggax_typ{becke,pw91x,debug97x,fila_1x,fila_2x,pbe_x,revpbe_x,rpbe_x,
       xpbe_x,brx89,brx2k,off}  */
  cpMini->cpopts.cp_gga = cp->cpopts.cp_gga;  
  cpMini->cpcoeffs_info.cp_laplacian_on = cp->cpcoeffs_info.cp_laplacian_on;
  cpMini->cpcoeffs_info.cp_ke_dens_on = cp->cpcoeffs_info.cp_ke_dens_on;
  cpMini->cpcoeffs_info.cp_tau_functional = cp->cpcoeffs_info.cp_tau_functional;
  sscanf(cp->pseudo.ggax_typ,"%s",cpMini->pseudo.ggax_typ);
  /* 3)\cp_ggac_typ{lyp,lypm1,pw91c,pbe_c,xpbe_c,tau1_c,off}               */
  sscanf(cp->pseudo.ggac_typ,"%s",cpMini->pseudo.ggac_typ);
  /* 4)\cp_sic{on,off} */
  cpMini->cpopts.cp_sic = cp->cpopts.cp_sic;
  /* 5)\cp_e_e_interact{on,off} */
  cpMini->cpopts.cp_nonint = cp->cpopts.cp_nonint;
  /* 6)\cp_norb{full_ortho,norm_only,no_constrnt,off} */
  cpMini->cpopts.cp_norb = cp->cpopts.cp_norb;
  /* 7)\cp_gauss{on,off} */
  cpMini->cpopts.cp_gauss = cp->cpopts.cp_gauss;
  /* 9)\cp_nl_list{on,off} */
  cpMini->pseudo.nl_cut_on = cp->pseudo.nl_cut_on;
  /* 10)\cp_mass_tau_def{#} some random number, anyway I'll not use this */
  cpParse->cp_mass_tau_def = 6.0;
  /* 11)\cp_mass_cut_def{#} some random number, anyway I'll not use this */
  cpParse->cp_mass_cut_def = 15.0;
  /* 12)\cp_energy_cut_def{#} */
  cpParse->cp_ecut_def = 2;
  /* 13)\cp_fict_KE{#} */
  cpMini->cpopts.te_ext = cp->cpopts.te_ext;
  /*  14)\cp_ptens{on,off} */
  cpMini->cpopts.cp_ptens_calc = cp->cpopts.cp_ptens_calc;
  /* 15)\cp_init_orthog{on,off} */
  cpMini->cpopts.cp_init_orthog = cp->cpopts.cp_init_orthog;
  /* 16)\cp_cg_line_min_len{#} */
  generalDataMini->minopts.cp_cg_line_min_len = general_data->minopts.cp_cg_line_min_len;
  /* 17)\cp_minimize_typ{min_std,min_cg,min_diis} */
  generalDataMini->minopts.cp_min_std = 0;
  generalDataMini->minopts.cp_min_cg = 1;
  generalDataMini->minopts.cp_min_diis = 0;
  cpMini->cpopts.cp_hess_calc = 0;
  /* 18)\cp_diis_hist_len{#} */
  generalDataMini->minopts.cp_diis_hist_len = general_data->minopts.cp_diis_hist_len;
  /* 19)\cp_orth_meth{gram_schmidt,lowdin,normalize,none} */
  cpMini->cpopts.cp_gs = 1;
  cp->cpopts.cp_low = 0;
  cp->cpopts.cp_normalize = 0;
  /* 20)\cp_restart_type{initial,restart_pos,restart_posvel,restart_all}*/
  cpParse->istart_cp = 0; //gen_wave
  /* 21)\diis_hist_len{#} I dont need it */
  generalDataMini->minopts.diis_hist_len = general_data->minopts.diis_hist_len;
  /* 22)\nlvps_list_skin{#}   */
  cpMini->pseudo.nlvps_skin = cp->pseudo.nlvps_skin;
  /* 23)\gradient_cutoff{#}   */
  cpMini->pseudo.gga_cut = cp->pseudo.gga_cut;
  /* 24)\zero_cp_vel{initial,periodic,no} I don't need this */
  cpMini->cpopts.zero_cp_vel = cp->cpopts.zero_cp_vel;
  /* 25)\cp_check_perd_size{#}   */
  cpMini->cpopts.icheck_perd_size = cp->cpopts.icheck_perd_size;
  /* 26)\cp_tol_edge_dist{#}   */
  cpMini->cpopts.tol_edge_dist = cp->cpopts.tol_edge_dist;
  /* 27)\cp_para_typ{#}   hybrid */
  cpMini->cpopts.cp_para_opt = 0;
  /* 28)\cp_dual_grid_opt{#} dual grid off */
  cpMini->cpopts.cp_dual_grid_opt = 0;
  /* 30)\cp_move_dual_box_opt{#}   */
  generalDataMini->cell.imov_cp_box = general_data->cell.imov_cp_box;
  /* 31)\cp_energy_cut_dual_grid_def{#} I don't need this*/
  cpParse->cp_ecut_dual_grid_def = 2;
  /* 32)\cp_alpha_conv_dual{#} */
  cpMini->pseudo.alpha_conv_dual = cp->pseudo.alpha_conv_dual;
  /* 33)\interp_pme_dual{#} */
  cpMini->pseudo.n_interp_pme_dual = cp->pseudo.n_interp_pme_dual;
  /* 34)\cp_elf_calc_frq{#} */
  cpMini->cpcoeffs_info.cp_elf_calc_frq = cp->cpcoeffs_info.cp_elf_calc_frq;
  /* 35)\cp_ngrid_skip{#} */
  cpMini->cpopts.cp_ngrid_skip = cp->cpopts.cp_ngrid_skip;
  /* 36)\cp_isok_opt{#} I don't need this */
  cpMini->cpopts.cp_isok_opt = cp->cpopts.cp_isok_opt;
  /* 37)\cp_hess_cut{#} */
  cpMini->cpcoeffs_info.cp_hess_cut = cp->cpcoeffs_info.cp_hess_cut;
  /* 38)\basis_set_opt{plane_wave,dvr} */
  cpMini->cpcoeffs_info.iopt_cp_pw  = 1;
  cpMini->cpcoeffs_info.iopt_cp_dvr = 0;
  /* 39)\dvr_grid_dens_def{#} I don't need this */
  cpParse->cp_dvrdens_def = 0.75;
  /* 40)\cp_nl_trunc_opt{#} */
  cpMini->cpopts.cp_nl_trunc_opt = cp->cpopts.cp_nl_trunc_opt;
  /* 41)\cp_wan_min_opt{#} I don't need it */
  cpMini->cpopts.cp_wan_min_opt = cp->cpopts.cp_wan_min_opt;
  /* 42)\cp_wan_opt{#} I dont need it */
  cpMini->cpopts.cp_wan_opt = cp->cpopts.cp_wan_opt;
  /* 43)\cp_wan_calc_frq{#} I dont need it */
  cpMini->cp_wannier.cp_wan_calc_frq = cp->cp_wannier.cp_wan_calc_frq;
  /* 44)\wan_func_typ{#} I dont need it */
  cpMini->cp_wannier.wan_func_typ = cp->cp_wannier.wan_func_typ;
  /* 45)\wan_diag_typ{#} I dont need it */
  cpMini->cp_wannier.wan_diag_typ = cp->cp_wannier.wan_diag_typ;
  /* 46)\cp_dip_calc_frq{#} */
  cpMini->cpcoeffs_info.cp_dip_calc_frq = cp->cpcoeffs_info.cp_dip_calc_frq;
  /* 47)\cp_nloc_wan_opt{#} I don't need it */
  cpMini->cpopts.cp_nloc_wan_opt = cp->cpopts.cp_nloc_wan_opt;
  /* 48)\cp_kinet_wan_opt{#} I don't need it */
  cpMini->cpopts.cp_kinet_wan_opt = cp->cpopts.cp_kinet_wan_opt;
  /* 49)\rcut_wan_orb{#} I dont need it */
  cpMini->cp_wannier.rcut_wan_orb = cp->cp_wannier.rcut_wan_orb;
  /* 50)\rcut_wan_nl{#} I dont need it */
  cpMini->cp_wannier.rcut_wan_nl = cp->cp_wannier.rcut_wan_nl;
  /* 51)\nmax_wan_orb{#} I dont need it */
  cpMini->cp_wannier.nmax_wan_orb = cp->cp_wannier.nmax_wan_orb;
  /* 52)\b3_cutoff{#} */
  cpMini->pseudo.b3_cut = cp->pseudo.b3_cut;
  /* 53)\b3_alpha{#}   */
  cpMini->pseudo.b3_alp = cp->pseudo.b3_alp;
  /* 54)\cp_wan_init_opt{#} I dont need it */
  cpMini->cpopts.cp_wan_init_opt = cp->cpopts.cp_wan_init_opt;
  /* 55)\cp_init_min_opt{#} */
  cpMini->cpopts.cp_init_min_opt =  cp->cpopts.cp_init_min_opt;
  /* 56)\iwrite_init_wcent{#} I dont need it */
  cpMini->cpopts.iwrite_init_wcent = cp->cpopts.iwrite_init_wcent;
  /* 57)\iwrite_init_worb{#} I dont need it */
  cpMini->cpopts.iwrite_init_worb = cp->cpopts.iwrite_init_worb;
  /* 58)\iwrite_init_state{#} I dont need it */
  cpMini->cpopts.iwrite_init_state = cp->cpopts.iwrite_init_state;
  /* 59)\dvr_clus_tmax{#}{#} I dont need it */
  cpMini->cp_dvr_clus.tmax = cp->cp_dvr_clus.tmax;
  /* 60)\dvr_clus_num_twindow{#} I dont need it */
  cpMini->cp_dvr_clus.num_twindow = cp->cp_dvr_clus.num_twindow ;
  /* 61)\dvr_clus_num_tquad{#} I dont need it */
  cpMini->cp_dvr_clus.num_tquad = cp->cp_dvr_clus.num_tquad;
  /* 62)\dvr_clus_grid_dens{#} I dont need it */
  cpMini->cp_dvr_clus.grid_dens = cp->cp_dvr_clus.grid_dens;
  /* 63)\dvr_clus_rmax{#} I dont need it */
  cpMini->cp_dvr_clus.rmax = cp->cp_dvr_clus.rmax;

/*=======================================================================*/
/*  IV) set_sim_params_vpot                                              */

  /* 1)\shift_inter_pe{on,off} cp_parse */
  classParse->ishift_pot = 2;
  classMini->energy_ctrl.iswit_vdw = class->energy_ctrl.iswit_vdw;
  classMini->interact.iswit_vdw = classMini->energy_ctrl.iswit_vdw;
  generalDataMini->stat_avg.iswit_vdw = classMini->energy_ctrl.iswit_vdw;
  /* 2)\inter_spline_pts{#} */
  classMini->interact.nsplin = class->interact.nsplin;
  bondedMini->ecor.nsplin = bonded->ecor.nsplin;
  bondedMini->ecor.nsplin_m2 = bonded->ecor.nsplin_m2;
  /* 3)\intra_block_min{} I dont need this */
  classMini->energy_ctrl.nblock_min = class->energy_ctrl.nblock_min;
  /* 4)\pten_inter_respa{#} I dont need this */
  classMini->interact.pten_inter_guess = class->interact.pten_inter_guess;
  /* 5)\pten_kin_respa{#} I dont need this */
  classMini->interact.pten_kin_guess = class->interact.pten_kin_guess;
  /* 6)\pseud_spline_pts{#} */
  cpMini->pseudo.nsplin_g = cp->pseudo.nsplin_g;
  /* 13)\ewald_interp_pme{#} */
  classMini->part_mesh.n_interp = class->part_mesh.n_interp;
  cpMini->cpscr.cpscr_atom_pme.n_interp = cp->cpscr.cpscr_atom_pme.n_interp;
  /* 7)\scratch_length{#} */
  bondedMini->intra_scr.nlen = bonded->intra_scr.nlen;
  classMini->for_scr.nlen = class->for_scr.nlen;
  classMini->part_mesh.nlen_pme = class->part_mesh.nlen_pme;
  cpMini->cpscr.cpscr_atom_pme.nlen_pme = cp->cpscr.cpscr_atom_pme.nlen_pme;
  /* 8)\ewald_alpha{#} I may need to change this */
  generalDataMini->ewald.alp_ewd = general_data->ewald.alp_ewd;
  /* 9)\ewald_kmax{#} class_parse ??? */
  classParse->kmax_ewd = 10;
  /* 10)\ewald_respa_kmax{#} class_parse ??? */
  classParse->kmax_res = 0;
  /* 11)\ewald_pme_opt{#} */
  classMini->part_mesh.pme_on = class->part_mesh.pme_on;
  /* 12)\ewald_kmax_pme{#} ??? */
  classMini->part_mesh.kmax_pme = class->part_mesh.kmax_pme;
  /* 14)\ewald_respa_pme_opt{#} I dont need this */
  classMini->part_mesh.pme_res_on = class->part_mesh.pme_res_on;
  /* 15)\ewald_respa_kmax_pme{#} I dont need this */
  classMini->part_mesh.kmax_pme_res = class->part_mesh.kmax_pme_res;
  /* 16)\ewald_respa_interp_pme{#} I dont need this */
  classMini->part_mesh.n_interp_res = class->part_mesh.n_interp_res;
  /* 17)\sep_VanderWaals{on,off} I dont need this */
  classMini->energy_ctrl.isep_vvdw = class->energy_ctrl.isep_vvdw;
  /* 18)\dielectric_opt{on,off} I dont need this */
  classMini->interact.dielectric_opt = class->interact.dielectric_opt;
  /* 19)\dielectric_rheal{#} I dont need this */
  classMini->interact.dielectric_rheal = class->interact.dielectric_rheal;
  /* 20)\dielectric_cut{#} I dont need this */
  classMini->interact.dielectric_cut = class->interact.dielectric_cut;
  /* 21)\dielectric_eps{#} I dont need this */
  classMini->interact.dielectric_eps = class->interact.dielectric_eps;
  /* 22)\std_intra_block{on,off} I dont need this */
  classMini->energy_ctrl.block_std_on = class->energy_ctrl.block_std_on;
  /* 23)\con_intra_block{on,off} I dont need this */
  classMini->energy_ctrl.block_con_on;
  /* 24)\inter_PE_calc_freq */
  generalDataMini->timeinfo.iget_pe_real_inter_freq = general_data->timeinfo.iget_pe_real_inter_freq;
  /*  25)\pme_paralell_opt{#} serial */
  classMini->part_mesh.pme_para_opt = 0;

/*=======================================================================*/
/*  V) set_sim_params_run                                                */
 
  /* 1)\init_resmpl_atm_vel{on,off} class_parse */
  classParse->ivx_smpl = 0;
  /* 2)\init_resmpl_cp_vel{on,off} cp_parse */
  cpParse->ivc_smpl = 0;
  /* 3)\init_resmpl_cp_nhc{on,off} cp_parse */
  cpParse->ivcnhc_smpl = 0;
  /* 4a)\resmpl_frq_atm_vel{#} */
  classMini->vel_samp_class.nvx_smpl = class->vel_samp_class.nvx_smpl;
  /* 4b)\rescale_frq_atm_vel{#} */
  classMini->vel_samp_class.nvx_scale = class->vel_samp_class.nvx_scale;
  /* 5)\respa_steps_lrf{#} */
  generalDataMini->timeinfo.nres_ter = general_data->timeinfo.nres_ter;
  generalDataMini->timeinfo.int_res_ter = general_data->timeinfo.int_res_ter;
  classMini->energy_ctrl.int_res_ter = class->energy_ctrl.int_res_ter;
  /* 6)\respa_steps_torsion{#} */
  generalDataMini->timeinfo.nres_tor = general_data->timeinfo.nres_tor;
  generalDataMini->timeinfo.int_res_tor = general_data->timeinfo.int_res_tor;
  /* 7)\respa_steps_intra{#} */
  generalDataMini->timeinfo.nres_tra = general_data->timeinfo.nres_tra;
  generalDataMini->timeinfo.int_res_tra = general_data->timeinfo.int_res_tra;
  classMini->energy_ctrl.int_res_tra = class->energy_ctrl.int_res_tra;
  /* 8)\respa_rheal{#} */
  classMini->interact.rheal_res = class->interact.rheal_res;
  /* 9)\shake_tol{#} */
  bondedMini->constrnt.tolshake = bonded->constrnt.tolshake;
  /* 10)\rattle_tol{#} */
  bondedMini->constrnt.tolratl = bonded->constrnt.tolratl;
  /* 11)\max_constrnt_iter{#} */
  bondedMini->constrnt.max_iter = bonded->constrnt.max_iter;
  /* 12)\init_rescale_atm_vel{on,off} class_parse */
  classParse->ivx_scale = 0;
  /* 13)\init_rescale_atm_nhc{on,off} class_parse */
  classParse->ivnhc_scale = 0;
  /* 14)\init_rescale_cp_vel{on,off} cp_parse */
  cpParse->ivc_scale = 0;
  /* 15)\resmpl_frq_cp_vel{ # } */
  cpMini->vel_samp_cp.nvc_smpl = cp->vel_samp_cp.nvc_smpl;
  /* 16)\group_con_tol{#} */
  bondedMini->grp_bond_con.tol = bonded->grp_bond_con.tol;
  /* 17)\cp_norb_tol{#} */
  cpMini->cpconstrnt.c_tolnorb = cp->cpconstrnt.c_tolnorb;
  /* 18)\cp_ks_rot_frq{#} */
  cpMini->cpcoeffs_info.ks_rot_on = cp->cpcoeffs_info.ks_rot_on;
  cpMini->cpcoeffs_info.n_ks_rot = cp->cpcoeffs_info.n_ks_rot;
  /* 19)\cp_shake_tol{#} */
  cpMini->cpconstrnt.c_tolshake = cp->cpconstrnt.c_tolshake;
  /* 20)\cp_rattle_tol{#} */
  cpMini->cpconstrnt.c_tolratl = cp->cpconstrnt.c_tolratl;
  /* 24)\cp_min_tol{#} */
  generalDataMin->minopts.tol_coef = general_data->minopts.tol_coef;
  /* 21)\cp_run_tol{#} */
  cpMini->cpopts.tol_coef = cp->cpopts.tol_coef;
  /* 22)\zero_com_vel{yes,no}*/
  classParse->zero_com_vel = 0;
  /* 23)\min_tol{#} */
  generalDataMini->minopts.tol_atom = general_data->minopts.tol_atom;
  /* 26)hess_opt{full_an,full_num,unit} */
  classMini->clatoms_info.hess_calc = 0;
  /* 27)\class_mass_scale_fact{#} */
  classMini->clatoms_info.mass_sc_fact = class->clatoms_info.mass_sc_fact;
  /* 28)\hmat_int_typ{normal,upper_triangle} */
  generalDataMini->cell.hmat_int_typ = general_data->cell.hmat_int_typ;
  /* 29)\hmat_cons_typ{none,ortho_rhom,mono_clin} */
  generalDataMini->cell.hmat_cons_typ = general_data->cell.hmat_cons_typ;
  /* 30)\hmat_cons_typ{none,ortho_rhom,mono_clin} */
  generalDataMini->minopts.min_atm_com_fix_opt = general_data->minopts.min_atm_com_fix_opt;
  /* 31)\rescale_frq_cp_vel{ # } */
  cpMini->vel_samp_cp.nvc_scal = cp->vel_samp_cp.nvc_scal;
  /* 32)\auto_rescale_cp_vel{ on/off } */
  cpMini->vel_samp_cp.iauto_vc_scal_opt = cp->vel_samp_cp.iauto_vc_scal_opt;
  /* 33)\auto_rescale_cp_vel_tol{ # } */
  cpMini->vel_samp_cp.vc_scal_tol = cp->vel_samp_cp.vc_scal_tol;
  /* 34)\wan_min_tol{#} */
  generalDataMini->minopts.tol_wan_coef = general_data->minopts.tol_wan_coef;

/*=======================================================================*/
/*  VI) set_sim_params_nhc                                               */

/*=======================================================================*/
/*  VII) set_sim_params_vol                                              */
  
  /* 1)\volume_tau{#} class_parse */
  classParse->tau_vol = 1000.0;
  /* 2)\volume_nhc_tau{#} class_parse */
  classParse->tau_vol_nhc = 1000.0
  /* 3)\periodicity{0,1,2,3} */
  generalDataMini->cell.iperd = general_data->cell.iperd;
  /* 4)\intra_perds{on,off} */
  generalDataMini->cell.intra_perds = general_data->cell.intra_perds;

/*=======================================================================*/
/*  VIII) set_sim_params_write                                           */

  /* 25)\mol_set_file */
  fileNameParse->molsetname = (char *)cmalloc(MAXWORD*sizeof(char));
  strcpy(filenameParse->molsetname,fragInfo->molSetName);

/*=======================================================================*/
/*  IX) set_sim_params_pimd                                              */

/*=======================================================================*/
/*  X) set_sim_params_velo                                               */
  
/*=======================================================================*/
/*  XI) set_sim_params_msqd                                              */

/*=======================================================================*/
/*  XII) set_sim_params_iikt_iso                                         */

/*=======================================================================*/
/*  XIII) set_sim_params_ickt_iso                                        */

/*=======================================================================*/
/*  XIV) set_sim_params_rdf                                              */

/*=======================================================================*/
/*  XV) set_sim_params_harmonic                                          */

/*=======================================================================*/
/*  XVI) set_sim_params_dafed                                            */

/*=======================================================================*/
/*  XVII) set_sim_params_stodft                                          */

/*=======================================================================*/
/*  XVII) set_sim_params_stodft                                          */




/*------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/




