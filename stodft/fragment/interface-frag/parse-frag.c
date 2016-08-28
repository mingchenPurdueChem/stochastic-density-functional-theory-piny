/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: parse-frag.c                                 */
/*                                                                          */
/* This subprogram generate all inputs for fragmentation                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_parse_entry.h"
#include "../proto_defs/proto_parse_local.h"
#include "../proto_defs/proto_sim_params_entry.h"
#include "../proto_defs/proto_mol_params_entry.h"
#include "../proto_defs/proto_intra_params_entry.h"
#include "../proto_defs/proto_coords_entry.h"
#include "../proto_defs/proto_coords_local.h"
#include "../proto_defs/proto_cp_ewald_entry.h"
#include "../proto_defs/proto_surf_params_entry.h"
#include "../proto_defs/proto_inter_params_entry.h"
#include "../proto_defs/proto_inter_params_local.h"
#include "../proto_defs/proto_vps_params_entry.h"
#include "../proto_defs/proto_lists_entry.h"
#include "../proto_defs/proto_scratch_entry.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_vel_sampl_class_entry.h"
#include "../proto_defs/proto_vel_sampl_cp_entry.h"
#include "../proto_defs/proto_vel_sampl_cp_local.h"
#include "../proto_defs/proto_coords_cp_entry.h"
#include "../proto_defs/proto_coords_cp_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_entry.h"
#include "../proto_defs/proto_communicate_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_integrate_cpmin_entry.h"
#include "../proto_defs/proto_energy_cp_entry.h"
#include "../proto_defs/proto_interface_frag_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Parse: note there is a noncommuting order of calls                      */
/*==========================================================================*/

void parseFrag(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp,
           ANALYSIS *analysis,CLASS *classMini,BONDED *bondedMini,
	   GENERAL_DATA *generalDataMini,CP *cpMini,ANALYSIS *analysisMini)

/*========================================================================*/
/*             Begin Routine                                              */
         {/*Begin subprogram: */
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int iii,cp_on,pimd_on,cp_md,i,pi_beads_proc_st;
  int ivx_flag,ivnhc_flag,ipos_flag;
  int nchrg,iperd,myid,np_states,np_forc,pi_beads,pme_on;
  int myid_state,myid_forc,myatm_start,myatm_end;
  int cp_dual_grid_opt_on;        /*dualed option flag for CP */
  int ncons,ncons_dn;

  CLASS_PARSE      classParse;
  CP_PARSE         cpParse;
  FILENAME_PARSE   fileNameParse;
  FREE_PARSE       freeParse;
  SPLINE_PARSE     splineParse;
  NULL_INTER_PARSE nullInterParse;

  int icontrol_proc  = general_data->error_check_on;
  int num_proc       = class->communicate.np;
  MPI_Comm world     = class->communicate.world;
  double *tot_memory = &(class->tot_memory);

  int iopt_cp_pw,iopt_cp_dvr;

/*========================================================================*/
/*   I) Zero the malloc size variables                                  */

  control_zero_mem(classMini,bondedMini,generalDataMini,cpMini,&classParse,
                   &nullInterParse);
  

/*========================================================================*/
/*   II) Set the sim parameters: Done first to get input file names       */
/*               (interface/sim_params/control_sim_params.c)              */

  //filename_parse.input_name  = (char *)cmalloc(MAXWORD*sizeof(char));
  //strcpy(filename_parse.input_name,input_name);

  copySimParam(general_data,bonded,class,cp,generalDataMini,bondedMini,
		classMini,cpMini,&classParse,&cpParse,&fileNameParse);
   
  /*
  controlSimParamsFrag(class,general_data,bonded,cp,analysis,classMini,generalDataMini,
		       bondedMini,cpMini,analysisMini,
		       &classParse,&cpParse);
  */

  
/*========================================================================*/
/*   III) Read in atom and CP parameters: Done second to get info         */
/*                                        needed for set_intra;           */
/*                                        some atom mallocing             */
/*               (interface/mol_params/control_mol_params.c)              */

  controlMolParamsFrag(class,general_data,bonded,cp,classMini,
		       generalDataMini,bondedMini,cpMini,&classParse,
                       &cpParse,&freeParse,&fileNameParse);
  
/*========================================================================*/
/*  IV) Read in atom, molecule connectivity data: Done before setting     */
/*                                                 therms;                */
/*                                                 majority atom mallocing*/
/*                                                 intramol mallocing     */
/*                                                 pressure mallocing     */
/*               (interface/intra_params/control_intra_params.c)          */

  
  control_intra_params(tot_memory,
		       &(classMini->clatoms_info),(classMini->clatoms_pos),
		       &(classMini->ghost_atoms),&(classMini->atommaps),
		       bondedMini,&fileNameParse,&freeParse,
		       &classParse,&null_inter_parse,
		       &(generalDataMini->simopts),&(classMini->communicate),
		       (classMini->surface.isurf_on));
  

/*========================================================================*/
/*    V) Communicate class interface: done before proceeding further      */

  /*
  if(num_proc>1){
    Barrier(world);
    communicate_interface(class,bonded,cp,general_data,&null_inter_parse,
                          &class_parse,&cp_parse);
  }//endif
  */

/*========================================================================*/
/*   VI) Assign Flags                                                     */

  
  cp_dual_grid_opt_on = cpMini->cpopts.cp_dual_grid_opt;

  pimd_on = generalDataMini->simopts.pimd
	   +generalDataMini->simopts.cp_pimd 
	   +generalDataMini->simopts.cp_wave_pimd
	   +generalDataMini->simopts.cp_wave_min_pimd 
	   +generalDataMini->simopts.debug_pimd
	   +generalDataMini->simopts.debug_cp_pimd;

  cp_on = generalDataMini->simopts.cp_min
	 +generalDataMini->simopts.cp_wave_min
         +generalDataMini->simopts.cp
	 +generalDataMini->simopts.cp_wave
         +generalDataMini->simopts.cp_pimd
	 +generalDataMini->simopts.cp_wave_pimd
         +generalDataMini->simopts.debug_cp
	 +generalDataMini->simopts.debug_cp_pimd
         +generalDataMini->simopts.cp_wave_min_pimd;


  cp_md = generalDataMini->simopts.cp   
         +generalDataMini->simopts.cp_wave
         +generalDataMini->simopts.debug_cp
         +generalDataMini->simopts.cp_pimd
         +generalDataMini->simopts.cp_wave_pimd;

  nchrg       = classMini->clatoms_info.nchrg;
  iperd       = generalDataMini->cell.iperd;
  myid        = classMini->communicate.myid;
  np_states   = classMini->communicate.np_states;
  np_forc     = classMini->communicate.np_forc;
  pi_beads    = classMini->clatoms_info.pi_beads;
  pme_on      = classMini->part_mesh.pme_on;
  iopt_cp_pw  = cpMini->cpcoeffs_info.iopt_cp_pw;
  iopt_cp_dvr = cpMini->cpcoeffs_info.iopt_cp_dvr;
  

/*========================================================================*/
/*    VII) Read in hmat. Do before set_cp_ewald                           */
/*                (interface/coords/read_coord.c)                         */
/*  Pass all atom positions into */

  /* Keep
  mall_coord(classMini,generalDataMini);
  mall_pressure(classMini,generalDataMini);  
  if(myid==0){//change
    read_hmat(class,general_data,&filename_parse,class_parse.istart,
              cp_dual_grid_opt_on,&(cp->cpewald.dbox_rat),
              &(cp->cpewald.box_rat));
  }//endif
  */
  /*
  if(num_proc>1){
    comm_cell_data(&(general_data->cell),
                   &(cp->cpewald.dbox_rat),
                   &(cp->cpewald.box_rat),world);

  }//endif
  */

/*========================================================================*/
/* VIII) Set up the ewald/cp: Done before setting intermol PE             */
/*                            CP/Ewald mallocing                          */
/*                (interface/cp_ewald/control_set_cp_ewald                */

/*--------------------------------------------------------------------------*/
  //Keep if((nchrg>0&&iperd>0)||cp_on==1){
/*--------------------------------------------------------------------------*/
/* Set up CP and Ewald stuff                                                */

    /* Keep
    if(iopt_cp_pw == 1){//change
      control_set_cp_ewald(&(general_data->simopts),&(general_data->cell),
                           &(cp->cpcoeffs_info),&(general_data->ewald),
                           &(cp->cpewald),&cp_parse,
                           &(cp->pseudo.gmin_true),
                           &(cp->pseudo.gmin_spl),
                           &(cp->pseudo.gmax_spl),
                           &(class->ewd_scr),(class_parse.kmax_ewd),
                           (class_parse.kmax_res),
                           tot_memory,general_data->timeinfo.int_res_ter,
                           &(class->part_mesh),&(bonded->ecor),myid,
                           cp->cpopts.cp_lsda,
                           general_data->minopts.cp_min_diis,
                           cp_dual_grid_opt_on); 
    }
    classMini->clatoms_info.alp_ewd = generalDataMini->ewald.alp_ewd;
    */
/*--------------------------------------------------------------------------*/
/*  Calculate Number of CP fictitious degrees of freedom                    */

    /* Keep
    if(cp_on){
      if(iopt_cp_pw==1){calculate_cp_nfree(cpMini);}
      cpMini->cpopts.te_ext /= (double)(cpMini->cpcoeffs_info.cp_nfree);
      cpMini->vel_samp_cp.div_scal = (double)(cpMini->cpcoeffs_info.cp_nfree);
    }//endif cp_on 
    */
/*--------------------------------------------------------------------------*/
  //Keep }/*endif nchrg > 0 or cp is on */
/*========================================================================*/
/*   IX) Build communication groups: must be done in serial and parallel  */
/*                                    after call to set_cp_ewald          */

  /*
  control_group_communicators(classMini,cpMini,cp_on); //redo
  //if(num_proc>1){ Barrier(world); }
  myid_state  = classMini->communicate.myid_state;
  myid_forc   = classMini->communicate.myid_forc;
  */

/*========================================================================*/
/*   IX) Create the FFT packages  */
  /* Keep
  //change
  if((nchrg>0&&iperd>0&&pme_on==1)||cp_on==1){
    if(myid_state<num_proc&&myid_state>=0){
      control_fft_pkg(&(cp->cp_sclr_fft_pkg3d_sm),&(cp->cp_para_fft_pkg3d_sm),
                      &(cp->cp_sclr_fft_pkg3d_dens_cp_box),&(cp->cp_para_fft_pkg3d_dens_cp_box),
                      &(cp->cp_sclr_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_lg),
                      &(general_data->pme_fft_pkg),&(general_data->pme_res_fft_pkg),
                      &(general_data->ewald),&(cp->cpewald),&(class->part_mesh),&(cp->cpcoeffs_info),
                      &(class->communicate),cp_on,cp->cpopts.cp_lsda,
		      tot_memory,general_data->timeinfo.int_res_ter,
                      cp->cpopts.cp_para_opt,cp_dual_grid_opt_on);
    }//endif
  }//endif
  */

/*========================================================================*/
/*   XI) Setup intermolecular potential stuff: interspline mallocing      */
/*                (interface/inter_params/control_inter_params.c)         */
  /* Keep
  //change
  control_inter_params(&(class->interact),&spline_parse,
                       &filename_parse,general_data->ewald.alp_ewd,
                       nchrg,class->clatoms_info.natm_tot,
                       class->atommaps.natm_typ,class->atommaps.atm_typ,
                       class->atommaps.iatm_atm_typ,iperd,
                       class_parse.ishift_pot,tot_memory,
                       general_data->timeinfo.int_res_ter,
                       myid,world,num_proc);
  */

/*========================================================================*/
/*    XII) Setup the surface potential if needed                          */

/*========================================================================*/
/*    XII) Setup pseudopotential stuff: pseudospline mallocing            */
/*                (interface/interparams/control_vps_params.c)            */

  //Keep if(cp_on==1){//change
  /*---------------------------------------------------------*/
  /* Create a list of ab initio atoms                        */
  /* List is used in gen_wave and also cp_dual_check routine */

    /* Keep
    make_cp_atom_list(&(class->atommaps),
                       class->clatoms_info.cp_atm_flag,
                       &(class->clatoms_info.nab_initio),
                       class->clatoms_info.natm_tot);
    if(myid_state<np_states){
      control_vps_params(&(cp->pseudo),&(general_data->cell),&filename_parse,
                         &spline_parse,class->atommaps.natm_typ,
                         class->atommaps.atm_typ,
                         tot_memory,class->clatoms_info.natm_tot,
                         class->clatoms_info.nab_initio,
                         cp->cpopts.cp_ptens_calc,cp_dual_grid_opt_on,
                         &(class->communicate),cp_parse.cp_ecut,
                         &(cp->cpcoeffs_info));
    }//endif
  }//endif
  */


/*========================================================================*/
/*    XIII) Set particle exclusions and ewald corrections                 */
/*                (interface/lists/set_exclude.c)                         */

  /*THIS HAS ALSO BEEN MODIFIED FOR CLASSICAL HCA EQUILIBRATION */
  /* Keep
  set_exclude(&(class->clatoms_info),&(class->ghost_atoms),bonded,
              &(bonded->excl),&null_inter_parse,
              iperd,tot_memory,
              general_data->ewald.alp_ewd,icontrol_proc);
  */

/*========================================================================*/
/*   XV) Set thermostats: Done before reading the coordinates;           */
/*                        atm NHC mallocing                              */
/*                (interface/coords/set_atm_NHC.c                        */

/*========================================================================*/
/*  XVI) Read in atm positions/velocities/NHCs:                           */
/*                (interface/coords/read_coord.c)                         */

   /* Keep
   //change
   read_coord(class,general_data,&filename_parse,
              class_parse.istart,cp_dual_grid_opt_on);
   */


/*========================================================================*/
/*  XVII) Spline the ewald corrections (needs particle positions)         */

  /* Keep
  if((nchrg > 0 && iperd > 0)){//change
    splin_ecor(&(bonded->ecor),&(general_data->ewald),(class->clatoms_pos),
               pi_beads,icontrol_proc,tot_memory);
  }else{
    general_data->ewald.self_erf = 1.0;
  }//endif
  */

/*========================================================================*/
/* XVIII) Set thermostats: Done before reading the coeffs                 */
/*                        CP NHC mallocing                                */
/*                (interface/coords/set_coef_NHC.c)                       */

  /* Keep
  if(cp_on==1){
   if(myid_state<np_states){//change
     mall_coef(cp,&(general_data->simopts),class->clatoms_info.pi_beads_proc);
     cp->cptherm_info.num_c_nhc      = 0;
     cp->cptherm_info.num_c_nhc_proc = 0;
     cp->cptherm_info.num_c_nhc_norm = 0;
     cp->cptherm_info.massiv_flag    = 0;
   }//endif
  }//endif
  */

/*========================================================================*/
/*   XIX)malloc scratch space                                            */
/*                (interface/scratch/mall_scratch.c)                     */

  //Keep control_mall_scratch(classMini,bondedMini,cpMini,generalDataMini);

/*========================================================================*/
/* XX) Read in coeffs/velocities:                                         */
/*                (interface/coords/read_coef.c)                          */
/*     And tidy up the dual option                                        */

  /* Keep
  if(cp_on==1){
    if(myid_state<np_states){
      read_coef(cp,general_data,class,&filename_parse,&cp_parse,tot_memory);
      if(myid == 0){cfree(&(filename_parse.vps_name[1]));} 
      if(cp->cpcoeffs_info.cp_elf_calc_frq >0 || 
        cp->cpcoeffs_info.cp_dip_calc_frq>0){
        mall_properties(cp);
      }
      if(cp->cpscr.cpscr_wannier.cp_wannier_on==1 ||
        cp->cpcoeffs_info.cp_dip_calc_frq > 0){
        mall_wannier(cp);
      }
      if(cp->cpopts.cp_nloc_wan_opt==1){
        read_wan_cent(cp,general_data);
      }
    }//endif np_state
  }//endif cp_on
  */
/*========================================================================*/
/* XXI) Set up branch root neighbor list data                             */

  /*Keep
  if(class->nbr_list.brnch_root_list_opt>0){//do I need this?
    control_brnch_root_list(class,bonded);
  }//endif
  */

/*========================================================================*/
/* XXI) Control Molecular Decomposition       */

  //Keep control_molec_decomp(class,bonded,general_data);//do I need this?

/*========================================================================*/
/*   X) Initialize path integral transformations: after group communicators*/

/*========================================================================*/
/*  XVII) Communicate classical coordinates to non-bead processors :      */
/*                                         after path_integral_init     */

/*========================================================================*/
/*   XXII) malloc neigbor list memory                                     */
/*                (interface/lists/mall_make_lists.c)                     */
 
  /* Keep
  //do I need this? 
  get_cut_skin(class->interact.cutskin,
               class->interact.cutskin_root,
               class->interact.cutoff,
               class->interact.cutskin_res,
               class->interact.cutskin_root_res,
               class->interact.cutoff_res,
               class->interact.skin,
               class->interact.spread,
               class->interact.brnch_root_skin,
               class->interact.nter_typ);

  pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;
  myatm_start = class->clatoms_info.myatm_start;
  myatm_end = class->clatoms_info.myatm_end;

  if(( (myid_state==0) || (np_forc==np_states!= 1)==0  )&&pimd_on==1){ //DY
    control_pimd_trans_mode(class,general_data);

    control_pimd_trans_pos(class,general_data);

  }//endif

  mall_make_lists(class,general_data,bonded,icontrol_proc);
  */

/*========================================================================*/
/* XXIII) Orthogonalize coefficients                                      */
  
  /*
  if(myid_state<np_states){
    control_init_cp_orthog(generalDataMini,cpMini,&cp_parse,cp_on,cp_md,myid);
  }//endif
  */

/*========================================================================*/
/*  XXIV) Assign/resample the initial class velocities                    */
/*                (interface/vel_sampl)                                   */

/*========================================================================*/
/*  XXV) Assign/resample the initial coef velocities                      */
/*                (interface/vel_sampl)                                   */

/*========================================================================*/
/*  XXV) Initialize UFED/dafed calculation                                */
/*                (dafed/)                                                */

/*========================================================================*/
/* XXVII) Flush the buffers                                               */
  
  fflush(stdout);
  fflush(stderr);

/*------------------------------------------------------------------------*/
}/*end routine*/ 
/*==========================================================================*/







