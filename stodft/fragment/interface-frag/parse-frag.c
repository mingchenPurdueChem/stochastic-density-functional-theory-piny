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
  int realSparseOpt;

  CLASS_PARSE      classParse;
  CP_PARSE         cpParse;
  FILENAME_PARSE   fileNameParse;
  FREE_PARSE       freeParse;
  SPLINE_PARSE     splineParse;
  NULL_INTER_PARSE nullInterParse;

  int icontrol_proc;
  int num_proc;
  MPI_Comm world     = classMini->communicate.world;

  int iopt_cp_pw,iopt_cp_dvr;
  double totMemFake;
  double *tot_memory = &totMemFake;
  *tot_memory = 0.0;

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

  icontrol_proc = generalDataMini->error_check_on;
  num_proc = classMini->communicate.np;
  
   
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


  controlIntraParamsFrag(tot_memory,classMini,generalDataMini,bondedMini,
			&fileNameParse,&freeParse,&classParse,&nullInterParse);

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
  
  //cpMini->cpopts.fftw3dFlag = 1;
  //cpMini->cpcoeffs_info.fftw3dFlag = 1;
  //cpMini->cpewald.fftw3dFlag = 1;
  // For non-cubic box, the fft in this package doesn't use all grid. Thus it's hard to 
  // predict fftw output
  cpMini->cpopts.fftw3dFlag = 1;
  cpMini->cpcoeffs_info.fftw3dFlag = 1;
  cpMini->cpewald.fftw3dFlag = 1;

  classMini->clatoms_info.ifirst_vps = 0;
  cpMini->cpcoeffs_info.itime_ks = 0;
  cpMini->cpewald.realSparseOpt = cpMini->cpopts.realSparseOpt;

/*========================================================================*/
/*    VII) Read in hmat. Do before set_cp_ewald                           */
/*                (interface/coords/read_coord.c)                         */
/*  Pass all atom positions into */

  mall_coord(classMini,generalDataMini);
  mall_pressure_frag(classMini,generalDataMini);  

  // Do this first
  initCoordHmatFFT(general_data,class,cp,generalDataMini,classMini,cpMini);

  if(myid==0){//change
    readHmatFrag(classMini,generalDataMini,cpMini,cp_dual_grid_opt_on,
		&(cpMini->cpewald.dbox_rat),&(cpMini->cpewald.box_rat));
  }//endif
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

  realSparseOpt = cpMini->cpewald.realSparseOpt;

/*--------------------------------------------------------------------------*/
  if((nchrg>0&&iperd>0)||cp_on==1){
/*--------------------------------------------------------------------------*/
/* Set up CP and Ewald stuff                                                */
    if(iopt_cp_pw == 1){//change
      if(realSparseOpt==0){
        controlSetCpEwaldFrag(generalDataMini,classMini,cpMini,bondedMini,
	  		      cp,class,general_data,bonded,&cpParse,&classParse);
      }
      else{
        controlSetCpEwaldFragSparse(generalDataMini,classMini,cpMini,bondedMini,
                              cp,class,general_data,bonded,&cpParse,&classParse);
      }
    }
    classMini->clatoms_info.alp_ewd = generalDataMini->ewald.alp_ewd;
/*--------------------------------------------------------------------------*/
/*  Calculate Number of CP fictitious degrees of freedom                    */

    if(cp_on){
      if(iopt_cp_pw==1){calculate_cp_nfree(cpMini);}
      cpMini->cpopts.te_ext /= (double)(cpMini->cpcoeffs_info.cp_nfree);
      cpMini->vel_samp_cp.div_scal = (double)(cpMini->cpcoeffs_info.cp_nfree);
    }//endif cp_on 

/*--------------------------------------------------------------------------*/
  }/*endif nchrg > 0 or cp is on */
/*========================================================================*/
/*   IX) Build communication groups: must be done in serial and parallel  */
/*                                    after call to set_cp_ewald          */

  controlGroupCommunicatorsFrag(classMini,cpMini,cp_on);
  //if(num_proc>1){ Barrier(world); }
  myid_state  = classMini->communicate.myid_state;
  myid_forc   = classMini->communicate.myid_forc;
  

/*========================================================================*/
/*   IX) Create the FFT packages  */
  //change
  //printf("cp_on %i myid_state %i num_proc %i\n",cp_on,myid_state,num_proc); 

  cpMini->cp_sclr_fft_pkg3d_sm.threadFlag = cpMini->cpopts.threadFlag;
  cpMini->cp_para_fft_pkg3d_sm.threadFlag = cpMini->cpopts.threadFlag;
  cpMini->cp_sclr_fft_pkg3d_dens_cp_box.threadFlag = cpMini->cpopts.threadFlag;
  cpMini->cp_para_fft_pkg3d_dens_cp_box.threadFlag = cpMini->cpopts.threadFlag;
  cpMini->cp_sclr_fft_pkg3d_lg.threadFlag = cpMini->cpopts.threadFlag;
  cpMini->cp_para_fft_pkg3d_lg.threadFlag = cpMini->cpopts.threadFlag;

  cpMini->cp_sclr_fft_pkg3d_sm.numThreads = cpMini->communicate.numThreads;
  cpMini->cp_para_fft_pkg3d_sm.numThreads = cpMini->communicate.numThreads;
  cpMini->cp_sclr_fft_pkg3d_dens_cp_box.numThreads = cpMini->communicate.numThreads;
  cpMini->cp_para_fft_pkg3d_dens_cp_box.numThreads = cpMini->communicate.numThreads;
  cpMini->cp_sclr_fft_pkg3d_lg.numThreads = cpMini->communicate.numThreads;
  cpMini->cp_para_fft_pkg3d_lg.numThreads = cpMini->communicate.numThreads;

  if((nchrg>0&&iperd>0&&pme_on==1)||cp_on==1){
    if(myid_state<num_proc&&myid_state>=0){
      controlFFTPkgFrag(generalDataMini,classMini,cpMini,cp);
    }//endif
  }//endif

/*========================================================================*/
/*   XI) Setup intermolecular potential stuff: interspline mallocing      */
/*                (interface/inter_params/control_inter_params.c)         */
  //change
  controlInterParamsFrag(generalDataMini,classMini,cpMini,bondedMini,cp,
			  &splineParse,&fileNameParse,&classParse);

/*========================================================================*/
/*    XII) Setup the surface potential if needed                          */

/*========================================================================*/
/*    XII) Setup pseudopotential stuff: pseudospline mallocing            */
/*                (interface/interparams/control_vps_params.c)            */

  if(cp_on==1){//change
  /*---------------------------------------------------------*/
  /* Create a list of ab initio atoms                        */
  /* List is used in gen_wave and also cp_dual_check routine */

    make_cp_atom_list(&(classMini->atommaps),
                       classMini->clatoms_info.cp_atm_flag,
                       &(classMini->clatoms_info.nab_initio),
                       classMini->clatoms_info.natm_tot);
    if(myid_state<np_states){
      controlVpsParamsFrag(generalDataMini,classMini,cpMini,&fileNameParse,
			    &splineParse,&cpParse);
    }//endif
    if(cpMini->pseudo.pseudoReal.pseudoRealFlag==1){
      copyNlppReal(general_data,bonded,class,cp,generalDataMini,bondedMini,
                   classMini,cpMini);
    }
  }//endif


/*========================================================================*/
/*    XIII) Set particle exclusions and ewald corrections                 */
/*                (interface/lists/set_exclude.c)                         */

  /*THIS HAS ALSO BEEN MODIFIED FOR CLASSICAL HCA EQUILIBRATION */
  set_exclude(&(classMini->clatoms_info),&(classMini->ghost_atoms),bondedMini,
              &(bondedMini->excl),&nullInterParse,
              iperd,tot_memory,
              generalDataMini->ewald.alp_ewd,/*icontrol_proc*/0);

/*========================================================================*/
/*   XV) Set thermostats: Done before reading the coordinates;           */
/*                        atm NHC mallocing                              */
/*                (interface/coords/set_atm_NHC.c                        */

/*========================================================================*/
/*  XVI) Read in atm positions/velocities/NHCs:                           */
/*                (interface/coords/read_coord.c)                         */

   //change
   /*
   read_coord(class,general_data,&filename_parse,
              class_parse.istart,cp_dual_grid_opt_on);
   */

/*========================================================================*/
/*  XVII) Spline the ewald corrections (needs particle positions)         */

  if((nchrg > 0 && iperd > 0)){//change
    splin_ecor(&(bondedMini->ecor),&(generalDataMini->ewald),(classMini->clatoms_pos),
               pi_beads,/*icontrol_proc*/0,tot_memory);
  }else{
    generalDataMini->ewald.self_erf = 1.0;
  }//endif

/*========================================================================*/
/* XVIII) Set thermostats: Done before reading the coeffs                 */
/*                        CP NHC mallocing                                */
/*                (interface/coords/set_coef_NHC.c)                       */

  if(cp_on==1){
   if(myid_state<np_states){//change
     mall_coef(cpMini,&(generalDataMini->simopts),classMini->clatoms_info.pi_beads_proc);
     cpMini->cptherm_info.num_c_nhc      = 0;
     cpMini->cptherm_info.num_c_nhc_proc = 0;
     cpMini->cptherm_info.num_c_nhc_norm = 0;
     cpMini->cptherm_info.massiv_flag    = 0;
   }//endif
  }//endif

/*========================================================================*/
/*   XIX)malloc scratch space                                            */
/*                (interface/scratch/mall_scratch.c)                     */

  //debug
  //cpMini->cpcoeffs_info.nstate_up = 16;
  //cpMini->cpcoeffs_info.nstate_dn = 16;
  controlMallScratchFrag(classMini,bondedMini,cpMini,generalDataMini);
  //cpMini->cpcoeffs_info.nstate_up = 4;
  //cpMini->cpcoeffs_info.nstate_dn = 4;


/*========================================================================*/
/* XX) Read in coeffs/velocities:                                         */
/*                (interface/coords/read_coef.c)                          */
/*     And tidy up the dual option                                        */

  if(cp_on==1){
    if(myid_state<np_states){
      readCoefFrag(cpMini,generalDataMini,classMini,&fileNameParse,&cpParse,tot_memory);

      printf("aaaaaaaaaa ccreal %lg\n",cpMini->cpcoeffs_pos[1].cre_up[448317]);
      cpMini->pseudo.pseudoReal.forceCalcFlag = 1;
      if(cpMini->pseudo.pseudoReal.pseudoRealFlag==1){
        initRealNlppWf(cpMini,classMini,generalDataMini);
      }
      cpMini->pseudo.pseudoReal.forceCalcFlag = 0;
      if(myid == 0){cfree(&(fileNameParse.vps_name[1]));} 
    }//endif np_state
  }//endif cp_on

  printf("bbbbbbbbbbbbbb ccreal %lg\n",cpMini->cpcoeffs_pos[1].cre_up[448317]);

/*========================================================================*/
/* XXI) Set up branch root neighbor list data                             */

  if(classMini->nbr_list.brnch_root_list_opt>0){//do I need this?
    control_brnch_root_list(classMini,bondedMini);
  }//endif

  printf("ccccccccccccc ccreal %lg\n",cpMini->cpcoeffs_pos[1].cre_up[448317]);

/*========================================================================*/
/* XXI) Control Molecular Decomposition       */

   control_molec_decomp(class,bonded,general_data);//do I need this?

  printf("dddddddddddd ccreal %lg\n",cpMini->cpcoeffs_pos[1].cre_up[448317]);

/*========================================================================*/
/*   X) Initialize path integral transformations: after group communicators*/

/*========================================================================*/
/*  XVII) Communicate classical coordinates to non-bead processors :      */
/*                                         after path_integral_init     */

/*========================================================================*/
/*   XXII) malloc neigbor list memory                                     */
/*                (interface/lists/mall_make_lists.c)                     */
 
  get_cut_skin(classMini->interact.cutskin,
               classMini->interact.cutskin_root,
               classMini->interact.cutoff,
               classMini->interact.cutskin_res,
               classMini->interact.cutskin_root_res,
               classMini->interact.cutoff_res,
               classMini->interact.skin,
               classMini->interact.spread,
               classMini->interact.brnch_root_skin,
               classMini->interact.nter_typ);

  pi_beads_proc_st = classMini->clatoms_info.pi_beads_proc_st;
  myatm_start = classMini->clatoms_info.myatm_start;
  myatm_end = classMini->clatoms_info.myatm_end;

  mallMakeListsFrag(classMini,generalDataMini,bondedMini,icontrol_proc);

  printf("eeeeeeeeee ccreal %lg\n",cpMini->cpcoeffs_pos[1].cre_up[448317]);

/*========================================================================*/
/* XXIII) Orthogonalize coefficients                                      */
  
  
  if(myid_state<np_states){
    control_init_cp_orthog(generalDataMini,cpMini,&cpParse,cp_on,cp_md,myid);
  }//endif
  

  printf("fffffffff ccreal %lg\n",cpMini->cpcoeffs_pos[1].cre_up[448317]);

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

/*------------------------------------------------------------------------*/
}/*end routine*/ 
/*==========================================================================*/







