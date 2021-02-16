/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                      Module: main                                        */
/*                                                                          */
/* This program performs MD on a classical potential energy surface (PES),  */
/* minimization on a classical PES,                                         */
/* MD on a mixed classical-density functional PES,                          */
/* minimization on a mixed classical-density functional PES,                */
/* PIMD on a classical PES,                                                 */
/* centroid minimization on a classical PES,                                */
/* PIMD on a mixed classical-density functional PES,                        */
/* centroid minimization on a mixed classical-density functional PES.       */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_main_entry.h"
#include "../proto_defs/proto_main_local.h"
#include "../proto_defs/proto_main_cp_local.h"
#include "../proto_defs/proto_parse_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_communicate_entry.h"
#include "../proto_defs/proto_stodft_entry.h"

#ifdef FFTW3
 #include "fftw3.h"
#endif
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 int main (int argc, char *argv[])
/*==========================================================================*/
   {/* begin routine */ 
/*==========================================================================*/
/*   Local Variables */

  int iii,is;
  CLASS class;
  BONDED bonded;
  GENERAL_DATA general_data;
  CP cp;
  ANALYSIS analysis;
  //debug ntask-------------------------
    //char processor_name[MPI_MAX_PROCESSOR_NAME];
    //int  namelen;
  //-----------------------------
#ifdef FAST_FILTER
  CLASS class2;
  BONDED bonded2;
  GENERAL_DATA general_data2;
  CP cp2;
  ANALYSIS analysis2;
  MPI_Group world_group;
  MPI_Group calc_group;
  int np_local;
  int myid_local;
  int i;
  int *index_local;
#endif

/*=======================================================================*/
/*  I)             Check for input file                                  */

  if(argc < 2) {
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("No input file specified\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

/*=======================================================================*/
/* II)            Initialize MPI                                         */

#ifdef FAST_FILTER
  /*
  // Here we wan to leave (np_local-1)'th process to do SVD and gemm ONLY
  // The (np_local-1)'th process will be on a different node

  Init(&argc,&argv,&class.communicate.world_ext);
  // Here we wan to leave (np_local-1)'th process to do SVD and gemm ONLY
  // The (np_local-1)'th process will be on a different node
  MPI_Comm_rank(class.communicate.world_ext,&class.communicate.myid_ext);
  MPI_Comm_size(class.communicate.world_ext,&class.communicate.np_ext);
  MPI_Comm_group(class.communicate.world_ext,&world_group);
  if(class.communicate.np_ext>1){
    //np_local = class.communicate.np_ext-1;
    np_local = class.communicate.np_ext;
  }
  else{
    // If I have only one process, nothing happens
    np_local = class.communicate.np_ext;
  }
  index_local = (int*)cmalloc(np_local*sizeof(int));
  for(i=0;i<np_local;i++)index_local[i] = i;
  MPI_Group_incl(world_group, np_local, index_local,&calc_group);
  MPI_Comm_create_group(class.communicate.world_ext,calc_group, 0, &class.communicate.world);
  class.communicate.np = -1;
  class.communicate.myid = -1;
  if(class.communicate.world!=MPI_COMM_NULL){
    Comm_size(class.communicate.world,&class.communicate.np);
    Comm_rank(class.communicate.world,&class.communicate.myid);
  }
  cfree(index_local);
  */
  Init(&argc,&argv,&class.communicate.world);
  Comm_size(class.communicate.world,&class.communicate.np);
  Comm_rank(class.communicate.world,&class.communicate.myid);
#else
  Init(&argc,&argv,&class.communicate.world);
  Comm_size(class.communicate.world,&class.communicate.np);
  Comm_rank(class.communicate.world,&class.communicate.myid);
#endif
  general_data.error_check_on = (class.communicate.myid==0?1:0);
  //debug ntask-------------------------
    //MPI_Get_processor_name(processor_name,&namelen);

    //fprintf(stderr,"Process %d of %d running on %s\n",
    //               class.communicate.myid, class.communicate.np, processor_name);
    //fflush(stderr);
  //-------------------------
#ifdef FAST_FILTER
  //MPI_Comm_dup(class.communicate.world_ext,&class2.communicate.world_ext);  
  MPI_Comm_dup(class.communicate.world,&class2.communicate.world);
  //class2.communicate.myid_ext = class.communicate.myid_ext;
  //class2.communicate.np_ext = class.communicate.np_ext;
  class2.communicate.myid = class.communicate.myid;
  class2.communicate.np = class.communicate.np;
  general_data2.error_check_on = (class.communicate.myid==0?1:0);
#endif
/*=======================================================================*/
/* III)            Invoke User Interface                                 */

  parse(&class,&bonded,&general_data,&cp,&analysis,argv[1]);

#ifdef FAST_FILTER
  parse(&class2,&bonded2,&general_data2,&cp2,&analysis2,argv[2]);
#endif
/*========================================================================*/
/* IV)              Perform Simulation                                    */

  /*----------------------------------------------------------------------*/
  /*  i) RUN MD/PIMD/CP/CPPIMD */
  if(general_data.simopts.md==1){
    control_md(&class,&bonded,&general_data,&analysis);
  }
  if(general_data.simopts.pimd==1){
    control_pimd(&class,&bonded,&general_data,&analysis);
  }
  if((general_data.simopts.cp_wave+general_data.simopts.cp) ==1){
      control_cp(&class,&bonded,&general_data,&cp,&analysis);
  }/*endif*/
  if(general_data.simopts.cp_pimd+general_data.simopts.cp_wave_pimd==1){
      control_cp_pimd(&class,&bonded,&general_data,&cp,&analysis);
  }/*endif*/

  /*----------------------------------------------------------------------*/
  /*  ii) DEBUG MD/PIMD/CP/CPPIMD */
  if(general_data.simopts.debug==1){
    control_debug(&class,&bonded,&general_data);
  }
  if(general_data.simopts.debug_pimd==1){
    control_debug_pimd(&class,&bonded,&general_data);
  }
  if(general_data.simopts.debug_cp==1){
    control_debug_cp(&class,&bonded,&general_data,&cp);
  }
  if(general_data.simopts.debug_cp_pimd==1){
    control_debug_cp_pimd(&class,&bonded,&general_data,&cp);
  }

  /*----------------------------------------------------------------------*/
  /*  iii) MINIMIZE MD/PIMD/CP/CPPIMD */
  if(general_data.simopts.minimize==1){
   control_min(&class,&bonded,&general_data,&analysis);  
  }
  if((general_data.simopts.cp_wave_min+general_data.simopts.cp_min) ==1 ){
    if(cp.cpopts.stodftOn==0){
      //Regular DFT scheme
      control_cp_min(&class,&bonded,&general_data,&cp,&analysis);
    }
    else{
      //Stochastic DFT scheme
#ifdef FAST_FILTER
      controlStodftMin(&class,&bonded,&general_data,&cp,&analysis,&class2,&bonded2,&general_data2,&cp2,&analysis2);
#else
      controlStodftMin(&class,&bonded,&general_data,&cp,&analysis);
#endif
    }
  }
  if((general_data.simopts.cp_wave_min_pimd) ==1 ){
    control_cp_pimd_min(&class,&bonded,&general_data,&cp,&analysis);
  }

/*==========================================================================*/
/* V)                Exit Program                                           */

  if(class.communicate.np>1){
   Barrier(class.communicate.world);
   Finalize();
  }/*endif*/
  fflush(stdout);
  exit(0); 
  return 0;

/*----------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/




