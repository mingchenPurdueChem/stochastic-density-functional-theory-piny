/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control-stodft.c                             */
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
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_stodft_local.h"

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

  MPI_Comm comm_states = class->communicate.comm_states;
  MPI_Comm world       = class->communicate.world;

  STAT_AVG *stat_avg = &(general_data->stat_avg);

  int myid = class->communicate.myid;
  int numProc = cp->communicate.np;
  int numTime = general_data->timeinfo.ntime;
  int iTime;
  int ip_now = 1;

  double elecEnergy,elecEnergyOld,elecEnergyOldTemp,elecEnergyTemp;
  double deltaEnergy;


/*======================================================================*/
/* I) Write to Screen                                                   */

 if(myid==0){
  PRINT_LINE_STAR;
  printf("Running Stochastic DFT.\n");
  PRINT_LINE_DASH;
 }/* endif */

/*======================================================================*/
/* II) Initialize the stochastic DFT                                    */

  initStodft(class,bonded,general_data,cp,ip_now);

/*======================================================================*/
/* III) Initial call to output_cp_min: need to open confp file          */

 


/*======================================================================*/
/* IV) Electronic Structure calculation for initial configuration       */
  
  scfStodft(class,bonded,general_data,cp,ip_now);


/*======================================================================*/
/* V) Loop over the specified number of time steps for geometric	*/
/*    optimization.							*/
 
  numTime = 0; //I don't want to see this loop in today's debug
  for(iTime=1;iTime<numTime;iTime++){

  /*---------------------------------------------------------------------*/
  /* 1) Update Atom Coordinate                                           */

//Let's temperaly test the filter by given the density and chemical potential, test how the filter looks like on KS eigenfunctions. 

  /*---------------------------------------------------------------------*/
  /* 2) Run SCF calculation                                              */

  elecEnergyOld     = stat_avg->cp_eke
	              + stat_avg->cp_enl
		      + stat_avg->cp_ehart
		      + stat_avg->cp_exc
		      + stat_avg->cp_eext;

  if(numProc>1){
   elecEnergyOldTemp = elecEnergyOld;
   Allreduce(&(elecEnergyOldTemp),&(elecEnergyOld),1,MPI_DOUBLE,MPI_SUM,0,
	     comm_states);
  }/*endif*/

  //Minimize with stochastic dft
  scfStodft(class,bonded,general_data,cp,ip_now);
  elecEnergy     = stat_avg->cp_eke
		 + stat_avg->cp_enl
		 + stat_avg->cp_ehart
		 + stat_avg->cp_exc
		 + stat_avg->cp_eext;
  if(numProc>1){
   elecEnergyTemp = elecEnergy;
   Allreduce(&(elecEnergyTemp),&(elecEnergy),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
  }/*endif*/
  deltaEnergy = fabs(elecEnergy-elecEnergyOld);

  /*----------------------------------------------------------------------*/
  /*   6)Calculate some simple averages                                   */

  /*-----------------------------------------------------------------------*/
  /*   7)Produce the output specified by the user                          */


  /*---------------------------------------------------------------------*/
  /*  8) Analysis Routine                                               */

 /*---------------------------------------------------------------------*/
 /*   Check for exit condition                                      */

  }//endfor iTime

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






