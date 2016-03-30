#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_dafed_entry.h"
#include "../proto_defs/proto_dafed_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void dafed_screen_io (DAFED_INFO *dinfo,DAFED *dafed,int itime)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */	
   
   int i;
   int type;
   int n_cv = dinfo->n_cv;

   double kin,kin_ave;
   double pot_harm;
   double therm_k,therm_k_ave;
   double therm_p;
   double temp,temp_ave;
   double pot_bias = dinfo->pot_bias;

   for(i=0;i<n_cv;i++){
      type =        dafed[i].type;
      kin =         dafed[i].kin;
      kin_ave =     dafed[i].kin_ave;
      pot_harm =    dafed[i].pot_harm;
      therm_k =     dafed[i].therm_k;
      therm_k_ave = dafed[i].therm_k_ave;
      therm_p     = dafed[i].therm_p;      
      temp =        dafed[i].temp;
      temp_ave =    dafed[i].temp_ave;
      switch(type){
        case 4: printf("-----dih------\n");
      }
      printf("kinetic energy %.6lg  %.6lg\n",kin,kin_ave/itime); 
      printf("harmonic potential energy %.6lg\n",pot_harm);
      printf("thermostat kinetic energy %.6lg  %.6lg\n",therm_k,therm_k_ave/itime);
      printf("thermostat potential %.6lg\n",therm_p);
      printf("temperature %.6lg  %.6lg\n",temp,temp_ave/itime); 
      printf("---------------\n");
   }
   printf("biasing potential %.6lg\n",pot_bias);  
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/   

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void dafed_trajectory_io(DAFED_INFO *dinfo,DAFED *dafed,int itime)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
   int n_cv = dinfo->n_cv;
   int i,j;
   int write_freq = dinfo->write_freq;
   int bias_write_freq = dinfo->bias_write_freq; 
   int n_lat_tot = dinfo->n_lat_tot;
   int bias_on = dinfo->bias_on;
   FILE *t_ext,*bias_out;
   BIAS_HIST *bhist = dinfo->bhist;
   char* traj_fn = dinfo->traj_fn;

   if(itime%write_freq==0){
     t_ext = fopen(traj_fn,"a");
     for(i=0;i<n_cv;i++){
        fprintf(t_ext,"%.6lg %.6lg ",dafed[i].s,dafed[i].q);	   
     }
     fprintf(t_ext,"\n");
     fclose(t_ext);
   }
   if(bias_on==1&&itime%bias_write_freq==0){
     bias_out = fopen("biasing_pot","a");
     for(i=0;i<n_lat_tot;i++){
	for(j=0;j<n_cv;j++)fprintf(bias_out,"%lg ",bhist[i].x[j]);     
	fprintf(bias_out,"%lg\n",bhist[i].gaussian_v);
     }
     fprintf(bias_out,"------------------\n");
     fclose(bias_out);
   }
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/
   
      
	  
