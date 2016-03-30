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
void check_boundaries(DAFED *dafed)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */	

   int bdy = dafed->bdy;
   double min = dafed->min;
   double max = dafed->max;
   double s = dafed->s;
   double diff = dafed->diff;
      
   if(bdy==1){ //hard boundary
     if(s<min){
       dafed->s = 2.0*min-s;
       dafed->vs *= -1.0;
     }
     if(s>max){
       dafed->s = 2.0*max-s;
       dafed->vs *= -1.0;
     }
   }
   if(bdy==2){ //periodic boundary
     if(s<min)dafed->s += diff;
     if(s>max)dafed->s -= diff;
   }
   /*The harmonic boundary is left to the force calculation part*/   
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/
     

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_dafed_ggmt(DAFED *daf)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */	
   GGMT *th = &(daf->therm);

   //printf("DAFED GGMT NOW!\n");
   int i,j,k;
   int nr = th->nr;
   int ny = th->ny;
   double *G = th->G;
   double *Q = th->Q;
   double *vt = th->vt;
   double *qt = th->qt;
   double *w  = th->w;
   double dT = th->dT;
   double tao = th->tau;
   double kTs = daf->kTs;
   double ms = daf->ms;
   double vs = daf->vs;
   double kt = vs*vs*ms;

   double aa,bb;
   /*for debug*/
   /*printf("-------\n");
   printf("nr %i\n",nr);
   printf("ny %i\n",ny);
   printf("Q1 %lg Q2 %lg\n",Q[1],Q[2]);
   printf("w1 %lg w2 %lg w3 %lg\n",w[1],w[2],w[3]);
   printf("dT %lg\n",dT);*/

   for(i=1;i<=nr;i++){
      for(j=1;j<ny;j++){
         G[1] = (kt-kTs)/Q[1];
         G[2] = (kt*kt/3.0-kTs*kTs)/Q[2];
         vt[1] += 0.25*w[j]*G[1];
         vt[2] += 0.25*w[j]*G[2];

         aa = exp(-0.125*w[j]*(vt[1]+kTs*vt[2]));
         //printf("aa1 %lg\n",aa);
         vs *= aa;
         kt = vs*vs*ms;
  
         bb = kt*vt[2]/3.0;
         vs *= sqrt(1.0/(1.0+0.5*w[j]*bb));

         aa = exp(-0.125*w[j]*(vt[1]+kTs*vt[2]));
         vs *= aa;
         kt = vs*vs*ms;

         qt[1] += 0.5*w[j]*vt[1];
         qt[2] += 0.5*w[j]*vt[2]*(kTs+kt);

         aa = exp(-0.125*w[j]*(vt[1]+kTs*vt[2]));
         vs *= aa;
         kt = vs*vs*ms;
  
         bb = kt*vt[2]/3.0;
         vs *= sqrt(1.0/(1.0+0.5*w[j]*bb));

         aa = exp(-0.125*w[j]*(vt[1]+kTs*vt[2]));
         vs *= aa;
         kt = vs*vs*ms;

         G[1] = (kt-kTs)/Q[1];
         G[2] = (kt*kt/3.0-kTs*kTs)/Q[2];
         vt[1] += 0.25*w[j]*G[1];
         vt[2] += 0.25*w[j]*G[2];
      }
   }
   //printf("vs ratio %lg\n",vs/daf->vs); 
   daf->vs = vs;

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void evolve_position(DAFED *daf)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */	
   daf->s += daf->dt*daf->vs;  
   check_boundaries(daf);
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void evolve_velocity(DAFED *daf)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */	
 //  printf(" d-AFED vel update\n");
  daf->vs += 0.5*daf->dt*daf->Fs/(daf->ms);  
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_dafed_0_to_dt2 (CLASS *class, GENERAL_DATA *general_data)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */	

   CLATOMS_INFO *clatoms_info = &(class->clatoms_info);

   DAFED *dafed = clatoms_info->dafed;
   DAFED_INFO *dinfo = &(clatoms_info->dinfo);
   int n_cv = dinfo->n_cv;
   int i;


   for(i=0;i<n_cv;i++){
      int_dafed_ggmt(&dafed[i]);
      evolve_velocity(&dafed[i]);
      evolve_position(&(dafed[i]));
   }      

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_dafed_dt2_to_dt (CLASS *class, GENERAL_DATA *general_data)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */	

   CLATOMS_INFO *clatoms_info = &(class->clatoms_info);

   DAFED *dafed = clatoms_info->dafed;
   DAFED_INFO *dinfo = &(clatoms_info->dinfo);
   int n_cv = dinfo->n_cv;
   int i;

   for(i=0;i<n_cv;i++){
      evolve_velocity(&(dafed[i]));
      int_dafed_ggmt(&(dafed[i]));
   }      
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_dafed_0_to_dt2_cl (CLASS *class, GENERAL_DATA *general_data,
                            int int_ggmt_flag,double dt)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

   CLATOMS_INFO *clatoms_info = &(class->clatoms_info);

   DAFED *dafed = clatoms_info->dafed;
   DAFED_INFO *dinfo = &(clatoms_info->dinfo);
   int n_cv = dinfo->n_cv;
   int i;


   for(i=0;i<n_cv;i++){
      dafed[i].dt = dt;
      if(int_ggmt_flag==1)int_dafed_ggmt(&dafed[i]);
      evolve_velocity(&dafed[i]);
      evolve_position(&(dafed[i]));
   }

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_dafed_dt2_to_dt_cl (CLASS *class, GENERAL_DATA *general_data,
                             int int_ggmt_flag,double dt)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

   CLATOMS_INFO *clatoms_info = &(class->clatoms_info);

   DAFED *dafed = clatoms_info->dafed;
   DAFED_INFO *dinfo = &(clatoms_info->dinfo);
   int n_cv = dinfo->n_cv;
   int i;

   for(i=0;i<n_cv;i++){
      dafed[i].dt = dt;
      evolve_velocity(&(dafed[i]));
      if(int_ggmt_flag==1)int_dafed_ggmt(&(dafed[i]));
   }
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_dafed_final (DAFED_INFO *dinfo,DAFED *dafed,TIMEINFO *timeinfo)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */	

   int n_cv              = dinfo->n_cv;
   int bias_on           = dinfo->bias_on;
   int steps_bias        = dinfo->steps_bias;
   int itime             = timeinfo->itime;
   int print_screen_freq = dinfo->print_screen_freq;
   int write_freq        = dinfo->write_freq;
   int i;

   double q,s,vs,ms;

   dinfo->ene_tot = 0.0;
   for(i=0;i<n_cv;i++){
      q = dafed[i].q;
      s = dafed[i].s;
      vs = dafed[i].vs;
      ms = dafed[i].ms;
      dafed[i].kin = 0.5*vs*vs*ms;      
      get_dafed_therm_energy(&(dafed[i].therm),&(dafed[i].therm_k),&(dafed[i].therm_p));
      dafed[i].temp = 2.0*dafed[i].kin*3.1577464e5;
      dafed[i].kin_ave += dafed[i].kin;
      dafed[i].temp_ave += dafed[i].temp;
      dafed[i].therm_k_ave += dafed[i].therm_k;
      dinfo->ene_tot += dafed[i].kin+dafed[i].pot_harm+dafed[i].therm_k+dafed[i].therm_p;
   }   
   if(bias_on==1&&itime%steps_bias==0)update_gau(dinfo,dafed);//update biasing potential 
  
   if(itime%print_screen_freq==0)dafed_screen_io(dinfo,dafed,itime);
   dafed_trajectory_io(dinfo,dafed,itime);

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/   

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void get_dafed_therm_energy(GGMT *th, double *ke,double *pe)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */	
  double *Q = th->Q;
  double *vt = th->vt;
  double *qt = th->qt;
  double kTs = th->kTs;
  *ke = 0.5*Q[1]*vt[1]*vt[1]+0.5*Q[2]*vt[2]*vt[2];
  *pe = kTs*(qt[1]+qt[2]);  
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/





