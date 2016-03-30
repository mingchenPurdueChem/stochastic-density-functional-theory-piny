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
void scale_params(DAFED *daf) 
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */	
  double dt = daf->dt;
  double tau = daf->therm.tau;
  printf("Scaling to Atomic Units (d-AFED)\n");
  printf(" - Tau -> (tau*dt)**2 \n");
 
  tau *= dt;
  daf->therm.tau = tau*tau; 
  daf->therm.dT   = dt;  
  printf(" - ms  -> amu to a.u.\n");
 
  daf->ms *= PROT_MASS;
  
  printf(" - kTs ->\n");
  
  daf->kTs /= BOLTZ;
  daf->therm.kTs = daf->kTs; 
  printf(" - ks ->\n");
  
  daf->ks *=(BOHR*BOHR)/BOLTZ; 
  daf->k_bdy *= (BOHR*BOHR)/BOLTZ;

  daf->vs = -sqrt(daf->kTs/daf->ms);
  daf->s  = -4.0;
  daf->Fs = 0;
 
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initialize_therm(GGMT *th)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */	
  double nr = (double)(th->nr);
  double dT = th->dT;
  double kTs = th->kTs;
  th->ny = 3;
  th->Q[1]  = th->kTs*th->tau;
  th->Q[2]  = 8.0*th->tau*pow(th->kTs,3)/3.0;
  th->w[1]  = 1.0/(2.0-pow(2.0,(1.0/3.0)));
  th->w[2]  = 1.0-(2.0*th->w[1]);
  th->w[3]  = th->w[1];
  th->w[1]  *= dT/nr;
  th->w[2]  *= dT/nr;
  th->w[3]  *= dT/nr;
  th->qt[1] = 1.0;   th->qt[2] = 1.0;
  th->vt[1] = -sqrt(th->kTs/th->Q[1]);
  th->vt[2] = -sqrt(th->kTs/th->Q[2]);
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void get_Phi_params(DAFED *dafed, CLATOMS_INFO *clatoms_info)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
  int i,j;
  int num_atm_list = dafed->num_atm_list;
  int *atm_list    = dafed->atm_list;
  int natm_tot     = clatoms_info->natm_tot;

  if(num_atm_list!=4){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("number of atoms in dihedral should be 4,\n");
    printf("instead of %i\n",num_atm_list);
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(0);
  }
  for(i=0;i<num_atm_list;i++){
     if(atm_list[i]<1||atm_list[i]>natm_tot){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Invalid atom index,should larger or equal to 1\n");
       printf("and smaller or equal to %i, but you have an index %i\n",natm_tot,atm_list[i]);
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(0);
     }
  }
  for(i=0;i<num_atm_list;i++){
     for(j=i+1;j<num_atm_list;j++){
        if(atm_list[i]==atm_list[j]){
          printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          printf("Invalid atom index,should not equal with each other.\n");
          printf("Index %i equals index %i\n",i,j);
          printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          exit(0);
        }
     }
  }
          
  dafed->min *= M_PI/180.0;
  dafed->max *= M_PI/180.0;
  printf("min %lg max %lg\n",dafed->min,dafed->max);

  printf("==================================\n");
  printf("d-AFED Radius of Gyration\n");
  printf("-----------------\n");  
  for(i=0;i<num_atm_list;i++)printf("Dih Atm:%i\n",atm_list[i]);  
  printf("==================================\n");
  fflush(stdout);
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initialize_dafed(CLASS *class, GENERAL_DATA *general_data)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */	
//#include "../typ_defs/typ_mask.h"
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  DAFED *dafed = clatoms_info->dafed;
  DAFED_INFO *dinfo = &(clatoms_info->dinfo); 
  //ATOMMAPS *am = &(class->atommaps);
  int i,j,n_bin_phi,steps_bias,num_atm_list;
  int div,res;
  int iflag;
  int num_phi = 0;
  int n_cv = dinfo->n_cv;
  int bias_on = dinfo->bias_on;
  int nres_tra = (general_data->timeinfo.nres_tra);
  int nres_tor = (general_data->timeinfo.nres_tor);
  int nres_ter = (general_data->timeinfo.nres_ter);
  int ggmt_nr  = dinfo->ggmt_nr;
  double dt       = general_data->timeinfo.dt;
  double ggmt_tau = dinfo->ggmt_tau;
  char* traj_fn   = dinfo->traj_fn;
  FILE *t1;

  int blockcounts[1];
  int ndinfo = 12;

  t1 = fopen(traj_fn,"w");
  fprintf(t1,"### ");

  //printf("n_cv %i\n",n_cv); //%%%%%%%%%%
  dinfo->print_screen_freq = general_data->filenames.iwrite_screen;
  for(i=0;i<n_cv;i++){
/*------------------- copy parameters if necessary-----------------------------*/
     //dafed[i].dt = dt/((double)(nres_tra*nres_tor*nres_ter)); //if no respa then dt/1.0
     dafed[i].dt = dt;
     //printf("dt %lg\n",dafed[i].dt);
     dafed[i].therm.tau = ggmt_tau;
     dafed[i].therm.nr  = ggmt_nr;
     dafed[i].therm.kTs = dafed[i].kTs;
     //printf("nr %i\n",dafed[i].therm.nr);
     num_atm_list = dafed[i].num_atm_list;
     //printf("num_atm_list %i %i\n",i,num_atm_list);
/*-----------------rescale parameters&initialize thermostat---------------------*/
     scale_params(&(dafed[i]));
     initialize_therm(&(dafed[i].therm));     
/*-----------------------set boundary condition----------------------------------*/
/*-----------------------malloc for force, values--------------------------------*/
     dafed[i].Fx = (double*)cmalloc(num_atm_list*sizeof(double));
     dafed[i].Fy = (double*)cmalloc(num_atm_list*sizeof(double));
     dafed[i].Fz = (double*)cmalloc(num_atm_list*sizeof(double));     
/*-----------set forces and value and statisticle values to zero-----------------*/
     //zero_forces(&(dafed[i]));
     for(j=0;j<num_atm_list;j++){
        dafed[i].Fx[j] = 0.0;
        dafed[i].Fy[j] = 0.0;
        dafed[i].Fz[j] = 0.0;
     }
/*---------set parameters for each CV---------------*/

     switch(dafed[i].type){
       case 4:
	  num_phi += 1;
	  get_Phi_params(&(dafed[i]),clatoms_info);
	  fprintf(t1,"s_phi%i Phi%i(r) ",num_phi,num_phi);
          break;
     }
/*-----------------------set boundary condition----------------------------------*/
     dafed[i].diff = dafed[i].max-dafed[i].min;
/*---------set parameters for energy---------------*/
     dafed[i].kin_ave = 0.0;
     dafed[i].therm_k_ave = 0.0;
     dafed[i].temp_ave = 0.0; 
/*---------set parameters for biasing potential---------------*/
     if(bias_on==1){
       if(dafed[i].bin_num==0){
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
         printf("number of bin is zero!!\n");
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
         exit(0);
       }
       dafed[i].w_bin = (dafed[i].max-dafed[i].min)/dafed[i].bin_num;
     }
  }
/*--------------------------------set output file--------------------------------*/     
  fprintf(t1,"\n");
  fclose(t1);
  if(dinfo->bias_on==1)bias_init(dinfo,dafed); 
  //for(i=0;i<n_cv;i++)printf("mas %lg\n",dafed[i].ms); 

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/
