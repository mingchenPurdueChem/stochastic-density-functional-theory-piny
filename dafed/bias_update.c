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

void bias_init(DAFED_INFO *dinfo,DAFED *dafed)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
  BIAS_HIST *bhist;         
  C_NBLT *nblt_array;
   
  double min,max;
  double w_bin;
  double *gau_v;
  int i,j,k,l;
  int n_cv = dinfo->n_cv;
  int type;
  int bin_num,lat_num;
  int *n_bin,*n_frag,*n_lat;
  int nblt_tot = 1;
  int ind,ind_s,sum_exp;
  int n_lat_tot;
  int n_lat_chunk;
  int *index; //need free
  int bias_readin = dinfo->bias_readin;
  int count;
  FILE *bias_in;

  n_lat_tot = 1;
  for(i=0;i<n_cv;i++){
     min = dafed[i].min;
     max = dafed[i].max;
     type = dafed[i].type;
     bin_num = dafed[i].bin_num;
     lat_num = bin_num+1;
     dafed[i].lat_num = lat_num;
     w_bin = (max-min)/bin_num;
     n_lat_tot *= lat_num;
     nblt_tot *= 2;
  }
  dinfo->n_lat_tot = n_lat_tot;
  n_lat = (int*)cmalloc(n_cv*sizeof(int));
  n_lat_chunk = n_lat_tot;
  for(i=0;i<n_cv;i++){
     lat_num = dafed[i].lat_num;
     n_lat_chunk /= lat_num;
     n_lat[i] = n_lat_chunk;
  }
  dinfo->n_lat = n_lat;
  dinfo->nblt_tot = nblt_tot;
  //gau_v = (double*)cmalloc(n_lat_tot*sizeof(double));
  //dinfo->gau_v = gau_v;
  bhist = (BIAS_HIST*)cmalloc(n_lat_tot*sizeof(BIAS_HIST));
  dinfo->bhist = bhist;
  for(i=0;i<n_lat_tot;i++)bhist[i].x = (double*)cmalloc(n_cv*sizeof(double));

/* initilize neibourhood lattice point */
/***************************************/
/*     (1,0) -------------- (1,1)      */
/*           |            |            */
/*           |            |            */
/*           |            |            */
/*           |            |            */
/*           |            |            */
/*           |            |            */
/*     (0,0) -------------- (0,1)      */
/***************************************/


  dinfo->nblt_array = (C_NBLT*)cmalloc(nblt_tot*sizeof(C_NBLT));
  nblt_array = dinfo->nblt_array;
  for(i=0;i<nblt_tot;i++){
     nblt_array[i].diff = (int*)cmalloc(n_cv*sizeof(int));
  }
  for(i=0;i<nblt_tot;i++){
     ind = i;
     sum_exp = 0;
     for(j=0;j<n_cv;j++){
        printf("i %i j %i\n",i,j);
        ind_s = (ind>>1)<<1;
        printf("ind %i ind_s %i\n",ind,ind_s);
        nblt_array[i].diff[j] = ind-ind_s;
        ind = ind>>1;
        sum_exp += nblt_array[i].diff[j];
     }
     if(sum_exp%2==0)nblt_array[i].sign = 1.0;
     else nblt_array[i].sign = -1.0;
  }

  /*get all lattice points' positions*/
  index = (int*)cmalloc(n_cv*sizeof(int));
  for(i=0;i<n_lat_tot;i++){
     get_ind_rev(index,n_lat,i,n_cv);
     for(j=0;j<n_cv;j++){
        min = dafed[j].min;
        max = dafed[j].max;
        w_bin = dafed[j].w_bin;
        bhist[i].x[j] = min+index[j]*w_bin;
     }
  }
  //readin the existing biasing potential
  if(bias_readin==1){
    bias_in = fopen("bias_readin","r");	
    if(bias_in==NULL){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");	    
      printf("The name of the input file of biasing potential\n");
      printf("must be bias_readin\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }
    else{
      count = 0;	    
      for(i=0;i<n_lat_tot;i++){
         fscanf(bias_in,"%lg",&(bhist[i].gaussian_v));
	 if(feof(bias_in)!=0)count += 1;
      }
      if(count!=n_lat_tot){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	printf("We need %i lattice point, but the file contains\n",n_lat_tot);
	printf("%i lattice point\n",count);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");	
	fflush(stdout);
	exit(0);
      }
    }
  }    
  
  free(index);

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void update_gau(DAFED_INFO *dinfo,DAFED *dafed)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  BIAS_HIST *bhist               = dinfo->bhist;

  int i,j;
  int type;
  int n_cv      = dinfo->n_cv;
  int n_lat_tot = dinfo->n_lat_tot;
  int *n_lat    = dinfo->n_lat;
  int *index;  //need free
  double diff;
  double s;
  double *x;   //lattice point
  double min,max;
  double dist;
  double A = dinfo->A;
  double sigma = dinfo->sigma;
  
  index = (int*)cmalloc(n_cv*sizeof(int));
  for(i=0;i<n_lat_tot;i++){
    x = bhist[i].x;
    dist = 0.0;
    for(j=0;j<n_cv;j++){
       min = dafed[j].min;
       max = dafed[j].max;
       s = dafed[j].s;
       type = dafed[j].type;
       diff = fabs(s-x[j]);
       if(type==4){
         if(diff>=M_PI)diff -= 2.0*M_PI;
       }
       dist += diff*diff;
       //printf("%lg ",x[j]);
    }
    //printf("\n");
    bhist[i].gaussian_v += A*exp(-sigma*dist);
  }
  //exit(0);

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

















