#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_dafed_energy.h"
#include "../proto_defs/proto_math.h"

/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

void force_bias_num(DAFED_INFO *dinfo,DAFED *dafed)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
  BIAS_HIST *bhist   = dinfo->bhist;
  C_NBLT *nblt_array = dinfo->nblt_array;
  int i,j,k;
  int lat_num;
  int sum_exp;
  int nblt_tot    = dinfo->nblt_tot;
  int n_lat_tot   = dinfo->n_lat_tot;
  int n_cv        = dinfo->n_cv;
  int *n_lat      = dinfo->n_lat;
  int *index      = (int*)malloc(n_cv*sizeof(int));
  int *index_temp = (int*)malloc(n_cv*sizeof(int));
  int *index_c    = (int*)malloc(n_cv*sizeof(int));
  int *type       = (int*)malloc(n_cv*sizeof(int));
  int ind,ind_c;
  int *diff,*diff_c;
  int bin_num;
  double min,max,w_bin;
  double pot_interp = 0.0;
  double *f_bias = (double*)malloc(n_cv*sizeof(double));
  double *s      = (double*)malloc(n_cv*sizeof(double));
  double *x,*x_c;
  double pre = 1.0;
  double prod;
  double pot_v;
  double lat_m,lat_p;

  /*printf("before\n");
  for(i=0;i<nblt_tot;i++){
     for(j=0;j<n_cv;j++)printf("%i ",nblt_array[i].diff[j]);
     printf("\n");
  }*/

  for(i=0;i<n_cv;i++){
     bin_num = dafed[i].bin_num;
     min = dafed[i].min;
     max = dafed[i].max;
     //printf("min %lg max %lg\n",min,max);
     w_bin = dafed[i].w_bin;
     s[i] = dafed[i].s;
     type[i] = dafed[i].type;
     index[i] = (int)((s[i]-min)/w_bin);
     if(index[i]>=bin_num)index[i] = bin_num-1; //prevent the round error
     lat_m = min+index[i]*w_bin;
     lat_p = min+(index[i]+1)*w_bin;
     pre *= w_bin;
     f_bias[i] = 0.0;
  }
  /*printf("after\n");
  for(i=0;i<nblt_tot;i++){
     for(j=0;j<n_cv;j++)printf("%i ",nblt_array[i].diff[j]);
     printf("\n");
  }*/
  for(i=0;i<nblt_tot;i++){
     diff = nblt_array[i].diff;
     diff_c = nblt_array[nblt_tot-1-i].diff;
     for(j=0;j<n_cv;j++){
        //printf("index %i diff %i diff_c %i\n",index[j],diff[j],diff_c[j]);
        index_temp[j] = index[j]+diff[j];
        index_c[j] = index[j]+diff_c[j];
        //printf("index_temp %i index_c %i lat_num %i\n",index_temp[j],index_c[j],dafed[j].lat_num);
     }
     ind   = get_ind(index_temp,n_lat,n_cv); 
     ind_c = get_ind(index_c,n_lat,n_cv);
     //trap
     if(ind>=n_lat_tot||ind_c>=n_lat_tot||ind<0||ind_c<0){
       for(j=0;j<n_cv;j++){
          printf("s %lg min %lg\n",dafed[j].s,dafed[j].min);
          printf("index %i diff %i diff_c %i\n",index[j],diff[j],diff_c[j]);
          printf("index_temp %i index_c %i lat_num %i\n",index_temp[j],index_c[j],dafed[j].lat_num);
       }
       printf("ind %i ind_c %i n_lat_tot %i\n",ind,ind_c,dinfo->n_lat_tot);
       exit(0);
     }
     pot_v = bhist[ind].gaussian_v;
     x_c = bhist[ind_c].x;
     prod = nblt_array[i].sign;
     for(j=0;j<n_cv;j++){
        if(type[j]==4)prod *= diff_perodic(x_c[j],s[j]);
        else prod *= x_c[j]-s[j];
     }
     pot_interp += pot_v*prod;
     for(j=0;j<n_cv;j++){
        prod = nblt_array[i].sign;
        for(k=0;k<n_cv;k++){
           if(k!=j){
              if(type[k]==4)prod *= diff_perodic(x_c[k],s[k]);
              else prod *= x_c[k]-s[k];
           }
        }
        f_bias[j] += pot_v*prod;
     }
  }
  pot_interp /= pre;
  dinfo->pot_bias = pot_interp;
  for(i=0;i<n_cv;i++){
     dafed[i].Fs += f_bias[i]/pre;
  } 
  //debug
  //printf("pot %lg\n",dinfo->pot_bias);
  //for(i=0;i<n_cv;i++)printf("%lg ",f_bias[i]); 
  //printf("\n");
  //
  free(index);
  free(index_temp);
  free(f_bias);
  free(s);
  free(index_c);
  free(type);

         
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/
        
/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
double diff_perodic(double x,double y)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
  double diff = x-y;
  if(diff>=M_PI)diff -= 2.0*M_PI;
  if(diff<-M_PI)diff += 2.0*M_PI;
  return diff;
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/






















 






