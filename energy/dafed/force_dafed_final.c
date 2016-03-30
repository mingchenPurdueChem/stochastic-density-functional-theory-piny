
#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_dafed_energy.h"
#include "../proto_defs/proto_math.h"


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

void force_dafed_final(CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos) { 


/*========================================================================*/
/*========================================================================*/

  int i,j;
  DAFED_INFO *dinfo = &(clatoms_info->dinfo);
  DAFED *dafed      = clatoms_info->dafed;
  int n_cv          = dinfo->n_cv;
  int num_atm_list;
  int *atm_list;
  double *clatoms_fx    = clatoms_pos->fx;
  double *clatoms_fy    = clatoms_pos->fy;
  double *clatoms_fz    = clatoms_pos->fz;
  double *Fx,*Fy,*Fz;
  double s,min,max,k_bdy;

  for(i=0;i<n_cv;i++){
     num_atm_list = dafed[i].num_atm_list;
     atm_list     = dafed[i].atm_list;
     Fx = dafed[i].Fx;
     Fy = dafed[i].Fy;
     Fz = dafed[i].Fz;
     for(j=0;j<num_atm_list;j++){
        clatoms_fx[atm_list[j]] += Fx[j];
        clatoms_fy[atm_list[j]] += Fy[j];
        clatoms_fz[atm_list[j]] += Fz[j];
     }
  }
  for(i=0;i<n_cv;i++){//force for bdy=3
     if(dafed[i].bdy==3){
       s = dafed[i].s;
       min = dafed[i].min;
       max = dafed[i].max;
       k_bdy = dafed[i].k_bdy;       
       if(s<min)dafed[i].Fs -= k_bdy*(s-min);
       if(s>max)dafed[i].Fs -= k_bdy*(s-max);
     }
  }     
  //for(i=0;i<n_cv;i++)printf("i %i force %lg ",i,dafed[i].Fs);
  //printf("\n");

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/
     

