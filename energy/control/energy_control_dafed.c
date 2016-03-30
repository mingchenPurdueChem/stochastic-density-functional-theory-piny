/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: energy_control.c                               */
/*                                                                          */
/* This routine calls the required force and PE routines                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_dafed_energy.h"
#include "../proto_defs/proto_math.h"


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

void energy_control_dafed(CLASS *class, BONDED *bonded, 
                    GENERAL_DATA *general_data)

/*========================================================================*/
  { /* Begin Routine */
/*========================================================================*/

  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS  *clatoms_pos  = &(class->clatoms_pos[1]);//no pimd
  
  DAFED_INFO *dinfo = &(clatoms_info->dinfo);
  DAFED *dafed      = clatoms_info->dafed;

  int i,j,k;
  int num_atm_list;
  int type;
  int n_cv    = dinfo->n_cv;
  int bias_on = dinfo->bias_on;
  
  for(i=0;i<n_cv;i++){
     //printf("min %lg max %lg\n",dafed[i].min,dafed[i].max);
     num_atm_list = dafed[i].num_atm_list;
     //printf("num_atm_list %i\n",num_atm_list);
     type = dafed[i].type;
     for(j=0;j<num_atm_list;j++){
        dafed[i].Fx[j] = 0.0;
        dafed[i].Fy[j] = 0.0;
        dafed[i].Fz[j] = 0.0;
     }
     switch(type){
       case 4: force_Phi(dinfo,&(dafed[i]),clatoms_pos);break;
     }
  }
/*add biasing potential*/
  if(bias_on==1){
    force_bias_num(dinfo,dafed);
  }

  force_dafed_final(clatoms_info,clatoms_pos);

/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/





