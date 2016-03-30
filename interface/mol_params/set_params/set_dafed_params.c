/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_mol_parms.c                              */
/*                                                                          */
/* This subprogram sets d-AFED parameters                                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_mol_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void set_dafed_params(FILENAME_PARSE *filename_parse,char fun_key[],
                      DICT_WORD dafed_dict[],int num_dafed_dict,DAFED *dafed,
                      int n_cv)
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/

  int i,index,list_length;
  char *cv_name,*bk_list;
  DAFED *dafed_temp;
  
  /*=======================================================================*/
  /* I) Check for missing key words*/
  
  for(i=1;i<num_dafed_dict;i++){
    if(dafed_dict[i].iuset==0 && dafed_dict[i].key_type==1){
      keyword_miss(dafed_dict,filename_parse->molsetname,fun_key,i);}
  }    /*endfor*/
  /*=======================================================================*/
  /* II) Set parameters */
  /*-----------------------------------------------------------------------*/ 
  /*  1)\atom_index{} */
  sscanf(dafed_dict[1].keyarg,"%d",&index);
  if(index>=n_cv||index<0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("invalid index %i n_cv %i\n",index,n_cv);
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }
  printf("dafed mol param %i\n",index);
  dafed_temp = &(dafed[index]);
  /*-----------------------------------------------------------------------*/ 
  /*  2)\name{} */  
  //if(strcasecmp(dafed_dict[2].keyarg,"Rgyr")==0){dafed_temp->type=2;}
  //if(strcasecmp(dafed_dict[2].keyarg,"NH")==0){dafed_temp->type=3;}
  //if(strcasecmp(dafed_dict[2].keyarg,"Ree")==0){dafed_temp->type=1;}
  if(strcasecmp(dafed_dict[2].keyarg,"Dih")==0){dafed_temp->type=4;}
  //if(strcasecmp(dafed_dict[2].keyarg,"Dih_cor")==0){dafed_temp->type=5;}
  //if(strcasecmp(dafed_dict[2].keyarg,"Nalpha")==0){dafed_temp->type=6;}
  //if(strcasecmp(dafed_dict[2].keyarg,"Nbeta")==0){dafed_temp->type=7;}
  //if(strcasecmp(dafed_dict[2].keyarg,"Alpha_cont")==0){dafed_temp->type=8;}
  //if(strcasecmp(dafed_dict[2].keyarg,"Gama_cont")==0){dafed_temp->type=9;}
  /*-----------------------------------------------------------------------*/ 
  /*  3)\atom_list{} */
  i=0;
  dafed_temp->atm_list = NULL;
  bk_list = strtok(dafed_dict[3].keyarg,",");
  while(bk_list!=NULL){
    i++;
    dafed_temp->atm_list = (int *)crealloc(dafed_temp->atm_list,i*sizeof(int));
    sscanf(bk_list,"%d",&(dafed_temp->atm_list[i-1]));
    bk_list = strtok(NULL,",");
    //printf("i %i\n",i);
  }
  dafed_temp->num_atm_list = i;
  /*-----------------------------------------------------------------------*/
  /* 4)\mass{} */
  sscanf(dafed_dict[4].keyarg,"%lg",&(dafed_temp->ms));
  /*-----------------------------------------------------------------------*/
  /* 5)\temperature{} */
  sscanf(dafed_dict[5].keyarg,"%lg",&(dafed_temp->kTs));
  /*-----------------------------------------------------------------------*/
  /* 6)\force_const{} */
  sscanf(dafed_dict[6].keyarg,"%lg",&(dafed_temp->ks));
  /*-----------------------------------------------------------------------*/
  /* 7)\cv_min{} */ //degree if dihedral
  sscanf(dafed_dict[7].keyarg,"%lg",&(dafed_temp->min));
  /*-----------------------------------------------------------------------*/
  /* 8)\cv_max{} */ //degree if dihedral
  sscanf(dafed_dict[8].keyarg,"%lg",&(dafed_temp->max));
  /*-----------------------------------------------------------------------*/
  /* 9)\bdry{} */
  if(strcasecmp(dafed_dict[9].keyarg,"hw")==0)dafed_temp->bdy = 1;
  if(strcasecmp(dafed_dict[9].keyarg,"per")==0)dafed_temp->bdy = 2;
  if(strcasecmp(dafed_dict[9].keyarg,"sw")==0)dafed_temp->bdy = 3;
  /*-----------------------------------------------------------------------*/
  /* 9)\bin_num{} */ //number of bins, lattice point = bin_num+1
  sscanf(dafed_dict[10].keyarg,"%i",&(dafed_temp->bin_num));


/*-----------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/




