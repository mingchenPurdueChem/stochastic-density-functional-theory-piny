/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_inter_params                         */
/*                                                                          */
/* This reads in and sets up the electron-atom interaction pseudopotential  */ 
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_inter_params_entry.h"
#include "../proto_defs/proto_inter_params_local.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_search_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void controlInterParamsFrag(GENERAL_DATA *generalDataMini,CLASS *classMini,
			    CP *cpMini,BONDED *bondedMini,CP *cp,
			    SPLINE_PARSE *spline_parse,FILENAME_PARSE *filename_parse
			    CLASS_PARSE *class_parse)

/*======================================================================*/
/*  Begin routine */
     {/*begin routine*/
/*======================================================================*/
/*          Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  INTERACT *interact = &(classMini->interact);
  DATA_BASE_ENTRIES *inter_base;          /* Lst: Database parameters    */
  CATM_LAB *cinter,*cinter_base;
  DICT_WORD *fun_dict;

  int i,j,iii;                           /* Num: Counters               */
  int ninter;                            /* Num: Number of interactions */
  int nsearch,natm_srch;
  int num_fun_dict;
  int ifirst,ityp;
  int nsplin_mall;
  int nsplin_mall_tot;
  int ninter_mall;
  int ninter_unique;                     /* Num: number of interactions with 
                                             unqiue paramters */
  int ninter_unique_mall;
  int nbase,nbase2,ibase_want;
  int ncharge = classMini->clatoms_info.nchrg;
  int natm_tot = classMini->clatoms_info.natm_tot;
  int natm_typ = classMini->atommaps.natm_typ;
  int iperd = generalDataMini->cell.iperd;
  int ishift_pot = class_parse->ishift_pot;
  int int_res_ter = generalDataMini->timeinfo.int_res_ter;
  int myid = classMini->communicate.myid;
  int num_proc = classMini->communicate.np;
  int iatm_typ[] = classMini->atommaps.iatm_atm_typ;
  NAME atm_typ[] = classMini->atommaps.atm_typ;
  MPI_Comm comm = classMini->communicate.world;
  
  int *inter_label;
  int *ifound,*isearch,*igood;           /* Lst: found,search goodness flags*/

  char typ[5];
  char *fun_key;


  double alp_ewd = generalDataMini->ewald.alp_ewd;
  double now_mem;                        /* Num: Memory allocated here  */

  double *eps,*sig;                      /* Lst: Lennard-Jones params   */
  double *awill,*bwill,*c6m,*c8m,*c10m;  /* Lst: Williams params        */
  double *cwill ,*rm_swit, *c9m;         /* Lst: Aziz-chen params       */
  double *temp_cutoff,*temp_cutoff_res,*temp_cutti;
 
/*======================================================================*/
/* 0) Write to screen                                                   */

/*======================================================================*/
/*  I) Malloc the memory                                                 */

  ninter = natm_typ*(natm_typ + 1)/2;
  ninter_mall = ninter;
  if((ninter_mall!=0)&&((ninter_mall %2)==0)){ninter_mall +=1;}
  inter_label = (int *) cmalloc(ninter_mall*sizeof(int))-1;
  eps        = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  sig        = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  awill      = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  bwill      = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  cwill      = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  rm_swit    = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  c6m        = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  c8m        = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  c9m        = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  c10m       = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  fun_key    = (char *)cmalloc(MAXWORD*sizeof(char));  
  cinter     = (CATM_LAB *)cmalloc(ninter*sizeof(CATM_LAB))-1;  
  ifound     = (int *)cmalloc(ninter*sizeof(int))-1;
  isearch    = (int *)cmalloc(ninter*sizeof(int))-1;
  igood      = (int *)cmalloc(ninter*sizeof(int))-1;

  temp_cutoff     = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  temp_cutoff_res = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  temp_cutti      = (double *) cmalloc(ninter_mall*sizeof(double))-1;

  interact->cutoff     = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  interact->cutoff_res = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  interact->cutti      = (double *) cmalloc(ninter_mall*sizeof(double))-1;

  interact->cutskin    = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  interact->cutskin_res= (double *) cmalloc(ninter_mall*sizeof(double))-1;
  interact->cutskin_root   = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  interact->cutskin_root_res= (double *) cmalloc(ninter_mall*sizeof(double))-1;

  interact->inter_map_index = (int *) cmalloc(ninter_mall*sizeof(int))-1;

/*======================================================================*/
/*  II) Set up the data structures                                      */

if(myid==0){
  ifirst =1;
  set_potfun_dict(&fun_dict,&num_fun_dict,ifirst);
  ityp = 0;
  for(i=1;i <= natm_typ; i++) {
    for(j=i;j <= natm_typ; j++) {
      ityp++;
      strcpy(cinter[ityp].atm1,atm_typ[i]);
      strcpy(cinter[ityp].atm2,atm_typ[j]);
    }/*endfor*/
  }/*endfor*/
  for(i=1;i<=ninter;i++){ifound[i]=0;}
  for(i=1;i<=ninter;i++){igood[i]=6;}
}/*endif*/

/*======================================================================*/
/*  III) Search the user defined data base                              */

if(myid==0){
  natm_srch = 2;
  if(strlen(filename_parse->user_inter_name) != 0){
    nsearch = 1;
    ibase_want = 1;
    count_data_base(filename_parse->user_inter_name,fun_dict,num_fun_dict,
                    &nbase,ibase_want);
    if(nbase>0){
      nbase2 = 2*nbase;
      inter_base  = (DATA_BASE_ENTRIES *)
                       cmalloc(nbase2*sizeof(DATA_BASE_ENTRIES))-1;
      cinter_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
      read_data_base(filename_parse->user_inter_name,fun_dict,num_fun_dict,
                     inter_base,cinter_base,ibase_want,nbase);
      search_base(nbase,nbase2,cinter_base,ninter,cinter,igood,ifound,
                  isearch,nsearch,natm_srch,filename_parse->user_inter_name);

      assign_base_inter(inter_base,nbase,ifound,ninter,
                        sig,eps,awill,bwill,cwill,rm_swit,c6m,c8m,c9m,c10m,
                        inter_label,interact->cutoff,interact->cutoff_res,
                        interact->cutti,
                        isearch,nsearch,cinter,cinter_base);
      cfree(&inter_base[1]);
      cfree(&cinter_base[1]);
    }/*endif*/
  }/*endif*/
}/*endif*/
/*======================================================================*/
/*  IV) Search the default defined data base                            */

if(myid==0){
  if(strlen(filename_parse->def_inter_name) != 0){
    nsearch = 2;
    ibase_want = 1;
    count_data_base(filename_parse->def_inter_name,fun_dict,num_fun_dict,
                    &nbase,ibase_want);
    if(nbase>0){
      nbase2 = 2*nbase;
      inter_base = (DATA_BASE_ENTRIES *)
                    cmalloc(nbase*sizeof(DATA_BASE_ENTRIES))-1;
      cinter_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
      read_data_base(filename_parse->def_inter_name,fun_dict,num_fun_dict,
                     inter_base,cinter_base,ibase_want,nbase);
      search_base(nbase,nbase2,cinter_base,ninter,cinter,igood,ifound,
                  isearch,nsearch,natm_srch,filename_parse->def_inter_name);

      assign_base_inter(inter_base,nbase,ifound,ninter,
                        sig,eps,awill,bwill,cwill,rm_swit,c6m,c8m,c9m,c10m,
                        inter_label,interact->cutoff,interact->cutoff_res,
                        interact->cutti,
                        isearch,nsearch,cinter,cinter_base);
      cfree(&inter_base[1]);
      cfree(&cinter_base[1]);
    }/*endif*/
  }/*endif*/
}/*endif*/

/*======================================================================*/
/* V) Check list for missing entries                                    */

if(myid==0){
  strcpy(typ,"inter");
  atmlst_not_found(ninter,cinter,ifound,natm_srch,typ);
}/*endif*/

/*======================================================================*/
/*Find unique values of epsilon,sigma,rcut for LJ and null interactions */
/* This involves rearranging all the intermolecular interactions.       */

if(myid==0){

 for(i=1; i<= ninter; i++){
   temp_cutoff[i]     = interact->cutoff[i];
   temp_cutoff_res[i] = interact->cutoff_res[i];
   temp_cutti[i]      = interact->cutti[i];
 }/*endfor*/

#define BYPASS_OFF
#ifdef BYPASS

 ninter_unique = ninter;
 for(i=1;i<=ninter;i++){
   interact->inter_map_index[i] = i;
 }/*endfor*/

#else
 sort_inter_params(eps,sig,
                   awill,bwill,cwill,
                   rm_swit,
                   c6m,c8m,c9m,c10m,
                   temp_cutoff,temp_cutoff_res,temp_cutti,
                   inter_label,interact->inter_map_index,
                   &ninter_unique,ninter);
#endif

}/*endif myid*/

if(num_proc > 1){ Bcast(&ninter_unique,1,MPI_INT,0,comm);}

 ninter_unique_mall = ninter_unique;
 if((ninter_unique_mall!=0)&&((ninter_unique_mall %2)==0))
    {ninter_unique_mall++;}

/*======================================================================*/
/* V) Broadcast the parameters                                          */

 if(num_proc>1){
   Bcast(&(sig[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(eps[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(awill[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(bwill[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(cwill[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(rm_swit[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(c6m[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(c8m[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(c9m[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(c10m[1]),ninter_unique,MPI_DOUBLE,0,comm);

   Bcast(&(inter_label[1]),ninter_unique,MPI_INT,0,comm);

   Bcast(&(temp_cutoff[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(temp_cutoff_res[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(temp_cutti[1]),ninter_unique,MPI_DOUBLE,0,comm);

   Bcast(&(interact->cutoff[1]),ninter,MPI_DOUBLE,0,comm);
   Bcast(&(interact->cutoff_res[1]),ninter,MPI_DOUBLE,0,comm);
   Bcast(&(interact->cutti[1]),ninter,MPI_DOUBLE,0,comm);
   Bcast(&(interact->inter_map_index[1]),ninter,MPI_INT,0,comm);
 }/*endif*/

/*======================================================================*/
/* VI) Allocate spline arrays                                           */

  interact->nter_typ = ninter;
  interact->nsplin_tot = interact->nsplin*ninter;

  nsplin_mall_tot = interact->nsplin*ninter_unique;
  if((nsplin_mall_tot!=0)&&((nsplin_mall_tot % 2)==0)){nsplin_mall_tot += 1;}

  nsplin_mall = interact->nsplin;
  if((nsplin_mall!=0)&&((nsplin_mall % 2)==0)){nsplin_mall += 1;}

  interact->cv0    = (double *) cmalloc(nsplin_mall_tot*sizeof(double))-1;
  interact->cdv0   = (double *) cmalloc(nsplin_mall_tot*sizeof(double))-1;
  interact->cv0_c  = (double *) cmalloc(nsplin_mall_tot*sizeof(double))-1;
  interact->cdv0_c = (double *) cmalloc(nsplin_mall_tot*sizeof(double))-1;

  interact->vcut_coul  = (double *) cmalloc(ninter_unique_mall*sizeof(double))-1;
  interact->rmin_spl   = (double *) cmalloc(ninter_unique_mall*sizeof(double))-1;
  interact->dr_spl     = (double *) cmalloc(ninter_unique_mall*sizeof(double))-1;
  interact->dri_spl    = (double *) cmalloc(ninter_unique_mall*sizeof(double))-1;


/*======================================================================*/
/*          Assign mall variables                                 */

  interact->ninter_mall        = ninter_mall;
  interact->nsplin_mall_tot    = nsplin_mall_tot;

/*=======================================================================*/
/* VII) Set up the splines for the real-space intermolecular potential   */
/*     energy and forces                                                 */

/* Spline the intermolecular interaction */

 set_inter_splin(sig,eps,awill,bwill,cwill,rm_swit,c6m,c8m,c9m,c10m,
                  alp_ewd,ishift_pot,inter_label,
                  spline_parse,interact,
                  temp_cutti,temp_cutoff,
                  ninter_unique,ncharge,iperd,myid);


/*=======================================================================*/
/*  VIII) Get the long range correction                                   */

  interact->clong = 0.0;
  interact->clong_res = 0.0;      
  if(iperd == 3) {

    get_clong(natm_tot,natm_typ,ninter,iatm_typ,c6m,
              interact->inter_map_index,
             &(interact->clong),
             &(interact->clong_res),temp_cutoff,temp_cutoff_res,
             interact->iswit_vdw,interact->rheal_res);

  } /* endif */
  
  if(myid==0){
   printf("Dispersion long range parameter %.15g \n",interact->clong);
   if(int_res_ter==1){
    printf("Dispersion long range parameter(RESPA) %g \n",
         interact->clong_res);
   }/*endif*/
  }/*endif*/

/*=======================================================================*/
/*  IX) Free temporary memory                                            */

  cfree(&inter_label[1]);
  cfree(&eps[1]);
  cfree(&sig[1]);
  cfree(&awill[1]);
  cfree(&bwill[1]);
  cfree(&cwill[1]);
  cfree(&rm_swit[1]);
  cfree(&c6m[1]);
  cfree(&c8m[1]);
  cfree(&c9m[1]);
  cfree(&c10m[1]);
  cfree(&temp_cutoff[1]);
  cfree(&temp_cutoff_res[1]);
  cfree(&temp_cutti[1]);

  if(myid==0){
    cfree(&fun_dict[1]);
  }/*endif*/
  cfree(fun_key);
  cfree(&cinter[1]);
  cfree(&ifound[1]);
  cfree(&isearch[1]);
  cfree(&igood[1]);

/*=======================================================================*/
/* X) Write to screen                                                    */

/*----------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/


