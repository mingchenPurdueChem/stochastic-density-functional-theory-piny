/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: all-control-frag.c                             */
/*                                                                          */
/* This file provide all modified control modules for fragmentation	    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"

#include "../proto_defs/proto_mol_params_entry.h"
#include "../proto_defs/proto_mol_params_local.h"
#include "../proto_defs/proto_inter_params_entry.h"
#include "../proto_defs/proto_inter_params_local.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_search_entry.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_communicate_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_communicate_local.h"
#include "../proto_defs/proto_cp_ewald_entry.h"
#include "../proto_defs/proto_cp_ewald_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_vps_params_entry.h"
#include "../proto_defs/proto_vps_params_local.h"

#include "../proto_defs/proto_interface_frag_local.h"

#define DEBUG_CLUS_CORR_OFF
#define CHECK_CLUS_OFF

#define MAX_INT 12.0
#define MIN3(A,B,C) (MIN(MIN((A),(B)),(C)))

#define ORIG_OFF
#define PME

#define DEBUG_DKNY_OFF

#define JUERG_FACTOR_ON
#ifdef  JUERG_FACTOR_ON
#define JUERG_FACTOR 0.72
#else
#define JUERG_FACTOR 1.0
#endif

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void controlInterParamsFrag(GENERAL_DATA *generalDataMini,CLASS *classMini,
	        CP *cpMini,BONDED *bondedMini,CP *cp,
	        SPLINE_PARSE *spline_parse,FILENAME_PARSE *filename_parse,
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
  int *iatm_typ = classMini->atommaps.iatm_atm_typ;
  NAME *atm_typ = classMini->atommaps.atm_typ;
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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void controlMolParamsFrag(CLASS *class,GENERAL_DATA *general_data,
                          BONDED *bonded,CP *cp,CLASS *classMini,
			  GENERAL_DATA *generalDataMini,BONDED *bondedMini,
			  CP *cpMini,CLASS_PARSE *classParse,
			  CP_PARSE *cpParse,FREE_PARSE *freeParse,
                          FILENAME_PARSE *filenameParse)

/*==========================================================================*/
    { /*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  CLATOMS_INFO *clatoms_info = &(classMini->clatoms_info);
  STODFTINFO *stodftInfo = cp->stodftInfo;
  CPCOEFFS_INFO *cpcoeffs_info = &(cpMini->cpcoeffs_info);
  DICT_MOL dict_mol;                        /* Dictionaries and sizes */
  DICT_MOL *dictMolAll;
  DICT_WORD *word;
  ATOMMAPS *atommaps = &(class->atommaps);

  int i,iii,jmol_typ,nmol_tot; 
  int nmol_typ,bond_free_num,bend_free_num,tors_free_num,rbar_sig_free_iopt;
  int ifirst,nline,nfun_key;
  int num_user,num_def,cp_on;
  int nhist,nhist_bar,nhist_sig,nfree,nsurf;
  int iperd          = generalDataMini->cell.iperd;
  int pi_beads       = classMini->clatoms_info.pi_beads;
  int np_forc        = classMini->communicate.np_forc;
  int numMolTyp	     = atommaps->nmol_typ;
  int iMolTyp;
  int iextend,ipress;
  int *mol_freeze_opt;
  int *imol_nhc_opt;
  int *nmol_jmol_typ;

  char *molsetname   = filenameParse->molsetname;
  char *fun_key;

  double now_memory;
  double *text_mol;
  double *text_nhc_mol;
  double *tot_memory = &(classMini->tot_memory);

  FILE *fp; 

  iextend = generalDataMini->ensopts.nvt 
       +generalDataMini->ensopts.npt_i  
       +generalDataMini->ensopts.npt_f
       +generalDataMini->ensopts.nst;
  ipress = generalDataMini->ensopts.npt_i
          +generalDataMini->ensopts.npt_f;

  cp_on = generalDataMini->simopts.cp_min 
         +generalDataMini->simopts.cp_wave_min 
         +generalDataMini->simopts.cp_wave_min_pimd
         +generalDataMini->simopts.cp
         +generalDataMini->simopts.cp_wave
         +generalDataMini->simopts.cp_pimd 
     +generalDataMini->simopts.cp_wave_pimd
         +generalDataMini->simopts.debug_cp
     +generalDataMini->simopts.debug_cp_pimd;

/*========================================================================*/
/* 0) Output to screen */

/*========================================================================*/
/* I) Set up dictionaries and malloc temporaries                          */

  ifirst = 1;

  dictMolAll = (DICT_MOL*)cmalloc(numMolTyp*sizeof(DICT_MOL));

  for(iMolTyp=0;iMolTyp<numMolTyp;iMolTyp++){
    set_mol_dict(&dictMolAll[iMolTyp].mol_dict,&dictMolAll[iMolTyp].num_mol_dict,
                 iextend,classParse->tau_nhc_def,
                 generalDataMini->statepoint.t_ext,ifirst);    
  }

  set_molset_fun_dict(&dict_mol.fun_dict,&dict_mol.num_fun_dict);
  /*
  set_mol_dict(&dict_mol.mol_dict,&dict_mol.num_mol_dict,
               iextend,classParse->tau_nhc_def,
               generalDataMini->statepoint.t_ext,ifirst);
  */
  set_wave_dict(&dict_mol.wave_dict,&dict_mol.num_wave_dict,cpParse);
  set_bond_free_dict(&dict_mol.bond_free_dict,&dict_mol.num_bond_free_dict);
  set_bend_free_dict(&dict_mol.bend_free_dict,&dict_mol.num_bend_free_dict);
  set_tors_free_dict(&dict_mol.tors_free_dict,&dict_mol.num_tors_free_dict); 
  set_rbar_free_dict(&dict_mol.rbar_free_dict,&dict_mol.num_rbar_free_dict); 
  set_user_base_dict(&dict_mol.user_base_dict,&dict_mol.num_user_base_dict);   
  set_def_base_dict(&dict_mol.def_base_dict,&dict_mol.num_def_base_dict);
  set_surf_dict(&dict_mol.surface_dict,&dict_mol.num_surface_dict);
  num_user = dict_mol.num_user_base_dict;
  num_def  = dict_mol.num_def_base_dict;

  word = (DICT_WORD *)cmalloc(sizeof(DICT_WORD));
  fun_key = (char *)cmalloc(MAXWORD*sizeof(char));

/*========================================================================*/
/* II) Read the moleset file and count the molecule types                 */
/*       and the free energy stuff                                        */
  
  nmol_typ            = 0;
  nline               = 0;
  nfun_key            = 0;
  bond_free_num       = 0;
  bend_free_num       = 0;
  tors_free_num       = 0;
  rbar_sig_free_iopt  = 0;
  nsurf               = 0;

  fp = cfopen(molsetname,"r");
  while(get_fun_key_cnt(fp,fun_key,&nline,&nfun_key,molsetname)){
    if(!strcasecmp(fun_key,"molecule_def")) {nmol_typ+=1;}
    if(!strcasecmp(fun_key,"bond_free_def")){bond_free_num+=1;}
    if(!strcasecmp(fun_key,"bend_free_def")){bend_free_num+=1;}
    if(!strcasecmp(fun_key,"tors_free_def")){tors_free_num+=1;}
    if(!strcasecmp(fun_key,"rbar_sig_free_def")){rbar_sig_free_iopt+=1;}
    if(!strcasecmp(fun_key,"surface_def")){nsurf+=1;}
  }/*endwhile*/
  fclose(fp);

  classMini->atommaps.nmol_typ    = nmol_typ;
  bondedMini->bond_free.num       = bond_free_num;
  bondedMini->bend_free.num       = bend_free_num;
  bondedMini->tors_free.num       = tors_free_num;
  bondedMini->rbar_sig_free.iopt  = rbar_sig_free_iopt ;
  bondedMini->rbar_sig_free.nfree = 0;

  classMini->surface.isurf_on = nsurf;

/*========================================================================*/
/* III) Malloc molecular data storage and initialize                      */
  
  now_memory = (double)((nmol_typ)*(sizeof(int)*2))*1.e-06;
  (*tot_memory) += now_memory;
  
  classMini->atommaps.nmol_jmol_typ  = (int *)cmalloc(nmol_typ*sizeof(int))-1;
  classMini->atommaps.mol_typ        = (NAME *)cmalloc(nmol_typ*sizeof(NAME))-1;
  classMini->atommaps.nres_1mol_jmol_typ  = (int *)cmalloc(nmol_typ*sizeof(int))-1;
  classMini->atommaps.jres_jmol_typ_strt  = (int *)cmalloc(nmol_typ*sizeof(int))-1;
  classMini->atommaps.icons_jres_jmol_typ = (int *)cmalloc(nmol_typ*sizeof(int))-1;
  classParse->ionfo_opt         = (int *)cmalloc(nmol_typ*sizeof(int))-1;
  classParse->ires_bond_conv    = (int *)cmalloc(nmol_typ*sizeof(int))-1;
  classParse->tau_nhc_mol       = (double*)cmalloc(nmol_typ*sizeof(double))-1;
  classParse->text_nhc_mol      = (double*)cmalloc(nmol_typ*sizeof(double))-1;
  classParse->imol_nhc_opt      = (int *)cmalloc(nmol_typ*sizeof(int))-1;  
  classParse->mol_freeze_opt    = (int *)cmalloc(nmol_typ*sizeof(int))-1;  
  classParse->mol_hydrog_mass_opt = (int *)cmalloc(nmol_typ*sizeof(int))-1;
  classParse->mol_hydrog_mass_val = 
                               (double *)cmalloc(nmol_typ*sizeof(double))-1;
  classParse->mol_hydrog_con_opt= (int *)cmalloc(nmol_typ*sizeof(int))-1;
  filenameParse->mol_param_name = (NAME *)cmalloc(nmol_typ*sizeof(NAME))-1;
  filenameParse->user_intra_name= (NAME *)cmalloc(num_user*sizeof(NAME))-1;
  filenameParse->def_intra_name = (NAME *)cmalloc(num_def*sizeof(NAME))-1;
  if(bond_free_num>0){
    bondedMini->bond_free.file       = (char *)cmalloc(MAXWORD*sizeof(char));
    freeParse->imoltyp_bond_free= (int *) cmalloc(2*sizeof(int))-1;
    freeParse->imol_bond_free   = (int *) cmalloc(2*sizeof(int))-1;
    freeParse->ires_bond_free   = (int *) cmalloc(2*sizeof(int))-1;
    freeParse->iatm_bond_free   = (int *) cmalloc(2*sizeof(int))-1;
  }
  if(bend_free_num>0){
    bondedMini->bend_free.file       = (char *)cmalloc(MAXWORD*sizeof(char));
    freeParse->imoltyp_bend_free= (int *) cmalloc(3*sizeof(int))-1;
    freeParse->imol_bend_free   = (int *) cmalloc(3*sizeof(int))-1;
    freeParse->ires_bend_free   = (int *) cmalloc(3*sizeof(int))-1;
    freeParse->iatm_bend_free   = (int *) cmalloc(3*sizeof(int))-1;
  }
  if(tors_free_num>0){
    bondedMini->tors_free.file       = (char *)cmalloc(MAXWORD*sizeof(char));
    freeParse->imoltyp_tors_free= (int *) cmalloc(8*sizeof(int))-1;
    freeParse->imol_tors_free   = (int *) cmalloc(8*sizeof(int))-1;
    freeParse->ires_tors_free   = (int *) cmalloc(8*sizeof(int))-1;
    freeParse->iatm_tors_free   = (int *) cmalloc(8*sizeof(int))-1;
  }

  if(rbar_sig_free_iopt>0){
    bondedMini->rbar_sig_free.file   = (char *)cmalloc(MAXWORD*sizeof(char));
  }
  classMini->clatoms_info.text_mol   = (double*)cmalloc(nmol_typ*sizeof(double))-1;

  text_mol       = classMini->clatoms_info.text_mol;
  text_nhc_mol   = classParse->text_nhc_mol;
  mol_freeze_opt = classParse->mol_freeze_opt;
  imol_nhc_opt   = classParse->imol_nhc_opt;

/*========================================================================*/
/* IV) Get molecular/CP setup information */

  controlSetMolParamsFrag(cpParse,classParse,filenameParse,freeParse,
                bondedMini,classMini,cpMini,generalDataMini,class,cp,
                general_data,&dict_mol,word,fun_key,&nfun_key,
                iextend,generalDataMini->statepoint.t_ext,ifirst,pi_beads,
                dictMolAll);

  /*
  controlSetMolParamsFrag(&(class->atommaps),&(cp->cpopts),
                         &(cp->cpcoeffs_info),cp_parse,class_parse,
                         bonded,&(class->surface),
                         filename_parse,free_parse,
	     clatoms_info,
                         &dict_mol,word,
                         fun_key,&nfun_key,iextend,
                         general_data->statepoint.t_ext,
                         ifirst,pi_beads,dictMolAll);
  */
  //exit(0);
  if(tors_free_num==1){
    tors_free_num = bondedMini->tors_free.num; /* changed in above routine */
  }/*endif*/

  for(i=1; i<= nmol_typ; i++){
    text_mol[i] = text_nhc_mol[i];   
  }/*endfor*/

/*========================================================================*/
/* V) Free some memory and malloc some other                              */

  if(bond_free_num>0){
    nhist = (bondedMini->bond_free.nhist);
    bonded->bond_free.hist = (double *)cmalloc(nhist*sizeof(double))-1;
  }/*endif*/

  if(bend_free_num>0){
    nhist = (bondedMini->bend_free.nhist);
    bonded->bend_free.hist = (double *)cmalloc(nhist*sizeof(double))-1;
  }/*endif*/

  if(tors_free_num==1){
    nhist = (bondedMini->tors_free.nhist);
    bonded->tors_free.hist = (double *)cmalloc(nhist*sizeof(double))-1;
  }/*endif*/
  if(tors_free_num==2){
    nhist = (bondedMini->tors_free.nhist);
    bonded->tors_free.hist_2d = cmall_mat(1,nhist,1,nhist);
  }/*endif*/

  if(rbar_sig_free_iopt>0){
    nhist_bar = (bondedMini->rbar_sig_free.nhist_bar);
    nhist_sig = (bondedMini->rbar_sig_free.nhist_sig);
    nfree     = (bondedMini->rbar_sig_free.nfree);
    bondedMini->rbar_sig_free.hist    = cmall_mat(1,nhist_bar,1,nhist_sig);
    bondedMini->rbar_sig_free.hist_rn = cmall_mat(1,nfree,1,nhist_bar);
  }/*endif*/

  cfree(fun_key);
  cfree(word);
  cfree(&dict_mol.fun_dict[1]);
  //cfree(&dict_mol.mol_dict[1]);
  cfree(&dict_mol.wave_dict[1]);
  cfree(&dict_mol.bond_free_dict[1]);
  cfree(&dict_mol.bend_free_dict[1]);
  cfree(&dict_mol.tors_free_dict[1]);
  cfree(&dict_mol.rbar_free_dict[1]);
  cfree(&dict_mol.user_base_dict[1]);
  cfree(&dict_mol.def_base_dict[1]);
  for(iMolTyp=0;iMolTyp<numMolTyp;iMolTyp++){
    cfree(&(dictMolAll[iMolTyp].mol_dict[1]));
  }
  cfree(dictMolAll);

  for(jmol_typ=1;jmol_typ<=nmol_typ;jmol_typ++){
    if((mol_freeze_opt[jmol_typ]>0)&&(ipress>0)){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No freezing allowed with NPT_I and NPT_F       \n");
      printf("Ensembles                                      \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }//endif
    if((mol_freeze_opt[jmol_typ]==1)&&(imol_nhc_opt[jmol_typ]>0)){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Do not use thermostating when freezing the     \n");
      printf("whole molecule.(Use none option)               \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
  }/*endfor*/

/*========================================================================*/
/* Error check */
 
  nmol_tot = 0;
  nmol_jmol_typ  = classMini->atommaps.nmol_jmol_typ;
  for(i=1;i<=nmol_typ;i++)nmol_tot += nmol_jmol_typ[i];
  if(nmol_tot<np_forc){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");   
    printf("Number of force level processors %d is greater than\n",np_forc);
    printf("the total number of molecules %d\n",nmol_tot);
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

  if(cp_on==1&&cpMini->cpopts.cp_dual_grid_opt>=1){
    if(cpParse->cp_ecut_dual_grid>cpParse->cp_ecut){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");   
      printf("The small dense grid cutoff is less than the large sparse");
      printf("grid cutoff. This might work, but I doubut it\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }
  }

/*------------------------------------------------------------------------*/
   }/*end routine*/
/*===========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void controlSetMolParamsFrag(CP_PARSE *cp_parse,CLASS_PARSE *class_parse,
			     FILENAME_PARSE *filename_parse,FREE_PARSE *free_parse,
			     BONDED *bondedMini,CLASS *classMini,CP *cpMini,
			     GENERAL_DATA *generalDataMini,CLASS *class,CP *cp,
			     GENERAL_DATA *general_data,DICT_MOL *dict_mol,
			     DICT_WORD *word,char *fun_key,int *nfun_key,
			     int iextend,double t_ext,int ifirst,int pi_beads,
			     DICT_MOL *dictMolAll)
/*=======================================================================*/
{/*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                               */

  ATOMMAPS *atommaps	    = &(classMini->atommaps);
  CPOPTS   *cpopts  = &(cpMini->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info	= &(cpMini->cpcoeffs_info);
  SURFACE *surface  = &(classMini->surface);
  CLATOMS_INFO *clatoms_info	= &(classMini->clatoms_info);
  ATOMMAPS *atommapsMacro   = &(class->atommaps);
  DICT_MOL *dictMolFrag;
  STODFTINFO *stodftInfo    = cp->stodftInfo;
  FRAGINFO *fragInfo	    = stodftInfo->fragInfo;

  int nline,nkey,i,num,ierr;
  int iMol,iType,iKey;
  int nmol_typ            = atommaps->nmol_typ;
  int num_fun_dict        = dict_mol->num_fun_dict;
  int numMolDict	  = dictMolAll[0].num_mol_dict;
  int bond_free_num       = bondedMini->bond_free.num;
  int bend_free_num       = bondedMini->bend_free.num;
  int tors_free_num       = bondedMini->tors_free.num;
  int rbar_sig_free_iopt  = bondedMini->rbar_sig_free.iopt;
  int nbar_bond;
  int countMol;
  int numMolTyp	      = atommapsMacro->nmol_typ;
  int iFrag   = fragInfo->iFrag;
  int numMolTypeNow   = fragInfo->numMolTypeFrag[iFrag];
  int typeInd;

  FILE *fp;
  int *mol_ind_chk;                                 /* Index check      */
  int *ifound;
  int *nmol_jmol_typ      = atommaps->nmol_jmol_typ;
  int *nres_1mol_jmol_typ = atommaps->nres_1mol_jmol_typ;
  int *imoltyp_rbar1_free;
  int *imoltyp_rbar2_free;
  int *imol_rbar1_free;
  int *imol_rbar2_free;
  int *ires_rbar1_free;
  int *ires_rbar2_free;
  int *imoltyp_bond_free = free_parse->imoltyp_bond_free;
  int *imol_bond_free     = free_parse->imol_bond_free;
  int *ires_bond_free    = free_parse->ires_bond_free;
  int *iatm_bond_free    = free_parse->iatm_bond_free;
  int *imoltyp_bend_free = free_parse->imoltyp_bend_free;
  int *imol_bend_free     = free_parse->imol_bend_free;
  int *ires_bend_free    = free_parse->ires_bend_free;
  int *iatm_bend_free    = free_parse->iatm_bend_free;
  int *imoltyp_tors_free = free_parse->imoltyp_tors_free;
  int *imol_tors_free     = free_parse->imol_tors_free;
  int *ires_tors_free    = free_parse->ires_tors_free; 
  int *iatm_tors_free    = free_parse->iatm_tors_free;
  int *molTypeFragNow	 = fragInfo->molTypeFrag[iFrag];
  int *molNumTypeFragNow = fragInfo->molNumTypeFrag[iFrag];
  int *numElecUpFragProc = fragInfo->numElecUpFragProc;
  int *numElecDnFragProc = fragInfo->numElecDnFragProc;
  
  char *molsetname        = filename_parse->molsetname;

/*=======================================================================*/
/* 0) Set up molecular index checking memory                             */

  mol_ind_chk     = (int *) cmalloc(nmol_typ*sizeof(int))-1;
  ifound          = (int *) cmalloc(num_fun_dict*sizeof(int))-1;
  for(i=1;i<=nmol_typ;i++){ mol_ind_chk[i]=0;}
  for(i=1;i<=num_fun_dict;i++){ ifound[i]=0;}

/*=======================================================================*/
/* I) Fetch a valid functional key word from molset file                 */

  fp = cfopen(molsetname,"r");

  nline          = 1;
  *nfun_key      = 0;
  word->iuset    = 0;
  word->key_type = 0;
  countMol = 0;
  while(get_fun_key(fp,fun_key,&nline,nfun_key,molsetname)){
    get_fun_key_index(fun_key,dict_mol->num_fun_dict,dict_mol->fun_dict,
                      nline,*nfun_key,molsetname,&num);

/*=======================================================================*/
/* II) Fetch the key words and key args of the functional key word       */
/*     and stick them into the correct dictionary                        */

    /*
    if(num==1){
      set_mol_dict(&(dict_mol->mol_dict),&(dict_mol->num_mol_dict),
                 iextend,class_parse->tau_nhc_def,t_ext,ifirst);
    }
    */
    nkey=0;
    while(get_word(fp,word,&nline,&nkey,*nfun_key,molsetname)){
      switch(num){
      case 1:
        put_word_dict(word,dictMolAll[countMol].mol_dict,dictMolAll[countMol].num_mol_dict,
                      fun_key,nline,nkey,*nfun_key,molsetname);
        break;
      case 2:
        put_word_dict(word,dict_mol->wave_dict,
                      dict_mol->num_wave_dict,fun_key,nline,nkey,
                      *nfun_key,molsetname);
        break;
      case 3:
        put_word_dict(word,dict_mol->bond_free_dict,
                     dict_mol->num_bond_free_dict,fun_key,nline,
                     nkey,*nfun_key,molsetname);
        break;
      case 4:
        put_word_dict(word,dict_mol->bend_free_dict,
                    dict_mol->num_bend_free_dict,fun_key,nline,
                    nkey,*nfun_key,molsetname);
        break;
      case 5:
        put_word_dict(word,dict_mol->tors_free_dict,
                    dict_mol->num_tors_free_dict,fun_key,nline,
                    nkey,*nfun_key,molsetname);
        break;
      case 6:
        put_word_dict(word,dict_mol->def_base_dict,
                    dict_mol->num_def_base_dict,fun_key,nline,
                    nkey,*nfun_key,molsetname);
        break;
      case 7:
        put_word_dict(word,dict_mol->user_base_dict,
                    dict_mol->num_user_base_dict,fun_key,nline,
                    nkey,*nfun_key,molsetname);
        break;
      case 8:
        put_word_dict(word,dict_mol->rbar_free_dict,
                    dict_mol->num_rbar_free_dict,fun_key,nline,
                    nkey,*nfun_key,molsetname);
        break;
      case 9:
        put_word_dict(word,dict_mol->surface_dict,
                    dict_mol->num_surface_dict,fun_key,nline,
                    nkey,*nfun_key,molsetname);
        break;
      } /* end switch dictionary fills*/
    } /* endwhile fetching keywords and keyargs*/
    if(num==1) countMol += 1;
  
/*=====================================================================*/
/* III) Assign the key args of the key words to the appropriate        */
/*      program variables                                              */

    ifound[num]=1;
    switch(num){
    /*
    case 1:
      set_mol_params(filename_parse,fun_key,dict_mol->mol_dict,
                   dict_mol->num_mol_dict,
                   class_parse,atommaps,mol_ind_chk,pi_beads);
      break;
    */
    case 2:
      sprintf(dict_mol->wave_dict[1].keyarg,"%i",numElecUpFragProc[iFrag]);
      sprintf(dict_mol->wave_dict[2].keyarg,"%i",numElecDnFragProc[iFrag]);
      set_wave_params(molsetname,fun_key,
                    dict_mol->wave_dict,dict_mol->num_wave_dict,
                    cpopts,cpcoeffs_info,cp_parse);
      break;      
    case 3:
      set_bond_free_params(molsetname,fun_key,
                        dict_mol->bond_free_dict,
                        dict_mol->num_bond_free_dict,
                        &(bondedMini->bond_free),free_parse,
                        atommaps->nmol_typ);
      break;
    case 4:
      set_bend_free_params(molsetname,fun_key,
                        dict_mol->bend_free_dict,
                        dict_mol->num_bend_free_dict,
                        &(bondedMini->bend_free),free_parse,
                        atommaps->nmol_typ);
      break;
      
    case 5:
      set_tors_free_params(molsetname,fun_key,
                        dict_mol->tors_free_dict,
                        dict_mol->num_tors_free_dict,
                        &(bondedMini->tors_free),free_parse,
                        atommaps->nmol_typ);
      break;
    case 6:
      set_def_base_params(filename_parse,dict_mol->def_base_dict,
                       dict_mol->num_def_base_dict);
      break;
    case 7:
      set_user_base_params(filename_parse,dict_mol->user_base_dict,
                        dict_mol->num_user_base_dict);
      break;
    case 8:
      set_rbar_free_params(molsetname,fun_key,
                        dict_mol->rbar_free_dict,
                        dict_mol->num_rbar_free_dict,
                        &(bondedMini->rbar_sig_free),free_parse,
                        atommaps->nmol_typ);
      break;
    case 9:
      set_surf_params(molsetname,fun_key,
                      dict_mol->surface_dict,
                      dict_mol->num_surface_dict,
                      surface); 
      break;
    } /* end switch assigning dict stuff to variables */
  } /*endwhile getting functional keywords data*/

/*=====================================================================*/
/* III) Assign mol_dict variables	       */

  dictMolFrag = (DICT_MOL*)cmalloc(numMolTypeNow*sizeof(DICT_MOL));

  for(iType=0;iType<numMolTypeNow;iType++){
    set_mol_dict(&dictMolFrag[iType].mol_dict,&dictMolFrag[iType].num_mol_dict,
                 iextend,class_parse->tau_nhc_def,
                 generalDataMini->statepoint.t_ext,1);
  }


  for(iType=0;iType<numMolTypeNow;iType++){
    typeInd = molTypeFragNow[iType];
    sprintf(dictMolFrag[iType].mol_dict[1].keyarg,"%i",iType+1);
    sprintf(dictMolFrag[iType].mol_dict[2].keyarg,"%i",molNumTypeFragNow[iType]);
    //printf("typeInd %i\n",typeInd);
    for(iKey=3;iKey<=numMolDict;iKey++){
      //printf("typeInd %i\n",typeInd);
      strcpy(dictMolFrag[iType].mol_dict[iKey].keyarg,dictMolAll[typeInd-1].mol_dict[iKey].keyarg);
    }
    for(iKey=1;iKey<=numMolDict;iKey++)dictMolFrag[iType].mol_dict[iKey].iuset = 1;
    //for(iKey=1;iKey<=numMolDict;iKey++)printf("%s %s\n",dictMolAll[typeInd-1].mol_dict[iKey].keyword,dictMolAll[typeInd-1].mol_dict[iKey].keyarg);
    set_mol_params(filename_parse,fun_key,dictMolFrag[iType].mol_dict,dictMolFrag[iType].num_mol_dict,
                   class_parse,atommaps,mol_ind_chk,pi_beads);
  }

/*=====================================================================*/
/* IV) Make sure everything that is required has been found            */ 
/*     if it hasn't either die or define defaults                      */ 

  if(ifound[1]==0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("  Required functional keyword %s not found     \n",
              dict_mol->fun_dict[1].keyword);
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/
  if(ifound[6]==0){
    set_def_base_params(filename_parse,dict_mol->def_base_dict,
                     dict_mol->num_def_base_dict);
  }/*endif*/
  if(ifound[7]==0){
    set_user_base_params(filename_parse,dict_mol->user_base_dict,
                      dict_mol->num_user_base_dict);
  }/*endif*/
  
/*=======================================================================*/
/* IV) Check  indices                                                    */

  for(i = 1;i<=nmol_typ;i++){
    if(mol_ind_chk[i]!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@  \n");
      printf("Molecule number %d specified %d times in file %s \n",
              i,mol_ind_chk[i],molsetname);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@  \n");
      fflush(stdout);
      exit(1);
    }/*endif*/
  } /*endfor*/

/*=======================================================================*/
/* V) Check  free energy indices                                         */

  ierr = 0;

  if(bond_free_num>0){
    for(i=1;i<=2;i++){
      if(imol_bond_free[i]>nmol_jmol_typ[imoltyp_bond_free[i]]){ierr=i;}
      if(ires_bond_free[i]>MAX(nres_1mol_jmol_typ[imoltyp_bond_free[i]],1))
                                                               {ierr=i;}
    }/*endfor*/
  }/*endif*/

  if(bend_free_num>0){
    for(i=1;i<=3;i++){
      if(imol_bend_free[i]>nmol_jmol_typ[imoltyp_bend_free[i]]){ierr=i;}
      if(ires_bend_free[i]>MAX(nres_1mol_jmol_typ[imoltyp_bend_free[i]],1))
                                                               {ierr=i;}
    }/*endfor*/
  }/*endif*/

  if(tors_free_num>0){tors_free_num=bondedMini->tors_free.num;}
  if(tors_free_num==1){
    for(i=1;i<=4;i++){
      if(imol_tors_free[i]>nmol_jmol_typ[imoltyp_tors_free[i]]){ierr=i;}
      if(ires_tors_free[i]>MAX(nres_1mol_jmol_typ[imoltyp_tors_free[i]],1))
                                                               {ierr=i;}
    }/*endfor*/
  }/*endif*/
  if(tors_free_num==2){
    for(i=1;i<=8;i++){
      if(imol_tors_free[i]>nmol_jmol_typ[imoltyp_tors_free[i]]){ierr=i;}
      if(ires_tors_free[i]>MAX(nres_1mol_jmol_typ[imoltyp_tors_free[i]],1))
                                                               {ierr=i;}
    }/*endfor*/
  }/*endif*/

  if(rbar_sig_free_iopt>0){
    imoltyp_rbar1_free = free_parse->imoltyp_rbar1_free;
    imoltyp_rbar2_free = free_parse->imoltyp_rbar2_free;
    imol_rbar1_free    = free_parse->imol_rbar1_free;
    imol_rbar2_free    = free_parse->imol_rbar2_free;
    ires_rbar1_free    = free_parse->ires_rbar1_free;
    ires_rbar2_free    = free_parse->ires_rbar2_free;
    nbar_bond          = free_parse->nbar_bond;

    for(i=1;i<=nbar_bond;i++){
     if(imol_rbar1_free[i]>nmol_jmol_typ[imoltyp_rbar1_free[i]]){ierr=1;}
     if(imol_rbar2_free[i]>nmol_jmol_typ[imoltyp_rbar2_free[i]]){ierr=2;}
     if(ires_rbar1_free[i]>
              MAX(nres_1mol_jmol_typ[imoltyp_rbar1_free[i]],1)){ierr=1;}
     if(ires_rbar2_free[i]>
              MAX(nres_1mol_jmol_typ[imoltyp_rbar2_free[i]],1)){ierr=2;}
    }/*endfor*/
  }/*endif*/

  if(ierr>0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Error specifing the %dth atom in a free energy def \n",ierr);
      printf("in set up file %s \n",molsetname);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*endif*/ 

/*=======================================================================*/
/* VI) Check for stochastic DFT                                          */
/*     (move to control_mol_params) */

/*=======================================================================*/
/* VII) Free memory                                                      */

  cfree(&mol_ind_chk[1]);
  cfree(&ifound[1]);
  cfree(dictMolFrag);

/*-----------------------------------------------------------------------*/
   }/*end routine*/
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void controlIntraParamsFrag(double *tot_memory,CLASS *classMini,
		    GENERAL_DATA *generalDataMini,BONDED *bonded,//actually bondedMini
                    FILENAME_PARSE *filename_parse,
                    FREE_PARSE *free_parse,CLASS_PARSE *class_parse,
                    NULL_INTER_PARSE *null_inter_parse)

/*========================================================================*/
/*     Begin routine                                                      */
{/*begin routine*/

/*========================================================================*/
/*     Local Variables                                                    */
  CLATOMS_INFO *clatoms_info = &(classMini->clatoms_info);
  CLATOMS_POS *clatoms_pos = classMini->clatoms_pos;
  GHOST_ATOMS *ghost_atoms = &(classMini->ghost_atoms);
  ATOMMAPS *atommaps = &(classMini->atommaps);
  SIMOPTS *simopts = &(generalDataMini->simopts);
  COMMUNICATE *communicate = &(classMini->communicate);

  DICT_INTRA dict_intra;
  BUILD_INTRA build_intra;
  RESBOND_PARSE resbond_parse;
  START_INDEX start_index;

  int i,iii;
  int jmol_typ;
  int ifirst,nresidue,mol_or_res;
  int nmol_tot,nres_bond;
  int nres_tot,ncon_tot,nmass_unphys;
  int myid = communicate->myid;
  int isurf_on = classMini->surface.isurf_on;

  double cpu1,cpu2;

  char *filename,*fun_key;
  FILE *fp;

  int pi_beads = clatoms_info->pi_beads;

/*=======================================================================*/
/*   I) Start the routine                                                */

  cputime(&cpu1);


/*=======================================================================*/
/* II) Initialize/malloc intra stuff                                     */
/*      (init_intra_params.c)                                            */

  filename               = (char *)cmalloc(MAXWORD*sizeof(char));  
  fun_key                 = (char *)cmalloc(MAXWORD*sizeof(char));  

  init_intra_params(clatoms_info,ghost_atoms,atommaps,&build_intra,
                    bonded,null_inter_parse,&resbond_parse,
	    filename_parse);

/*=======================================================================*/
/* III) Set up the dictionaries                                          */
/*      (set_intra_dict.c)                                               */

  ifirst = 1;
  dict_intra.word         = (DICT_WORD *)cmalloc(sizeof(DICT_WORD));
  set_intra_fun_dict(&dict_intra.fun_dict,&dict_intra.num_fun_dict,ifirst);
  set_atm_dict(&dict_intra.atm_dict,&dict_intra.num_atm_dict,ifirst);
  set_intra_dict(&dict_intra.intra_dict,&dict_intra.num_intra_dict,ifirst);
  set_mol_name_dict(&dict_intra.mol_name_dict,&dict_intra.num_mol_name_dict,
                ifirst);
  set_res_name_dict(&dict_intra.res_name_dict,&dict_intra.num_res_name_dict,
                ifirst);
  set_res_def_dict(&dict_intra.res_def_dict,&dict_intra.num_res_def_dict,
                ifirst);
  set_res_bond_dict(&dict_intra.res_bond_dict,&dict_intra.num_res_bond_dict,
                ifirst);

  strcpy(atommaps->atm_typ[1],"HELP");

/*=======================================================================*/
/* V) Zero the list counters                                             */
 
   bonded->bond.npow             = 0;
   bonded->bond.ncon             = 0;
   null_inter_parse->nbond_nul   = 0;
   bonded->bend.npow             = 0;
   bonded->bend.ncon             = 0;
   null_inter_parse->nbend_nul   = 0;
   bonded->tors.npow             = 0;
   bonded->tors.ncon             = 0;
   bonded->tors.nimpr            = 0;
   null_inter_parse->ntors_nul   = 0;
   bonded->onfo.num              = 0;
   null_inter_parse->nonfo_nul   = 0;
   bonded->bend_bnd.num          = 0;
   clatoms_info->natm_tot        = 0;
   atommaps->nres_typ            = 0;
   atommaps->nfreeze             = 0;
   atommaps->natm_typ            = 0;
   ghost_atoms->nghost_tot       = 0;
   ghost_atoms->natm_comp_max    = 0;
   bonded->grp_bond_con.num_21   = 0;
   bonded->grp_bond_con.num_33   = 0;
   bonded->grp_bond_con.num_43   = 0;
   bonded->grp_bond_con.num_23   = 0;
   bonded->grp_bond_con.num_46   = 0;
   bonded->grp_bond_con.ntyp_21  = 0;
   bonded->grp_bond_con.ntyp_33  = 0;
   bonded->grp_bond_con.ntyp_43  = 0;
   bonded->grp_bond_con.ntyp_23  = 0;
   bonded->grp_bond_con.ntyp_46  = 0;
   bonded->grp_bond_watts.num_33 = 0;
   bonded->grp_bond_watts.ntyp_33= 0;

/*=======================================================================*/
/* IV) Loop over molecular parameter files                               */

  for(jmol_typ=1;jmol_typ<=atommaps->nmol_typ;jmol_typ++){

/*-----------------------------------------------------------------------*/
/*  0) Write to the screen                                               */

/*-----------------------------------------------------------------------*/
/* 1) Store the present list counter values                              */

    start_index.nbond_pow    = bonded->bond.npow;
    start_index.nbond_con    = bonded->bond.ncon;
    start_index.nbond_nul    = null_inter_parse->nbond_nul;
    start_index.nbend_pow    = bonded->bend.npow;
    start_index.nbend_con    = bonded->bend.ncon;
    start_index.nbend_nul    = null_inter_parse->nbend_nul;
    start_index.ntors_pow    = bonded->tors.npow;
    start_index.ntors_con    = bonded->tors.ncon;
    start_index.ntors_nul    = null_inter_parse->ntors_nul;
    start_index.nonfo        = bonded->onfo.num;
    start_index.nonfo_nul    = null_inter_parse->nonfo_nul;
    start_index.nbend_bnd    = bonded->bend_bnd.num;
    start_index.natm         = clatoms_info->natm_tot;
    start_index.nfreeze      = atommaps->nfreeze;
    start_index.nghost_tot   = ghost_atoms->nghost_tot;
    start_index.ngrp_21      = bonded->grp_bond_con.num_21;
    start_index.ngrp_33      = bonded->grp_bond_con.num_33;
    start_index.ngrp_43      = bonded->grp_bond_con.num_43;
    start_index.ngrp_23      = bonded->grp_bond_con.num_23;
    start_index.ngrp_46      = bonded->grp_bond_con.num_46;
    start_index.ngrp_watt_33 = bonded->grp_bond_watts.num_33;
   
/*------------------------------------------------------------------------*/
/*  2) Count and error check the molecular parm file:                     */
/*     (fetch_residue.c)                                                  */

    strcpy(filename,filename_parse->mol_param_name[jmol_typ]);
    nresidue = atommaps->nres_1mol_jmol_typ[jmol_typ];/* spec in mol_set_file*/

    mol_or_res = 1;
    check_parmfile(filename,&(dict_intra.num_fun_dict),&(dict_intra.fun_dict),
                    fun_key,nresidue,&nres_bond,mol_or_res);
                                          /* in fetch_residue.c */
    resbond_parse.nres_bond = nres_bond; 
    resbond_parse.nresidue  = nresidue;
    resbond_parse_realloc(&resbond_parse);
                                          /* in manipulate_res_bonds.c */

/*-----------------------------------------------------------------------*/
/* 3) Read in the molecule name functional keyword: molecule type        */
/*     number of atoms and/or number of residues                         */
/*     (fetch_residue.c)                                                 */


    fetch_molname(filename,&dict_intra,atommaps,
                  fun_key,jmol_typ,nresidue);
/*-----------------------------------------------------------------------*/
/* 4) Read in the molecule residue definitions:                          */
/*       index,natom,parm_file,name                                      */
/*     (fetch_residue.c)                                                 */

    if(nresidue>0){
       fetch_residue_defs(filename,&dict_intra,
                          atommaps,filename_parse,
                          fun_key,nresidue,
                          jmol_typ,&build_intra);

    }else{
      fetch_residue_def0(atommaps,&build_intra,jmol_typ);  
      strcpy(filename_parse->res_param_name[1],
         filename_parse->mol_param_name[jmol_typ]);
        /* assigns the molecule type to the residue type 
       when there are no explicit residues           */
    }/*endif*/

/*-----------------------------------------------------------------------*/
/* 5) Read in the molecule resbonds then map them to the residues:       */
/*     Residue types(pairs),residue indices (pairs), bond sites(pairs),  */
/*     bond files(pairs), bond modifier and bond labels read in          */
/*     (fetch_residue.c and manipulate_res_bonds.c)                      */

    if(nres_bond>0){
       fetch_residue_bonds(filename,&dict_intra,fun_key,
                           resbond_parse.resbond_prm,
                           atommaps,nresidue,
                           nres_bond,jmol_typ,pi_beads); 
                    /* in fetch_residue.c */
       map_residue_bonds(&resbond_parse); /* in manipulate_res_bonds.c */
    }else{
       for(i=1;i<=MAX(nresidue,1);i++){
         resbond_parse.nres_bond1_jres[i]=0;resbond_parse.res_bond1_off[i]=0;
         resbond_parse.nres_bond2_jres[i]=0;resbond_parse.res_bond2_off[i]=0;
       }/*endfor*/  
    }/*endif*/
/*-----------------------------------------------------------------------*/
/* 6) Read in the residues: Set the atoms,bonds,bends,torsions,etc.      */
/*                          of each residue. Modificiations of residues  */
/*                          involved in bonds performed.                 */
/*                          Molecules with no residues treated here also */
/*     (control_res_params.c)                                           */


    control_res_params(tot_memory,clatoms_info,ghost_atoms,
                        atommaps,bonded,&resbond_parse,
	    &build_intra,
                        filename_parse,free_parse,class_parse,
                        null_inter_parse,filename,
	    &dict_intra,fun_key,jmol_typ);

/*------------------------------------------------------------------------*/
/* 7) Bond the residues together and create intramolecular interactions   */
/*    based on the molecular connectivity                                 */
/*    (residue_bond.c)                                                    */

    if(nres_bond>0){
       resbond_parse.ionfo  = class_parse->ionfo_opt[jmol_typ];
       resbond_parse.iconv  = class_parse->ires_bond_conv[jmol_typ];
       residue_bond(clatoms_info,clatoms_pos,
                    atommaps,bonded,&resbond_parse,
	    &build_intra,jmol_typ,
                    class_parse->mol_hydrog_con_opt[jmol_typ]);
    }/*endif*/

/*------------------------------------------------------------------------*/
/* 7.5) Print out progress in intramolecular connectivity lists           */

#ifdef DEBUG
    printf("End Fetch: nbonds=%d,nbends=%d,ntors=%d,nimpr=%d,nonfo=%d\n",
     bonded->bond.npow,bonded->bend_bnd.num,
     bonded->tors.npow,bonded->tors.nimpr,bonded->onfo.num);
#endif

/*------------------------------------------------------------------------*/
/* 7.7) Implement the freeze option                                       */

 fetch_freeze(class_parse,atommaps,&build_intra,&start_index,clatoms_info,
              jmol_typ);

/*------------------------------------------------------------------------*/
/* 7.7.5) Check for consistency with zero_com_vel                          */

  if(atommaps->nfreeze > 0 && class_parse->zero_com_vel==1){
    printf("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("Frozen atoms not compatible with zeroing the center of mass\n");
    printf("If you want to freeze atoms, set the zero_com_vel to `no'\n");
    printf("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

/*------------------------------------------------------------------------*/
/* 7.8) Implement the hydrog_mass option                                  */

 fetch_hydrog_mass(class_parse,atommaps,&build_intra,&start_index,
                   clatoms_info,jmol_typ);

/*-----------------------------------------------------------------------*/
/* 8) Replicate the molecule now that it has been constructed            */
/*    (replicate_mol.c)                                                    */

   if(atommaps->nmol_jmol_typ[jmol_typ] > 1){
     replicate_mol(clatoms_info,ghost_atoms,atommaps,&build_intra,bonded,
                   null_inter_parse,&start_index,jmol_typ);
   }/*endif*/

/*-----------------------------------------------------------------------*/

 }/*endfor:jmol_typ*/

/*=======================================================================*/
/*  Check for frozen atoms involved in constraints */

  freeze_con_check(atommaps,&(bonded->bond),&(bonded->bend),&(bonded->tors),
                   &(bonded->grp_bond_con));

/*=======================================================================*/
/*  DEBUG */

#ifdef DEBUG
  printf("Vomitting \n");mal_verify(1);
  vomit_intra_list(clatoms_info,ghost_atoms,atommaps,bonded,null_inter_parse);
#endif

/*=======================================================================*/
/*  VI)Tidy up                                                           */

 close_intra_params_frag(clatoms_info,clatoms_pos,ghost_atoms,atommaps,
                    &build_intra,bonded,null_inter_parse,tot_memory,
                    simopts,communicate->np_forc);

/*=======================================================================*/
/*  VIII) Output to screen                                                */

  cputime(&cpu2);

/*========================================================================*/
/* Check for physical masses */

  nmass_unphys = 0;
  for(i=1;i<=clatoms_info->natm_tot;i++){
    if(clatoms_info->mass[i]<900.0){nmass_unphys++;}
  }/*endfor*/

#define NONEXPERT
#ifdef NONEXPERT
  if(nmass_unphys>0){
    printf("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("There were %d atoms with masses less than 1/2 AMU\n",nmass_unphys);
    printf("This might be OK for experts, but not in general.\n");
    printf("If you are performing cp minimization, use the option\n");
    printf("class_mass_scale_fact in sim_run_def, instead.\n");
    printf("Expert users can disable this error by undefining the\n");
    printf("NONEXPERT ifdef in control_intra_params.c on line 358\n");
    printf("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/
#endif

/*========================================================================*/
/*  IX) Set the intramolecular potential params                          */

  setIntraPotentFrag(bonded,&build_intra,
               filename_parse->def_intra_name,
               filename_parse->user_intra_name);
  
#ifdef DEBUG
   printf("Vomitting \n");mal_verify(1);
   vomit_intra_potent(bonded,&build_intra);
#endif

/*========================================================================*/
/*  X) Free memory                */                                     
  
  cfree(fun_key);
  cfree(filename);

  cfree(&dict_intra.atm_dict[1]);
  cfree(&dict_intra.intra_dict[1]);
  cfree(&dict_intra.mol_name_dict[1]);
  cfree(&dict_intra.word[0]);
  cfree(&dict_intra.fun_dict[1]);

  cfree(&build_intra.mask_atm[1]);
  cfree(&build_intra.bond_site[1]);
  cfree(&build_intra.index_atm[1]);
  cfree(&build_intra.iatm_ind_chk[1]);
  cfree(&build_intra.cbond_typ_pow[1]);
  cfree(&build_intra.cbond_typ_con[1]);
  cfree(build_intra.cbond_typ_now);
  cfree(&build_intra.cbend_typ_pow[1]);
  cfree(&build_intra.cbend_typ_con[1]);
  cfree(build_intra.cbend_typ_now);
  cfree(&build_intra.ctors_typ_pow[1]);
  cfree(&build_intra.ctors_typ_con[1]);
  cfree(build_intra.ctors_typ_now);
  cfree(&build_intra.confo_typ[1]);
  cfree(build_intra.confo_typ_now); 
  cfree(&build_intra.cbend_bnd_typ[1]);
  cfree(build_intra.cbend_bnd_typ_now);

/*=======================================================================*/
/*  X) Write out synopsis                                               */

  nmol_tot=0;
  nres_tot=0;
  for(i=1;i<=atommaps->nmol_typ;i++){
    nmol_tot+=(atommaps->nmol_jmol_typ)[i];
    nres_tot+=(atommaps->nres_1mol_jmol_typ)[i]*
              (atommaps->nmol_jmol_typ)[i];
  }/*endfor*/
  atommaps->nres_tot = nres_tot;   

  /*
  printf("There are %d molecular units      \n",nmol_tot);
  printf("There are %d molecular unit types \n",atommaps->nmol_typ);
  printf("There are %d residue units        \n",nres_tot);
  printf("There are %d residue unit types   \n",atommaps->nres_typ);
  printf("There are %d atoms                \n",clatoms_info->natm_tot);
  printf("There are %d degrees of freedom   \n",clatoms_info->nfree);
  printf("There are %d ghost atoms          \n",ghost_atoms->nghost_tot);
  printf("There are %d atom types           \n",atommaps->natm_typ);
  printf("There are %d charged atoms        \n",clatoms_info->nchrg);
  printf("There are %d power series bonds   \n",bonded->bond.npow);
  printf("There are %d constrained bonds    \n",bonded->bond.ncon);
  printf("There are %d null bonds           \n",null_inter_parse->nbond_nul);
  printf("There are %d power series bends   \n",bonded->bend.npow);
  printf("There are %d constrained bends    \n",bonded->bend.ncon);
  printf("There are %d null bends           \n",null_inter_parse->nbend_nul);
  printf("There are %d Urey-Bradley bends   \n",bonded->bend_bnd.num);
  printf("There are %d power series torsions\n",bonded->tors.npow
                                               -bonded->tors.nimpr);
  printf("There are %d improper torsions    \n",bonded->tors.nimpr);
  printf("There are %d constrained torsions \n",bonded->tors.ncon);
  printf("There are %d null torsions        \n",null_inter_parse->ntors_nul);
  printf("There are %d lj onefours          \n",bonded->onfo.num);
  printf("There are %d null onefours        \n",null_inter_parse->nonfo_nul);
  printf("There are %d free energy bonds    \n",bonded->bond_free.num);
  printf("There are %d free energy bends    \n",bonded->bend_free.num);
  printf("There are %d free energy torsions \n",bonded->tors_free.num);
  printf("There are %d free energy rbar-sig \n",bonded->rbar_sig_free.nfree);
  printf("There are %d 21 group constraints \n",bonded->grp_bond_con.num_21);
  printf("There are %d 23 group constraints \n",bonded->grp_bond_con.num_23);
  printf("There are %d 33 group constraints \n",bonded->grp_bond_con.num_33);
  printf("There are %d 43 group constraints \n",bonded->grp_bond_con.num_43);
  printf("There are %d 46 group constraints \n",bonded->grp_bond_con.num_46);
  printf("There are %d 33 group Watts       \n",bonded->grp_bond_watts.num_33);
  printf("There is  %d surface              \n",isurf_on);
  printf("\n");
  */
  
  if((simopts->debug+simopts->debug_cp+simopts->debug_pimd)==1){
    printf("Enter an integer ");scanf("%d",&iii);
  }/*endif*/

/*=======================================================================*/
/*   V) If np_forc > 1 and there are non-group constraints, die.         */

#ifdef DEVELOP
  if(((communicate->np_forc)>1)&&
      (bonded->bond.ncon+bonded->bend.ncon+bonded->tors.ncon>0)){
    printf("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("Classical force parallel routine implemented for group\n");
    printf("constraints only.\n");
    printf("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/
#endif
  
/*=======================================================================*/
/*   V) Assign mall variables  */

/*=======================================================================*/
/*  XI) No atm minimization with constraints                             */

  ncon_tot = bonded->grp_bond_con.num_21
           + bonded->grp_bond_con.num_23
           + bonded->grp_bond_con.num_33
           + bonded->grp_bond_con.num_43
           + bonded->grp_bond_con.num_46
           + bonded->bond.ncon
           + bonded->bend.ncon
           + bonded->tors.ncon;
  if(((simopts->cp_min)==1)&&(ncon_tot>0)){
    printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("Atomic position minimization with constraints under CP \n");
    printf("not implemented\n");
    printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/ 

/*--------------------------------------------------------------------------*/
  }/*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  control_set_intra_potent:                                               */
/*==========================================================================*/

void setIntraPotentFrag(BONDED *bonded,BUILD_INTRA *build_intra,
                      NAME def_intra_name[],NAME user_intra_name[])

/*========================================================================*/
    { /* begin routine */
/*========================================================================*/
  /*            Local Variables                                             */
  int i,iii,ifirst,isum,n;
  int num_fun_dict;
  DICT_WORD *fun_dict;
  /*=======================================================================*/
  /*-----------------------------------------------------------------------*/
  /*=======================================================================*/
  /*     Control set_intra_potent:                                         */
  /*=======================================================================*/
  /* 0) Set up                                                             */

  ifirst =1;
  set_potfun_dict(&fun_dict,&num_fun_dict,ifirst);

  /*=======================================================================*/
  /* I) set the bond parameters                                            */
  isum = bonded->bond.ntyp_pow
       + bonded->bond.ntyp_con
       + bonded->grp_bond_con.ntyp_21
       + bonded->grp_bond_con.ntyp_23
       + bonded->grp_bond_con.ntyp_33
       + bonded->grp_bond_watts.ntyp_33
       + bonded->grp_bond_con.ntyp_43
       + bonded->grp_bond_con.ntyp_46;
  if(isum > 0){
    set_bond_potent(&(bonded->bond),&(bonded->grp_bond_con),
                        &(bonded->grp_bond_watts),build_intra,
                        def_intra_name[1],user_intra_name[1],
                        fun_dict,num_fun_dict);
  }/*endif*/

  /*=======================================================================*/
  /* II) set the bend parameters                                           */
  if((bonded->bend.ntyp_con) > 0){
    set_bend_potent(&(bonded->bend),build_intra,
                        def_intra_name[2],user_intra_name[2],
                        fun_dict,num_fun_dict);
  }/*endif*/

  /*=======================================================================*/
  /* IV) set the onfo parameters                                            */
  if((bonded->onfo.ntyp) > 0){
    set_onfo_potent(&(bonded->onfo),build_intra,
                        def_intra_name[4],user_intra_name[4],
                        fun_dict,num_fun_dict);
  }/*endif*/

  /*=======================================================================*/
  /* V) set the bend_bnd parameters                                        */
  if((bonded->bend_bnd.ntyp) > 0){
    n = bonded->bend_bnd.ntyp;
    build_intra->ibend_bnd_typ_pure = (int *)cmalloc(n*sizeof(int))-1;
    build_intra->ibend_bnd_typ_map = (int *)cmalloc(n*sizeof(int))-1;
    set_bend_bnd_potent(&(bonded->bend_bnd),build_intra,
                        def_intra_name[2],user_intra_name[2],
                        fun_dict,num_fun_dict);
    extract_pure_bends(&(bonded->bend_bnd),&(bonded->bend),build_intra);
    cfree(&(build_intra->ibend_bnd_typ_pure[1]));
    cfree(&(build_intra->ibend_bnd_typ_map[1]));
  }/*endif*/

  /*=======================================================================*/
  /* III) set the tors parameters                                          */
  if((bonded->tors.ntyp_pow+bonded->tors.ntyp_con) > 0){
    set_tors_potent(&(bonded->tors),build_intra,
                        def_intra_name[3],user_intra_name[3],
                        fun_dict,num_fun_dict);
  }/*endif*/

/*=======================================================================*/
  cfree(&fun_dict[1]);

/*-----------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Control ewald/cp g-space initialization */
/*==========================================================================*/

void controlSetCpEwaldFrag(GENERAL_DATA *generalDataMini,CLASS *classMini,
			   CP *cpMini,BONDED *bondedMini,CP *cp,CLASS *class,
			   GENERAL_DATA *general_data,BONDED *bonded,
                           CP_PARSE *cp_parse,CLASS_PARSE *class_parse)
/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*************************************************************************/
/* Our FFT grid is numGridDim[2](c)*numGridDim[1](b)*numGridDim[0](a)	 */
/* The k space coef number should be 2n+1 for each dimension, from -n to */
/* n. and n=numGridDim[i]-1. The total k space coef number is ((2na+1)*	 */
/* (2nb+1)*(2nc+1)-1)/2. Wf and density share the same k space	     */
/*************************************************************************/
/*=======================================================================*/
/*            Local variable declarations:                               */

  CPOPTS *cpopts	= &(cpMini->cpopts);
  SIMOPTS *simopts	= &(generalDataMini->simopts);
  CELL *cell		= &(generalDataMini->cell);
  CPCOEFFS_INFO *cpcoeffs_info	= &(cpMini->cpcoeffs_info);
  EWALD *ewald		= &(generalDataMini->ewald);
  CPEWALD *cpewald	= &(cpMini->cpewald);
  PSEUDO *pseudo	= &(cpMini->pseudo);
  EWD_SCR *ewd_scr	= &(classMini->ewd_scr);
  ECOR *ecor		= &(bondedMini->ecor);
  STODFTINFO *stodftInfo    = cp->stodftInfo;
  FRAGINFO *fragInfo	    = stodftInfo->fragInfo;
  
  int iFrag = fragInfo->iFrag;
  int idum1=1,idum2=1,idum3=1;
  int i,j,k,iii;                            /* Num: Debug tool                    */
  int nmall;
  int cp_on;                          /* Num: CP on flag                    */
  int nk1,nk2,nk3;                    /* Num: Number of k vec on small grid */
  int nkf1_dens_cp_box,nkf2_dens_cp_box,nkf3_dens_cp_box;               
		      /* Num: Number of k vec on large  grid*/
  int nkf1,nkf2,nkf3;                 /* Num: Number of k vec on mega  grid */

  int nktot,nktot_res,nktot_sm;       /* Num: Spherically cutoff k vecs     */
  int nktot_dens_cp_box;
  int ncoef,ncoef_dens_cp_box,ncoef_l;
  int ngrid_tot,nlen_pme;             /* Num: PME sizes                     */
  int ngrid_a_res,ngrid_b_res,ngrid_c_res;
  int ngrid_a, ngrid_b, ngrid_c;
  int pme_b_opt;
  int iperd	= cell->iperd;
  int box_rat	    = cpewald->box_rat;
  int kmax_ewd	    = class_parse->kmax_ewd;
  int kmax_res	    = class_parse->kmax_res;
  int cp_lsda	    = cpopts->cp_lsda;
  int cp_dual_grid_opt_on = cpopts->cp_dual_grid_opt;
  int numGridFragProc = fragInfo->numGridFragProc[iFrag];

  int *kmaxv;                         /* Lst: K-vector ranges               */
  int *kmax_cp_tmp,*kmaxv_res,cp_on_tmp;
  int *kmax_cp;
  int *kmaxv_dens_cp_box;
  int *kmax_cp_dens_cp_box;
  int *numGridDim = fragInfo->numGridFragDim[iFrag];

  double ecut_now;                    /* Num: Energy cutoff                 */
  double ecut_dens_cp_box_now;        /* Num: Energy cutoff                 */
  double ecut_res,ecut_tmp;
  double ecut_lg;                     /* Num: Energy cutoff for dens        */
  double ecut_sm;                     /* Num: Energy cutoff for wf          */
  double deth,deth_cp,side;           /* Num: Volumes and sizes             */
  double now_mem;                     /* Num: Memory used here              */
  double gmin_spl_tmp,gmin_true_tmp;  /* Num : Min/Max g-vectors            */
  double gmax_spl_tmp;
  double dbox_rat  =   cpewald->dbox_rat;

  double *gmin_true = &(pseudo->gmin_true);
  double *gmin_spl  = &(pseudo->gmin_spl);
  double *gmax_spl  = &(pseudo->gmax_spl);
  double *bfact_r, *bfact_i;
  double *hmati_ewd;                  /* Num: Inverse ewald h matrix        */
  double *hmati_ewd_cp;               /* Num: Inverse cp ewald h matrix     */
  double *hmat_ewd    = cell->hmat_ewd;
  double *hmat_ewd_cp = cell->hmat_ewd_cp;

/*=======================================================================*/
/* 0) Output to screen                                                   */

/*=======================================================================*/
/* I) Set cp switch and initialize respa kvectors                        */
 
  cp_on = 1;
  // int_res_ter==0
  ewald->nktot_res=0;
  ecor->nktot_res=0;

/*=======================================================================*/
/* II) Allocate simple memory                                            */

  hmati_ewd      = (double *) cmalloc((size_t)9*sizeof(double))-1;
  hmati_ewd_cp   = (double *) cmalloc((size_t)9*sizeof(double))-1;
  kmaxv          =    (int *) cmalloc((size_t)3*sizeof(int))-1;
  kmaxv_res      =    (int *) cmalloc((size_t)3*sizeof(int))-1;
  kmax_cp_tmp    =    (int *) cmalloc((size_t)3*sizeof(int))-1;

  cpewald->kmax_cp = (int *) cmalloc((size_t)3*sizeof(int))-1;
  kmax_cp          = cpewald->kmax_cp;
  cpewald->kmax_cp_dens_cp_box = (int *) cmalloc((size_t)3*sizeof(int))-1;
  kmax_cp_dens_cp_box          = cpewald->kmax_cp_dens_cp_box;
  if(cp_dual_grid_opt_on >= 1){ 
    kmaxv_dens_cp_box = (int *) cmalloc((size_t)3*sizeof(int))-1;
  }/*endif*/     


/*==========================================================================*/
/* III) Get inverse cell matrix and convert ewald_alpha                     */

  gethinv(hmat_ewd,hmati_ewd,&deth,iperd);    
  gethinv(hmat_ewd_cp,hmati_ewd_cp,&deth_cp,iperd);    

  side  = pow(deth,(1.0/3.0));  
  (ewald->alp_ewd) /= side;
  
/*==========================================================================*/
/* IV) Calculate cutoff, count number k vectors, Malloc and Fill            */

/*----------------------------------------------------------------------*/
/* A) Calculate cutoff, count number k vectors, malloc and fill        */
/*    Without the dual box this is standard grid for cp/ewald          */
/*    With the dual box this is the small box calculation              */
   
  calc_cutoff(kmax_ewd,&ecut_now,&(cp_parse->cp_ecut),cp_on,
                kmax_cp,kmaxv,hmati_ewd_cp,deth_cp);  

  kmaxv[1] = numGridDim[0]/2-1;
  kmaxv[2] = numGridDim[1]/2-1;
  kmaxv[3] = numGridDim[2]/2-1;
  //printf("kmaxv %i %i %i\n",kmaxv[1],kmaxv[2],kmaxv[3]);

  countkvec3d(&(ewald->nktot),ecut_now,kmaxv,hmati_ewd_cp);

  nktot                   = ewald->nktot;
  cpcoeffs_info->ecut     = ecut_now;
  cpcoeffs_info->ncoef_l  = nktot+1;
  ncoef_l                 = nktot+1;
  ecor->ecut              = 4.0*ecut_now;
  ewald->ecut             = 4.0*ecut_now;
  ewald->nkc_max          = kmaxv[3];

  //ecor->ecut              = 4.0*ecut_now;
  //ewald->ecut             = 4.0*ecut_now;
  //ecor->ecut = bonded->ecor.ecut; // I don't think I need ecor but just to prevent segfault
  //ewald->ecut = general_data->ewald.ecut;
  //ecut_now = 1.0e30; // A big number so that all grids included

  //kmaxv[1] = numGridDim[2]/2-1; 
  //kmaxv[2] = numGridDim[1]/2-1;
  //kmaxv[3] = numGridDim[0]/2-1;
  //old cubic
  //ewald->nkc_max = kmaxv[3];
  kmax_cp[1] = kmaxv[1];
  kmax_cp[2] = kmaxv[2];
  kmax_cp[3] = kmaxv[3];
  //nktot = ((kmaxv[1]*2+1)*(kmaxv[2]*2+1)*(kmaxv[3]*2+1)-1)/2;
  //ewald->nktot = nktot;
  //cpcoeffs_info->ncoef_l = nktot+1;
  kmax_cp_dens_cp_box[1] = kmax_cp[1];
  kmax_cp_dens_cp_box[2] = kmax_cp[2];
  kmax_cp_dens_cp_box[3] = kmax_cp[3];

  //printf("kmax_cp %i %i %i\n",kmax_cp[1],kmax_cp[2],kmax_cp[3]);
  //printf("nktot %i\n",nktot);



/*----------------------------------------------------------------------*/
/* A.1) For dualing : Calculate cutoff and count kvectors for the large */
/*      box and save the small box.                                     */
 
/*----------------------------------------------------------------------*/
/* B) Malloc                                                            */

  nmall =  nktot+1;   if((nmall % 2)==0){nmall++;}
  ewald->nktot_mall = nmall;

  ewald->kastr = (int *) cmalloc(nmall*sizeof(int))-1;
  ewald->kbstr = (int *) cmalloc(nmall*sizeof(int))-1;
  ewald->kcstr = (int *) cmalloc(nmall*sizeof(int))-1;
  ewald->ibrk1 = (int *) cmalloc(nmall*sizeof(int))-1;
  ewald->ibrk2 = (int *) cmalloc(nmall*sizeof(int))-1;

/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/* C) Fill                                                                */
  setkvec3d(nktot,ecut_now,kmaxv,hmati_ewd,
            ewald->kastr,ewald->kbstr,ewald->kcstr,
            ewald->ibrk1,ewald->ibrk2,cp_on,
            gmin_spl,gmin_true,gmax_spl);

  /* 
  for(i=1;i<=nktot;i++){
    printf("ka %i kb %i kc %i\n",ewald->kastr[i],ewald->kbstr[i],ewald->kcstr[i]);
  }

  exit(0);
  */


/*------------------------------------------------------------------------*/
/* C) Fill DENS_CP_BOX                                                    */

/*=======================================================================*/
/* V) Setup PME                                                          */

/*=======================================================================*/
/* VI) Set up RESPA k vectors and RESPA PME                              */

/*=======================================================================*/
/* VII) Setup CP coefs and k vectors                                     */

  if(cp_on == 1) {

/*--------------------------------------------------------------------*/
/*  A)  Count the k-vectors                                           */
    
    ecut_sm = ecut_now;
    countkvec3d_sm(&(cpewald->nktot_sm),ecut_sm,
		   kmax_cp_dens_cp_box,hmati_ewd_cp);
    nktot_sm = cpewald->nktot_sm;
    cpcoeffs_info->ncoef   = nktot_sm+1;
    ncoef                  = nktot_sm+1;
    printf("ecut_sm %lg ncoef %i cpcoeffs_info->ncoef_l %i\n",
            ecut_sm,ncoef,cpcoeffs_info->ncoef_l);
    
    /*    
    ecut_sm = ecut_now;  
    nktot_sm = nktot;
    cpewald->nktot_sm = nktot;
    cpcoeffs_info->ncoef = nktot_sm+1;
    ncoef = nktot_sm+1;
    */
/*--------------------------------------------------------------------*/
/*  B)  Malloc                                                       */
    nmall =  (nktot_sm+1);if((nmall % 2)==0){nmall++;}
    cpewald->nktot_cp_sm_mall = nmall;
    cpewald->kastr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
    cpewald->kbstr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
    cpewald->kcstr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
    cpewald->ibrk1_sm = (int *) cmalloc(nmall*sizeof(int))-1;
    cpewald->ibrk2_sm = (int *) cmalloc(nmall*sizeof(int))-1;

    nmall =  ncoef; if((nmall % 2)==0){nmall++;}
    cpcoeffs_info->cmass = (double *)cmalloc(nmall*sizeof(double))-1;
/*--------------------------------------------------------------------*/
/*  C)  Fill and check                                                */

    setkvec3d_sm(nktot_sm,ecut_sm,kmax_cp_dens_cp_box,hmati_ewd_cp,
                 cpewald->kastr_sm,cpewald->kbstr_sm,cpewald->kcstr_sm,
                 cpewald->ibrk1_sm,cpewald->ibrk2_sm,
                 &(cpewald->gw_gmin),&(cpewald->gw_gmax));
    cpewald->gCutoffKe = cp->cpewald.gCutoffKe;

    if(cp_dual_grid_opt_on == 0 && cp_on == 1){
      check_kvec(ewald->nktot,ewald->kastr,ewald->kbstr,ewald->kcstr,nktot_sm,
                  cpewald->kastr_sm,cpewald->kbstr_sm,cpewald->kcstr_sm);
    }
/*--------------------------------------------------------------------*/
/*  D) Set up the cp masses                                           */

    set_cpmass(ncoef,cpewald->kastr_sm,
               cpewald->kbstr_sm,cpewald->kcstr_sm,
               cpcoeffs_info->cmass,hmati_ewd_cp,
               &(cp_parse->cp_mass_tau_def),cp_parse->cp_mass_cut_def,
               &(cpcoeffs_info->icmass_unif));
   }/*endif:cpon*/

/*=======================================================================*/
/* VIII) Output time has arrived                                         */

/*=======================================================================*/
/* IX) Free excess memory                                                */

   cfree(&(hmati_ewd)[1]);
   cfree(&(hmati_ewd_cp)[1]);
   cfree(&(kmaxv)[1]);
   cfree(&(kmaxv_res)[1]);
   cfree(&(kmax_cp_tmp)[1]);

/*-----------------------------------------------------------------------*/ 
  }/* end routine */
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Control ewald/cp g-space initialization */
/*==========================================================================*/

void controlSetCpEwaldFragSparse(GENERAL_DATA *generalDataMini,CLASS *classMini,
			   CP *cpMini,BONDED *bondedMini,CP *cp,CLASS *class,
			   GENERAL_DATA *general_data,BONDED *bonded,
                           CP_PARSE *cp_parse,CLASS_PARSE *class_parse)
/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*************************************************************************/
/* Our FFT grid is numGridDim[2](c)*numGridDim[1](b)*numGridDim[0](a)	 */
/* The k space coef number should be 2n+1 for each dimension, from -n to */
/* n. and n=numGridDim[i]-1. The total k space coef number is ((2na+1)*	 */
/* (2nb+1)*(2nc+1)-1)/2. Wf and density share the same k space	     */
/*************************************************************************/
/*=======================================================================*/
/*            Local variable declarations:                               */

  CPOPTS *cpopts	= &(cpMini->cpopts);
  SIMOPTS *simopts	= &(generalDataMini->simopts);
  CELL *cell		= &(generalDataMini->cell);
  CPCOEFFS_INFO *cpcoeffs_info	= &(cpMini->cpcoeffs_info);
  EWALD *ewald		= &(generalDataMini->ewald);
  CPEWALD *cpewald	= &(cpMini->cpewald);
  PSEUDO *pseudo	= &(cpMini->pseudo);
  EWD_SCR *ewd_scr	= &(classMini->ewd_scr);
  ECOR *ecor		= &(bondedMini->ecor);
  STODFTINFO *stodftInfo    = cp->stodftInfo;
  FRAGINFO *fragInfo	    = stodftInfo->fragInfo;
  
  int iFrag = fragInfo->iFrag;
  int idum1=1,idum2=1,idum3=1;
  int i,j,k,iii;                            /* Num: Debug tool                    */
  int nmall;
  int cp_on;                          /* Num: CP on flag                    */
  int nk1,nk2,nk3;                    /* Num: Number of k vec on small grid */
  int nkf1_dens_cp_box,nkf2_dens_cp_box,nkf3_dens_cp_box;               
		      /* Num: Number of k vec on large  grid*/
  int nkf1,nkf2,nkf3;                 /* Num: Number of k vec on mega  grid */

  int nktot,nktot_res,nktot_sm;       /* Num: Spherically cutoff k vecs     */
  int nktot_dens_cp_box;
  int ncoef,ncoef_dens_cp_box,ncoef_l;
  int ngrid_tot,nlen_pme;             /* Num: PME sizes                     */
  int ngrid_a_res,ngrid_b_res,ngrid_c_res;
  int ngrid_a, ngrid_b, ngrid_c;
  int pme_b_opt;
  int iperd	= cell->iperd;
  int box_rat	    = cpewald->box_rat;
  int kmax_ewd	    = class_parse->kmax_ewd;
  int kmax_res	    = class_parse->kmax_res;
  int cp_lsda	    = cpopts->cp_lsda;
  int cp_dual_grid_opt_on = cpopts->cp_dual_grid_opt;
  int numGridFragProc = fragInfo->numGridFragProc[iFrag];

  int *kmaxv;                         /* Lst: K-vector ranges               */
  int *kmax_cp_tmp,*kmaxv_res,cp_on_tmp;
  int *kmax_cp,*kmax_rho;
  int *kmaxv_dens_cp_box;
  int *kmax_cp_dens_cp_box;
  int *kmax_dummy;
  int *numGridDim = fragInfo->numGridFragDim[iFrag];

  double ecut_now;                    /* Num: Energy cutoff                 */
  double ecut_rho;
  double ecut_dens_cp_box_now;        /* Num: Energy cutoff                 */
  double ecut_res,ecut_tmp;
  double ecut_lg;                     /* Num: Energy cutoff for dens        */
  double ecut_sm;                     /* Num: Energy cutoff for wf          */
  double deth,deth_cp,side;           /* Num: Volumes and sizes             */
  double now_mem;                     /* Num: Memory used here              */
  double gmin_spl_tmp,gmin_true_tmp;  /* Num : Min/Max g-vectors            */
  double gmax_spl_tmp;
  double dbox_rat  =   cpewald->dbox_rat;

  double *gmin_true = &(pseudo->gmin_true);
  double *gmin_spl  = &(pseudo->gmin_spl);
  double *gmax_spl  = &(pseudo->gmax_spl);
  double *bfact_r, *bfact_i;
  double *hmati_ewd;                  /* Num: Inverse ewald h matrix        */
  double *hmati_ewd_cp;               /* Num: Inverse cp ewald h matrix     */
  double *hmat_ewd    = cell->hmat_ewd;
  double *hmat_ewd_cp = cell->hmat_ewd_cp;

/*=======================================================================*/
/* 0) Output to screen                                                   */

/*=======================================================================*/
/* I) Set cp switch and initialize respa kvectors                        */
 
  cp_on = 1;
  // int_res_ter==0
  ewald->nktot_res=0;
  ecor->nktot_res=0;

/*=======================================================================*/
/* II) Allocate simple memory                                            */

  hmati_ewd      = (double *) cmalloc((size_t)9*sizeof(double))-1;
  hmati_ewd_cp   = (double *) cmalloc((size_t)9*sizeof(double))-1;
  kmaxv          =    (int *) cmalloc((size_t)3*sizeof(int))-1;
  kmaxv_res      =    (int *) cmalloc((size_t)3*sizeof(int))-1;
  kmax_cp_tmp    =    (int *) cmalloc((size_t)3*sizeof(int))-1;
  kmax_dummy     =    (int *) cmalloc((size_t)3*sizeof(int))-1;

  cpewald->kmax_cp = (int *) cmalloc((size_t)3*sizeof(int))-1;
  cpewald->kmax_rho = (int *) cmalloc((size_t)3*sizeof(int))-1;
  kmax_cp          = cpewald->kmax_cp;
  kmax_rho         = cpewald->kmax_rho;
  cpewald->kmax_cp_dens_cp_box = (int *) cmalloc((size_t)3*sizeof(int))-1;
  kmax_cp_dens_cp_box          = cpewald->kmax_cp_dens_cp_box;
  if(cp_dual_grid_opt_on >= 1){ 
    kmaxv_dens_cp_box = (int *) cmalloc((size_t)3*sizeof(int))-1;
  }/*endif*/     


/*==========================================================================*/
/* III) Get inverse cell matrix and convert ewald_alpha                     */

  /*
  printf("hmat_ewd %lg %lg %lg hmat_ewd_cp %lg %lg %lg\n",
         hmat_ewd[1],hmat_ewd[5],hmat_ewd[9],
         hmat_ewd_cp[1],hmat_ewd_cp[5],hmat_ewd_cp[9]);
  */
  gethinv(hmat_ewd,hmati_ewd,&deth,iperd);    
  gethinv(hmat_ewd_cp,hmati_ewd_cp,&deth_cp,iperd);    

  side  = pow(deth,(1.0/3.0));  
  (ewald->alp_ewd) /= side;
  
/*==========================================================================*/
/* IV) Calculate cutoff, count number k vectors, Malloc and Fill            */

/*----------------------------------------------------------------------*/
/* A) Calculate cutoff, count number k vectors, malloc and fill        */
/*    Without the dual box this is standard grid for cp/ewald          */
/*    With the dual box this is the small box calculation              */
   
  calc_cutoff(kmax_ewd,&ecut_now,&(cp_parse->cp_ecut),cp_on,
              kmax_cp,kmax_dummy,hmati_ewd_cp,deth_cp);  
  // rho cutoff
  //printf("1111111 ecut_rho %lg\n",cpewald->eCutoffRho);
  calc_cutoff(kmax_ewd,&ecut_rho,&(cpewald->eCutoffRho),cp_on,
                kmax_rho,kmax_dummy,hmati_ewd_cp,deth_cp);
  //printf("ecut_now %lg ecut_rho %lg\n",ecut_now,ecut_rho);
  /*
  printf("ecut_now %lg kmax_cp %i %i %i kmaxv %i %i %i hmati_ewd_cp %lg %lg %lg\n",ecut_now,kmax_cp[1],kmax_cp[2],kmax_cp[3],kmaxv[1],kmaxv[2],kmaxv[3],hmati_ewd_cp[1],hmati_ewd_cp[5],hmati_ewd_cp[9]);
  */
  
  //printf("ecut %lg\n",cp_parse->cp_ecut);
  if(kmax_rho[1]>numGridDim[0]/2-1||kmax_rho[2]>numGridDim[1]/2-1||
     kmax_rho[3]>numGridDim[2]/2-1){

    kmax_rho[1] = numGridDim[0]/2-1;
    kmax_rho[2] = numGridDim[1]/2-1;
    kmax_rho[3] = numGridDim[2]/2-1;
  }
  kmaxv[1] = kmax_rho[1];
  kmaxv[2] = kmax_rho[2];
  kmaxv[3] = kmax_rho[3];
  // Just in case we have larger kmax_cp due to rounding error
  if(kmax_cp[1]>kmax_rho[1])kmax_cp[1] = kmax_rho[1];
  if(kmax_cp[2]>kmax_rho[2])kmax_cp[2] = kmax_rho[2];
  if(kmax_cp[3]>kmax_rho[3])kmax_cp[3] = kmax_rho[3];
  //printf("kmaxv %i %i %i\n",kmaxv[1],kmaxv[2],kmaxv[3]);

  countkvec3d_sm(&(ewald->nktot),ecut_rho,kmaxv,hmati_ewd_cp);

  nktot                   = ewald->nktot;
  cpcoeffs_info->ecut     = ecut_now;
  cpcoeffs_info->ncoef_l  = nktot+1;
  ncoef_l                 = nktot+1;
  ecor->ecut              = ecut_rho;
  ewald->ecut             = ecut_rho;
  ewald->nkc_max          = kmaxv[3];

  //ecor->ecut              = 4.0*ecut_now;
  //ewald->ecut             = 4.0*ecut_now;
  //ecor->ecut = bonded->ecor.ecut; // I don't think I need ecor but just to prevent segfault
  //ewald->ecut = general_data->ewald.ecut;
  //ecut_now = 1.0e30; // A big number so that all grids included

  //kmaxv[1] = numGridDim[2]/2-1; 
  //kmaxv[2] = numGridDim[1]/2-1;
  //kmaxv[3] = numGridDim[0]/2-1;
  //old cubic
  //ewald->nkc_max = kmaxv[3];
  //nktot = ((kmaxv[1]*2+1)*(kmaxv[2]*2+1)*(kmaxv[3]*2+1)-1)/2;
  //ewald->nktot = nktot;
  //cpcoeffs_info->ncoef_l = nktot+1;
  kmax_cp_dens_cp_box[1] = kmax_cp[1];
  kmax_cp_dens_cp_box[2] = kmax_cp[2];
  kmax_cp_dens_cp_box[3] = kmax_cp[3];

  //printf("kmax_cp %i %i %i\n",kmax_cp[1],kmax_cp[2],kmax_cp[3]);
  //printf("nktot %i\n",nktot);



/*----------------------------------------------------------------------*/
/* A.1) For dualing : Calculate cutoff and count kvectors for the large */
/*      box and save the small box.                                     */
 
/*----------------------------------------------------------------------*/
/* B) Malloc                                                            */

  nmall =  nktot+1;   if((nmall % 2)==0){nmall++;}
  ewald->nktot_mall = nmall;

  ewald->kastr = (int *) cmalloc(nmall*sizeof(int))-1;
  ewald->kbstr = (int *) cmalloc(nmall*sizeof(int))-1;
  ewald->kcstr = (int *) cmalloc(nmall*sizeof(int))-1;
  ewald->ibrk1 = (int *) cmalloc(nmall*sizeof(int))-1;
  ewald->ibrk2 = (int *) cmalloc(nmall*sizeof(int))-1;

/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/* C) Fill                                                                */
  setkvec3d_sm(nktot,ecut_rho,kmaxv,hmati_ewd,
            ewald->kastr,ewald->kbstr,ewald->kcstr,
            ewald->ibrk1,ewald->ibrk2,
            gmin_spl,gmax_spl);
  //printf("1111111111111111\n");

  /*
   for(i=1;i<=nktot;i++){
     printf("nkabcccccccc %i %i %i %i\n",i,
            ewald->kastr[i],ewald->kbstr[i],ewald->kcstr[i]);
   }
   fflush(stdout);
   exit(0);
  */

  ewald->kastr[nktot+1] = 0.0;
  ewald->kbstr[nktot+1] = 0.0;
  ewald->kcstr[nktot+1] = 0.0;

  *gmin_true = *gmin_spl;
  *gmin_spl *=0.75;
  *gmax_spl *= 4.0/3.0;

  /* 
  for(i=1;i<=nktot;i++){
    printf("ka %i kb %i kc %i\n",ewald->kastr[i],ewald->kbstr[i],ewald->kcstr[i]);
  }

  exit(0);
  */

/*------------------------------------------------------------------------*/
/* C) Fill DENS_CP_BOX                                                    */

/*=======================================================================*/
/* V) Setup PME                                                          */

/*=======================================================================*/
/* VI) Set up RESPA k vectors and RESPA PME                              */

/*=======================================================================*/
/* VII) Setup CP coefs and k vectors                                     */

  if(cp_on == 1) {

/*--------------------------------------------------------------------*/
/*  A)  Count the k-vectors                                           */
    
    ecut_sm = ecut_now;

    countkvec3d_sm(&(cpewald->nktot_sm),ecut_sm,
                   kmax_cp_dens_cp_box,hmati_ewd_cp);


    //cpewald->nktot_sm = nktot;
    //nktot_sm = nktot;
    nktot_sm = cpewald->nktot_sm;
    cpcoeffs_info->ncoef   = nktot_sm+1;
    ncoef                  = nktot_sm+1;
    
    
    //printf("ecut_sm %lg ncoef %i cpcoeffs_info->ncoef_l %i\n",
    //        ecut_sm,ncoef,cpcoeffs_info->ncoef_l);
    
    /*    
    ecut_sm = ecut_now;  
    nktot_sm = nktot;
    cpewald->nktot_sm = nktot;
    cpcoeffs_info->ncoef = nktot_sm+1;
    ncoef = nktot_sm+1;
    */
/*--------------------------------------------------------------------*/
/*  B)  Malloc                                                       */
    nmall =  (nktot_sm+1);if((nmall % 2)==0){nmall++;}
    cpewald->nktot_cp_sm_mall = nmall;
    cpewald->kastr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
    cpewald->kbstr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
    cpewald->kcstr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
    cpewald->ibrk1_sm = (int *) cmalloc(nmall*sizeof(int))-1;
    cpewald->ibrk2_sm = (int *) cmalloc(nmall*sizeof(int))-1;

    nmall =  ncoef; if((nmall % 2)==0){nmall++;}
    cpcoeffs_info->cmass = (double *)cmalloc(nmall*sizeof(double))-1;
/*--------------------------------------------------------------------*/
/*  C)  Fill and check                                                */

    setkvec3d_sm(nktot_sm,ecut_sm,kmax_cp_dens_cp_box,hmati_ewd_cp,
                 cpewald->kastr_sm,cpewald->kbstr_sm,cpewald->kcstr_sm,
                 cpewald->ibrk1_sm,cpewald->ibrk2_sm,
                 &(cpewald->gw_gmin),&(cpewald->gw_gmax));
    cpewald->gCutoffKe = cp->cpewald.gCutoffKe;

    if(cp_dual_grid_opt_on == 0 && cp_on == 1){
      check_kvec(ewald->nktot,ewald->kastr,ewald->kbstr,ewald->kcstr,nktot_sm,
                  cpewald->kastr_sm,cpewald->kbstr_sm,cpewald->kcstr_sm);
    }
/*--------------------------------------------------------------------*/
/*  D) Set up the cp masses                                           */

    set_cpmass(ncoef,cpewald->kastr_sm,
               cpewald->kbstr_sm,cpewald->kcstr_sm,
               cpcoeffs_info->cmass,hmati_ewd_cp,
               &(cp_parse->cp_mass_tau_def),cp_parse->cp_mass_cut_def,
               &(cpcoeffs_info->icmass_unif));
   }/*endif:cpon*/

/*=======================================================================*/
/* VIII) Output time has arrived                                         */

/*=======================================================================*/
/* IX) Free excess memory                                                */

   cfree(&(hmati_ewd)[1]);
   cfree(&(hmati_ewd_cp)[1]);
   cfree(&(kmaxv)[1]);
   cfree(&(kmaxv_res)[1]);
   cfree(&(kmax_cp_tmp)[1]);

   //fflush(stdout);
   //exit(0);
/*-----------------------------------------------------------------------*/ 
  }/* end routine */
/*=======================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void controlGroupCommunicatorsFrag(CLASS *class,CP *cp,int cp_on)

/*==========================================================================*/
{/* begin routine */
/*==========================================================================*/
/*          Local variable declarations                                     */

 int iii;
 int num_proc = class->communicate.np;
 int myid     = class->communicate.myid;
 int comm_compare_result;
 MPI_Comm world;

/*==========================================================================*/
/* I) Build different communicators                                        */

 if(num_proc>1){

   build_communicate_groups(&(class->communicate),cp_on);

 }else{

   class->communicate.myid_bead        = 0;
   class->communicate.myid_bead_prime  = 0;
   class->communicate.myid_state       = 0;
   class->communicate.myid_forc        = 0;
   class->communicate.myid_forc_source = 0;
   class->communicate.myid_forc_target = 0;

   /*
   Comm_dup(class->communicate.world,&class->communicate.comm_forc);
   Comm_dup(class->communicate.world,&class->communicate.comm_forc_source);
   Comm_dup(class->communicate.world,&class->communicate.comm_forc_target);
   Comm_dup(class->communicate.world,&class->communicate.comm_states);
   Comm_dup(class->communicate.world,&class->communicate.comm_beads);
   Comm_dup(class->communicate.world,&class->communicate.comm_faux);
   */

 }/*endif*/

/*==========================================================================*/
/* II) Duplicate communicate into force package */

 /*-----------------------------------------------------------------------*/
 /* i) Force level communicator */

  class->class_comm_forc_pkg.myid     = class->communicate.myid_forc;
  class->class_comm_forc_pkg.num_proc = class->communicate.np_forc;
  //Comm_dup(class->communicate.comm_forc,&class->class_comm_forc_pkg.comm);

  class->class_comm_forc_pkg.dbl_num_proc =
                            (double)(class->communicate.np_forc);

 /*-----------------------------------------------------------------------*/
 /* ii) Source-Force level communicator */

  class->class_comm_forc_pkg.plimpton_ind.myid_source
                          = class->communicate.myid_forc_source;
  class->class_comm_forc_pkg.plimpton_ind.num_proc_source
                          = class->communicate.np_forc_src;
  /*
  Comm_dup(class->communicate.comm_forc_source,
          &(class->class_comm_forc_pkg.plimpton_ind.source));
  */
 /*-----------------------------------------------------------------------*/
 /* iii) Target-Force level communicator */

  class->class_comm_forc_pkg.plimpton_ind.myid_target
                          = class->communicate.myid_forc_target;
  class->class_comm_forc_pkg.plimpton_ind.num_proc_target
                          = class->communicate.np_forc_trg;
  /*
  Comm_dup(class->communicate.comm_forc_target,
          &(class->class_comm_forc_pkg.plimpton_ind.target));
  */

/*==========================================================================*/
/* II) Set up some Bead level parallel stuff                                */

  class->clatoms_info.pi_beads_proc_st    =
                     (class->communicate.myid_bead_prime)*
                     (class->clatoms_info.pi_beads_proc)+1;

  class->clatoms_info.pi_beads_proc_end   =
                     (class->communicate.myid_bead_prime+1)*
                     (class->clatoms_info.pi_beads_proc);

  cp->cpcoeffs_info.pi_beads_proc_st  = class->clatoms_info.pi_beads_proc_st;
  cp->cpcoeffs_info.pi_beads_proc_end = class->clatoms_info.pi_beads_proc_end;

/*==========================================================================*/
/* III) Duplicate class->communicate into cp->communicate                   */
/*      then set up the cp_comm_package                                     */

 /*-------------------------------------------------------------------------*/
 /* i) Duplicate */

  cp->communicate.np             = class->communicate.np;
  cp->communicate.np_beads       = class->communicate.np_beads;
  cp->communicate.np_states      = class->communicate.np_states;
  cp->communicate.np_forc        = class->communicate.np_forc;
  cp->communicate.np_forc_src    = class->communicate.np_forc_src;
  cp->communicate.np_forc_trg    = class->communicate.np_forc_trg;
  cp->communicate.numThreads     = class->communicate.numThreads;

  cp->communicate.myid             = class->communicate.myid;
  cp->communicate.myid_bead        = class->communicate.myid_bead;
  cp->communicate.myid_bead_forc   = class->communicate.myid_bead_forc;
  cp->communicate.myid_state       = class->communicate.myid_state;
  cp->communicate.myid_forc        = class->communicate.myid_forc;
  cp->communicate.myid_forc_source = class->communicate.myid_forc_source;
  cp->communicate.myid_forc_target = class->communicate.myid_forc_target;


  if(cp->communicate.np>1){

#ifdef DUPLICATE_NOT_CORRECT
    Comm_dup(class->communicate.world,&(cp->communicate.world));
    Comm_dup(class->communicate.comm_beads,&(cp->communicate.comm_beads));
    Comm_dup(class->communicate.comm_beads_forc,
              &(cp->communicate.comm_beads_forc));
    Comm_dup(class->communicate.comm_states,&(cp->communicate.comm_states));
    Comm_dup(class->communicate.comm_forc,&(cp->communicate.comm_forc));
    Comm_dup(class->communicate.comm_forc_source,
              &(cp->communicate.comm_forc_source));
    Comm_dup(class->communicate.comm_forc_target,
              &(cp->communicate.comm_forc_target));
#endif

    cp->communicate.world            = class->communicate.world;
    cp->communicate.comm_beads       = class->communicate.comm_beads;
    cp->communicate.comm_beads_forc  = class->communicate.comm_beads_forc;
    cp->communicate.comm_states      = class->communicate.comm_states;
    cp->communicate.comm_forc        = class->communicate.comm_forc;
    cp->communicate.comm_forc_source = class->communicate.comm_forc_source;
    cp->communicate.comm_forc_target = class->communicate.comm_forc_target;


  }else{

    cp->communicate.world            = class->communicate.world;
    cp->communicate.comm_beads       = class->communicate.comm_beads;
    cp->communicate.comm_beads_forc  = class->communicate.comm_beads_forc;
    cp->communicate.comm_states      = class->communicate.comm_states;
    cp->communicate.comm_forc        = class->communicate.comm_forc;
    cp->communicate.comm_forc_source = class->communicate.comm_forc_source;
    cp->communicate.comm_forc_target = class->communicate.comm_forc_target;
    cp->communicate.comm_faux        = class->communicate.comm_faux;

  }/*endif : cp->communicate.np>1*/

 /*-------------------------------------------------------------------------*/
 /* ii) Build the cp package */

  //Comm_dup(class->communicate.world,&world);


  if(cp->cpcoeffs_info.iopt_cp_pw==1){
    build_cp_comm_pkg(cp,world); //Contains no src/target stuff
  }else if (cp->cpcoeffs_info.iopt_cp_dvr==1){
    build_cp_comm_pkg_dvr(cp,world);
  }else{
    printf("@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("Wrong basis set in building communication package\n");
    printf("@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }

/*==========================================================================*/
   }/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Control the setup of the FFT packages                                    */
/*==========================================================================*/

void controlFFTPkgFrag(GENERAL_DATA *generalDataMini,CLASS *classMini,CP *cpMini,
	    CP *cp)

/*=========================================================================*/
     {/*Begin Routine*/
/*=========================================================================*/
/*    Local Variables */

  PARA_FFT_PKG3D *cp_sclr_fft_pkg_sm = &(cpMini->cp_sclr_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_para_fft_pkg_sm = &(cpMini->cp_para_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg_dens_cp_box = &(cpMini->cp_sclr_fft_pkg3d_dens_cp_box);
  PARA_FFT_PKG3D *cp_para_fft_pkg_dens_cp_box = &(cpMini->cp_para_fft_pkg3d_dens_cp_box);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg_lg = &(cpMini->cp_sclr_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_para_fft_pkg_lg = &(cpMini->cp_para_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg_sparse = &(cpMini->cp_sclr_fft_pkg3d_sparse);
  PARA_FFT_PKG3D *cp_para_fft_pkg_sparse = &(cpMini->cp_para_fft_pkg3d_sparse);
  PARA_FFT_PKG3D *pme_fft_pkg = &(generalDataMini->pme_fft_pkg);
  PARA_FFT_PKG3D *pme_res_fft_pkg = &(generalDataMini->pme_res_fft_pkg);
  EWALD *ewald = &(generalDataMini->ewald);
  CPEWALD *cpewald = &(cpMini->cpewald);
  PART_MESH *part_mesh = &(classMini->part_mesh);
  CPCOEFFS_INFO *cpcoeffs_info = &(cpMini->cpcoeffs_info);
  COMMUNICATE *communicate = &(classMini->communicate);
  CPOPTS *cpopts = &(cpMini->cpopts);
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  

  int iFrag = fragInfo->iFrag;
  int nkf1,nkf2,nkf3,nfft_ext,iii;
  int nkf1_dens_cp_box,nkf2_dens_cp_box,nkf3_dens_cp_box;
          /* used for denisty on cp grid when have 2 boxes*/
  int cp_on	  = 1;  // Always do cp
  int cp_lsda	      = cpopts->cp_lsda;
  int cp_para_opt     = cpopts->cp_para_opt;
  int cp_dual_grid_opt_on = cpopts->cp_dual_grid_opt;
  int nstate_up          = cpcoeffs_info->nstate_up;
  int nstate_dn          = cpcoeffs_info->nstate_dn;
  int ncoef              = cpcoeffs_info->ncoef;
  int ncoef_dens_cp_box  = cpcoeffs_info->ncoef_l_dens_cp_box;
  int ncoef_l            = cpcoeffs_info->ncoef_l;
  int nktot              = ewald->nktot;
  int nktot_dens_cp_box  = cpewald->nktot_dens_cp_box;
  int box_rat            = cpewald->box_rat;
  int nktot_res          = ewald->nktot_res;
  int *kmax_cp           = cpewald->kmax_cp;
  int *kmax_cp_dens_cp_box = cpewald->kmax_cp_dens_cp_box;
  int myid               = communicate->myid;
  int np_states          = communicate->np_states;
  int myid_state         = communicate->myid_state;
  int pme_on             = part_mesh->pme_on;
  int pme_res_on         = part_mesh->pme_res_on;
  int myid_forc          = communicate->myid_forc;
  int np_forc            = communicate->np_forc;
  int ngrid_a            = part_mesh->ngrid_a;
  int ngrid_b            = part_mesh->ngrid_b;
  int n_interp           = part_mesh->n_interp;
  int n_interp_res       = part_mesh->n_interp_res;
  int pme_para_opt       = part_mesh->pme_para_opt;
  int *numGridDim    = fragInfo->numGridFragDim[iFrag];
  int realSparseOpt      = cpewald->realSparseOpt;

/*=========================================================================*/
/* 0) Print to screen and check for nproc > nstate error */

/*=========================================================================*/
/* 0.1) Set CP FFT Size  */

  nkf1 = numGridDim[0];
  nkf2 = numGridDim[1];
  nkf3 = numGridDim[2];
  //printf("nkf1 %i nkf2 %i nkf3 %i\n",nkf1,nkf2,nkf3);
  //nkf1 = 4*(kmax_cp_dens_cp_box[1]+1);
  //nkf2 = 4*(kmax_cp_dens_cp_box[2]+1);
  //nkf3 = 4*(kmax_cp_dens_cp_box[3]+1);

/*=========================================================================*/
/* I) DENS_CP_BOX CP scalar package                                        */

/*=========================================================================*/
/* II) DENSITY_CP_BOX  parallel package                                    */
/*       This package must be made for both hybrid and full_g options      */

/*=========================================================================*/
/* I) Large CP scalar package                                              */


 if(cp_on == 1 && cp_para_opt == 0){/* hybrid option */
    cp_sclr_fft_pkg_lg->nkf1       = nkf1;
    cp_sclr_fft_pkg_lg->nkf2       = nkf2;
    cp_sclr_fft_pkg_lg->nkf3       = nkf3;
      
    cp_sclr_fft_pkg_lg->nktot      = nktot;
    cp_sclr_fft_pkg_lg->ncoef      = ncoef_l;

    cp_sclr_fft_pkg_lg->myid       = 0;
    cp_sclr_fft_pkg_lg->myidp1     = 1;
    cp_sclr_fft_pkg_lg->num_proc   = 1;
    cp_sclr_fft_pkg_lg->comm       = communicate->comm_faux;

    cp_sclr_fft_pkg_lg->numThreads = 1;

    create_para_fft_pkg3d(cp_sclr_fft_pkg_lg,
                          ewald->kastr,ewald->kbstr,
                          ewald->kcstr,cp_dual_grid_opt_on);
  }/*endif*/


/*=========================================================================*/
/* II) Large CP parallel package                                           */

  if(cp_on == 1){
    /*
    nkf1 = 4*(kmax_cp[1]+1);
    nkf2 = 4*(kmax_cp[2]+1);
    nkf3 = 4*(kmax_cp[3]+1);
    */
    cp_para_fft_pkg_lg->nkf1       = nkf1;
    cp_para_fft_pkg_lg->nkf2       = nkf2;
    cp_para_fft_pkg_lg->nkf3       = nkf3;
      
    cp_para_fft_pkg_lg->nktot      = nktot;
    cp_para_fft_pkg_lg->ncoef      = ncoef_l;
   
    cp_para_fft_pkg_lg->myid       = myid_state;
    cp_para_fft_pkg_lg->myidp1     = myid_state+1;
    cp_para_fft_pkg_lg->num_proc   = np_states;
    cp_para_fft_pkg_lg->comm       = communicate->comm_states;

    create_para_fft_pkg3d(cp_para_fft_pkg_lg,
                          ewald->kastr,ewald->kbstr,
                          ewald->kcstr,cp_dual_grid_opt_on);
  }/*endif*/

/*=========================================================================*/
/* III) Small  CP scalar package                                           */

  if(cp_on == 1 && cp_para_opt == 0){/* hybrid option */
    /*
    nkf1 = 4*(kmax_cp_dens_cp_box[1]+1);
    nkf2 = 4*(kmax_cp_dens_cp_box[2]+1);
    nkf3 = 4*(kmax_cp_dens_cp_box[3]+1);
    */
    cp_sclr_fft_pkg_sm->nkf1       = nkf1;
    cp_sclr_fft_pkg_sm->nkf2       = nkf2;
    cp_sclr_fft_pkg_sm->nkf3       = nkf3;
      
    cp_sclr_fft_pkg_sm->nktot      = ncoef-1;
    cp_sclr_fft_pkg_sm->ncoef      = ncoef;

    cp_sclr_fft_pkg_sm->myid       = 0;
    cp_sclr_fft_pkg_sm->myidp1     = 1;
    cp_sclr_fft_pkg_sm->num_proc   = 1;
    cp_sclr_fft_pkg_sm->comm       = communicate->comm_faux;

    create_para_fft_pkg3d(cp_sclr_fft_pkg_sm,
                          cpewald->kastr_sm,cpewald->kbstr_sm,
                          cpewald->kcstr_sm,cp_dual_grid_opt_on);

  }/*endif*/

  if(cp_on==1&&cp_para_opt==0&&realSparseOpt==1){
    //printf("11111111111111111111111111111 initial fft");
    cp_sclr_fft_pkg_sparse->nkf1 = nkf1;
    cp_sclr_fft_pkg_sparse->nkf2 = nkf2;
    cp_sclr_fft_pkg_sparse->nkf3 = nkf3;
    cp_sclr_fft_pkg_sparse->nktot = ncoef-1;
    cp_sclr_fft_pkg_sparse->ncoef = ncoef;
    cp_sclr_fft_pkg_sparse->myid = 0;
    cp_sclr_fft_pkg_sparse->myidp1 = 1;
    cp_sclr_fft_pkg_sparse->num_proc = 1;
    cp_sclr_fft_pkg_sparse->comm = communicate->comm_faux;

    create_para_fft_pkg3d(cp_sclr_fft_pkg_sparse,
                          cpewald->kastr_sm,cpewald->kbstr_sm,
                          cpewald->kcstr_sm,cp_dual_grid_opt_on);

    cp_para_fft_pkg_sparse->nkf1 = nkf1;
    cp_para_fft_pkg_sparse->nkf2 = nkf2;
    cp_para_fft_pkg_sparse->nkf3 = nkf3;
    cp_para_fft_pkg_sparse->nktot = ncoef-1;
    cp_para_fft_pkg_sparse->ncoef = ncoef;
    cp_para_fft_pkg_sparse->myid = myid_state;
    cp_para_fft_pkg_sparse->myidp1 = myid_state+1;
    cp_para_fft_pkg_sparse->num_proc = np_states;
    cp_para_fft_pkg_sparse->comm = communicate->comm_states;

    create_para_fft_pkg3d(cp_para_fft_pkg_sparse,
                          cpewald->kastr_sm,cpewald->kbstr_sm,
                          cpewald->kcstr_sm,cp_dual_grid_opt_on);
  }
/*=========================================================================*/
/* IV) Small  CP parallel package                                         */

  if(cp_on == 1 && cp_para_opt == 1){/* full g option */
    /*
    nkf1 = 4*(kmax_cp_dens_cp_box[1]+1);
    nkf2 = 4*(kmax_cp_dens_cp_box[2]+1);
    nkf3 = 4*(kmax_cp_dens_cp_box[3]+1);
    */

    cp_para_fft_pkg_sm->nkf1       = nkf1;
    cp_para_fft_pkg_sm->nkf2       = nkf2;
    cp_para_fft_pkg_sm->nkf3       = nkf3;
      
    cp_para_fft_pkg_sm->nktot      = ncoef-1;
    cp_para_fft_pkg_sm->ncoef      = ncoef;

    cp_para_fft_pkg_sm->myid       = myid_state;
    cp_para_fft_pkg_sm->myidp1     = myid_state+1;
    cp_para_fft_pkg_sm->num_proc   = np_states;
    cp_para_fft_pkg_sm->comm       = communicate->comm_states;

    create_para_fft_pkg3d(cp_para_fft_pkg_sm,
                         cpewald->kastr_sm,cpewald->kbstr_sm,
                          cpewald->kcstr_sm,cp_dual_grid_opt_on);

  }/*endif*/

/*=========================================================================*/
/* V) PME package                                                         */

/*=========================================================================*/
/* VI) PME_RES package                                                      */

/*=========================================================================*/
/* VI) Output */

/*-------------------------------------------------------------------------*/
}/*end routine */
/*=========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void controlVpsParamsFrag(GENERAL_DATA *generalDataMini,CLASS *classMini,
	      CP *cpMini,FILENAME_PARSE *filename_parse,
	      SPLINE_PARSE *spline_parse,CP_PARSE *cp_parse)

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  PSEUDO *pseudo = &(cpMini->pseudo);
  CELL *cell = &(generalDataMini->cell);
  ATOMMAPS *atommaps = &(classMini->atommaps);
  CLATOMS_INFO *clatoms_info = &(classMini->clatoms_info);
  CPOPTS *cpopts = &(cpMini->cpopts);
  COMMUNICATE *communicate = &(classMini->communicate);
  CPCOEFFS_INFO *cpcoeffs_info = &(cpMini->cpcoeffs_info);
  
  CVPS *cvps_typ;
  VPS_FILE *vps_file;          /* Fle:  Pseudopotential file          */
  DICT_WORD *word;
  DICT_WORD *fun_dict;
  DICT_WORD *vps_dict,*vps_dict_tmp;
  NAME *atm_typ = atommaps->atm_typ;
  MPI_Comm comm = communicate->world;
 
  int natm_typ = atommaps->natm_typ;
  int natm_tot = clatoms_info->natm_tot;
  int natm_ab_init = clatoms_info->nab_initio;
  int cp_ptens_calc = cpopts->cp_ptens_calc;
  int cp_dual_grid_opt = cpopts->cp_dual_grid_opt;
  int cp_lsda = cpopts->cp_lsda;
  
  int i;                       /* Num:  For loop counters             */
  int ifound;                  /* Num:  Data base match flag          */
  int ishift,ishift2,ishift3;  /* Num:  Angular momentum shifts       */
  int ngh_now;                 /* Num:  Number of ngh points read in  */
  int ngh_max;                 /* Num:  Max number of ngh points      */
  int num_fun_dict;
  int num_vps_dict,ifirst,iii;
  int natm_typ_mall,natm_mall,nsplin_mall,norm_mall,nlist_mall,nsplin_mall_dvr;
  int nmall_gh,natm_typ_gh = 0;   /* starting malloc value for gauss-hermite */
  int iopt_cp_dvr = cpcoeffs_info->iopt_cp_dvr;
  int myid      = communicate->myid;
  int num_proc  = communicate->np;

  char *filename;              /* Char: temp file name                */
  char *fun_key;

  double now_mem;              /* Num:  Current memory usage          */
  double ecut_cp = cp_parse->cp_ecut;
  double dummy0,dummy1,dummy2,dummy3;
  double alpha_conv_dual = pseudo->alpha_conv_dual;
  double vol_cp   = cell->vol_cp;
  double *hmat     = cell->hmat;
  double grid_h;

/*==========================================================================*/
/* 0) Output                                                                */

/*==========================================================================*/
/* 0.5) Toast bad spline points : (note conversion back to Ry)              */

  if((pseudo->nsplin_g<4000)&&(2.0*ecut_cp>60.0)&&(iopt_cp_dvr== 0)){
     if(myid==0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Now, dude, lets clean up those files and use a \n");
       printf("reasonable number of psuedo spline points at large \n");
       printf("cutoffs. Its clear, %d points at %g Ry, won't do!\n",
               pseudo->nsplin_g,2.0*ecut_cp);
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
     }//endif
     exit(1);
  }//endif

/*==========================================================================*/
/* I) Convert alpha_conv_dual and Allocate some temporary character arrays  */

  alpha_conv_dual        /= (pow(vol_cp,1.0/3.0));
  pseudo->alpha_conv_dual = alpha_conv_dual;

  if(myid==0){
    filename_parse->vps_name = (NAME *) cmalloc(natm_typ*sizeof(NAME))-1;
    vps_file  = (VPS_FILE *) cmalloc(natm_typ*sizeof(VPS_FILE))-1;
    pseudo->vps_file = (VPS_FILE *) cmalloc(natm_typ*sizeof(VPS_FILE))-1;
    fun_key   = (char *)cmalloc(MAXWORD*sizeof(char));  
    filename  = (char *)cmalloc(MAXWORD*sizeof(char));  
    word      = (DICT_WORD *)cmalloc(sizeof(DICT_WORD))-1;  
    cvps_typ  = (CVPS *)cmalloc(sizeof(CVPS));  
    ifirst    = 1;
    set_potfun_dict(&fun_dict,&num_fun_dict,ifirst);
    set_potvps_dict(&vps_dict,&num_vps_dict,ifirst);      
    set_potvps_dict(&vps_dict_tmp,&num_vps_dict,ifirst);      
  }//endif


/*==========================================================================*/
/* III) Malloc up the vps stuff                                             */ 

  natm_typ_mall      = natm_typ;

  if( (natm_typ_mall % 2)==0){natm_typ_mall++;}

  pseudo->n_ang      = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
  pseudo->loc_opt    = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
  pseudo->ivps_label = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
  pseudo->iformat    = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
  pseudo->rcut_nl    = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
  pseudo->q_pseud    = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;

  pseudo->nrad_0 = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
  pseudo->nrad_1 = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
  pseudo->nrad_2 = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
  pseudo->nrad_3 = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;

  pseudo->nl_alp = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
  pseudo->nl_beta = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
  pseudo->nl_filter = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;

  pseudo->phi0_0 = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
  pseudo->phi0_1 = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
  pseudo->phi0_2 = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;

/*==========================================================================*/
/* III) Loop over all unique atom types and get the vps stuff               */ 

  if(myid==0){
    pseudo->n_rad_max    = 0;
    pseudo->n_ang_max    = 0;
    pseudo->n_ang_max_kb = 0;
    pseudo->n_ang_max_gh = 0;
    ngh_max = 0;
    for(i=1;i<=natm_typ;i++) {
  /*--------------------------------------------------------------------------*/
  /*     A) First search the user defined data base                           */

      ifound = 0;
      strcpy(cvps_typ->atm1,atm_typ[i]);
      if(strcasecmp(filename_parse->user_vps_name,"")!=0) {
	search_base_vps(filename_parse->user_vps_name,
			 cvps_typ,fun_dict,num_fun_dict,
			 &vps_dict_tmp,vps_dict,num_vps_dict,&ifound);
	if(ifound==1){
	  set_vps_params(vps_dict,
		   filename_parse->user_vps_name,fun_key,
		   &(pseudo->ivps_label[i]),&(pseudo->iformat[i]),filename,
		   &(pseudo->loc_opt[i]),&(pseudo->n_ang[i]),
		   &(pseudo->rcut_nl[i]),&ngh_now,
		   &(pseudo->nrad_0[i]),
		   &(pseudo->nrad_1[i]),&(pseudo->nrad_2[i]),
		   &(pseudo->nrad_3[i]),
		   &(pseudo->nl_alp[i]),&(pseudo->nl_beta[i]),
		   &(pseudo->nl_filter[i]),
		   &(pseudo->phi0_0[i]),
		   &(pseudo->phi0_1[i]),&(pseudo->phi0_2[i]));
	}//endif
      }//endif
  /*--------------------------------------------------------------------------*/
  /*     B) If you haven't found it search the default data base              */

      if(ifound == 0) {
	search_base_vps(filename_parse->def_vps_name,
			 cvps_typ,fun_dict,num_fun_dict,
			 &vps_dict_tmp,vps_dict,num_vps_dict,&ifound);
	if(ifound==1){
	  set_vps_params(vps_dict,
		   filename_parse->def_vps_name,fun_key,
		   &(pseudo->ivps_label[i]),&(pseudo->iformat[i]),filename,
		   &(pseudo->loc_opt[i]),&(pseudo->n_ang[i]),
		   &(pseudo->rcut_nl[i]),&ngh_now,
		   &(pseudo->nrad_0[i]),
		   &(pseudo->nrad_1[i]),&(pseudo->nrad_2[i]),
		   &(pseudo->nrad_3[i]),
		   &(pseudo->nl_alp[i]),&(pseudo->nl_beta[i]),
		   &(pseudo->nl_filter[i]),
		   &(pseudo->phi0_0[i]),
		   &(pseudo->phi0_1[i]),&(pseudo->phi0_2[i]));

	}//endif
      }//endif
      //printf("vnl_kb_flag %i\n",pseudo->vnl_kb_flag);
  /*--------------------------------------------------------------------------*/
  /*     C) Make sure you have now found this puppy, if not exit              */

      if(ifound==0){
	 printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	 printf("Electron pseudopotential interaction with\n"); 
	 printf("%s\n",atm_typ[i]);
	 printf("not found in default interaction data base\n");
	 printf("pi_md.vps\n");
	 if(strlen(filename_parse->user_vps_name) > 0)  {
	     printf("or in user defined pseudopot data base\n");
	     printf("%s\n",filename_parse->user_vps_name);
	 /*endif*/}
	putchar('\n');
	printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	fflush(stdout);
	exit(1);
      }//endif

  /*--------------------------------------------------------------------------*/
  /*     D) Find maximum angular momentum component                           */
      pseudo->n_ang_max = (pseudo->n_ang_max>pseudo->n_ang[i]? 
			   pseudo->n_ang_max:pseudo->n_ang[i]);
      if(pseudo->ivps_label[i]!=2){
        pseudo->n_ang_max_kb = (pseudo->n_ang_max_kb>pseudo->n_ang[i]? 
			        pseudo->n_ang_max_kb:pseudo->n_ang[i]);
      }else{
        pseudo->n_ang_max_gh = (pseudo->n_ang_max_gh>pseudo->n_ang[i]? 
			        pseudo->n_ang_max_gh:pseudo->n_ang[i]);
        natm_typ_gh++;
      }//endif
      pseudo->n_rad_max = (pseudo->n_rad_max>pseudo->nrad_0[i]? 
			   pseudo->n_rad_max:pseudo->nrad_0[i]);
      pseudo->n_rad_max = (pseudo->n_rad_max>pseudo->nrad_1[i]?
			   pseudo->n_rad_max:pseudo->nrad_1[i]);
      pseudo->n_rad_max = (pseudo->n_rad_max>pseudo->nrad_2[i]?
			   pseudo->n_rad_max:pseudo->nrad_2[i]);
      pseudo->n_rad_max = (pseudo->n_rad_max>pseudo->nrad_3[i]?
			   pseudo->n_rad_max:pseudo->nrad_3[i]);
      ngh_max           = MAX(ngh_max,ngh_now);
      strcpy(vps_file[i].name,filename);
      strcpy(filename_parse->vps_name[i],filename);
      strcpy(pseudo->vps_file[i].name,filename);
    }//endfor natm_typ
    pseudo->ngh = ngh_max;
    pseudo->natm_typ_gh  = natm_typ_gh;
  }//endif : myid==0

  if(num_proc>1){
    Bcast(&(pseudo->ivps_label[1]),natm_typ,MPI_INT,0,comm);
    Bcast(&(pseudo->loc_opt[1]),natm_typ,MPI_INT,0,comm);
    Bcast(&(pseudo->n_ang[1]),natm_typ,MPI_INT,0,comm);
    Bcast(&(pseudo->rcut_nl[1]),natm_typ,MPI_DOUBLE,0,comm);
    Bcast(&(pseudo->nrad_0[1]),natm_typ,MPI_INT,0,comm);
    Bcast(&(pseudo->nrad_1[1]),natm_typ,MPI_INT,0,comm);
    Bcast(&(pseudo->nrad_2[1]),natm_typ,MPI_INT,0,comm);
    Bcast(&(pseudo->nrad_3[1]),natm_typ,MPI_INT,0,comm);

    Bcast(&(pseudo->nl_alp[1]),natm_typ,MPI_DOUBLE,0,comm);
    Bcast(&(pseudo->nl_beta[1]),natm_typ,MPI_DOUBLE,0,comm);
    Bcast(&(pseudo->nl_filter[1]),natm_typ,MPI_INT,0,comm);

    Bcast(&(pseudo->phi0_0[1]),natm_typ,MPI_DOUBLE,0,comm);
    Bcast(&(pseudo->phi0_1[1]),natm_typ,MPI_DOUBLE,0,comm);
    Bcast(&(pseudo->phi0_2[1]),natm_typ,MPI_DOUBLE,0,comm);


    Bcast(&(pseudo->n_rad_max),1,MPI_INT,0,comm);
    Bcast(&(pseudo->n_ang_max),1,MPI_INT,0,comm);
    Bcast(&(pseudo->n_ang_max_kb),1,MPI_INT,0,comm);
    Bcast(&(pseudo->n_ang_max_gh),1,MPI_INT,0,comm);

    Bcast(&(pseudo->ngh),1,MPI_INT,0,comm);
    Bcast(&(pseudo->natm_typ_gh),1,MPI_INT,0,comm);
    Bcast(&(ngh_max),1,MPI_INT,0,comm);

    natm_typ_gh = pseudo->natm_typ_gh;
  }//endif

/*==========================================================================*/
/*  III) Allocate more memory for pseudopotentials                          */

/*--------------------------------------------------------------------------*/
/* i) Malloc Pseudo spline and other stuff */
  pseudo->nsplin_g_tot = (pseudo->n_ang_max+1)*(pseudo->n_rad_max)
                         *(pseudo->nsplin_g)*natm_typ;
  nsplin_mall = pseudo->nsplin_g_tot;
  norm_mall   = (pseudo->n_ang_max+1)*natm_typ*
                 (pseudo->n_rad_max)*(pseudo->n_rad_max);
  if((nsplin_mall % 2)==0){nsplin_mall++;}
  if((norm_mall % 2)==0){norm_mall++;}
  pseudo->nsplin_g_mall = nsplin_mall;
  pseudo->nvpsnorm_mall = norm_mall;

  pseudo->vps0 = (double*)cmalloc(nsplin_mall*sizeof(double))-1;
  pseudo->vps1 = (double*)cmalloc(nsplin_mall*sizeof(double))-1;
  pseudo->vps2 = (double*)cmalloc(nsplin_mall*sizeof(double))-1;
  pseudo->vps3 = (double*)cmalloc(nsplin_mall*sizeof(double))-1;
 
  now_mem = (nsplin_mall*4*sizeof(double))*1.e-06;

  if(cp_ptens_calc==1||iopt_cp_dvr){
     pseudo->dvps0 = (double*)cmalloc(nsplin_mall*sizeof(double))-1;
     pseudo->dvps1 = (double*)cmalloc(nsplin_mall*sizeof(double))-1;
     pseudo->dvps2 = (double*)cmalloc(nsplin_mall*sizeof(double))-1;
     pseudo->dvps3 = (double*)cmalloc(nsplin_mall*sizeof(double))-1;
  }//endif
  pseudo->vpsnorm = (double*)cmalloc(norm_mall*sizeof(double))-1;
  pseudo->gzvps   = (double*)cmalloc(natm_typ_mall*sizeof(double))-1;
  pseudo->gzvps0  = (double*)cmalloc(natm_typ_mall*pseudo->n_rad_max*
                                       sizeof(double))-1;
  pseudo->nrad_max_l = (int*)cmalloc(natm_typ_mall*5*sizeof(int))-1;

/*--------------------------------------------------------------------------*/
/* iii) Pseudo list : MAJOR HACKET JOB */

  natm_mall = natm_tot;
  if((natm_mall%2)==0)natm_mall++;
  nlist_mall = (pseudo->n_ang_max+1)*natm_tot;
  if((nlist_mall%2)==0){nlist_mall++;}
  pseudo->n_ang_mall    = pseudo->n_ang_max;
  pseudo->nlist_mall    = nlist_mall;

  pseudo->x0w   = (double*)cmalloc(natm_mall*sizeof(double))-1;
  pseudo->y0w   = (double*)cmalloc(natm_mall*sizeof(double))-1;
  pseudo->z0w   = (double*)cmalloc(natm_mall*sizeof(double))-1;
  pseudo->np_nl        = (int*)cmalloc((pseudo->n_ang_max+1)*sizeof(int))-1;
  pseudo->np_nl_gh     = (int*)cmalloc((pseudo->n_ang_max+1)*sizeof(int))-1;
  pseudo->ip_nl        = (int*)cmalloc(nlist_mall*sizeof(int))-1;
  pseudo->ip_nl_gh     = (int*)cmalloc(nlist_mall*sizeof(int))-1;
  pseudo->ip_nl_rev    = (int*)cmalloc(nlist_mall*sizeof(int))-1;
  pseudo->ip_nl_rev_gh = (int*)cmalloc(nlist_mall*sizeof(int))-1;

  pseudo->map_nl = (int*)cmalloc(natm_mall*sizeof(int))-1;
  pseudo->ip_loc_cp_box = (int*)cmalloc(natm_mall *sizeof(int))-1;

  pseudo->np_nl_rad_str = cmall_int_mat(1,pseudo->n_ang_max+1,
                                         1,pseudo->n_rad_max);
  pseudo->np_nl_rad_end = cmall_int_mat(1,pseudo->n_ang_max+1,
                                         1,pseudo->n_rad_max);
  pseudo->rgh = (double*)cmalloc(ngh_max*sizeof(double))-1;
  nmall_gh =  ngh_max*(pseudo->natm_typ_gh)*(pseudo->n_ang_max+1);
  pseudo->wgh = (double*)cmalloc((nmall_gh)*sizeof(double))-1;

  for(i=1;i<=nmall_gh;i++)pseudo->wgh[i] = 0.0;

  if(iopt_cp_dvr==1){
    pseudo->nsplin_r = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
    pseudo->dr_spl   = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
    pseudo->rmin     = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
    /*
    pseudo->rmax     = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
    pseudo->z_1      = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
    pseudo->z_2      = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
    pseudo->alp_1    = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
    pseudo->alp_2    = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;

    now_mem   = ( natm_typ_mall*sizeof(int)
	  +   natm_typ_mall*7*sizeof(double))*1.e-06;
    *tot_memory += now_mem;
    */
  }

/*--------------------------------------------------------------------------*/
/* iv) Output */

/*==========================================================================*/
/*  IV) Spline up the stuff                                                 */


   for(i=1;i<=natm_typ;i++){
     pseudo->nrad_max_l[1] = MAX(pseudo->nrad_max_l[1],pseudo->nrad_0[i]);
     pseudo->nrad_max_l[2] = MAX(pseudo->nrad_max_l[2],pseudo->nrad_1[i]);
     pseudo->nrad_max_l[3] = MAX(pseudo->nrad_max_l[3],pseudo->nrad_2[i]);
     pseudo->nrad_max_l[4] = MAX(pseudo->nrad_max_l[4],pseudo->nrad_3[i]);
   }/*endfor*/

   pseudo->dg_spl = ((pseudo->gmax_spl)-(pseudo->gmin_spl))
                    /((double)(pseudo->nsplin_g));

   for(i=1;i<=natm_typ*(pseudo->n_rad_max);i++){
     (pseudo->gzvps0)[i] = 0;
   }/*endfor*/

   if(iopt_cp_dvr== 0){ /*PW basis */
     for(i=1;i<=natm_typ;i++) {
       ishift  = (i-1)*(pseudo->n_ang_max+1)*(pseudo->nsplin_g)
                      *(pseudo->n_rad_max);
       ishift2 = (i-1)*(pseudo->n_ang_max+1)*(pseudo->n_rad_max)
                      *(pseudo->n_rad_max);
       ishift3 = (i-1)*(pseudo->n_rad_max);
       if(myid==0){strcpy(filename,vps_file[i].name);}
       if(cp_ptens_calc == 1){
         make_vps_splin(filename,pseudo->loc_opt[i],pseudo->n_ang[i],
                         pseudo->ivps_label[i],pseudo->iformat[i],
                         pseudo->nsplin_g,pseudo->dg_spl,
                         pseudo->gmin_spl,
                         pseudo->gmax_spl,pseudo->gmin_true,
                         &(pseudo->vps0)[ishift],&(pseudo->vps1)[ishift],
                         &(pseudo->vps2)[ishift],&(pseudo->vps3)[ishift],
                         &(pseudo->dvps0)[ishift],&(pseudo->dvps1)[ishift],
                         &(pseudo->dvps2)[ishift],&(pseudo->dvps3)[ishift],
                         &(pseudo->gzvps)[i],&(pseudo->gzvps0)[ishift3],
                         &(pseudo->q_pseud[i]),
                         &(pseudo->vpsnorm)[ishift2],
                         (pseudo->nrad_0[i]),
                         (pseudo->nrad_1[i]),(pseudo->nrad_2[i]),
                         (pseudo->nrad_3[i]),
                         cp_ptens_calc,myid,comm,num_proc,
                         cp_dual_grid_opt,alpha_conv_dual,pseudo->n_rad_max,
                         pseudo->n_ang_max_gh,
                         &(pseudo->ngh),pseudo->rgh,pseudo->wgh);
       }else{
         make_vps_splin(filename,pseudo->loc_opt[i],pseudo->n_ang[i],
                         pseudo->ivps_label[i],pseudo->iformat[i],
                         pseudo->nsplin_g,pseudo->dg_spl,
                         pseudo->gmin_spl,
                         pseudo->gmax_spl,pseudo->gmin_true,
                         &(pseudo->vps0)[ishift],&(pseudo->vps1)[ishift],
                         &(pseudo->vps2)[ishift],&(pseudo->vps3)[ishift],
                         &dummy0,&dummy1,
                         &dummy2,&dummy3,
                         &(pseudo->gzvps)[i],&(pseudo->gzvps0)[ishift3],
                         &(pseudo->q_pseud[i]),
                         &(pseudo->vpsnorm)[ishift2],
                         (pseudo->nrad_0[i]),
                         (pseudo->nrad_1[i]),(pseudo->nrad_2[i]),
                         (pseudo->nrad_3[i]),
                         cp_ptens_calc,myid,comm,num_proc,
                         cp_dual_grid_opt,alpha_conv_dual,pseudo->n_rad_max,
                         pseudo->n_ang_max_gh,
                         &(pseudo->ngh),pseudo->rgh,pseudo->wgh);
       }/* endif */
     }/*endfor*/
   }else{     
     grid_h = hmat[1]/(double)(cpcoeffs_info->grid_nx);
     if(grid_h < hmat[5]/(double)(cpcoeffs_info->grid_ny)){
       grid_h = hmat[5]/(double)(cpcoeffs_info->grid_ny);
     }
     if(grid_h < hmat[9]/(double)(cpcoeffs_info->grid_nz)){
       grid_h = hmat[9]/(double)(cpcoeffs_info->grid_nz);
     }

     if( cp_ptens_calc == 1){
       if(myid==0){
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("The pressure tensor is not implemented for DVR CP \n");
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
       }/*endif*/
       exit(1);
     }/*endif*/
     if(pseudo->natm_typ_gh > 0){  /* NO GAUSS-HERMITE PSEUDOPOTENTIALS */
       if(myid==0){
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("Gauss-Hermite integration for pseudopotentials has not \n");
         printf("been implemented for DVR CP \n");
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
       }/*endif*/
       exit(1);
     }/*endif*/
     if( pseudo->n_rad_max > 1 ){  /* NO GOEDECKER PSEUDOPOTENTIALS */
       if(myid==0){
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("Goedecker pseudopotentials have not ");
         printf("been implemented for DVR CP \n");
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
       }/*endif*/
       exit(1);
     }/*endif*/

     pseudo->nsplin_r_tot = 0;
     nsplin_mall_dvr = nsplin_mall;

     for(i=1;i<=natm_typ;i++) {
       if(myid==0){strcpy(filename,vps_file[i].name);}
       make_vps_splin_dvr(filename,pseudo->loc_opt[i],natm_typ,i,
                          &nsplin_mall_dvr,pseudo,
                          &(pseudo->nsplin_r_tot),
                          &(pseudo->gzvps)[i],&(pseudo->gzvps0)[i],
                          (pseudo->nl_filter[i]), (pseudo->phi0_0[i]),
                          (pseudo->phi0_1[i]),(pseudo->phi0_2[i]),
                          (pseudo->nl_alp[i]),(pseudo->nl_beta[i]),
                          myid,comm,num_proc,grid_h,cell->iperd);

     }

   }

/*==========================================================================*/
/*  V) Asign VNL Flags                                                      */

/*==========================================================================*/
/*  VI) Free                                                 */

  if(myid==0){
    cfree(&vps_file[1]);
    cfree(fun_key);
    cfree(filename);
    cfree(&word[1]);
    cfree(&fun_dict[1]);
    cfree(&vps_dict[1]);
    cfree(&vps_dict_tmp[1]);
    cfree(cvps_typ);
  }/*endif*/

/*==========================================================================*/
/*  VI) Output                                                              */

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/



