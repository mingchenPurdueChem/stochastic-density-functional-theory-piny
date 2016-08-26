/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control-mol-params-frag.c                    */
/*                                                                          */
/* This subprogram passes molecular setup to fragmentation                  */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_mol_params_entry.h"
#include "../proto_defs/proto_mol_params_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_interface_frag_local.h"

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
                general_data,dict_mol,word,fun_key,&nfun_key,
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
  exit(0);
  if(tors_free_num==1){
    tors_free_num = bondedMini->tors_free.num; /* changed in above routine */
  }/*endif*/

  for(i=1; i<= nmol_typ; i++){
    text_mol[i] = text_nhc_mol[i];   
  }/*endfor*/

  if(stodftOn==1&&readCoeffFlag==1){
    if((numStateStoUp!=cpcoeffs_info->nstate_up)||
        (numStateStoDn!=cpcoeffs_info->nstate_dn)){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("You want to readin stochastic orbitals but the number of \n");
      printf("orbitals is inconsistent! Please check the parameters!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }
  }

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

  if(cp_on==1&&cp->cpopts.cp_dual_grid_opt>=1){
    if(cp_parse->cp_ecut_dual_grid>cp_parse->cp_ecut){
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
		GENERAL_DATA *general_data,DICT_MOL *dict_mol,DICT_WORD *word,
                char *fun_key,int *nfun_key,
                int iextend,double t_ext,int ifirst,int pi_beads,
		DICT_MOL *dictMolAll)
/*=======================================================================*/
{/*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                               */

  ATOMMAPS *atommaps	    = &(classMini->atommaps);
  CPOPTS   *cpopts	= &(cpMini->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info	= &(cpMini->cpcoeffs_info);
  SURFACE *surface	= &(classMini->surface);
  CLATOMS_INFO *clatoms_info	= &(classMini->clatoms_info);
  ATOMMAPS *atommapsMacro   = &(class->atommaps);
  DICT_MOL *dictMolFrag;
  STODFTINFO *stodftInfo    = &(cp->stodftInfo);
  FRAGINFO *fragInfo	    = &(stodftInfo->fragInfo);

  int nline,nkey,i,num,ierr;
  int iMol,iType,iKey;
  int nmol_typ            = atommaps->nmol_typ;
  int num_fun_dict        = dict_mol->num_fun_dict;
  int bond_free_num       = bondedMini->bond_free.num;
  int bend_free_num       = bondedMini->bend_free.num;
  int tors_free_num       = bondedMini->tors_free.num;
  int rbar_sig_free_iopt  = bondedMini->rbar_sig_free.iopt;
  int nbar_bond;
  int countMol;
  int numMolTyp	      = atommapsMacro->nmol_typ;
  int iFrag	  = fragInfo->iFrag;
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
      set_wave_params(molsetname,fun_key,
                    dict_mol->wave_dict,dict_mol->num_wave_dict,
                    cpopts,cpcoeffs_info,cp_parse);
      break;      
    case 3:
      set_bond_free_params(molsetname,fun_key,
                        dict_mol->bond_free_dict,
                        dict_mol->num_bond_free_dict,
                        &(bonded->bond_free),free_parse,
                        atommaps->nmol_typ);
      break;
    case 4:
      set_bend_free_params(molsetname,fun_key,
                        dict_mol->bend_free_dict,
                        dict_mol->num_bend_free_dict,
                        &(bonded->bend_free),free_parse,
                        atommaps->nmol_typ);
      break;
      
    case 5:
      set_tors_free_params(molsetname,fun_key,
                        dict_mol->tors_free_dict,
                        dict_mol->num_tors_free_dict,
                        &(bonded->tors_free),free_parse,
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
                        &(bonded->rbar_sig_free),free_parse,
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
/* III) Assign mol_dict variables		       */

  for(iType=0;iType<numMolTypeNow;iType++){
    typeInd = molTypeFragNow[iType];
    sprintf(dictMolFrag[iType].mol_dict[1].keyarg,"%i",iType);
    sprintf(dictMolFrag[iType].mol_dict[2].keyarg,"%i",molNumTypeFragNow[iType]);
    for(iKey=3;iKey<14;iKey++){
      strcpy(dictMolFrag[iType].mol_dict[iKey].keyarg,dictMolAll[typeInd].mol_dict[iKey].keyarg);
    }
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

  if(tors_free_num>0){tors_free_num=bonded->tors_free.num;}
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








