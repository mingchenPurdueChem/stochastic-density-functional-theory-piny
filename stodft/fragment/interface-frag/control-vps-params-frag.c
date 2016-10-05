/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_vps_params                           */
/*                                                                          */
/* This reads in and sets up the electron-atom interaction pseudopotential  */ 
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_vps_params_entry.h"
#include "../proto_defs/proto_search_entry.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_vps_params_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"

#define DEBUG_DKNY_OFF

#define JUERG_FACTOR_ON
#ifdef  JUERG_FACTOR_ON
#define JUERG_FACTOR 0.72
#else
#define JUERG_FACTOR 1.0
#endif

typedef struct vps_file{
 char name[MAXWORD];
}VPS_FILE;

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

  if( (pseudo->nsplin_g< 4000) && (2.0*ecut_cp > 60.0) && (iopt_cp_dvr == 0)){
     if(myid==0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Now, dude, lets clean up those files and use a \n");
       printf("reasonable number of psuedo spline points at large \n");
       printf("cutoffs. Its clear, %d points at %g Ry, won't do!\n",
               pseudo->nsplin_g,2.0*ecut_cp);
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
     }/*endif*/
     exit(1);
  }/*endif*/

/*==========================================================================*/
/* I) Convert alpha_conv_dual and Allocate some temporary character arrays  */

   alpha_conv_dual        /= (pow(vol_cp,1.0/3.0));
   pseudo->alpha_conv_dual = alpha_conv_dual;

  if(myid==0){
   filename_parse->vps_name = (NAME *) cmalloc(natm_typ*sizeof(NAME))-1;
   vps_file  = (VPS_FILE *) cmalloc(natm_typ*sizeof(VPS_FILE))-1;
   fun_key   = (char *)cmalloc(MAXWORD*sizeof(char));  
   filename  = (char *)cmalloc(MAXWORD*sizeof(char));  
   word      = (DICT_WORD *)cmalloc(sizeof(DICT_WORD))-1;  
   cvps_typ  = (CVPS *)cmalloc(sizeof(CVPS));  
   ifirst    = 1;
   set_potfun_dict(&fun_dict,&num_fun_dict,ifirst);
   set_potvps_dict(&vps_dict,&num_vps_dict,ifirst);      
   set_potvps_dict(&vps_dict_tmp,&num_vps_dict,ifirst);      
 }/*endif*/


/*==========================================================================*/
/* III) Malloc up the vps stuff                                             */ 

   natm_typ_mall      = natm_typ;

   if( (natm_typ_mall % 2)==0){natm_typ_mall++;}

   pseudo->n_ang      = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
   pseudo->loc_opt    = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
   pseudo->ivps_label = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
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
                           &(pseudo->ivps_label[i]),filename,
                           &(pseudo->loc_opt[i]),&(pseudo->n_ang[i]),
                           &(pseudo->rcut_nl[i]),&ngh_now,
                           &(pseudo->nrad_0[i]),
                           &(pseudo->nrad_1[i]),&(pseudo->nrad_2[i]),
                           &(pseudo->nrad_3[i]),
                           &(pseudo->nl_alp[i]),&(pseudo->nl_beta[i]),
                           &(pseudo->nl_filter[i]),
                           &(pseudo->phi0_0[i]),
                           &(pseudo->phi0_1[i]),&(pseudo->phi0_2[i]));
         }/*endif*/
     }/*endif*/
/*--------------------------------------------------------------------------*/
/*     B) If you haven't found it search the default data base              */

      if(ifound == 0) {
         search_base_vps(filename_parse->def_vps_name,
                         cvps_typ,fun_dict,num_fun_dict,
                         &vps_dict_tmp,vps_dict,num_vps_dict,&ifound);
         if(ifound==1){
            set_vps_params(vps_dict,
                           filename_parse->def_vps_name,fun_key,
                           &(pseudo->ivps_label[i]),filename,
                           &(pseudo->loc_opt[i]),&(pseudo->n_ang[i]),
                           &(pseudo->rcut_nl[i]),&ngh_now,
                           &(pseudo->nrad_0[i]),
                           &(pseudo->nrad_1[i]),&(pseudo->nrad_2[i]),
                           &(pseudo->nrad_3[i]),
                           &(pseudo->nl_alp[i]),&(pseudo->nl_beta[i]),
                           &(pseudo->nl_filter[i]),
                           &(pseudo->phi0_0[i]),
                           &(pseudo->phi0_1[i]),&(pseudo->phi0_2[i]));

         }/*endif*/
      }/*endif*/
/*--------------------------------------------------------------------------*/
/*     C) Make sure you have now found this puppy, if not exit              */

      if(ifound == 0) {
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
     }/*endif*/

/*--------------------------------------------------------------------------*/
/*     D) Find maximum angular momentum component                           */
      pseudo->n_ang_max = (pseudo->n_ang_max > pseudo->n_ang[i] ? 
                           pseudo->n_ang_max : pseudo->n_ang[i]);

     if(pseudo->ivps_label[i] != 2){
      pseudo->n_ang_max_kb = (pseudo->n_ang_max_kb > pseudo->n_ang[i] ? 
                              pseudo->n_ang_max_kb : pseudo->n_ang[i]);
     }else{
      pseudo->n_ang_max_gh = (pseudo->n_ang_max_gh > pseudo->n_ang[i] ? 
                              pseudo->n_ang_max_gh : pseudo->n_ang[i]);
      natm_typ_gh++;
     }/*endif*/

      pseudo->n_rad_max = (pseudo->n_rad_max > pseudo->nrad_0[i] ? 
                           pseudo->n_rad_max : pseudo->nrad_0[i]);
      pseudo->n_rad_max = (pseudo->n_rad_max > pseudo->nrad_1[i] ?
                           pseudo->n_rad_max : pseudo->nrad_1[i]);
      pseudo->n_rad_max = (pseudo->n_rad_max > pseudo->nrad_2[i] ?
                           pseudo->n_rad_max : pseudo->nrad_2[i]);
      pseudo->n_rad_max = (pseudo->n_rad_max > pseudo->nrad_3[i] ?
                           pseudo->n_rad_max : pseudo->nrad_3[i]);
      ngh_max           = MAX(ngh_max,ngh_now);

      strcpy(vps_file[i].name,filename);
      strcpy(filename_parse->vps_name[i],filename);
   }/*endfor natm_typ*/
      pseudo->ngh = ngh_max;
      pseudo->natm_typ_gh  = natm_typ_gh;
 }/*endif : myid==0*/

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

 }/*endif*/

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

   pseudo->vps0    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
   pseudo->vps1    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
   pseudo->vps2    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
   pseudo->vps3    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
 
   now_mem   = ( nsplin_mall*4 *sizeof(double))*1.e-06;
  *tot_memory += now_mem;

   if(cp_ptens_calc == 1 || iopt_cp_dvr){

     pseudo->dvps0    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
     pseudo->dvps1    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
     pseudo->dvps2    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
     pseudo->dvps3    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
   }/* endif */
   pseudo->vpsnorm = (double *) cmalloc(norm_mall*sizeof(double))-1;
   pseudo->gzvps   = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;

   pseudo->gzvps0  = (double *) cmalloc(natm_typ_mall*pseudo->n_rad_max*
                                        sizeof(double))-1;
   pseudo->nrad_max_l = (int *)cmalloc(natm_typ_mall*5*sizeof(int))-1;

/*--------------------------------------------------------------------------*/
/* iii) Pseudo list : MAJOR HACKET JOB */

   natm_mall = natm_tot;
   if((natm_mall % 2)==0){natm_mall++;}
   nlist_mall = (pseudo->n_ang_max+1)*natm_tot;
   if((nlist_mall % 2)==0){nlist_mall++;}
   pseudo->n_ang_mall    = pseudo->n_ang_max;
   pseudo->nlist_mall    = nlist_mall;

   pseudo->x0w   = (double *) cmalloc(natm_mall*sizeof(double))-1;
   pseudo->y0w   = (double *) cmalloc(natm_mall*sizeof(double))-1;
   pseudo->z0w   = (double *) cmalloc(natm_mall*sizeof(double))-1;

   pseudo->np_nl        = (int *) cmalloc((pseudo->n_ang_max+1)*sizeof(int))-1;
   pseudo->np_nl_gh     = (int *) cmalloc((pseudo->n_ang_max+1)*sizeof(int))-1;
   pseudo->ip_nl        = (int *) cmalloc(nlist_mall*sizeof(int))-1;
   pseudo->ip_nl_gh     = (int *) cmalloc(nlist_mall*sizeof(int))-1;
   pseudo->ip_nl_rev    = (int *) cmalloc(nlist_mall*sizeof(int))-1;
   pseudo->ip_nl_rev_gh = (int *) cmalloc(nlist_mall*sizeof(int))-1;

   pseudo->map_nl = (int *) cmalloc(natm_mall*sizeof(int))-1;
   pseudo->ip_loc_cp_box = (int *) cmalloc(natm_mall *sizeof(int))-1;

   pseudo->np_nl_rad_str  = cmall_int_mat(1,pseudo->n_ang_max+1,
                                          1,pseudo->n_rad_max);
   pseudo->np_nl_rad_end  = cmall_int_mat(1,pseudo->n_ang_max+1,
                                          1,pseudo->n_rad_max);

   pseudo->rgh = (double *)cmalloc(ngh_max*sizeof(double))-1;

   nmall_gh =  ngh_max*(pseudo->natm_typ_gh)*(pseudo->n_ang_max+1);

   pseudo->wgh = (double *)cmalloc((nmall_gh)*sizeof(double))-1;

  for(i=1; i<= nmall_gh; i++){
   pseudo->wgh[i] = 0.0;
  }/*endfor*/ 

  if(iopt_cp_dvr == 1){
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
                         pseudo->ivps_label[i],
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
                         pseudo->ivps_label[i],
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
/*  V) Free                                                 */

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



