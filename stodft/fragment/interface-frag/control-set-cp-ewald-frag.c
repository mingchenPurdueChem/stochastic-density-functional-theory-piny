/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_ewald                                    */
/*                                                                          */
/* This reads in and sets up the k-space                                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_cp_ewald_entry.h"
#include "../proto_defs/proto_cp_ewald_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define DEBUG_CLUS_CORR_OFF
#define CHECK_CLUS_OFF

#define MAX_INT 12.0
#define MIN3(A,B,C) (MIN(MIN((A),(B)),(C)))

#define ORIG_OFF
#define PME

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Control ewald/cp g-space initialization */
/*==========================================================================*/

void controlSetCpEwaldFrag(GENERAL_DATA *generalDataMini,CLASS *classMini,
			   CP *cpMini,BONDED *bondedMini,CP *cp,CLASS *class,
			    GENERAL_DATA *general_data,BONDED *bonded,
                          CP_PARSE *cp_parse)
/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*************************************************************************/
/* Our FFT grid is numGridDim[0](c)*numGridDim[1](b)*numGridDim[2](a)	 */
/* The k space coef number should be 2n+1 for each dimension, from -n to */
/* n. and n=numGridDim[i]-1. The total k space coef number is ((2na+1)*	 */
/* (2nb+1)*(2nc+1)-1)/2. Wf and density share the same k space		 */
/*************************************************************************/
/*=======================================================================*/
/*            Local variable declarations:                               */

  CPOPTS *cpopts		= &(cpMini->cpopts);
  SIMOPTS *simopts		= &(generalDataMini->simopts);
  CELL *cell			= &(generalDataMini->cell);
  CPCOEFFS_INFO *cpcoeffs_info	= &(cpMini->cpcoeffs_info);
  EWALD *ewald			= &(generalDataMini->ewald);
  CPEWALD *cpewald		= &(cpMini->cpewald);
  PSEUDO *pseudo		= &(cpMini->pseudo);
  EWD_SCR *ewd_scr		= &(classMini->ewd_scr);
  ECOR *ecor			= &(bondedMini->ecor);
  STODFTINFO *stodftInfo	= cp->stodftInfo;
  FRAGINFO *fragInfo		= stodftInfo->fragInfo;
  
  int iFrag = fragInfo->iFrag;
  int idum1=1,idum2=1,idum3=1;
  int iii;                            /* Num: Debug tool                    */
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
  int iperd		= cell->iperd;
  int box_rat		= cpewald->box_rat;
  int kmax_ewd		= cp_parse->kmax_ewd;
  int kmax_res		= cp_parse->kmax_res;
  int cp_lsda		= cpopts->cp_lsda;
  int cp_dual_grid_opt_on = cpopts->cp_dual_grid_opt_on;
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
  double *gmax_spl  = &(pseodu->gmax_spl);
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
   
   /*
   calc_cutoff(kmax_ewd,&ecut_now,&(cp_parse->cp_ecut),cp_on,
                kmax_cp,kmaxv,hmati_ewd_cp,deth_cp);  
   printf("!!!!!!!!!!!!!! ecut_now %lg\n",ecut_now);

   countkvec3d(&(ewald->nktot),ecut_now,kmaxv,hmati_ewd_cp); 


   nktot                   = ewald->nktot;
   cpcoeffs_info->ecut     = ecut_now;
   cpcoeffs_info->ncoef_l  = nktot+1;
   ncoef_l                 = nktot+1;
   ecor->ecut              = 4.0*ecut_now;
   ewald->ecut             = 4.0*ecut_now;
   ewald->nkc_max          = kmaxv[3];
   */
  ecor->ecut = bonded->ecor.ecut; // I don't think I need ecor but just to prevent segfault
  ewald->ecut = general_data->ewald.ecut;
  ecut_now = 1.0e30; // A big number so that all grids included

  kmaxv[1] = numGridDim[2]/2-1; 
  kmaxv[2] = numGridDim[1]/2-1;
  kmaxv[3] = numGridDim[0]/2-1;
  ewald->nkc_max = kmaxv[3];
  kmax_cp[1] = kmaxv[1];
  kmax_cp[2] = kmaxv[2];
  kmax_cp[3] = kmaxv[3];
  nktot = ((kmaxv[1]*2+1)*(kmaxv[2]*2+1)*(kmaxv[3]*2+1)-1)/2;
  ewald->nktot = nktot;
  cpcoeffs_info->ncoef_l = nktot+1;
  kmax_cp_dens_cp_box[1] = kmax_cp[1];
  kmax_cp_dens_cp_box[2] = kmax_cp[2];
  kmax_cp_dens_cp_box[3] = kmax_cp[3];


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
      /*
      ecut_sm = ecut_dens_cp_box_now;
      countkvec3d_sm(&(cpewald->nktot_sm),ecut_sm,
                     kmax_cp_dens_cp_box,hmati_ewd_cp);
      nktot_sm = cpewald->nktot_sm;
      cpcoeffs_info->ncoef   = nktot_sm+1;
      ncoef                  = nktot_sm+1;
      */    
    ecut_now = 1.0e30;  
    nktot_sm = nktot;
    cpewald->nktot_sm = nktot;
    cpcoeffs_info->ncoef = nktot_sm+1;
    ncoef = nktot_sm+1;
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
/* Control the setup of the FFT packages                                    */
/*==========================================================================*/

void controlFFTPkgFrag(GENERAL_DATA *generalDataMini,CLASS *classMini,CP *cpMini,
			CP *cp)

/*=========================================================================*/
     {/*Begin Routine*/
/*=========================================================================*/
/*    Local Variables */

  PARA_FFT_PKG3D *cp_sclr_fft_pkg_sm = &(cpMini->cp_sclr_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_para_fft_pkg_sm = &(cpMini->cp_para_fft_pkg_sm);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg_dens_cp_box = &(cpMini->cp_sclr_fft_pkg_dens_cp_box);
  PARA_FFT_PKG3D *cp_para_fft_pkg_dens_cp_box = &(cpMini->cp_para_fft_pkg_dens_cp_box);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg_lg = &(cpMini->cp_sclr_fft_pkg_lg);
  PARA_FFT_PKG3D *cp_para_fft_pkg_lg = &(cpMini->cp_para_fft_pkg_lg);
  PARA_FFT_PKG3D *pme_fft_pkg = &(cpMini->pme_fft_pkg);
  PARA_FFT_PKG3D *pme_res_fft_pkg = &(cpMini->pme_res_fft_pkg);
  EWALD *ewald = &(generalDataMini->ewald);
  CPEWALD *cpewald = &(cpMini->cpewald);
  PART_MESH *part_mesh = &(classMini->part_mesh);
  CPCOEFFS_INFO *cpcoeffs_info = &(cpMini->cpcoeffs_info);
  COMMUNICATE *communicate = &(classMini->communicate);
  CPOPTS *cpopts = &(cpMini->cpopts);
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo->stodftInfo->fragInfo;
  

  int iFrag = fragInfo->iFrag;
  int nkf1,nkf2,nkf3,nfft_ext,iii;
  int nkf1_dens_cp_box,nkf2_dens_cp_box,nkf3_dens_cp_box;
          /* used for denisty on cp grid when have 2 boxes*/
  int cp_on		  = 1;  // Always do cp
  int cp_lsda		  = cpopts->cp_lsda;
  int cp_para_opt	  = cpopts->cp_para_opt;
  int cp_dual_grid_opt_on = cpopts->cp_dual_grid_opt_on;
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
  int *numGridDim	 = fragInfo->numGridFragDim[iFrag];

/*=========================================================================*/
/* 0) Print to screen and check for nproc > nstate error */

/*=========================================================================*/
/* 0.1) Set CP FFT Size  */

  nkf1 = numGridDim[2];
  nkf2 = numGridDim[1];
  nkf3 = numGridDim[0];
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


