/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_ewald (DVR)                              */
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

void control_set_cp_ewald_dvr(SIMOPTS *simopts,CELL *cell,
                          CPCOEFFS_INFO *cpcoeffs_info,
                          EWALD *ewald, CPEWALD *cpewald,CP_PARSE *cp_parse,
                          double *gmin_spl,double *gmax_spl,
                          int kmax_ewd, double *tot_memory,int int_res_ter,
                          ECOR *ecor, int myid, int cp_dual_grid_opt_on, int pme_on)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations:                               */

   double ecut_now;                    /* Num: Energy cutoff                 */
   double ecut_dens_cp_box_now;        /* Num: Energy cutoff                 */
   double ecut_sm;                     /* Num: Energy cutoff for wf          */
   double *hmati_ewd;                  /* Num: Inverse ewald h matrix        */
   double *hmati_ewd_cp;               /* Num: Inverse cp ewald h matrix     */
   double deth,deth_cp,side;           /* Num: Volumes and sizes             */
   double now_mem;                     /* Num: Memory used here              */
   int nmall;
   int cp_on;                          /* Num: CP on flag                    */

   int nktot,nktot_res,nktot_sm;       /* Num: Spherically cutoff k vecs     */
   int nktot_dens_cp_box;
   int ncoef,ncoef_dens_cp_box;

   int *kmaxv;                         /* Lst: K-vector ranges               */
   int *kmax_cp;
   int *kmaxv_dens_cp_box;
   int *kmax_cp_dens_cp_box;

   /*--------------------------------------*/
   /* Local pointers                       */

   int iperd           = cell->iperd;

   double *hmat_ewd    = cell->hmat_ewd;
   double *hmat_ewd_cp = cell->hmat_ewd_cp;

   double cp_dvrdens     = cp_parse->cp_dvrdens_def;
   double la             = cell->hmat[1];
   double lb             = cell->hmat[5];
   double lc             = cell->hmat[9];
   double nx2,ny2,nz2;


/*=======================================================================*/
/* 0) Output to screen                                                   */

   if(myid==0){
     printf("\n");PRINT_LINE_STAR;
     printf("Setting up DVR grids\n");
     PRINT_LINE_DASH;printf("\n");
   }/*endif*/

/*=======================================================================*/
/* I) Set cp switch and initialize respa kvectors                        */

   cp_on = simopts->cp_min 
          +simopts->cp_wave_min
          +simopts->cp
          +simopts->cp_wave
          +simopts->debug_cp
          +simopts->cp_pimd
          +simopts->debug_cp_pimd+simopts->cp_wave_pimd
          +simopts->cp_wave_min_pimd;

   if(int_res_ter==0){ewald->nktot_res=0;ecor->nktot_res=0;}

   if(cp_on != 1) {
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("CP has to be on!! DVR is only for CPMD part of PINY.\n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
   }


/*=======================================================================*/
/* II) Allocate simple memory                                            */

   hmati_ewd      = (double *) cmalloc((size_t)9*sizeof(double))-1;
   hmati_ewd_cp   = (double *) cmalloc((size_t)9*sizeof(double))-1;
   kmaxv          =    (int *) cmalloc((size_t)3*sizeof(int))-1;

   cpewald->kmax_cp = (int *) cmalloc((size_t)3*sizeof(int))-1;
   kmax_cp          = cpewald->kmax_cp;
   cpewald->kmax_cp_dens_cp_box = (int *) cmalloc((size_t)3*sizeof(int))-1;
   kmax_cp_dens_cp_box          = cpewald->kmax_cp_dens_cp_box;

   if(cp_dual_grid_opt_on != 0){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("         No dual-grid option for DVR.\n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
   }/*endif*/


/*==========================================================================*/
/* III) Get inverse cell matrix and convert ewald_alpha                     */

   gethinv(hmat_ewd,hmati_ewd,&deth,iperd);    
   gethinv(hmat_ewd_cp,hmati_ewd_cp,&deth_cp,iperd);    

   side  = pow(deth,(1.0/3.0));  
   (ewald->alp_ewd) /= side;

/*==========================================================================*/
/* IV) Set-up DVR grids  */

   cpcoeffs_info->dvr_grid_dens = cp_dvrdens;

   nx2 = (int)(la/(cpcoeffs_info->dvr_grid_dens));
   ny2 = (int)(lb/(cpcoeffs_info->dvr_grid_dens));
   nz2 = (int)(lc/(cpcoeffs_info->dvr_grid_dens));

   cpcoeffs_info->grid_nx = 2*nx2+1;
   cpcoeffs_info->grid_ny = 2*ny2+1;
   cpcoeffs_info->grid_nz = 2*nz2+1;

/*==========================================================================*/
/* V) Count number k vectors, Malloc and Fill                               */

  /* A) count number k vectors */

   calc_kmaxv(kmaxv,kmax_cp,&(cpcoeffs_info->grid_nx), &(cpcoeffs_info->grid_ny),
             &(cpcoeffs_info->grid_nz),&ecut_now,deth_cp,kmax_ewd,myid);

   ewald->nktot = (cpcoeffs_info->grid_nx)*(cpcoeffs_info->grid_ny)
                 *(cpcoeffs_info->grid_nz) -1;

   nktot                   = ewald->nktot;
   cpcoeffs_info->ecut     = ecut_now;
   cpcoeffs_info->ncoef_l  = nktot+1;
   cpcoeffs_info->ncoef    = nktot+1;
   ecor->ecut              = 4.0*ecut_now;
   ewald->ecut             = 4.0*ecut_now;
   ewald->nkc_max          = kmaxv[3];

   cpcoeffs_info->scale_fact = la*lb*lc/((double)(nktot+1));

   /* B) For dualing : No dual-grid option for DVR */

   ecut_dens_cp_box_now   = ecut_now;
   kmax_cp_dens_cp_box[1] = kmax_cp[1];
   kmax_cp_dens_cp_box[2] = kmax_cp[2];
   kmax_cp_dens_cp_box[3] = kmax_cp[3];

   /* C) Malloc  */

   nmall =  nktot+1;   
   if((nmall % 2)==0){nmall++;}

   ewald->nktot_mall = nmall;

   now_mem    = (nmall*(sizeof(double)*0 + sizeof(int)*5))*1.e-06;
   if(int_res_ter!=0){
      now_mem  += (nmall*(sizeof(double)*0 + sizeof(int)*6 ))*1.e-06;
   }/*endif*/
   (*tot_memory) += now_mem;

   ewald->kastr = (int *) cmalloc(nmall*sizeof(int))-1;
   ewald->kbstr = (int *) cmalloc(nmall*sizeof(int))-1;
   ewald->kcstr = (int *) cmalloc(nmall*sizeof(int))-1;
   ewald->ibrk1 = (int *) cmalloc(nmall*sizeof(int))-1;
   ewald->ibrk2 = (int *) cmalloc(nmall*sizeof(int))-1;
   if(int_res_ter != 0) {
      ewald->nktot_res_mall = nmall;
      ewald->kastr_res      = (int *) cmalloc(nmall*sizeof(int))-1;
      ewald->kbstr_res      = (int *) cmalloc(nmall*sizeof(int))-1;
      ewald->kcstr_res      = (int *) cmalloc(nmall*sizeof(int))-1;
      ewald->ibrk1_res      = (int *) cmalloc(nmall*sizeof(int))-1;
      ewald->ibrk2_res      = (int *) cmalloc(nmall*sizeof(int))-1;
      ewald->ibrk3          = (int *) cmalloc(nmall*sizeof(int))-1;
   }/*endif*/

   if(myid==0){
      printf("Ewald allocation: %g Mbytes; Total memory %g Mbytes\n",
              now_mem,*tot_memory);
   }/*endif*/

   /* D) Fill */

   setkvec3d_all(nktot,kmaxv,hmati_ewd, ewald->kastr,ewald->kbstr,ewald->kcstr,
                 ewald->ibrk1,ewald->ibrk2,gmin_spl,gmax_spl);

/*=======================================================================*/
/* V) Setup PME                                                          */

   if(pme_on==1){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Particle-Mesh-Ewald (PME) part has been removed\n"),
     printf("for the DVR part of CPMD implementation\n");
     printf("There is no spherical truncation in reciprocal \n");
     printf("space vectors in DVR, and, therefore, current\n");
     printf("PME implementation may not work for DVR\n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
   }/*endif*/

/*=======================================================================*/
/* VI) Set up RESPA k vectors and RESPA PME                              */

   if(int_res_ter != 0) {
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("    RESPA option in DVR???? You're kidding!\n"),
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
   }/*endif*/

/*=======================================================================*/
/* VII) Setup CP coefs and k vectors                                     */

    /* A)  Count the k-vectors */

    /* Spherical truncation is not used in DVR. However, it is */
    /* necessary to set up FFT package. */

    ecut_sm = ecut_dens_cp_box_now;
    countkvec3d_sm(&(cpewald->nktot_sm),ecut_sm,
                    kmax_cp_dens_cp_box,hmati_ewd_cp);

    nktot_sm = cpewald->nktot_sm;
    ncoef                  = nktot_sm+1;
    cpcoeffs_info->ncoef_pw = ncoef;

    /*  B)  Malloc  */

    nmall =  (nktot_sm+1);
    if((nmall % 2)==0){nmall++;}
    cpewald->nktot_cp_sm_mall = nmall;
    now_mem = (nmall*(sizeof(double)*0 + sizeof(int)*5 ))*1.e-06;
    *tot_memory += now_mem;
    cpewald->kastr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
    cpewald->kbstr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
    cpewald->kcstr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
    cpewald->ibrk1_sm = (int *) cmalloc(nmall*sizeof(int))-1;
    cpewald->ibrk2_sm = (int *) cmalloc(nmall*sizeof(int))-1;

    /*  C)  Fill and check */

    setkvec3d_sm(nktot_sm,ecut_sm,kmax_cp_dens_cp_box,hmati_ewd_cp,
                 cpewald->kastr_sm,cpewald->kbstr_sm,cpewald->kcstr_sm,
                 cpewald->ibrk1_sm,cpewald->ibrk2_sm,
                 &(cpewald->gw_gmin),&(cpewald->gw_gmax));

/*--------------------------------------------------------------------*/
/*  VIII) Set up the cp masses                                           */

    nmall =  cpcoeffs_info->ncoef;
    now_mem += (nmall*sizeof(double))*1.e-06;
    *tot_memory += now_mem;

    cpcoeffs_info->cmass = (double *)cmalloc(nmall*sizeof(double))-1;

    if(myid==0){
      printf("CP allocation: %g Mbytes; Total memory %g Mbytes\n",
                 now_mem,*tot_memory);
    }/*endif*/


    set_cpmass_dvr(cpcoeffs_info,&(cp_parse->cp_mass_tau_def),
                   &(cpcoeffs_info->icmass_unif));

    if(myid==0){
      printf("\nYour cp_mass is %g a.u.\n",cpcoeffs_info->cmass[1]/cpcoeffs_info->scale_fact);
    }

/*=======================================================================*/
/* VIII) Output time has arrived                                         */

   if(myid ==0){

     printf("Your DVR grid for density and orbitals is %d by %d by %d\n",
             cpcoeffs_info->grid_nx,cpcoeffs_info->grid_ny,cpcoeffs_info->grid_nz);
     printf("Total number of grid points is %d\n",cpcoeffs_info->ncoef); 
     if(iperd==3){
       printf("Same number of k-vectors are used for the FFT of density\n");
     }
     putchar('\n');

   }/*endif:myid==0*/

/*=======================================================================*/
/* IX) Free excess memory                                                */

   cfree(&(hmati_ewd)[1]);
   cfree(&(hmati_ewd_cp)[1]);
   cfree(&(kmaxv)[1]);

/*-----------------------------------------------------------------------*/ 
  }/* end routine */
/*=======================================================================*/

/*==========================================================================*/
/* Control the setup of the FFT packages                                    */
/*==========================================================================*/

void control_fft_pkg_dvr(PARA_FFT_PKG3D *cp_sclr_fft_pkg_sm,
                         PARA_FFT_PKG3D *cp_para_fft_pkg_sm,
                         PARA_FFT_PKG3D *cp_sclr_fft_pkg_lg,
                         PARA_FFT_PKG3D *cp_para_fft_pkg_lg,
                         EWALD* ewald,CPEWALD *cpewald,
                         CPCOEFFS_INFO *cpcoeffs_info,
                         COMMUNICATE *communicate,int cp_on, double *tot_memory,
                         int cp_para_opt,int cp_dual_grid_opt_on,
                         CPSCR *cpscr)

/*=========================================================================*/
     {/*Begin Routine*/
/*=========================================================================*/
/*    Local Variables */

  int nkf1 = cpcoeffs_info->grid_nx;
  int nkf2 = cpcoeffs_info->grid_ny;
  int nkf3 = cpcoeffs_info->grid_nz;
  int nmall;

  int nstate_up          = cpcoeffs_info->nstate_up;
  int ncoef              = cpcoeffs_info->ncoef;
  int ncoef_pw           = cpcoeffs_info->ncoef_pw;
  int nktot              = ewald->nktot;
  int myid               = communicate->myid;
  int np_states          = communicate->np_states;
  int myid_state         = communicate->myid_state;

  double now_mem;

/*=========================================================================*/
/* 0) Print to screen and check for nproc > nstate error */

  if(myid == 0){
    printf("\n");PRINT_LINE_STAR;
    printf("Setting up FFTs\n");
    PRINT_LINE_DASH;printf("\n");
  }/* endif myid */

  if(cp_on != 1){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("CP has to be on (cp_on==1) to use DVR as your basis\n");
    printf("Fail to set-up DVR/FFT package.\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }
  /*================MA==================*/
  /*  if(nstate_up < np_states){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("Number of states less than number of processors\n");
    printf("If possible, reduce number of processors to be\n");
    printf("less than the number of states or run a bigger system.\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
    }*/  /* endif */

/*=========================================================================*/
/* 1) Set CP FFT Size  */

  nkf1 = cpcoeffs_info->grid_nx;
  nkf2 = cpcoeffs_info->grid_ny;
  nkf3 = cpcoeffs_info->grid_nz;

  if(cp_dual_grid_opt_on >= 1){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf(" No dual-grid option for DVR basis. \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

/*=========================================================================*/
/* II) scalar package                                                       */


  if(cp_para_opt == 0){/* hybrid option */

    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf(" No hybrid option for DVR basis. \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);

    cp_sclr_fft_pkg_lg->nkf1       = nkf1;
    cp_sclr_fft_pkg_lg->nkf2       = nkf2;
    cp_sclr_fft_pkg_lg->nkf3       = nkf3;

    cp_sclr_fft_pkg_lg->nktot      = nktot;
    cp_sclr_fft_pkg_lg->ncoef      = ncoef;

    cp_sclr_fft_pkg_lg->myid       = 0;
    cp_sclr_fft_pkg_lg->myidp1     = 1;
    cp_sclr_fft_pkg_lg->num_proc   = 1;
    cp_sclr_fft_pkg_lg->comm       = communicate->comm_faux;

    create_para_fft_pkg3d(cp_sclr_fft_pkg_lg,
                          ewald->kastr,ewald->kbstr,
                          ewald->kcstr,cp_dual_grid_opt_on);

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
/* III) parallel package                                                 */

  if(cp_para_opt == 1){/* full g option */


    cp_para_fft_pkg_sm->nkf1       = nkf1;
    cp_para_fft_pkg_sm->nkf2       = nkf2;
    cp_para_fft_pkg_sm->nkf3       = nkf3;

    cp_para_fft_pkg_sm->nktot      = ncoef_pw-1;
    cp_para_fft_pkg_sm->ncoef      = ncoef_pw;

    cp_para_fft_pkg_sm->myid       = myid_state;
    cp_para_fft_pkg_sm->myidp1     = myid_state+1;
    cp_para_fft_pkg_sm->num_proc   = np_states;
    cp_para_fft_pkg_sm->comm       = communicate->comm_states;

    create_para_fft_pkg3d(cp_para_fft_pkg_sm,
                          cpewald->kastr_sm,cpewald->kbstr_sm,
                          cpewald->kcstr_sm,cp_dual_grid_opt_on);

    /*large FFT grid: There is only one FFT grid in DVR */
    /* This is redundant */

    /*
    cp_para_fft_pkg_lg->nkf1       = nkf1;
    cp_para_fft_pkg_lg->nkf2       = nkf2;
    cp_para_fft_pkg_lg->nkf3       = nkf3;

    cp_para_fft_pkg_lg->nktot      = ncoef-1;
    cp_para_fft_pkg_lg->ncoef      = ncoef;

    cp_para_fft_pkg_lg->myid       = myid_state;
    cp_para_fft_pkg_lg->myidp1     = myid_state+1;
    cp_para_fft_pkg_lg->num_proc   = np_states;
    cp_para_fft_pkg_lg->comm       = communicate->comm_states;


    create_para_fft_pkg3d(cp_para_fft_pkg_lg,
                          ewald->kastr,ewald->kbstr,
                          ewald->kcstr,cp_dual_grid_opt_on);
    */


     /* Do we need scalar package? */


  }/*endif*/


/*=========================================================================*/
/* V) Build DVR FFT map: SCALAR !!! */


  build_dvr_fft_map(cpscr,cp_para_fft_pkg_sm);

  nmall = 2.0*nkf1*nkf2*nkf2;
  now_mem    = (nmall*(sizeof(int)))*1.e-06;
  (*tot_memory) += now_mem;

  if(myid_state==0){
    printf("\nDVR FFT-map allocation: %g Mbytes; Total memory %g Mbytes\n",
              now_mem,*tot_memory);
  }/*endif*/


/*=========================================================================*/
/* V) Output */

  if(myid == 0){
    printf("\n");PRINT_LINE_DASH;
    printf("Finished setting up FFTs\n");
    PRINT_LINE_STAR;printf("\n");
  }/* endif myid */

/*-------------------------------------------------------------------------*/
     }/*end routine */
/*=========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void build_dvr_fft_map(CPSCR *cpscr,PARA_FFT_PKG3D *para_fft_pkg3d)

/*=========================================================================*/
     {/*Begin Routine*/
/*=========================================================================*/

   int nkf1       = para_fft_pkg3d->nkf1;
   int nkf2       = para_fft_pkg3d->nkf2;
   int nkf3       = para_fft_pkg3d->nkf3;

   int nfft_proc  = para_fft_pkg3d->nfft_proc;
   int nfft2_proc = nfft_proc/2;

   int i,j,k,icount,itemp,iii,ioff;
   int ip,jp,kp;
   int kb_str,kb_end;
   int *ifft,*jfft,*kfft;
   int *fft_map,*fft_map_rev;


/*============================================================================*/
/* 0) memory allocation */

   cpscr->cpscr_wave.ifft = (int *) cmalloc(nkf1*sizeof(int))-1;
   cpscr->cpscr_wave.jfft = (int *) cmalloc(nkf2*sizeof(int))-1;
   cpscr->cpscr_wave.kfft = (int *) cmalloc(nkf3*sizeof(int))-1;

   cpscr->cpscr_wave.fft_map = (int *) cmalloc(nkf1*nkf2*nkf3*sizeof(int))-1;
   cpscr->cpscr_wave.fft_map_rev = (int *) cmalloc(nkf1*nkf2*nkf3*sizeof(int))-1;

   ifft = cpscr->cpscr_wave.ifft;
   jfft = cpscr->cpscr_wave.jfft;
   kfft = cpscr->cpscr_wave.kfft;

   fft_map = cpscr->cpscr_wave.fft_map;
   fft_map_rev = cpscr->cpscr_wave.fft_map_rev;

/*===========================================================================*/
/* I) These 1-D maps are used to rearrange the FFT for DVR-CP */

   for(i=1; i <= nkf1; i++){
     ifft[i] = (i-1) + (nkf1+1)/2;
     ifft[i] = (ifft[i] >= nkf1 ? ifft[i]-nkf1 : ifft[i]);
   }

   for(j=1; j <= nkf2; j++){
     jfft[j] = (j-1) + (nkf2+1)/2;
     jfft[j] = (jfft[j] >= nkf2 ? jfft[j]-nkf2 : jfft[j]);
   }

   for(k=1; k <= nkf3; k++){
     kfft[k] = (k-1) + (nkf3+1)/2;
     kfft[k] = (kfft[k] >= nkf3 ? kfft[k]-nkf3 : kfft[k]);
   }

/*===========================================================================*/
/* II) This map is used to rearrange arrays rhocr rhoci fcre and fcim        */

    /* First Build Map */

   icount = 0;
   for(k=1; k <= nkf3; k++){
     kp = kfft[k];
     for(j=1; j <= nkf2; j++){
       jp = jfft[j];
       for(i=1; i <= nkf1; i++){
         ip = ifft[i];
         /* icount++; */
         icount=(k-1)*nkf2*nkf1+(j-1)*nkf1+i;
         fft_map[icount] = kp*nkf2*nkf1 + jp*nkf1 + ip + 1;
       }
     }
   } 

/*==========================================================================*/
/* Second rearrange order such that g=o term located at end of map */
/* Recall that nkf1, nkf2, nkf3 are all odd numbers */

   ioff = (nkf1*nkf2*nkf3)/2;
   itemp = fft_map[ioff+1];


   for(i=1; i <= (nkf1*nkf2*nkf3)/2 ; i++){
     fft_map[ioff+i] = fft_map[ioff+i+1];
   }

   fft_map[nkf1*nkf2*nkf3] = itemp;

/*-------------------------------------------------------------*/
/*  map for unpacking fft for 3-1D FFTS */

   icount = 0;
   for(k = nkf3; k >= 1; k--){
     for(j = nkf2; j >= 1; j--){
       for(i = nkf1; i >= 1; i--){
         icount++;
         fft_map_rev[icount] = (k-1)*nkf2*nkf1 + (j-1)*nkf1 + i;
       }
     }
   }

/*-------------------------------------------------------------------------*/
}/*end routine */
/*=========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Allocate memory and set up F-matrix for DVR cluster boundary condition   */
/*==========================================================================*/

void get_dvr_cp_clus_fmat(CP *cp, CELL *cell, double *tot_memory)

/*=========================================================================*/
{/*Begin Routine*/
/*=========================================================================*/

  int i,j,k,id,it,ngrid;

  int ntwind  = cp->cp_dvr_clus.num_twindow;
  int ntquad  = cp->cp_dvr_clus.num_tquad;
  double tmax = cp->cp_dvr_clus.tmax;
  double grid_dens = cp->cp_dvr_clus.grid_dens;
  double rmax = cp->cp_dvr_clus.rmax;

  int grid_nx  = cp->cpcoeffs_info.grid_nx;
  int grid_ny  = cp->cpcoeffs_info.grid_ny;
  int grid_nz  = cp->cpcoeffs_info.grid_nz;

  double *tq_all, *wq_all;
  double *fmat_x, *fmat_y, *fmat_z;

  double *hmat   = cell->hmat;

  double *c, *t, *w, *x, *y, *z, *x_dvr, *y_dvr, *z_dvr;
  double dx,dy,dz, dxx, dyy, dzz,dxs,dys,dzs, dxxs,dyys,dzzs, sqrt_tq;
  double sa,sb,sc,r,pi,tq,wq,plimit,xbar,ybar,zbar;
  double alpha,beta,pwr,high,low,a,b,factor,sum,delta,xi,xj,yi,yj,zi,zj;
  double endpts[2],arg, memory_now;
  int kind = 1;
  int kpts = 0;
  int nx,ny,nz;

  int idiv,irem, istart, iend, my_num_t;

  int np_states = cp->communicate.np_states;
  int myid_state = cp->communicate.myid_state;

  if(myid_state==0 ){
    printf("\n");PRINT_LINE_STAR;
    printf("Set up DVR cluster boundary condition\n");
    PRINT_LINE_DASH;
    printf("\n");
  }

/*=========================================================================*/
/* I) Distribute t quadrature potins over processors */

  idiv =  ntwind*ntquad/np_states;
  irem = (ntwind*ntquad % np_states);
  cp->cp_dvr_clus.my_num_t = idiv;
  if(myid_state < irem) {
     cp->cp_dvr_clus.my_num_t = idiv+1;
  }/*endif*/
  if(myid_state <= irem) {
    istart = myid_state*(idiv+1)+1;
  }else {
    istart = irem*(idiv+1) + (myid_state-irem)*idiv+1;
  }/*endif*/
  iend = istart + cp->cp_dvr_clus.my_num_t -1;
  my_num_t = cp->cp_dvr_clus.my_num_t;

/*=========================================================================*/
/* II) set up quadrature points for off-diagonal F-matrix */

  ngrid = rmax/grid_dens;
  if (ngrid%2 == 0) ngrid += 1;
  nx = ngrid;
  ny = ngrid;
  nz = ngrid;

  if(myid_state==0){
    printf("Number of grid points for off-diagonal F-matrix evaluation : %d\n",ngrid);
  }

  dxx = rmax/((double)(nx));
  dyy = rmax/((double)(ny));
  dzz = rmax/((double)(nz));

  x = (double *)cmalloc(nx*sizeof(double))-1;
  y = (double *)cmalloc(ny*sizeof(double))-1;
  z = (double *)cmalloc(nz*sizeof(double))-1;

  r=-0.5*rmax-dxx/2.0;
  for (k=1;k<=nx;k++){
    r=r+dxx;
    x[k]= r;
  }
  r=-0.5*rmax-dyy/2.0;
  for (k=1;k<=ny;k++){
    r=r+dyy;
    y[k]= r;
  }
  r=-0.5*rmax-dzz/2.0;
  for (k=1;k<=nz;k++){
    r=r+dzz;
    z[k]= r;
  }

/*========================================================================*/
/* III) Set up DVR points */

  x_dvr = (double *)cmalloc(grid_nx*sizeof(double))-1;
  y_dvr = (double *)cmalloc(grid_ny*sizeof(double))-1;
  z_dvr = (double *)cmalloc(grid_nz*sizeof(double))-1;


  dx = hmat[1]/((double)(grid_nx));
  sa = 1.0/((double) grid_nx);
  for(k=1;k<=grid_nx;k++){
    r = sa*((double)(k-1)) - 0.5;
    x_dvr[k] = r*hmat[1]+dx/2.0;
  }
  dy = hmat[5]/((double)(grid_ny));
  sb = 1.0/((double) grid_ny);
  for(k=1;k<=grid_ny;k++){
    r = sb*((double)(k-1)) - 0.5;
    y_dvr[k] = r*hmat[5]+dy/2.0;
  }
  dz = hmat[9]/((double)(grid_nz));
  sc = 1.0/((double) grid_nz);
  for(k=1;k<=grid_nz;k++){
    r = sc*((double)(k-1)) - 0.5;
    z_dvr[k] = r*hmat[9]+dz/2.0;
  }

/*========================================================================*/
/* IV) Set up Quadrature points for t */

  pi = M_PI;
  c = (double *)cmalloc(ntquad*sizeof(double))-1;
  t = (double *)cmalloc(ntquad*sizeof(double))-1;
  w = (double *)cmalloc(ntquad*sizeof(double))-1;
  cp->cp_dvr_clus.tq_all = (double *)cmalloc(my_num_t*sizeof(double))-1;
  cp->cp_dvr_clus.wq_all = (double *)cmalloc(my_num_t*sizeof(double))-1;

  tq_all = cp->cp_dvr_clus.tq_all;
  wq_all = cp->cp_dvr_clus.wq_all;

  GAUSSQ(&kind, &ntquad, &alpha, &beta, &kpts, &endpts[0], &c[1], &t[1], &w[1]);

  pwr = log10(tmax);
  low = 0.0;
  a = 0.0;
  id=0;
  for (i=1;i<=ntwind;i++){
    high = pwr-ntwind+i;
    b=pow(10.0,high);
    factor = (b-a)/2.0;
    for (j=1;j<=ntquad;j++){
      id +=1;
      if(id >= istart && id <= iend){
        tq_all[id-istart+1] = (b+a)/2.0+(b-a)/2.0*t[j];
        wq_all[id-istart+1] = w[j]*factor;
      }
    }
    low = high;
    a=pow(10.0,low);
  }

/*============================================================================*/
/* V) Allocate memory for F-matrix */


  cp->cp_dvr_clus.fmat_x =
           (double *)cmalloc(my_num_t*grid_nx*(grid_nx+1)/2*sizeof(double))-1;
  cp->cp_dvr_clus.fmat_y =
           (double *)cmalloc(my_num_t*grid_ny*(grid_ny+1)/2*sizeof(double))-1;
  cp->cp_dvr_clus.fmat_z =
           (double *)cmalloc(my_num_t*grid_nz*(grid_nz+1)/2*sizeof(double))-1;

  fmat_x = cp->cp_dvr_clus.fmat_x;
  fmat_y = cp->cp_dvr_clus.fmat_y;
  fmat_z = cp->cp_dvr_clus.fmat_z;

  memory_now = (grid_nx*(grid_nx+1) + grid_ny*(grid_ny+1) + grid_nz*(grid_nz+1))/2
               *my_num_t*(sizeof(double))*1.0e-6;

  cp->cp_dvr_clus.Fx =
           (double *)cmalloc(grid_nx*grid_nx*sizeof(double))-1;
  cp->cp_dvr_clus.Fy =
           (double *)cmalloc(grid_ny*grid_ny*sizeof(double))-1;
  cp->cp_dvr_clus.Fz =
           (double *)cmalloc(grid_nz*grid_nz*sizeof(double))-1;

  memory_now += (grid_nx*grid_nx + grid_ny*grid_ny + grid_nz*grid_nz)*(sizeof(double))
                *1.0e-6;

  *tot_memory += (memory_now);

  if(myid_state==0){
    printf("Allocation for DVR cluster boundary condition: %g Mbytes\n", memory_now);
    printf("Total memory: %g Mbytes\n",*tot_memory);
  }/*endif for myid=0*/

/*=========================================================================*/
/* VI) Compute F-matrix */

  /* X-coordinate*/

  plimit = 2.0*dxx/dx;
  plimit = plimit*plimit;

  id=0;
  for(it=1;it<=my_num_t;it++){

    tq = tq_all[it];

    dxs = dx*tq;
    sqrt_tq = sqrt(tq);

    for(i=1;i<=grid_nx;i++){
      xi = x_dvr[i];
      for(j=1;j<=i;j++){
        id += 1;
        xj = x_dvr[j];
        xbar = xi - xj;

        if(i == j){
          arg = pi/2.0/dx/tq;
          fmat_x[id] = sqrt(dx)*gerf(arg);
        }else{

          if(tq*tq >= plimit) {
            sum=0.0;
            for(k=1;k<=nx;k++){
              delta = x[k]*x[k];
              sum = sum + exp(-delta)*
                          sin(pi*(x[k]+xbar*tq)/dxs)/pi/(x[k]+xbar*tq)*sqrt(dxs)*dxx;
            }
          }else{
            ngrid = rmax/dx/tq*2.0;
            if(ngrid%2 == 0) ngrid += 1;
            dxxs = rmax/((double) ngrid);
            sum=0.0;
            r= -rmax/2.0-dxxs/2.0;
            for(k=1;k<=ngrid;k++){
              r=r+dxxs;
              delta = r*r;
              sum = sum + exp(-delta)*
                          sin(pi*(r+xbar*tq)/dxs)/pi/(r+xbar*tq)*sqrt(dxs)*dxxs;
            }
          } /* endif plimit */
          fmat_x[id] = sum/sqrt_tq;
        } /* endif diagonal */
      }
    }

  }/*end for t*/
  /* Y-Coordinate */

  if(fabs(hmat[1]-hmat[5]) < 1.0e-12 && grid_nx == grid_ny){
    if(myid_state==0) printf("Copying F-matrix for X to Y\n");
    id=0;
    for(it=1;it<=my_num_t;it++){
      for(i=1;i<=grid_nx;i++){
        for(j=1;j<=i;j++){
          id += 1;
          fmat_y[id] = fmat_x[id];
        }
      }
    }
  }else{

    plimit = 2.0*dyy/dy;
    plimit = plimit*plimit;

    id=0;
    for(it=1;it<=my_num_t;it++){

      tq = tq_all[it];

      dys = dy*tq;
      sqrt_tq = sqrt(tq);

      for(i=1;i<=grid_ny;i++){
        yi = y_dvr[i];
        for(j=1;j<=i;j++){
          id += 1;
          yj = y_dvr[j];
          ybar = yi - yj;

          if(i == j){
            arg = pi/2.0/dy/tq;
            fmat_y[id] = sqrt(dy)*gerf(arg);

          }else{

            if(tq*tq >= plimit) {
              sum=0.0;
              for(k=1;k<=ny;k++){
                delta = y[k]*y[k];
                sum = sum + exp(-delta)*
                          sin(pi*(y[k]+ybar*tq)/dys)/pi/(y[k]+ybar*tq)*sqrt(dys)*dyy;
              }
            }else{
              ngrid = rmax/dy/tq*2.0;
              if(ngrid%2 == 0) ngrid += 1;
              dyys = rmax/((double) ngrid);
              sum=0.0;
              r= -rmax/2.0-dyys/2.0;
              for(k=1;k<=ngrid;k++){
                r=r+dyys;
                delta = r*r;
                sum = sum + exp(-delta)*
                          sin(pi*(r+ybar*tq)/dys)/pi/(r+ybar*tq)*sqrt(dys)*dyys;
              }
            } /* endif plimit */
            fmat_y[id] = sum/sqrt_tq;
          } /* endif diagonal */
        }
      }

    }/*end for it*/
  }/*endif y==x */

  /* Z-coordinate */

  if(fabs(hmat[1]-hmat[9]) < 1.0e-12 && grid_nx == grid_nz){
    if(myid_state==0) printf("Copying F-matrix for X to Z\n");
    id=0;
    for(it=1;it<=my_num_t;it++){
      for(i=1;i<=grid_nz;i++){
        for(j=1;j<=i;j++){
          id += 1;
          fmat_z[id] = fmat_x[id];
        }
      }
    }
  }else{

    plimit = 2.0*dzz/dz;
    plimit = plimit*plimit;

    id=0;
    for(it=1;it<=my_num_t;it++){

      tq = tq_all[it];

      dzs = dz*tq;
      sqrt_tq = sqrt(tq);

      for(i=1;i<=grid_nz;i++){
        zi = z_dvr[i];
        for(j=1;j<=i;j++){
          id += 1;
          zj = z_dvr[j];
          zbar = zi - zj;

          if(i == j){
            arg = pi/2.0/dz/tq;
            fmat_z[id] = sqrt(dz)*gerf(arg);
          }else{

            if(tq*tq >= plimit) {
              sum=0.0;
              for(k=1;k<=nz;k++){
                delta = z[k]*z[k];
                sum = sum + exp(-delta)*
                          sin(pi*(z[k]+zbar*tq)/dzs)/pi/(z[k]+zbar*tq)*sqrt(dzs)*dzz;
              }
            }else{
              ngrid = rmax/dz/tq*2.0;
              if(ngrid%2 == 0) ngrid += 1;
              dzzs = rmax/((double) ngrid);
              sum=0.0;
              r= -rmax/2.0-dzzs/2.0;
              for(k=1;k<=ngrid;k++){
                r=r+dzzs;
                delta = r*r;
                sum = sum + exp(-delta)*
                          sin(pi*(r+zbar*tq)/dzs)/pi/(r+zbar*tq)*sqrt(dzs)*dzzs;
              }
            } /* endif plimit */
            fmat_z[id] = sum/sqrt_tq;
          } /* endif diagonal */
        }
      }

    }/*endif it*/

  }/*endif x==z */

  if(myid_state==0) printf("Done with F-matrix calculation\n\n");

  cfree(&x[1]);
  cfree(&y[1]);
  cfree(&z[1]);
  cfree(&c[1]);
  cfree(&t[1]);
  cfree(&w[1]);
  cfree(&x_dvr[1]);
  cfree(&y_dvr[1]);
  cfree(&z_dvr[1]);


/*-------------------------------------------------------------------------*/
}/* end routine */
/*=========================================================================*/



