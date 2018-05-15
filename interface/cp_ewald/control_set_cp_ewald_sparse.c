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
/* Control ewald/cp g-space initialization. In sparse case we only have one */
/* cutoff, that is the large cutoff in k space		        */
/*==========================================================================*/

void control_set_cp_ewald_sparse(SIMOPTS *simopts,CELL *cell,
                          CPCOEFFS_INFO *cpcoeffs_info,
                          EWALD *ewald, CPEWALD *cpewald,CP_PARSE *cp_parse,
                          double *gmin_true,double *gmin_spl,double *gmax_spl,
                          EWD_SCR *ewd_scr,int kmax_ewd,int kmax_res,
                          double *tot_memory,int int_res_ter,
                          PART_MESH *part_mesh,ECOR *ecor, int myid,
                          int cp_lsda,int cp_min_diis,int cp_dual_grid_opt_on) 

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations:                               */

   int realSparseOpt = cpewald->realSparseOpt;
   double ecut_now;                    /* Num: Energy cutoff                 */
   double ecut_dens_cp_box_now;        /* Num: Energy cutoff                 */
   double ecut_res,ecut_tmp;
   double ecut_lg;                     /* Num: Energy cutoff for dens        */
   double ecut_sm;                     /* Num: Energy cutoff for wf          */
   double *hmati_ewd;                  /* Num: Inverse ewald h matrix        */
   double *hmati_ewd_cp;               /* Num: Inverse cp ewald h matrix     */
   double deth,deth_cp,side;           /* Num: Volumes and sizes             */
   double now_mem;                     /* Num: Memory used here              */
   int idum1=1,idum2=1,idum3=1;
   int iii;                            /* Num: Debug tool                    */
   int nmall;
   int cp_on;                          /* Num: CP on flag                    */
   int nk1,nk2,nk3;                    /* Num: Number of k vec on small grid */
   int nkf1_dens_cp_box,
       nkf2_dens_cp_box,
       nkf3_dens_cp_box;               /* Num: Number of k vec on large  grid*/
   int nkf1,nkf2,nkf3;                 /* Num: Number of k vec on mega  grid */

   int nktot,nktot_res,nktot_sm;       /* Num: Spherically cutoff k vecs     */
   int nktot_dens_cp_box;
   int ncoef,ncoef_dens_cp_box,ncoef_l;

   int ngrid_tot,nlen_pme;             /* Num: PME sizes                     */
   int ngrid_a_res,ngrid_b_res,ngrid_c_res;
   int ngrid_a, ngrid_b, ngrid_c;

   int *kmaxv;                         /* Lst: K-vector ranges               */
   int *kmax_cp_tmp,*kmaxv_res,cp_on_tmp;
   int *kmax_cp;
   int *kmaxv_dens_cp_box;
   int *kmax_cp_dens_cp_box;

   int pme_b_opt;
   double *bfact_r, *bfact_i;

   double gmin_spl_tmp,gmin_true_tmp;  /* Num : Min/Max g-vectors            */
   double gmax_spl_tmp;

   /*--------------------------------------*/
   /* Local pointers                       */

   int pme_on          = part_mesh->pme_on;
   int pme_res_on      = part_mesh->pme_res_on;
   int n_interp        = part_mesh->n_interp;
   int n_interp_res    = part_mesh->n_interp_res;
   int iperd           = cell->iperd;

   double *hmat_ewd    = cell->hmat_ewd;
   double *hmat_ewd_cp = cell->hmat_ewd_cp;

   int box_rat      =   cpewald->box_rat;
   double dbox_rat  =   cpewald->dbox_rat;

/*=======================================================================*/
/* 0) Output to screen                                                   */

   if(myid==0){
     printf("\n");PRINT_LINE_STAR;
     printf("Setting up reciprocal space\n");
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

/*=======================================================================*/
/* II) Allocate simple memory                                            */

   hmati_ewd      = (double *) cmalloc((size_t)9*sizeof(double))-1;
   hmati_ewd_cp   = (double *) cmalloc((size_t)9*sizeof(double))-1;
   kmaxv          =    (int *) cmalloc((size_t)3*sizeof(int))-1;
   kmaxv_res      =    (int *) cmalloc((size_t)3*sizeof(int))-1;
   kmax_cp_tmp    =    (int *) cmalloc((size_t)3*sizeof(int))-1;

   if(cp_on == 1){
     cpewald->kmax_cp = (int *) cmalloc((size_t)3*sizeof(int))-1;
     kmax_cp          = cpewald->kmax_cp;
     cpewald->kmax_cp_dens_cp_box = (int *) cmalloc((size_t)3*sizeof(int))-1;
     kmax_cp_dens_cp_box          = cpewald->kmax_cp_dens_cp_box;
     if(cp_dual_grid_opt_on >= 1){
       kmaxv_dens_cp_box = (int *) cmalloc((size_t)3*sizeof(int))-1;
     }/*endif*/
   }/* endif cp_on */


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

   //In the sparse case every cutoff is large cutoff. 
   //cp_pase->cp_ecut: large
   
   calc_cutoff(kmax_ewd,&ecut_now,&(cp_parse->cp_ecut),cp_on,
                kmax_cp,kmaxv,hmati_ewd_cp,deth_cp);  
   // kmaxv=2*kmax_cp. I don't need this. Recall that this function is only 
   // used in cp_on==1 case, therefore I don't have to consider pure Ewald.
   if(cp_on==1){
     kmaxv[1] = kmax_cp[1];
     kmaxv[2] = kmax_cp[2];
     kmaxv[3] = kmax_cp[3];
   }
     
   //Now I need to count the total k point. However I have *4 in 
   //ecut_now in countkvec3d. Therefore I need to use countkvec3d_sm 
   //k range: -kmax_cp->kmax_cp

   printf("!!!!!!!!!!!!!! ecut_now %.16lg\n",ecut_now);

   countkvec3d_sm(&(ewald->nktot),ecut_now,kmax_cp,hmati_ewd_cp); 

   //Now I have everything for large cutoff, I don't need 4.0* anymore

   nktot                   = ewald->nktot;
   cpcoeffs_info->ecut     = ecut_now;
   cpcoeffs_info->ncoef_l  = nktot+1;
   ncoef_l                 = nktot+1;
   ecor->ecut              = ecut_now;
   ewald->ecut             = ecut_now;
   ewald->nkc_max          = kmaxv[3];
   printf("kmaxv[1] %i %i %i kmax_cp[1] %i %i %i\n",kmaxv[1],kmaxv[2],kmaxv[3],
	  kmax_cp[1],kmax_cp[2],kmax_cp[3]);

/*----------------------------------------------------------------------*/
/* A.1) For dualing : Calculate cutoff and count kvectors for the large */
/*      box and save the small box.                                     */
if(cp_on == 1){
   switch(cp_dual_grid_opt_on){
    case 0:  
      ecut_dens_cp_box_now   = ecut_now;
      kmax_cp_dens_cp_box[1] = kmax_cp[1];
      kmax_cp_dens_cp_box[2] = kmax_cp[2];
      kmax_cp_dens_cp_box[3] = kmax_cp[3];
    break;
    case 1: //Careful, not debugged for sparse  
      ecut_dens_cp_box_now               = ecut_now;
      nktot_dens_cp_box                  = nktot;
      cpewald->nktot_dens_cp_box         = nktot;
      cpcoeffs_info->ncoef_l_dens_cp_box = nktot + 1;
      ncoef_dens_cp_box                  = nktot + 1;
      cpcoeffs_info->ecut_dens_cp_box    = ecut_now;

      kmaxv_dens_cp_box[1]   = 2*kmaxv[1];
      kmaxv_dens_cp_box[2]   = 2*kmaxv[2];
      kmaxv_dens_cp_box[3]   = 2*kmaxv[3];

      kmax_cp_dens_cp_box[1] = kmax_cp[1];
      kmax_cp_dens_cp_box[2] = kmax_cp[2];
      kmax_cp_dens_cp_box[3] = kmax_cp[3];

      nk1 = 4*(kmax_cp[1]+1);
      nk2 = 4*(kmax_cp[2]+1);
      nk3 = 4*(kmax_cp[3]+1);
/*!!!!  DERIVE THIS EXPRESSION !!!! */
/* How mamy g-vectors do you need to get the density on the large grid */
/* correctly if you fft to g-space with a spherical cut-off and fft back */
/* Glenn's analysis indicates thou shalt not truncate g space  at all   */
/*    for this option , great since this makes the pme a more important */
/*    method  */

      kmaxv[1] = 2*(box_rat*(kmax_cp[1]+1)-1);
      kmaxv[2] = 2*(box_rat*(kmax_cp[2]+1)-1);
      kmaxv[3] = 2*(box_rat*(kmax_cp[3]+1)-1);

      countkvec3d(&(ewald->nktot),ecut_now,kmaxv,hmati_ewd);

      nktot                   = ewald->nktot;
      cpcoeffs_info->ecut     = ecut_now;
      cpcoeffs_info->ncoef_l  = nktot+1;
      ncoef_l                 = nktot+1;
      ecor->ecut              = 4.0*ecut_now;
      ewald->ecut             = 4.0*ecut_now;
      ewald->nkc_max          = kmaxv[3];
    break;
#ifdef  ORIG
    case 2://Careful, not debugged for sparse
      ecut_dens_cp_box_now               = ecut_now;
      nktot_dens_cp_box                  = nktot;
      cpewald->nktot_dens_cp_box         = nktot;
      cpcoeffs_info->ncoef_l_dens_cp_box = nktot + 1;
      ncoef_dens_cp_box                  = nktot + 1;
      cpcoeffs_info->ecut_dens_cp_box    = ecut_now;

      kmaxv_dens_cp_box[1]   = 2*kmaxv[1];
      kmaxv_dens_cp_box[2]   = 2*kmaxv[2];
      kmaxv_dens_cp_box[3]   = 2*kmaxv[3];

      kmax_cp_dens_cp_box[1] = kmax_cp[1];
      kmax_cp_dens_cp_box[2] = kmax_cp[2];
      kmax_cp_dens_cp_box[3] = kmax_cp[3];

      nk1 = 4*(kmax_cp[1]+1);
      nk2 = 4*(kmax_cp[2]+1);
      nk3 = 4*(kmax_cp[3]+1);
      kmaxv[1] = 2*(box_rat*(kmax_cp[1]+1)-1);
      kmaxv[2] = 2*(box_rat*(kmax_cp[2]+1)-1);
      kmaxv[3] = 2*(box_rat*(kmax_cp[3]+1)-1);

      countkvec3d(&(ewald->nktot),ecut_now,kmaxv,hmati_ewd);

      nktot                   = ewald->nktot;
      cpcoeffs_info->ecut     = ecut_now;
      cpcoeffs_info->ncoef_l  = nktot+1;
      ncoef_l                 = nktot+1;
      ecor->ecut              = 4.0*ecut_now;
      ewald->ecut             = 4.0*ecut_now;
      ewald->nkc_max          = kmaxv[3];
    break;
#endif
#ifdef  PME
    case 2://Careful, not debugged for sparse
      ecut_dens_cp_box_now               = ecut_now;
      nktot_dens_cp_box                  = nktot;
      cpewald->nktot_dens_cp_box         = nktot;
      cpcoeffs_info->ncoef_l_dens_cp_box = nktot + 1;
      ncoef_dens_cp_box                  = nktot + 1;
      cpcoeffs_info->ecut_dens_cp_box    = ecut_now;

      kmaxv_dens_cp_box[1]   = 2*kmaxv[1];
      kmaxv_dens_cp_box[2]   = 2*kmaxv[2];
      kmaxv_dens_cp_box[3]   = 2*kmaxv[3];

      kmax_cp_dens_cp_box[1] = kmax_cp[1];
      kmax_cp_dens_cp_box[2] = kmax_cp[2];
      kmax_cp_dens_cp_box[3] = kmax_cp[3];

      calc_cutoff(kmax_ewd,&ecut_now,&(cp_parse->cp_ecut_dual_grid),cp_on,
                  kmax_cp,kmaxv,hmati_ewd,deth);  
      countkvec3d(&(ewald->nktot),ecut_now,kmaxv,hmati_ewd); 

      nktot                   = ewald->nktot;
      cpcoeffs_info->ecut     = ecut_now;
      cpcoeffs_info->ncoef_l  = nktot+1;
      ncoef_l                 = nktot+1;
      ecor->ecut              = 4.0*ecut_now;
      ewald->ecut             = 4.0*ecut_now;
      ewald->nkc_max          = kmaxv[3];
    break;
#endif
   }/*end switch*/
}/*endif cp_on*/
 
/*----------------------------------------------------------------------*/
/* B) Malloc                                                            */

   nmall =  nktot+1;   if((nmall % 2)==0){nmall++;}
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

/*------------------------------------------------------------------------*/

   if( cp_dual_grid_opt_on >= 1 && cp_on == 1){
     nmall      =  nktot_dens_cp_box+1;   if((nmall % 2)==0){nmall++;}
     now_mem    = (nmall*(sizeof(double)*0 + sizeof(int)*5))*1.e-06;
     (*tot_memory) += now_mem;

     cpewald->kastr_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int))-1;
     cpewald->kbstr_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int))-1;
     cpewald->kcstr_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int))-1;
     cpewald->ibrk1_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int))-1;
     cpewald->ibrk2_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int))-1;

     if(myid==0){
       printf("cp grid allocation: %g Mbytes; Total memory %g Mbytes\n",
               now_mem,*tot_memory);
     }/*endif*/
   }/*endif cp_dual_grid_opt_on*/

/*------------------------------------------------------------------------*/
/* C) Fill                                                                */
   //For sparse option, we don't have 4.0*ecut. So we need to use setkvec3d_sm
   //gmin_spl, gmax_spl is scaled in setkvec3d, we need to get the same scaling 
   //here. For real space nonlocal pp, we need to have gmaxTrueLg double.

   printf("nktot %i ecut_now %.16lg kmaxv %i %i %i\n",nktot,ecut_now,kmaxv[1],
	  kmaxv[2],kmaxv[3]);
   printf("boxinv %.16lg %.16lg %.16lg\n",hmati_ewd[1],hmati_ewd[2],hmati_ewd[3]);
   setkvec3d_sm(nktot,ecut_now,kmaxv,hmati_ewd,
             ewald->kastr,ewald->kbstr,ewald->kcstr,
             ewald->ibrk1,ewald->ibrk2,
             gmin_spl,gmax_spl);
   if(cp_on==1){
     ewald->kastr[nktot+1] = 0.0;
     ewald->kbstr[nktot+1] = 0.0;
     ewald->kcstr[nktot+1] = 0.0;

     *gmin_true = *gmin_spl;
     *gmin_spl *=0.75;
     *gmax_spl *= 4.0/3.0;
   }
   //printf("1111111111 gmax_spl %lg\n",*gmax_spl);
   cpewald->gmaxTrueLg = *gmax_spl*1.499;
/*------------------------------------------------------------------------*/
/* C) Fill DENS_CP_BOX                                                    */

   //We have cp_dual_grid_opt_on=0 so we dont touch this

   if( cp_dual_grid_opt_on >= 1 && cp_on == 1){
     setkvec3d(nktot_dens_cp_box,ecut_dens_cp_box_now,
               kmaxv_dens_cp_box,hmati_ewd_cp,
               cpewald->kastr_dens_cp_box,
               cpewald->kbstr_dens_cp_box,cpewald->kcstr_dens_cp_box,
               cpewald->ibrk1_dens_cp_box,cpewald->ibrk2_dens_cp_box,cp_on,
               &gmin_spl_tmp,&gmin_true_tmp,&gmax_spl_tmp);
   }/*endif cp_dual_grid_opt_on*/

   if( cp_dual_grid_opt_on == 2 && cp_on == 1){
     *gmin_spl  = gmin_spl_tmp;
     *gmin_true = gmin_true_tmp;
     *gmax_spl  = gmax_spl_tmp;
   }/*endif*/

/*=======================================================================*/
/* V) Setup PME                                                          */

   //set_pme_grid will assume we are using a small ecut_now and build a 
   //box that can contain 4*ecut_now. Let's hack this part by scale down 
   //ecut_now. This need to be tested and be careful.

   if(pme_on==1){
/*-----------------------------------------------------------------------*/
/*  A) Get grid size */

      part_mesh->nktot_pme = nktot;
      set_pme_grid_sm(ecut_now,deth,hmati_ewd,kmaxv,
                  &(part_mesh->ngrid_a),&(part_mesh->ngrid_b),
                  &(part_mesh->ngrid_c),n_interp,
                  part_mesh->kmax_pme);
      part_mesh->ecut     = ecut_now;
      ngrid_a = part_mesh->ngrid_a;
      ngrid_b = part_mesh->ngrid_b;
      ngrid_c = part_mesh->ngrid_c;
      if((ngrid_a<n_interp) || (ngrid_b<n_interp) || (ngrid_c<n_interp)){
       if(myid==0){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        printf("The PME n_interp parameter > number of grid points \n");
        printf("This is not allowed\n");      
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       }/*endif*/
       fflush(stdout);
       exit(1);      
      }/*endif*/
/*-----------------------------------------------------------------------*/
/*  B) Malloc */
      now_mem = ( (nktot)*(sizeof(double))
                +3*(n_interp) *(sizeof(double))  )*1.e-06;
      if(pme_res_on == 1) {
        now_mem += ((nktot) *(sizeof(double)) )*1.e-06;
      }/*endif*/
     *tot_memory += now_mem;

   /*------------------------*/
   /*  Malloc the bweights   */
     nmall                     =  nktot;
     part_mesh->bweight_tot    = (double *) cmalloc(nmall*sizeof(double))-1;
     bfact_r =  part_mesh->bweight_tot;
     bfact_i =  part_mesh->bweight_tot;
     if(pme_res_on == 1) {
      part_mesh->bweight_tot_res = (double *) cmalloc(nmall*sizeof(double))-1;
     }/*endif*/

   /*------------------------*/
   /*  Malloc some scratch  */
      nmall                          = n_interp;
      part_mesh->ninterp_tot_mall    = n_interp;

      part_mesh->aj    = (double *)cmalloc(nmall*sizeof(double))-1;
      part_mesh->rn    = (double *)cmalloc(nmall*sizeof(double))-1;
      part_mesh->rn1   = (double *)cmalloc(nmall*sizeof(double))-1;

   /*-----------------------*/
   /*  Output               */
     if(myid==0){
      printf("PME allocation: %g Mbytes; Total memory %g Mbytes\n",
                                           now_mem,*tot_memory);
     }/*endif*/
/*-----------------------------------------------------------------------*/
/*  C) Create PME Bweight */
      pme_b_opt = 1;
      set_pme_wght(nktot,ewald->kastr,ewald->kbstr,ewald->kcstr,
                   ngrid_a,ngrid_b,ngrid_c,
                   idum1,idum2,idum3,
                   pme_b_opt,bfact_r,bfact_i,
                   part_mesh->bweight_tot,n_interp,
                   part_mesh->aj,part_mesh->rn,part_mesh->rn1);
   }/*endif*/

/*=======================================================================*/
/* VI) Set up RESPA k vectors and RESPA PME                              */

   // We dont have RESPA so forget about this

   if(int_res_ter != 0) {
/*-----------------------------------------------------------------------*/
/*  A) Repsa K vectors */
      ecut_tmp = 0.0;cp_on_tmp=0;

      calc_cutoff(kmax_res,&ecut_res,&ecut_tmp,cp_on_tmp,
                  kmax_cp_tmp,kmaxv_res,hmati_ewd,deth);

      ecor->ecut_res   = 4.0*ecut_res;
      ewald->ecut_res  = 4.0*ecut_res;
      countkvec3d(&(ewald->nktot_res),ecut_res,kmaxv_res,hmati_ewd);
      nktot_res         = ewald->nktot_res;
      ecor->nktot_res   = nktot_res;
      setkvec3d(nktot_res,ecut_res,kmaxv_res,hmati_ewd,
                ewald->kastr_res,ewald->kbstr_res,ewald->kcstr_res,
                ewald->ibrk1_res,ewald->ibrk2_res,cp_on_tmp,
                &gmin_spl_tmp,&gmin_true_tmp,&gmax_spl_tmp);
      setkvec3d_res(kmax_res,hmati_ewd,
                    ewald->kastr,ewald->kbstr,ewald->kcstr,ewald->ibrk3,
                    nktot,nktot_res);
      if(pme_res_on==1){
/*-----------------------------------------------------------------------*/
/*  B) Set the PME GRID */
       part_mesh->nktot_pme_res = nktot_res;
       set_pme_grid(ecut_res,deth,hmati_ewd,kmaxv_res,
                  &(part_mesh->ngrid_a_res),&(part_mesh->ngrid_b_res),
                  &(part_mesh->ngrid_c_res),n_interp_res,
                  part_mesh->kmax_pme_res);
       part_mesh->ecut_res     = 4.0*ecut_res;
       ngrid_a_res = part_mesh->ngrid_a_res;
       ngrid_b_res = part_mesh->ngrid_b_res;
       ngrid_c_res = part_mesh->ngrid_c_res;
       if((ngrid_a_res<n_interp_res) || (ngrid_b_res<n_interp_res) || 
          (ngrid_c_res<n_interp_res)){
        if(myid==0){
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         printf("The RESPA PME n_interp parameter > number of grid points \n");
         printf("This is not allowed\n");      
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        }/*endif*/
        fflush(stdout);
        exit(1);      
       }/*endif*/
/*-----------------------------------------------------------------------*/
/*  C) Create PME Bweight */
       pme_b_opt = 1;
       set_pme_wght(nktot_res,ewald->kastr_res,ewald->kbstr_res,
                  ewald->kcstr_res,ngrid_a_res,ngrid_b_res,ngrid_c_res,
                  idum1,idum2,idum3,
                  pme_b_opt,bfact_r,bfact_i,
                  part_mesh->bweight_tot_res,n_interp_res,
                  part_mesh->aj,part_mesh->rn,part_mesh->rn1);
      }/*endif*/
    }/*endif*/

/*=======================================================================*/
/* VII) Setup CP coefs and k vectors                                     */

   if(cp_on == 1) {

/*--------------------------------------------------------------------*/
/*  A)  Count the k-vectors                                           */

      //I don't need to redo the count thing since it would be the 
      //same as large cutoff

      ecut_sm = ecut_dens_cp_box_now;
   
      //countkvec3d_sm(&(cpewald->nktot_sm),ecut_sm,
      //               kmax_cp_dens_cp_box,hmati_ewd_cp);
      cpewald->nktot_sm = nktot;
      nktot_sm = nktot;
      //nktot_sm = cpewald->nktot_sm;
      cpcoeffs_info->ncoef   = nktot_sm+1;
      ncoef                  = nktot_sm+1;
/*--------------------------------------------------------------------*/
/*  B)  Malloc                                                       */
      nmall =  (nktot_sm+1);if((nmall % 2)==0){nmall++;}
      cpewald->nktot_cp_sm_mall = nmall;
      now_mem = (nmall*(sizeof(double)*0 + sizeof(int)*5 ))*1.e-06;
      *tot_memory += now_mem;
      cpewald->kastr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
      cpewald->kbstr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
      cpewald->kcstr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
      cpewald->ibrk1_sm = (int *) cmalloc(nmall*sizeof(int))-1;
      cpewald->ibrk2_sm = (int *) cmalloc(nmall*sizeof(int))-1;

      nmall =  ncoef; if((nmall % 2)==0){nmall++;}
      cpcoeffs_info->cmass = (double *)cmalloc(nmall*sizeof(double))-1;

      if(myid==0){
        printf("CP allocation: %g Mbytes; Total memory %g Mbytes\n",
                 now_mem,*tot_memory);
      }/*endif*/
/*--------------------------------------------------------------------*/
/*  C)  Fill and check                                                */

      printf("ecut_sm %.16lg nktot_sm %i\n",ecut_sm,nktot_sm);
      printf("%i %i %i\n",kmax_cp_dens_cp_box[1],kmax_cp_dens_cp_box[2],kmax_cp_dens_cp_box[3]);
      printf("boxinv %.16lg %.16lg %.16lg\n",hmati_ewd_cp[1],
	     hmati_ewd_cp[2],hmati_ewd_cp[3]);
      setkvec3d_sm(nktot_sm,ecut_sm,kmax_cp_dens_cp_box,hmati_ewd_cp,
                   cpewald->kastr_sm,cpewald->kbstr_sm,cpewald->kcstr_sm,
                   cpewald->ibrk1_sm,cpewald->ibrk2_sm,
                   &(cpewald->gw_gmin),&(cpewald->gw_gmax));
      //cpewald->gmaxTrueSm = cpewald->gw_gmax*0.74;
      cpewald->gmaxTrueSm = cpewald->gw_gmax;

    if(cp_dual_grid_opt_on == 0 && cp_on == 1){
      check_kvec(ewald->nktot,ewald->kastr,ewald->kbstr,ewald->kcstr,nktot_sm,
                  cpewald->kastr_sm,cpewald->kbstr_sm,cpewald->kcstr_sm);
    }

      if(cp_dual_grid_opt_on >= 1 && cp_on == 1){
        check_kvec(cpewald->nktot_dens_cp_box,cpewald->kastr_dens_cp_box,
                   cpewald->kbstr_dens_cp_box,
                   cpewald->kcstr_dens_cp_box,nktot_sm,
                   cpewald->kastr_sm,cpewald->kbstr_sm,cpewald->kcstr_sm);
      }/*endif cp_dual_grid_opt_on */   
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

   if(myid ==0){
 /*-----------------------------------------------------------------------*/
 /* A) CP output */

     // we will have a fft prepare routine that only used in this sparse case

     if(cp_on == 1) {
       switch(cp_dual_grid_opt_on){
         case 0:
           nkf1 = 2*(kmax_cp[1]+1);
           nkf2 = 2*(kmax_cp[2]+1);
           nkf3 = 2*(kmax_cp[3]+1);
         break;
         case 1:
           nkf1 = 2*box_rat*(kmax_cp_dens_cp_box[1] + 1);
           nkf2 = 2*box_rat*(kmax_cp_dens_cp_box[2] + 1);
           nkf3 = 2*box_rat*(kmax_cp_dens_cp_box[3] + 1);

           nkf1_dens_cp_box = 2*(kmax_cp_dens_cp_box[1] + 1);
           nkf2_dens_cp_box = 2*(kmax_cp_dens_cp_box[2] + 1);
           nkf3_dens_cp_box = 2*(kmax_cp_dens_cp_box[3] + 1);
          break;
#ifdef ORIG
         case 2:
           nkf1 = 2*box_rat*(kmax_cp_dens_cp_box[1] + 1);
           nkf2 = 2*box_rat*(kmax_cp_dens_cp_box[2] + 1);
           nkf3 = 2*box_rat*(kmax_cp_dens_cp_box[3] + 1);

           nkf1_dens_cp_box = 2*(kmax_cp_dens_cp_box[1] + 1);
           nkf2_dens_cp_box = 2*(kmax_cp_dens_cp_box[2] + 1);
           nkf3_dens_cp_box = 2*(kmax_cp_dens_cp_box[3] + 1);
         break;
#endif
#ifdef  PME
         case 2:
           nkf1 = 2*(kmax_cp[1] + 1);
           nkf2 = 2*(kmax_cp[2] + 1);
           nkf3 = 2*(kmax_cp[3] + 1);

           nkf1_dens_cp_box = 2*(kmax_cp_dens_cp_box[1] + 1);
           nkf2_dens_cp_box = 2*(kmax_cp_dens_cp_box[2] + 1);
           nkf3_dens_cp_box = 2*(kmax_cp_dens_cp_box[3] + 1);
         break;
#endif
       }/*end switch*/

       nk1  = 2*(kmax_cp_dens_cp_box[1]+1);
       nk2  = 2*(kmax_cp_dens_cp_box[2]+1);
       nk3  = 2*(kmax_cp_dens_cp_box[3]+1);

       ecut_lg = ecut_now;
       printf("Your large cp-fft grid is  %d by %d by %d\n",nkf1,nkf2,nkf3);
       printf("There are  %d total k-vectors ",ncoef_l); 
       printf("upon spherical truncation. \n");
       printf("The large energy cutoff is 4*Ecut= %f Ryd\n",2.0*ecut_lg);
       putchar('\n');

       if(cp_dual_grid_opt_on >= 1){
         printf("Your large density grid for the cp-fft box is %d by %d by %d",
                 nkf1_dens_cp_box,nkf2_dens_cp_box,nkf3_dens_cp_box);
         printf("\nThere are  %d total k-vectors ",ncoef_dens_cp_box); 
         printf("upon spherical truncation. \n");
         printf("The large energy cutoff is 4*Ecut= %f Ryd\n",8.0*ecut_sm);
         putchar('\n');
       }/*endif cp_dual_grid_opt_on */

       printf("Your small cp-fft grid is  %d by %d by %d\n",nk1,nk2,nk3);
       printf("There are %d total k-vectors ",ncoef); 
       printf("upon spherical truncation. \n");
       printf("The small energy cutoff is Ecut= %f Ryd\n",2.0*ecut_sm);

     }else {
 /*-----------------------------------------------------------------------*/
 /* B) Non-CP output */
       printf("You are using %d k-vectors in your ewald sum\n",nktot);
       printf("Your reciprocal space shape: (-%d,%d) by (-%d,%d) by (-%d,%d)",
              kmaxv[1],kmaxv[1],kmaxv[2],kmaxv[2],kmaxv[3],kmaxv[3]);
       printf("\n");
     } /* endif:cp_on */

 /*-----------------------------------------------------------------------*/
 /* C) PME output  */
     if(pme_on==1){
       printf("Your particle mesh grid shape: %d by %d by %d\n",
                                                   ngrid_a,ngrid_b,ngrid_c);
       if((pme_res_on==1)&&(int_res_ter != 0)&&(kmax_res > 0)){
         printf("Your respa particle mesh grid shape: %d by %d by %d\n",
                                       ngrid_a_res,ngrid_b_res,ngrid_c_res);
       }/*endif*/
     }/*endif*/

     putchar('\n');PRINT_LINE_DASH;
     printf("Completed reciprocal space set up\n");
     PRINT_LINE_STAR;printf("\n");

   }/*endif:myid==0*/

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

void control_fft_pkg_sparse(PARA_FFT_PKG3D *cp_sclr_fft_pkg_sm,
                     PARA_FFT_PKG3D *cp_para_fft_pkg_sm,
                     PARA_FFT_PKG3D *cp_sclr_fft_pkg_dens_cp_box,
                     PARA_FFT_PKG3D *cp_para_fft_pkg_dens_cp_box,
                     PARA_FFT_PKG3D *cp_sclr_fft_pkg_lg,
                     PARA_FFT_PKG3D *cp_para_fft_pkg_lg,
		     PARA_FFT_PKG3D *cp_sclr_fft_pkg_sparse,
		     PARA_FFT_PKG3D *cp_para_fft_pkg_sparse,
                     PARA_FFT_PKG3D *pme_fft_pkg,
                     PARA_FFT_PKG3D *pme_res_fft_pkg,
                     EWALD* ewald,CPEWALD *cpewald,
                     PART_MESH *part_mesh, CPCOEFFS_INFO *cpcoeffs_info,
                     COMMUNICATE *communicate,int cp_on,int cp_lsda,
                     double *tot_memory,int int_res_ter,
                     int cp_para_opt,int cp_dual_grid_opt_on)

/*=========================================================================*/
     {/*Begin Routine*/
/*=========================================================================*/
/*    Local Variables */

  int nkf1,nkf2,nkf3,nfft_ext,iii;
  int nkf1Sm,nkf2Sm,nkf3Sm;
  int nkf1_dens_cp_box,nkf2_dens_cp_box,nkf3_dens_cp_box;
          /* used for denisty on cp grid when have 2 boxes*/

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
  int realSparseOpt      = cpewald->realSparseOpt;

/*=========================================================================*/
/* 0) Print to screen and check for nproc > nstate error */

  if(myid == 0){
    printf("\n");PRINT_LINE_STAR;
    printf("Setting up FFTs\n");
    PRINT_LINE_DASH;printf("\n");
  }/* endif myid */

  if(cp_on == 1){
    if(cp_lsda == 0){ 
     if(nstate_up < np_states){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       printf("Number of states less than number of processors\n");
       printf("If possible, reduce number of processors to be\n");
       printf("less than the number of states or run a bigger system.\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/* endif */
    } else {
     if(nstate_up + nstate_dn < np_states){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       printf("Number of states less than number of processors\n");
       printf("If possible, reduce number of processors to be\n");
       printf("less than the number of states or run a bigger system.\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/* endif */
    }/* endif lsda */
  }/* endif */

/*=========================================================================*/
/* 0.1) Set CP FFT Size  */

  if(cp_on==1){
    nkf1 = 2*(kmax_cp_dens_cp_box[1]+1);
    nkf2 = 2*(kmax_cp_dens_cp_box[2]+1);
    nkf3 = 2*(kmax_cp_dens_cp_box[3]+1);
  }/*endif*/

/*=========================================================================*/
/* I) DENS_CP_BOX CP scalar package                                        */

 if(cp_dual_grid_opt_on >= 1 && cp_para_opt == 0){/* hybrid option */

    cp_sclr_fft_pkg_dens_cp_box->nkf1       = nkf1;
    cp_sclr_fft_pkg_dens_cp_box->nkf2       = nkf2;
    cp_sclr_fft_pkg_dens_cp_box->nkf3       = nkf3;
      
    cp_sclr_fft_pkg_dens_cp_box->nktot      = nktot_dens_cp_box;
    cp_sclr_fft_pkg_dens_cp_box->ncoef      = ncoef_dens_cp_box;

    cp_sclr_fft_pkg_dens_cp_box->myid       = 0;
    cp_sclr_fft_pkg_dens_cp_box->myidp1     = 1;
    cp_sclr_fft_pkg_dens_cp_box->num_proc   = 1;
    cp_sclr_fft_pkg_dens_cp_box->comm       = communicate->comm_faux;


    create_para_fft_pkg3d(cp_sclr_fft_pkg_dens_cp_box,
                          cpewald->kastr_dens_cp_box,
                          cpewald->kbstr_dens_cp_box,
                          cpewald->kcstr_dens_cp_box,cp_dual_grid_opt_on);
  }/*endif*/


/*=========================================================================*/
/* II) DENSITY_CP_BOX  parallel package                                    */
/*       This package must be made for both hybrid and full_g options      */

  if(cp_dual_grid_opt_on >= 1){

    cp_para_fft_pkg_dens_cp_box->nkf1       = nkf1;
    cp_para_fft_pkg_dens_cp_box->nkf2       = nkf2;
    cp_para_fft_pkg_dens_cp_box->nkf3       = nkf3;
      
    cp_para_fft_pkg_dens_cp_box->nktot      = nktot_dens_cp_box;
    cp_para_fft_pkg_dens_cp_box->ncoef      = ncoef_dens_cp_box;
   
    cp_para_fft_pkg_dens_cp_box->myid       = myid_state;
    cp_para_fft_pkg_dens_cp_box->myidp1     = myid_state+1;
    cp_para_fft_pkg_dens_cp_box->num_proc   = np_states;
    cp_para_fft_pkg_dens_cp_box->comm       = communicate->comm_states;

    create_para_fft_pkg3d(cp_para_fft_pkg_dens_cp_box,
                          cpewald->kastr_dens_cp_box,
                          cpewald->kbstr_dens_cp_box,
                          cpewald->kcstr_dens_cp_box,cp_dual_grid_opt_on);
  }/*endif*/

/*=========================================================================*/
/* I) Large CP scalar package                                              */


 if(cp_on == 1 && cp_para_opt == 0){/* hybrid option */

   switch(cp_dual_grid_opt_on){
    case 0:
     nkf1 = 2*(kmax_cp[1]+1);
     nkf2 = 2*(kmax_cp[2]+1);
     nkf3 = 2*(kmax_cp[3]+1);
    break;
    case 1:
     nkf1 = 2*box_rat*(kmax_cp_dens_cp_box[1]+1);
     nkf2 = 2*box_rat*(kmax_cp_dens_cp_box[2]+1);
     nkf3 = 2*box_rat*(kmax_cp_dens_cp_box[3]+1);
    break;

#ifdef ORIG
    case 2:
     nkf1 = 2*box_rat*(kmax_cp_dens_cp_box[1]+1);
     nkf2 = 2*box_rat*(kmax_cp_dens_cp_box[2]+1);
     nkf3 = 2*box_rat*(kmax_cp_dens_cp_box[3]+1);
    break;
#endif
#ifdef  PME
    case 2:
     nkf1 = 2*(kmax_cp[1]+1);
     nkf2 = 2*(kmax_cp[2]+1);
     nkf3 = 2*(kmax_cp[3]+1);
    break;
#endif
   }/*end switch */

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

   switch(cp_dual_grid_opt_on){
    case 0:
     nkf1 = 2*(kmax_cp[1]+1);
     nkf2 = 2*(kmax_cp[2]+1);
     nkf3 = 2*(kmax_cp[3]+1);
    break;
    case 1:
     nkf1 = 2*box_rat*(kmax_cp_dens_cp_box[1]+1);
     nkf2 = 2*box_rat*(kmax_cp_dens_cp_box[2]+1);
     nkf3 = 2*box_rat*(kmax_cp_dens_cp_box[3]+1);
    break;
#ifdef ORIG
    case 2:
     nkf1 = 2*box_rat*(kmax_cp_dens_cp_box[1]+1);
     nkf2 = 2*box_rat*(kmax_cp_dens_cp_box[2]+1);
     nkf3 = 2*box_rat*(kmax_cp_dens_cp_box[3]+1);
    break;
#endif
#ifdef  PME
    case 2:
     nkf1 = 2*(kmax_cp[1]+1);
     nkf2 = 2*(kmax_cp[2]+1);
     nkf3 = 2*(kmax_cp[3]+1);
    break;
#endif
   }/*end switch */

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

    nkf1 = 2*(kmax_cp_dens_cp_box[1]+1);
    nkf2 = 2*(kmax_cp_dens_cp_box[2]+1);
    nkf3 = 2*(kmax_cp_dens_cp_box[3]+1);

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
/* III) Small  CP sparse scalar package                                    */

  if(cp_on==1&&cp_para_opt==0&&realSparseOpt==1){
    nkf1 = 2*(kmax_cp_dens_cp_box[1]+1);
    nkf2 = 2*(kmax_cp_dens_cp_box[2]+1);
    nkf3 = 2*(kmax_cp_dens_cp_box[3]+1);
    //nkf1 = 4*(kmax_cp_dens_cp_box[1]+1);
    //nkf2 = 4*(kmax_cp_dens_cp_box[2]+1);
    //nkf3 = 4*(kmax_cp_dens_cp_box[3]+1);


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

    nkf1 = 2*(kmax_cp_dens_cp_box[1]+1);
    nkf2 = 2*(kmax_cp_dens_cp_box[2]+1);
    nkf3 = 2*(kmax_cp_dens_cp_box[3]+1);

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

  if(cp_on==0 && pme_on ==1 ){

    pme_fft_pkg->nkf1       = ngrid_a;
    pme_fft_pkg->nkf2       = ngrid_b;
    pme_fft_pkg->nkf3       = part_mesh->ngrid_c;
      
    pme_fft_pkg->nktot      = ewald->nktot;
    pme_fft_pkg->ncoef      = ewald->nktot+1;

    if(pme_para_opt==2){
     pme_fft_pkg->myid       = myid_forc;
     pme_fft_pkg->myidp1     = myid_forc+1;
     pme_fft_pkg->num_proc   = np_forc;
     pme_fft_pkg->comm       = communicate->comm_forc;
   }else{
     pme_fft_pkg->myid       = 0;
     pme_fft_pkg->myidp1     = 1;
     pme_fft_pkg->num_proc   = 1;
     pme_fft_pkg->comm       = communicate->comm_faux;
   }/*endif*/

    create_para_fft_pkg3d(pme_fft_pkg,
                          ewald->kastr,ewald->kbstr,
                          ewald->kcstr,cp_dual_grid_opt_on);

    pme_fft_pkg->scale_opt = 0;
#ifdef HP_VECLIB
    if(pme_fft_pkg->igeneric_opt==0){pme_fft_pkg->scale_opt = -1;}
#endif

    if(pme_para_opt==2){
      nfft_ext  = 2*ngrid_a*(  (pme_fft_pkg->nfft_ka_proc)
                             + (pme_fft_pkg->skb_fft_ka_proc) - 1
                             + (pme_fft_pkg->ekb_fft_ka_proc) - ngrid_b
                             + (n_interp-1)*ngrid_b );
      pme_fft_pkg->nfft_size = MAX(pme_fft_pkg->nfft_size,nfft_ext);

      if(np_forc>1){create_pme_comm_full_g(n_interp,pme_fft_pkg);}
    }/*endif*/

    if((pme_para_opt==1)&&(np_forc>1)){
       create_pme_comm_hybr(np_forc,myid_forc,communicate->comm_forc,
                             pme_fft_pkg);
    }/*endif*/

  }/*endif*/


/*=========================================================================*/
/* VI) PME_RES package                                                      */

  if(cp_on==0 && int_res_ter == 1 && pme_res_on==1 && nktot_res > 0){

    pme_res_fft_pkg->nkf1       = part_mesh->ngrid_a_res;
    pme_res_fft_pkg->nkf2       = part_mesh->ngrid_b_res;
    pme_res_fft_pkg->nkf3       = part_mesh->ngrid_c_res;
      
    pme_res_fft_pkg->nktot      = nktot_res;
    pme_res_fft_pkg->ncoef      = nktot_res+1;

    if(pme_para_opt==2){
      pme_res_fft_pkg->myid       = myid_forc;
      pme_res_fft_pkg->myidp1     = myid_forc+1;
      pme_res_fft_pkg->num_proc   = np_forc;
      pme_res_fft_pkg->comm       = communicate->comm_forc;
    }else{
      pme_res_fft_pkg->myid       = 0;
      pme_res_fft_pkg->myidp1     = 1;
      pme_res_fft_pkg->num_proc   = 1;
      pme_res_fft_pkg->comm       = communicate->comm_faux;
    }/*endif*/

    create_para_fft_pkg3d(pme_res_fft_pkg,
                          ewald->kastr_res,ewald->kbstr_res,
                          ewald->kcstr_res,cp_dual_grid_opt_on);

    pme_res_fft_pkg->scale_opt = 0;
#ifdef HP_VECLIB
    if(pme_res_fft_pkg->igeneric_opt==0){pme_res_fft_pkg->scale_opt = -1;}
#endif

    if(np_forc>1 && pme_para_opt==2){
      create_pme_comm_full_g(n_interp_res,pme_res_fft_pkg);
    }/*endif*/

    pme_fft_pkg->nfft_size = MAX(pme_fft_pkg->nfft_size,
                                 pme_res_fft_pkg->nfft_size);

    if((pme_para_opt==1)&&(np_forc>1)){ 
      create_pme_comm_hybr(np_forc,myid_forc,communicate->comm_forc,
                           pme_res_fft_pkg);
    }/*endif*/

  }/*endif*/

/*=========================================================================*/
/* VI) Output */

  if(myid == 0){
    printf("\n");PRINT_LINE_DASH;
    printf("Finished setting up FFTs\n");
    PRINT_LINE_STAR;printf("\n");
  }/* endif myid */

/*-------------------------------------------------------------------------*/
     }/*end routine */
/*=========================================================================*/

