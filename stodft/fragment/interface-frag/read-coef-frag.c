/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: read_coef.c                                  */
/*                                                                          */
/* This subprogram provides output for a MD on a                            */
/* LD-classical potential energy surface (LD-PES)                           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_coords_cp_entry.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

void read_coef_alloc_init(CP *, int ,double *);
void read_coef_init_nhc(CP *);

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void readCoefFrag(CP *cp,GENERAL_DATA *general_data,CLASS *class,
               FILENAME_PARSE *filename_parse,
               CP_PARSE *cp_parse,double *tot_memory)

/*==========================================================================*/
/*               Begin subprogram:                                          */
  {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  FILE *fp_dnameci;   

/* Local pointers */
  SIMOPTS *simopts       = &(general_data->simopts);
  int istart             = cp_parse->istart_cp;
  int myid               = cp->communicate.myid;
  MPI_Comm world         = cp->communicate.world;
  int cp_wave_min        = simopts->cp_wave_min;
  int cp_wave_min_pimd   = simopts->cp_wave_min_pimd;
  int cp_min             = simopts->cp_min;
  int cp_min_on;
  int initial_spread_opt = simopts->initial_spread_opt;
  int np_states          = cp->communicate.np_states;
  int num_proc           = cp->communicate.np;
  char *dnameci          = filename_parse->dnameci;

  int ibinary            = cp->cpopts.iread_coef_binary;
  int iii;

  int iopt_cp_pw         = cp->cpcoeffs_info.iopt_cp_pw;
  int iopt_cp_dvr        = cp->cpcoeffs_info.iopt_cp_dvr;

  cp_min_on = cp_wave_min + cp_wave_min_pimd + cp_min;

/*========================================================================*/
/*  V) Allocate and initialize coefficient arrays                         */

  if(iopt_cp_pw){
    read_coef_alloc_init(cp,cp_min_on,tot_memory);
    if(num_proc>1){Barrier(world);}
  }

/*========================================================================*/
/*  VII) Read/Spread the coefficients  */

  if(iopt_cp_pw==1){
    gen_wave_frag(class,general_data,cp,cp_parse,filename_parse->vps_name);
  }

  
  read_coef_init_nhc(cp);

/*-----------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coef_alloc_init(CP *cp,int cp_min_on,double *tot_memory)
              
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*    Local Variables   */
#include "../typ_defs/typ_mask.h"
  int myid = cp->communicate.myid;
  int i,ip,nread,is;
  int par_size_up,par_size_dn,ncoef_up_tot;
  int nstate,nstate2,ncoef_dn_tot;
  double *cre,*cim,*vcre,*vcim,*fcre,*fcim;
  int *ioff_up,*ioff_upt,*ioff_dn,*ioff_dnt;
  double mem_test;

/*  Local Pointers */
  int pi_beads       = cp->cpcoeffs_info.pi_beads;
  int pi_beads_proc  = cp->cpcoeffs_info.pi_beads_proc;
  int nstate_up_proc = cp->cpcoeffs_info.nstate_up_proc;
  int nstate_dn_proc = cp->cpcoeffs_info.nstate_dn_proc;
  int nstate_up      = cp->cpcoeffs_info.nstate_up;
  int nstate_dn      = cp->cpcoeffs_info.nstate_dn;
  int cp_lda         = cp->cpopts.cp_lda;
  int cp_lsda        = cp->cpopts.cp_lsda;
  int ncoef_up       = cp->cpcoeffs_info.ncoef;
  int ncoef_dn       = cp->cpcoeffs_info.ncoef;
  int cp_norb        = cp->cpopts.cp_norb;

  int np_states               = cp->communicate.np_states;
  int ncoef_up_proc           = cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
  int ncoef_dn_proc           = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;
  int nstate_max_up           =cp->cp_comm_state_pkg_up.nstate_max;
  int nstate_ncoef_proc_max_up=cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
  int nstate_max_dn           =cp->cp_comm_state_pkg_dn.nstate_max;
  int nstate_ncoef_proc_max_dn=cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;

  if(cp_lda == 1){ncoef_dn = 0;}



/*==========================================================================*/
/* 0) Calculate the sizes */

  par_size_up = nstate_max_up*(nstate_ncoef_proc_max_up);
  ncoef_up_tot = nstate_up_proc*ncoef_up;
  ncoef_up_tot = MAX(ncoef_up_tot,par_size_up);

  par_size_dn = nstate_max_dn*(nstate_ncoef_proc_max_dn);
  if(cp_lda ==1){par_size_dn = 0;}
  ncoef_dn_tot = MAX(nstate_dn_proc*ncoef_dn,1);
  ncoef_dn_tot = MAX(ncoef_dn_tot,par_size_dn);

  nstate  = MAX(nstate_up,nstate_dn);
  nstate2 = nstate*nstate;

/*==========================================================================*/
/* I) Malloc the variables */

 for(i=1;i<=pi_beads_proc;i++){
  cp->cpcoeffs_pos[i].cre_up =(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].cim_up =(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].cre_dn =(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].cim_dn =(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].fcre_up=(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].fcim_up=(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].fcre_dn=(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].fcim_dn=(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].vcre_up=(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].vcim_up=(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].vcre_dn=(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].vcim_dn=(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  // Test only
  cp->cpcoeffs_pos[i].kfcre_up = (double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].kfcim_up = (double *)cmalloc(ncoef_up_tot*sizeof(double))-1;


  cp->cpcoeffs_pos[i].ksmat_up    =(double *)cmalloc(nstate2*sizeof(double))-1;
  cp->cpcoeffs_pos[i].ksmat_dn    =(double *)cmalloc(nstate2*sizeof(double))-1;
  cp->cpcoeffs_pos[i].ksmat_eig_up=(double *)cmalloc(nstate*sizeof(double))-1;
  cp->cpcoeffs_pos[i].ksmat_eig_dn=(double *)cmalloc(nstate*sizeof(double))-1;
  cp->cpcoeffs_pos[i].norbmat_up  =(double *)cmalloc(nstate2*sizeof(double))-1;
  cp->cpcoeffs_pos[i].norbmat_dn  =(double *)cmalloc(nstate2*sizeof(double))-1;
  cp->cpcoeffs_pos[i].norbmati_up =(double *)cmalloc(nstate2*sizeof(double))-1;
  cp->cpcoeffs_pos[i].norbmati_dn =(double *)cmalloc(nstate2*sizeof(double))-1;
  cp->cpcoeffs_pos[i].ovmat_eigv_up=(double *)
                                  cmalloc(nstate2*sizeof(double))-1;
  cp->cpcoeffs_pos[i].ovmat_eigv_dn=(double *)
                                  cmalloc(nstate2*sizeof(double))-1;
  if(cp_min_on > 0){  /* Then we need to allocate the diagonal Hessian */
    cp->cpcoeffs_pos[i].cp_hess_re_up = (double *) cmalloc(ncoef_up*sizeof(double))-1;
    cp->cpcoeffs_pos[i].cp_hess_im_up = (double *) cmalloc(ncoef_up*sizeof(double))-1;
    cp->cpcoeffs_pos[i].cp_hess_re_dn = (double *) cmalloc(ncoef_dn*sizeof(double))-1;
    cp->cpcoeffs_pos[i].cp_hess_im_dn = (double *) cmalloc(ncoef_dn*sizeof(double))-1;
  }/* endif cp_min_on */
 }/*endfor*/

  cp->cpcoeffs_info.ioff_up  = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_dn  = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_upt = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_dnt = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpopts.occ_up          = (double *) cmalloc(nstate*sizeof(double))-1;
  cp->cpopts.occ_dn          = (double *) cmalloc(nstate*sizeof(double))-1;
  cp->cpopts.rocc_sum_up = (double *) cmalloc(nstate*nstate*sizeof(double))-1;
  cp->cpopts.rocc_sum_dn = (double *) cmalloc(nstate*nstate*sizeof(double))-1;

/*==========================================================================*/
/* II) Assign the offsets */

  ioff_up  = cp->cpcoeffs_info.ioff_up;
  ioff_upt = cp->cpcoeffs_info.ioff_upt;
  ioff_dn  = cp->cpcoeffs_info.ioff_dn;
  ioff_dnt = cp->cpcoeffs_info.ioff_dnt;

  for(is=1;is<=nstate;is++){ioff_up[is]=(is-1)*ncoef_up;}
  for(is=1;is<=nstate;is++){ioff_dn[is]=(is-1)*ncoef_dn;}

  if(np_states==1){
   for(is=1;is<=nstate;is++){ioff_upt[is]=(is-1)*ncoef_up;}
   for(is=1;is<=nstate;is++){ioff_dnt[is]=(is-1)*ncoef_dn;}
  }else{
   for(is=1;is<=nstate;is++){ioff_upt[is]=(is-1)*ncoef_up_proc;}
   for(is=1;is<=nstate;is++){ioff_dnt[is]=(is-1)*ncoef_dn_proc;}
  }/*endif*/

/*========================================================================*/
/* III) Initialize coeficient arrays and form flags                       */

 for(ip=1;ip<=pi_beads_proc;ip++){

   cp->cpcoeffs_pos[ip].icoef_form_up  = 0;
   cp->cpcoeffs_pos[ip].ivcoef_form_up = 0;
   cp->cpcoeffs_pos[ip].ifcoef_form_up = 1;
   cp->cpcoeffs_pos[ip].icoef_orth_up  = 1;
   if(cp_norb>0){cp->cpcoeffs_pos[ip].icoef_orth_up = 0;}
   cp->cpcoeffs_pos[ip].ivcoef_orth_up  = 1;
   if(cp_norb>0){cp->cpcoeffs_pos[ip].ivcoef_orth_up = 0;}

   cre = cp->cpcoeffs_pos[ip].cre_up;
   cim = cp->cpcoeffs_pos[ip].cim_up;
   vcre = cp->cpcoeffs_pos[ip].vcre_up;
   vcim = cp->cpcoeffs_pos[ip].vcim_up;
   fcre = cp->cpcoeffs_pos[ip].fcre_up;
   fcim = cp->cpcoeffs_pos[ip].fcim_up;
   nread = ncoef_up_tot;
   for(i=1;i<=nread;i++){
    cre[i] = 0.0;
    cim[i] = 0.0;
    vcre[i] = 0.0;
    vcim[i] = 0.0;
    fcre[i] = 0.0;
    fcim[i] = 0.0;
   }/*endfor*/

   if(cp_lsda==1){

    cp->cpcoeffs_pos[ip].icoef_form_dn  = 0; 
    cp->cpcoeffs_pos[ip].ivcoef_form_dn = 0;
    cp->cpcoeffs_pos[ip].ifcoef_form_dn = 1;
    cp->cpcoeffs_pos[ip].icoef_orth_dn  = 1;
    if(cp_norb>0){cp->cpcoeffs_pos[ip].icoef_orth_dn = 0;}
    cp->cpcoeffs_pos[ip].ivcoef_orth_dn  = 1;
    if(cp_norb>0){cp->cpcoeffs_pos[ip].ivcoef_orth_dn = 0;}

    cre = cp->cpcoeffs_pos[ip].cre_dn;
    cim = cp->cpcoeffs_pos[ip].cim_dn;
    vcre = cp->cpcoeffs_pos[ip].vcre_dn;
    vcim = cp->cpcoeffs_pos[ip].vcim_dn;
    fcre = cp->cpcoeffs_pos[ip].fcre_dn;
    fcim = cp->cpcoeffs_pos[ip].fcim_dn;
    nread = ncoef_dn_tot;
    for(i=1;i<=nread;i++){
     cre[i] = 0.0;
     cim[i] = 0.0;
     vcre[i] = 0.0;
     vcim[i] = 0.0;
     fcre[i] = 0.0;
     fcim[i] = 0.0;
    }/*endfor*/ 

   }/*endif:lsda*/

 }/*endfor:pi_bead*/

/*-----------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coef_init_nhc(CP *cp)
              
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

 int ichain,inhc,ip;
 double **c_nhc;

/* Local pointers */
 int pi_beads         = cp->cpcoeffs_info.pi_beads;
 int pi_beads_proc    = cp->cpcoeffs_info.pi_beads_proc;
 int num_c_nhc        = cp->cptherm_info.num_c_nhc;
 int num_c_nhc_proc   = cp->cptherm_info.num_c_nhc_proc;
 int massiv_flag      = cp->cptherm_info.massiv_flag;
 int len_c_nhc        = cp->cptherm_info.len_c_nhc;

/*==========================================================================*/
/* I) standard thermos */

    if(massiv_flag == 0 && num_c_nhc >0){

     for(ip=1;ip<=pi_beads_proc;ip++){
      c_nhc = cp->cptherm_pos[ip].c_nhc;
      for(ichain=1;ichain<=len_c_nhc;ichain++){
       for(inhc=1;inhc<=num_c_nhc_proc;inhc++){
         c_nhc[ichain][inhc] = 0.0;
       } /*endfor*/
      } /*endfor*/
     } /*endfor*/

    }/*endif*/

/*==========================================================================*/
/* II) Massive thermos */

   if(massiv_flag != 0 && num_c_nhc >0){

     for(ip=1;ip<=pi_beads_proc;ip++){
        cp->cptherm_pos[ip].c_nhc_massiv = 0.0;
     }/*endfor*/

    }/*endif*/

/*-----------------------------------------------------------------------*/
  } /* end routine */
/*==========================================================================*/

