/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_utilities                                */
/*                                                                          */
/* This subprogram provides some integrator utility routines                */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_integrate_cp_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_cpke_dvr(CPCOEFFS_INFO *cpcoeffs_info,CPCOEFFS_POS_DVR *cpcoeffs_pos_dvr,
                  STAT_AVG *stat_avg, int cp_lsda,int np_states)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    int is,icoef,i,j,iii;

    double kinet_cp,kinet_cp_tmp,kinet_old,istart,iend;
    double *dvrvc_up  = cpcoeffs_pos_dvr->dvrvc_up;
    double *dvrvc_dn  = cpcoeffs_pos_dvr->dvrvc_dn;

    int ivcoef_form_up        = cpcoeffs_pos_dvr->ivcoef_form_up;
    int ivcoef_form_dn        = cpcoeffs_pos_dvr->ivcoef_form_dn;
    double *cmass             = cpcoeffs_info->cmass;
    int *cpcoeffs_ioff_up     = cpcoeffs_info->ioff_upt;
    int *cpcoeffs_ioff_dn     = cpcoeffs_info->ioff_dnt;

    int nstate_up             = cpcoeffs_info->nstate_up;
    int nstate_dn             = cpcoeffs_info->nstate_dn;
    int icmoff_up             = cpcoeffs_info->icoef_start_up-1;
    int icmoff_dn             = cpcoeffs_info->icoef_start_dn-1;

    int ncoef_up, ncoef_dn;


/*========================================================================*/
/* 0) Checks and Assigns */

  if(np_states>1){
    if((ivcoef_form_up)!=1){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP vectors are not in transposed form \n");
     printf("in get_cpke_dvr \n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((ivcoef_form_dn)!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("in get_cpke_dvr \n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/

  if(np_states==1){
     ncoef_up     = cpcoeffs_info->ncoef;
     ncoef_dn     = cpcoeffs_info->ncoef;
  }else{
     ncoef_up     = cpcoeffs_info->nstate_ncoef_proc_up;
     ncoef_dn     = cpcoeffs_info->nstate_ncoef_proc_dn;
  }/*endif*/

/*========================================================================*/
/* I) Get the CP kinetic energy                                           */

   kinet_cp   = 0.0;
   for(is=1;is<=nstate_up;is++) {
     for(i=1;i<=ncoef_up;i++) {
       icoef = i+cpcoeffs_ioff_up[is];
       kinet_cp += (dvrvc_up[icoef]*dvrvc_up[icoef]*cmass[(i+icmoff_up)]);
     }/*endfor*/
   }/*endfor*/
   if( (cp_lsda == 1) && (nstate_dn != 0) ){
     for(is=1;is<=nstate_dn;is++) {
       for(i=1;i<=ncoef_dn;i++) {
         icoef = i+cpcoeffs_ioff_dn[is];
         kinet_cp += (dvrvc_dn[icoef]*dvrvc_dn[icoef]*cmass[(i+icmoff_dn)]);
       }/*endfor*/
     }/*endfor*/
   }/* endif */

   kinet_cp /= 2.0;
   stat_avg->kinet_cp   = kinet_cp;

/*------------------------------------------------------------------------*/
   /*end routine*/}
/*========================================================================*/

