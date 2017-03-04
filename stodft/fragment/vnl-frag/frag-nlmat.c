/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: frag-nlmat.c                                 */
/*                                                                          */
/* This routine calculate the kinectic energy, non-local pseudo-potential   */
/* energy and nuclei force correction.                                      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_frag_local.h"
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
calcNonLocalMatrix(CP *cp, CP *cpMini)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
/**************************************************************************/
/* This function calculate non-local pseudopotential matrix w.r.t. frag-  */
/* -ment MO, as well as the force component.				  */
/**************************************************************************/

/*======================================================================*/
/* 0) Check the forms                                                   */

  if(icoef_orth_up!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The Up coefficients must be in orthogonal form    \n");
    printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/
  if(cp_lsda==1 && nstate_dn != 0){
   if(icoef_orth_dn!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The dn coefficients must be in orthogonal form    \n");
    printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
  }/*endif*/

  if(np_states>1){
   if((icoef_form_up+ifcoef_form_up)!=0){
    //printf("icoef_form_up %i ifcoef_form_up %i\n",icoef_form_up,ifcoef_form_up);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The up coefs and coef forces must not be in transposed form\n");
    printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
   if(cp_lsda==1 && nstate_dn != 0){
    if((icoef_form_dn+ifcoef_form_dn)!=0){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("The dn coefs and coef forces must not be in transposed form\n");
     printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
    }/*endif*/
   }/*endif*/
  }/*endif*/

/*======================================================================*/
/* III) Determine the maximum open non-local angular momentum channel   */
/*      for Kleinman-Bylander and Goedecker pseudo potentials           */

  nl_max_kb = -1;
  for(i=1;i<=(n_ang_max_kb+1);i++){
   if(np_nl[i]>0){nl_max_kb=i-1;}
  }/*endfor*/

  nl_max_gh = -1;
  for(i=1;i<=(n_ang_max_gh+1);i++){
   if(np_nl_gh[i]>0){nl_max_gh=i-1;}
  }/*endfor*/

/*======================================================================*/
/* IV) Determine the maximum number of atoms in any                     */
/*       open angular momentum channel                                  */

  np_nlmax_kb = 1;
  for(i = 1;i<=(nl_max_kb+1);i++){
   np_nlmax_kb = MAX(np_nlmax_kb,np_nl[i]);
  }/*endfor*/

  np_nlmax_gh = 1;
  for(i = 1;i<=(nl_max_gh+1);i++){
   np_nlmax_gh = MAX(np_nlmax_gh,np_nl_gh[i]);
  }/*endfor*/

  np_nlmax_all = (np_nlmax_gh > np_nlmax_kb ? np_nlmax_gh : np_nlmax_kb);




/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

