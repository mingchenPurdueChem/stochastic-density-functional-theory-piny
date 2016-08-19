/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: cp_con_utils.c                               */
/*                                                                          */
/* Contains utilities for wave function S matrix evals and constraint solvers*/
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define HAND_MULT_OFF
#define MAX_ITER 1000

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_ovlap_mat_same_dvr(double *cr,int icoef_form,
                           double *omat,double *omat_tmp,double alpha_in,
                           double beta_in,int *ioff,
                           CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int ncoefm1,iii;
  int nstate2;
  int is,js;
  int ind_ij,ind_nc_i,ind_nc_j;
  int np_states  = cp_comm_state_pkg->num_proc;
  int myid_state = cp_comm_state_pkg->myid;
  int nstate     = cp_comm_state_pkg->nstate;
  int ncoef      = cp_comm_state_pkg->ncoef;
  int itransp = 0;
  int inormal = 1;
  double alpha=alpha_in;
  double beta = beta_in;
  MPI_Comm comm_states = cp_comm_state_pkg->comm;
  MPI_Comm world       = cp_comm_state_pkg->world;

/*========================================================================*/
/* I) Serial matrix multiplies                                            */

  if(np_states == 1){

/*------------------------------------------------------------------------*/
   /* A) Zero matrix                                    */

    nstate2=nstate*nstate;
    for(is=1;is <= nstate2;is++){
      omat[is] = 0.0;
    }/* endfor */

/*------------------------------------------------------------------------*/
   /* B) Sum over r-space                              */

    GEN_MATMUL(&(cr[1]),&ncoef,&itransp,&(cr[1]),&ncoef,&inormal,
               &(omat[1]),&nstate,&nstate,&nstate,&ncoef,
               &alpha,&beta);

  }/*endif*/
/*========================================================================*/
/* II) Parallel matrix multiplies with g=0 correction                    */

  if(np_states != 1){

    if(icoef_form !=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_ovlmat_mat_same \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/

    nstate2=nstate*nstate;
    for(is=1;is <= nstate2;is++){omat_tmp[is] = 0.0;}

    cp_par_ovlap_same_dvr(cr,icoef_form,omat_tmp,alpha,beta,cp_comm_state_pkg);

    Allreduce(&(omat_tmp[1]),&(omat[1]),nstate2,MPI_DOUBLE,MPI_SUM,0,comm_states);

  }/* endif */


/*========================================================================*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_par_ovlap_same_dvr(double *c,int icoef_form,
                           double *omat,double alpha,double beta,
                           CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  double wght_use,alpha_loc,beta_loc;
  int ioff,joff,is,js,ig,iii;
  int ind,jnd;
  int  itransp = 0;
  int  inormal = 1;

  int nstate                = cp_comm_state_pkg->nstate;
  int nstate_max            = cp_comm_state_pkg->nstate_max;
  int ncoef                 = cp_comm_state_pkg->ncoef;
  int nstate_proc           = cp_comm_state_pkg->nstate_proc;
  int nstate_proc_max       = cp_comm_state_pkg->nstate_proc_max;
  int nstate_ncoef_proc     = cp_comm_state_pkg->nstate_ncoef_proc;
  int num_proc              = cp_comm_state_pkg->num_proc;
  int myid                  = cp_comm_state_pkg->myid;
  MPI_Comm comm  = cp_comm_state_pkg->comm;
  MPI_Comm world = cp_comm_state_pkg->world;

/*========================================================================*/
/* I) Check the forms */

  if(icoef_form !=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed form\n");
    printf("on state processor %d in cp_par_ovlap_mat_same \n",myid);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/*========================================================================*/
/* II) Matrix multiply */


#ifdef HAND_MULT
  ioff = 0;
  for(is=1;is<=nstate;is++){
    joff = 0;
    for(js=1;js<=is;js++){
      ind  = (is-1)*nstate + js;
      for(ig=1;ig<=nstate_ncoef_proc;ig++){
        omat[ind] += (c[(joff+ig)]*c[(ioff+ig)]);
      }/*endfor*/
      jnd = (js-1)*nstate + is;
      omat[jnd] = omat[ind];
      joff += nstate_ncoef_proc;
    }/*endfor*/
    ioff += nstate_ncoef_proc;
  }/*endfor*/
#else
  alpha_loc = 1.0;
  beta_loc  = beta;
  GEN_MATMUL(&(c[1]),&nstate_ncoef_proc,&itransp,&(c[1]),
             &nstate_ncoef_proc,&inormal,
             &(omat[1]),&nstate,&nstate,&nstate,&nstate_ncoef_proc,
             &alpha_loc,&beta_loc);
#endif

/*========================================================================*/
 }/* end routine */
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_ovlap_mat_diff_dvr(double *cr1,int i1coef_form,
                           double *cr2,int i2coef_form,
                           double *omat,double *omat_tmp,
                           double alpha_in,double beta_in,int *ioff,
                           CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int ncoefm1,iproc,iii;
  int nstate2;
  int is,js;
  int ind_ij,ind_nc_i,ind_nc_j;
  int np_states  = cp_comm_state_pkg->num_proc;
  int myid_state = cp_comm_state_pkg->myid;
  int nstate     = cp_comm_state_pkg->nstate;
  int ncoef      = cp_comm_state_pkg->ncoef;
  int itransp    = 0;
  int inormal    = 1;
  double alpha=alpha_in;
  double beta = beta_in;
  MPI_Comm comm_states = cp_comm_state_pkg->comm;
  MPI_Comm world       = cp_comm_state_pkg->world;

  double tot;

/*========================================================================*/
/*========================================================================*/
/* I) Serial matrix multiplies                                            */


 if(np_states == 1){

/*------------------------------------------------------------------------*/
/* A) Zero matrix                                    */

   nstate2=nstate*nstate;
   for(is=1;is <= nstate2;is++){
      omat[is] = 0.0;
   }/* endfor */

/*------------------------------------------------------------------------*/
/* B) Sum over r-space                              */

   GEN_MATMUL(&(cr1[1]),&ncoef,&itransp,&(cr2[1]),&ncoef,&inormal,
              &(omat[1]),&nstate,&nstate,&nstate,&ncoef,
              &alpha,&beta);

 }/*endif*/

/*========================================================================*/
/* II) Parallel matrix multiplies                   */

 if(np_states != 1){

   if((i1coef_form+i2coef_form) !=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_ovlmat_mat_diff \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
   }/*endif*/

  nstate2=nstate*nstate;
  for(is=1;is <= nstate2;is++){omat_tmp[is] = 0.0;}

  cp_par_ovlap_diff_dvr(cr1,i1coef_form,cr2,i2coef_form,omat_tmp,
                        alpha,beta,cp_comm_state_pkg);
  Allreduce(&(omat_tmp[1]),&(omat[1]),nstate2,MPI_DOUBLE,MPI_SUM,0,
                                                       comm_states);

 }/* endif */


/*========================================================================*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_par_ovlap_diff_dvr(double *c1,int i1coef_form,
                           double *c2,int i2coef_form,
                           double *omat,double alpha,double beta,
                           CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

 double alpha_loc,beta_loc;
 int ioff,joff,is,js,ig,iproc,iii;
 int ind,jnd;
 int itransp = 0;
 int inormal = 1;

 int nstate                = cp_comm_state_pkg->nstate;
 int nstate_max            = cp_comm_state_pkg->nstate_max;
 int ncoef                 = cp_comm_state_pkg->ncoef;
 int nstate_proc           = cp_comm_state_pkg->nstate_proc;
 int nstate_proc_max       = cp_comm_state_pkg->nstate_proc_max;
 int nstate_ncoef_proc_max = cp_comm_state_pkg->nstate_ncoef_proc_max;
 int nstate_ncoef_proc     = cp_comm_state_pkg->nstate_ncoef_proc;
 int num_proc              = cp_comm_state_pkg->num_proc;
 int myid                  = cp_comm_state_pkg->myid;
 MPI_Comm comm   = cp_comm_state_pkg->comm;
 MPI_Comm world  = cp_comm_state_pkg->world;
/*========================================================================*/
/* I) Check the forms */

    if((i1coef_form+i2coef_form) !=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_par_ovlap_mat_diff \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/

/*========================================================================*/
/* II) Matrix multiply with correct treatment for g=0                     */


#ifdef HAND_MULT
  ioff = 0;
  for(is=1;is<=nstate;is++){
   joff = 0;
   for(js=1;js<=nstate;js++){
     ind  = (is-1)*nstate + js;
     for(ig=1;ig<=nstate_ncoef_proc;ig++){
       omat[ind] += alpha*(c1[(joff+ig)]*c2[(ioff+ig)]);
     }
     joff += nstate_ncoef_proc;
   }
   ioff += nstate_ncoef_proc;
  }
#else


  alpha_loc = alpha;
  beta_loc  = beta;

  GEN_MATMUL(&(c1[1]),&nstate_ncoef_proc,&itransp,&(c2[1]),
             &nstate_ncoef_proc,&inormal,
             &(omat[1]),&nstate,&nstate,&nstate,&nstate_ncoef_proc,
             &alpha_loc,&beta_loc);

#endif

/*========================================================================*/
    }/* end routine */
/*========================================================================*/

