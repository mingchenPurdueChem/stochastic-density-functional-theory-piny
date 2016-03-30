/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: cp_con.c                                     */
/*                                                                          */
/* Function performs orthogonalization of the wave functions                */
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

#define DEBUG_GS_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_shake_dvr(double *creal,int icoef_form, int icoef_orth,
                  double *vcreal,int ivcoef_form,double *creal_old,
                  int icoef_form_old, double c_tolshake,double dt,
                  CPSCR_OVMAT *cpscr_ovmat,int *ioff,double *occ,
                  double *rocc_sum,int cp_norb,int *iter_shake_cp,
                  CP_COMM_STATE_PKG *cp_comm_state_pkg,double dfact)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int i,j,is,js,nstate2,ind;
   double alpha,beta;
   double scale;

   /* Local pointer declarations */
   int np_states  = cp_comm_state_pkg->num_proc;
   int myid_state = cp_comm_state_pkg->myid;
   int nstate     = cp_comm_state_pkg->nstate;
   double *a      = cpscr_ovmat->ovlap1;
   double *b      = cpscr_ovmat->ovlap2;
   double *p1     = cpscr_ovmat->ovlap3;
   double *p2     = cpscr_ovmat->ovlap4;
   double *pscr   = cpscr_ovmat->ovlap5;  /* don't panic */
   double *xl2    = cpscr_ovmat->ovlap5;
   double *xlamb  = cpscr_ovmat->ovlap6;
   double *xlold  = cpscr_ovmat->ovlap7;

   double *rocc_sum_loc = cpscr_ovmat->ovlap8;
   double *occ_loc      = cpscr_ovmat->state_vec1;

   double occ_fact;
   int ncoef  = cp_comm_state_pkg->nstate_ncoef_proc;

/*========================================================================*/
/* 0) Checks                                                              */

  if(cp_norb>0){
    if(icoef_orth !=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in non-orthogonal form\n");
      printf("on state processor %d in cp_shake \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
  }/*endif*/

  if(np_states>1){
    if((icoef_form+ivcoef_form+icoef_form_old) !=3){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_shake \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
  }/*endif*/

/*========================================================================*/
/* I) Calculate overlap matrices between guess coeffs and old coeffs      */

  nstate2 = (nstate)*(nstate);

  for(is=1;is <= nstate2;is++){
    a[is] = b[is] = 0.0;
  }/* endfor */

  /* a) Assign diagonal elements according to occupation numbers */
  for(is=1;is <= nstate;is++){
    ind  = (is-1)*(nstate) + is;
    a[ind] = b[ind] = 1.0;
  }

  /* b) Calculate overlap integrals and construct matrices in equation */
  alpha = 1.0;
  beta  = 0.0;
  cp_ovlap_mat_same_dvr(creal,icoef_form,p1,pscr,alpha,beta,ioff,cp_comm_state_pkg);
  cp_ovlap_mat_diff_dvr(creal,icoef_form,creal_old,icoef_form_old,
                        p2,pscr,alpha,beta,ioff,cp_comm_state_pkg);

  occ_fact = ( occ[1] == 2 ? 2 : 1); /* HACK JOB for HYDROGEN */

  for(is=1;is <= nstate2; is++){
    a[is] -= p1[is]*(1.0/(occ_fact*dfact));
    b[is] -= p2[is]*(1.0/(occ_fact*dfact));
  }/*endfor*/

/*========================================================================*/
/* II) Iterative solution of matrix equation                              */

  /* a) Initialization of Lagrange multiplier matrix    */
  for(is=1; is <= nstate2; is++){
    rocc_sum_loc[is] = 1.0/2.0;
  }
  for(is=1; is <= nstate; is++){
    occ_loc[is] = 1.0;
  }
  for(is=1;is <=nstate2; is++){
    xlold[is] = xlamb[is] = rocc_sum_loc[is]*a[is];
  }/* endfor */

  /* b) Set up matrices and iteratively solve */

  cp_iter_mat_solve_shake(xlamb,xl2,xlold,a,b,p1,p2,occ_loc,rocc_sum_loc,
                          &nstate,iter_shake_cp,c_tolshake);

/*========================================================================*/
/* III) Apply constraint force to coefficients and coefficient velocities */

  alpha = 1.0;  beta = 1.0;

  cp_rotation_prim_dvr(creal_old,icoef_form_old,creal,icoef_form,xlamb,
                       alpha,beta,ioff,cp_comm_state_pkg);

  alpha = 1.0/dt;
  cp_rotation_prim_dvr(creal_old,icoef_form_old,vcreal,ivcoef_form,xlamb,
                       alpha,beta,ioff,cp_comm_state_pkg);

/*========================================================================*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_rattle_dvr(double *creal,int icoef_form, int icoef_orth,
                   double *vcreal,int ivcoef_form,CPSCR_OVMAT *cpscr_ovmat,
                   int *ioff,double *rocc_sum,
                   CP_COMM_STATE_PKG *cp_comm_state_pkg,int cp_norb,double dfact)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int is,js;
  int nstate2;
  int ind1,ind2;
  double alpha,beta;

  /*  Local pointer declarations   */
  int np_states   = cp_comm_state_pkg->num_proc;
  int myid_state  = cp_comm_state_pkg->myid;
  int nstate      = cp_comm_state_pkg->nstate;
  double *a       = cpscr_ovmat->ovlap1;
  double *b       = cpscr_ovmat->ovlap2;
  double *ylamb   = cpscr_ovmat->ovlap3;
  double *a_scr   = cpscr_ovmat->ovlap4;

/*========================================================================*/
/* 0) Checks                                                              */
  if(cp_norb>0){
    if(icoef_orth !=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in non-orthogonal form\n");
      printf("on state processor %d in cp_rattle \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
     }/*endif*/
  }/*endif*/

  if(np_states>1){
    if((icoef_form+ivcoef_form) !=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_rattle \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
  }/*endif*/

/*========================================================================*/
/* I) Calculate overlap matrices between coeffs and coeff velociites      */

   nstate2 = nstate*nstate;

   for(is=1;is <= nstate2;is++){a[is] = ylamb[is] = 0.0;}

   alpha = 1.0;
   beta  = 0.0;
   cp_ovlap_mat_diff_dvr(creal,icoef_form,vcreal,ivcoef_form,
                        a,a_scr,alpha,beta,ioff,cp_comm_state_pkg);

/*========================================================================*/
/* II) Calculate Lagrange multiplier matrix                               */

   for(is=1;is <=nstate; is++){
     for(js=1;js <=nstate; js++){
       ind1  = (is-1)*nstate + js;
       ind2 =  (js-1)*nstate + is;
       b[ind1] = a[ind2];
     }/* endfor js*/
   }/* endfor is */

   for(is=1;is <= nstate2;is++){
     ylamb[is] = -(a[is] + b[is])/(4.0*dfact);
   }

/*========================================================================*/
/* III) Apply constraint force to coefficient velocities                  */

  alpha = 1.0;  beta = 1.0;

  cp_rotation_prim_dvr(creal,icoef_form,vcreal,ivcoef_form,ylamb,
                       alpha,beta,ioff,cp_comm_state_pkg);

/*========================================================================*/
}/*end routine*/
/*========================================================================*/


