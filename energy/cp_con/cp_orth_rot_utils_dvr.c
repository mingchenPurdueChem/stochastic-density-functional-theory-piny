/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: cp_rotbasis_utils.c                          */
/*                                                                          */
/* Contains utilities rotating into and out of various basises              */
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

#define HAND_ROT_OFF

/*==========================================================================*/
/* Control gram_schmidt orthog */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_cp_gram_schmidt_dvr(double *creal,int icoef_form,
                                 double *occ,double *omat, double *omat_tmp,
                                 int *ioff, double scale_fact,
                                 CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nstate   = cp_comm_state_pkg->nstate;
  int ncoef    = cp_comm_state_pkg->ncoef;
  int num_proc = cp_comm_state_pkg->num_proc;
  int myid_state= cp_comm_state_pkg->myid;

/*========================================================================*/
/* 0) Do nothing if nstate = 0                                            */

  if(nstate == 0) return;

/*========================================================================*/
/* I) Scalar  GS                                                          */

  if(num_proc ==1){
    cp_gram_schmidt_scalar_dvr(creal,occ,omat,ioff,scale_fact,nstate,ncoef);
  }/*endif*/

/*========================================================================*/
/* II) Parallel GS                                                        */

  if(num_proc > 1){
    if(icoef_form !=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in control_cp_gram_schmidt \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
    cp_gram_schmidt_par_dvr(creal,icoef_form,occ,ioff,omat,omat_tmp,
                            scale_fact,cp_comm_state_pkg);
  }/*endif*/

/*------------------------------------------------------------------------*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/* Gram_schmidt orthog in scalar*/
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_gram_schmidt_scalar_dvr(double *cre,double *occ,
                                double *ovlap,int *ioff,double scale_fact,
                                int nstate,int ncoef)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */

  double cmix;
  double anorm,scale;
  double dvr_fact;
  int len,ig,is,js,joff;
  int i1,i2;
  int ind,jnd,nstate2;
  double tot;
/*==========================================================================*/
/* I) Normalize the states                                                 */

/*--------------------------------------------------------------------------*/
/* i) Find overlap of the each state with iself                            */

  for(is=1;is<=nstate;is++){
    ovlap[is]    = 0.0;
    i1 = ioff[is]+1; i2=ioff[is]+ncoef;
    for(ig=i1;ig<=i2;ig++){
      ovlap[is] += (cre[ig]*cre[ig]);
    }/*endfor*/
  }/*endfor:is*/

/*-------------------------------------------------------------------------*/
/* ii) Normalize each state to 1 for convenience                           */

  for(is=1;is<=nstate;is++){
    scale = sqrt(1.0/ovlap[is]);
    i1 = ioff[is]+1; i2=ioff[is]+ncoef;
    for(ig=i1;ig<=i2;ig++){
      cre[ig] *= scale;
    }/*endfor:ig*/
  }/*endfor:is*/

/*========================================================================*/
/* III) Gram schmidt orthogonalize                                        */

  for(is=1;is<=nstate-1;is++){

/*-----------------------------------------------------------------------*/
/* i) Find overlap of the is^th state with all higher states             */
/*    Note the joff value                                                */

   len = nstate-is;
   for(js=1;js<=len;js++){
     joff = ioff[(is+js)];
     ovlap[js]     = 0.0;
     for(ig=1;ig<=ncoef;ig++){
      ovlap[js]     += (cre[(ig+ioff[is])]*cre[(ig+joff)]);
     }/*endfor:ig*/
   }/*endfor:js*/

/*-----------------------------------------------------------*/
/* ii) Orthogonalize all higher states to the is^th state    */
/*     keeping each state properly normalized                */

   for(js=1;js<=len;js++){
     joff = ioff[(is+js)];
     cmix  = -ovlap[js];
     anorm = 1.0 + 2.0*cmix*ovlap[js] + cmix*cmix;
     scale = sqrt(1.0/anorm);
     for(ig=1;ig<=ncoef;ig++){
       cre[(ig+joff)] = (cre[(ig+ioff[is])]*cmix + cre[(ig+joff)])*scale;
     }/*endfor:ig*/
   }/*endfor:js*/

  }/*endfor:is loop*/


/*========================================================================*/
/* III) Adjust norms to occupation numbers                                */

  for(is=1;is<=nstate;is++){
    scale = sqrt(occ[is]/scale_fact);
    i1 = ioff[is]+1; i2=ioff[is]+ncoef;
    for(ig=i1;ig<=i2;ig++){
      cre[ig] *= scale;
    }/*endfor:ig*/
  }/*endfor:is*/

/*========================================================================*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/* Gram_schmidt orthog in parallel */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_gram_schmidt_par_dvr(double *cre,int icoef_form,
                             double *occ,int *ioff_vec,
                             double *ovlap, double *ovlap_tmp,
                             double scale_fact,
                             CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*======================================================================*/
/*                Begin Routine */
   {/*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */
#include "../typ_defs/typ_mask.h"

  double wght;
  double cmix;
  double anorm,scale;
  int len,ig,is,ioff,js,joff;
  int iii,jnd,nstate2,ind,i1,i2;

  int nstate                = cp_comm_state_pkg->nstate;
  int nstate_ncoef_proc     = cp_comm_state_pkg->nstate_ncoef_proc;
  int num_proc              = cp_comm_state_pkg->num_proc;
  int myid                  = cp_comm_state_pkg->myid;
  MPI_Comm comm             = cp_comm_state_pkg->comm;

/*========================================================================*/
/* 0) Checks                                                              */

  if(icoef_form !=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed form\n");
    printf("on state processor %d in cp_gram_schmidt_par \n",myid);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/*==========================================================================*/
/* I) Define useful constant                                                */

  wght = (myid != num_proc-1 ? 2.0 : 1.0);

/*==========================================================================*/
/* II) Normalize the states                                                 */

/*--------------------------------------------------------------------------*/
/* i) Find overlap of the each state with iself                            */

  ioff = 0;
  for(is=1;is<=nstate;is++){
    ovlap_tmp[is]    = 0.0;
    i1 = 1+ioff;i2=ioff+nstate_ncoef_proc;

    for(ig=i1;ig<=i2;ig++){
      ovlap_tmp[is] += (cre[ig]*cre[ig]);
    }/*endfor:ig*/
    ioff            += nstate_ncoef_proc;
  }/*endfor:is*/

  Allreduce(&(ovlap_tmp[1]),&(ovlap[1]),nstate,MPI_DOUBLE,MPI_SUM,0,comm);

/*-------------------------------------------------------------------------*/
/* ii) Normalize each state to 1                                           */

  ioff = 0;
  for(is=1;is<=nstate;is++){
    scale = sqrt(1.0/ovlap[is]);
    i1 = 1+ioff;i2=ioff+nstate_ncoef_proc;
    for(ig=i1;ig<=i2;ig++){
      cre[ig] *= scale;
    }/*endfor:ig*/
    ioff += nstate_ncoef_proc;
  }/*endfor:is*/

/*========================================================================*/
/* III) Gram schmidt orthogonalize                                         */

  ioff = 0;
  for(is=1;is<=nstate-1;is++){

/*-----------------------------------------------------------------------*/
/* i) Find overlap of the is^th state with all higher states             */
/*    Note the joff starting value                                       */

   len = nstate-is;
   joff = ioff + nstate_ncoef_proc;
   for(js=1;js<=len;js++){
     ovlap_tmp[js]     = 0.0;
     for(ig=1;ig<=nstate_ncoef_proc ;ig++){
      ovlap_tmp[js]     += (cre[(ig+ioff)]*cre[(ig+joff)]);
     }/*endfor:ig*/
     joff += nstate_ncoef_proc;
   }/*endfor:js*/

   Allreduce(&(ovlap_tmp[1]),&(ovlap[1]),len,MPI_DOUBLE,MPI_SUM,0,comm);

/*-----------------------------------------------------------*/
/* ii) Orthogonalize all higher states to the is^th state    */
/*     keeping each state properly normalized                */

   joff = ioff + nstate_ncoef_proc;
   for(js=1;js<=len;js++){
     cmix  = -ovlap[js];
     anorm = 1.0 + 2.0*cmix*ovlap[js] + cmix*cmix;
     scale = sqrt(1.0/anorm);
     for(ig=1;ig<=nstate_ncoef_proc;ig++){
       cre[(ig+joff)] = (cre[(ig+ioff)]*cmix+cre[(ig+joff)])*scale;
     }/*endfor:ig*/
     joff += nstate_ncoef_proc;
   }/*endfor:js*/

/*-----------------------------------------------------------*/
/* iv) increment ioff, the is guys offset                    */

   ioff += nstate_ncoef_proc;
  }/*endfor:is*/

/*========================================================================*/
/* VI) Adjust the norms to the occupation numbers                         */

  for(is=1;is<=nstate;is++){
    scale = sqrt(occ[is]/scale_fact);
    i1 = ioff_vec[is]+1; 
    i2=ioff_vec[is]+nstate_ncoef_proc;
    for(ig=i1;ig<=i2;ig++){
      cre[ig] *= scale;
    }/*endfor:ig*/
  }/*endfor:is*/

/*--------------------------------------------------------------------------*/
}/*end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Rotate into a basis in which the overlap matrix is the diagonal matrix   */
/* of occupation numbers                                                    */
/*==========================================================================*/

void cp_rotate_coef_ortho_dvr(double *creal,
                              int icoef_form,int *icoef_orth,
                              double *norbmat,double *norbmati,
                              double *ovmat_eigv,
                              double *c1_temp,
                              double *occ,int *ioff,
                              double *max_off_diag,double *max_diag,
                              double scale_fact,
                              CPSCR_OVMAT *cpscr_ovmat,
                              CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*         Local Variables           */

  int icoef_orth_now;
  int myid_state = cp_comm_state_pkg->myid;
  int np_states  = cp_comm_state_pkg->num_proc;
  double scale_fact_loc;

  int nstate     = cp_comm_state_pkg->nstate;
  int ncoef      = cp_comm_state_pkg->ncoef;
  int is,js,ig;
  double tot;
/*========================================================================*/
/*  0) Checks                                                             */

  if(*icoef_orth!=0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs must be in nonorthogonal form    \n");
    printf("on state processor %d in cp_rotate_coef_ortho   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1){
    if(icoef_form!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form \n");
      printf("on state processor %d in cp_rotate_coef_ortho  \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
  }/*endif*/

/*========================================================================*/
/*  I) Construct the transformation matrix, its inverse and a             */
/*     transformation matrix to a basis in which S is diagonal but not    */
/*     necessarily equal to the occupation numbers                        */

  icoef_orth_now = 0;

  cp_construct_orth_rotmat_dvr(creal,icoef_form,icoef_orth_now,
                               norbmat,norbmati,ovmat_eigv,
                               occ,ioff,max_off_diag,max_diag,scale_fact,
                               cpscr_ovmat,
                               cp_comm_state_pkg);

/*========================================================================*/
/*  II) Rotate using the transformation matrix and flip the orth flag     */

  cp_rotate_vector_dvr(creal,icoef_form,norbmat,ioff,c1_temp,cp_comm_state_pkg);

  *icoef_orth = 1;

/*========================================================================*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*         Construct St^{1/2} and St^{-1/2} such that                       */
/*     St^{-1/2}*S*St^{-1/2} = diag matrix of occ num                       */
/*==========================================================================*/

void cp_construct_orth_rotmat_dvr(double *creal, int icoef_form,int icoef_orth,
                                  double *norbmat,double *norbmati,
                                  double *ovmat_eigv, double *occ,int *ioff,
                                  double *max_off_diag,double *max_diag,
                                  double scale_fact, CPSCR_OVMAT *cpscr_ovmat,
                                  CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int ierr,ind,is,js,nstate2;
  int job=1;
  double alpha,beta,occ_rt;
  int itransp = 0;
  int inormal = 1;

/*             Local pointer declarations                                */
  int nstate       = cp_comm_state_pkg->nstate;
  int np_states    = cp_comm_state_pkg->num_proc;
  int myid_state   = cp_comm_state_pkg->myid;
  double *ovmat    = cpscr_ovmat->ovlap1;
  double *ov_scr   = cpscr_ovmat->ovlap2;
  double *oveigs   = cpscr_ovmat->state_vec1;
  double *rs_scr1  = cpscr_ovmat->state_vec2;
  double *rs_scr2  = cpscr_ovmat->state_vec3;
  double *occ_tmp  = cpscr_ovmat->state_vec4;
  double tot;
  int ncoef        = cp_comm_state_pkg->ncoef;
/*========================================================================*/
/* I) Check the form of the coefs                                         */

  if(icoef_orth!=0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in nonorthogonal form \n");
    printf("on state processor %d in cp_construct_orth_rotmat  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1 && icoef_form != 1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed form \n");
    printf("on state processor %d in cp_construct_orth_rotmat \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/*========================================================================*/
/* II) Calculate overlap matrix, S, of nonorthogonal orbitals             */

  alpha= scale_fact;
  beta = 1.0;

  cp_ovlap_mat_same_dvr(creal,icoef_form,ovmat,ov_scr,alpha,beta,ioff,cp_comm_state_pkg);

/*========================================================================*/
/* III) Calculate maximum off-diagonal element of S                       */

  for(is=1;is <= nstate; is++){
    for(js=is+1;js <=nstate; js++){
      ind = (is-1)*nstate + js;
      *max_off_diag = MAX(ovmat[ind],*max_off_diag);
    }/*endfor*/
  }/*endfor*/

/*========================================================================*/
/* IV) Calculate maximum diagonal element of S                            */

  for(is=1;is <= nstate; is++){
    ind = (is-1)*nstate + is;
    *max_diag = MAX(ovmat[ind],*max_diag);
  }/*endfor*/

/*========================================================================*/
/* V) Diagonalize the overlap matrix, S                                   */

  nstate2 = nstate*nstate;
  RS(&nstate,&nstate,&(ovmat[1]),&(oveigs[1]),&job,&(ovmat_eigv[1]),
     &(rs_scr1[1]),&(rs_scr2[1]),&ierr);

/*========================================================================*/
/* Shuffle me                                                             */

  rotate_occ_orth_shuffle(occ,occ_tmp,nstate);

/*========================================================================*/
/* VI) Construct the norb transformation matrix (St^{-1/2} and St^{1/2})  */
     /* i) Get inverse square root of eigenvalues                         */

  for(is=1;is <= nstate; is++){
    oveigs[is]  = 1.0/sqrt(oveigs[is]);
  }/*endfor*/

/*------------------------------------------------------------------------*/
/* ii) Construct the diagonal representation of St^{1/2} and St^{-1/2}*/

  for(is=1;is <= nstate2; is++){
    norbmat[is] = 0.0;
    norbmati[is] = 0.0;
  }/*endfor*/
  for(is=1;is <= nstate; is++){
    ind = (is-1)*nstate + is;
    occ_rt = sqrt(occ_tmp[is]);
    norbmat[ind]  = oveigs[is]*occ_rt;
    norbmati[ind] = 1.0/norbmat[ind];
  }/*endfor*/

/*------------------------------------------------------------------------*/
    /* iii) Rotate diag reps to true reps   */

  alpha = 1.0;  beta = 0.0;
  GEN_MATMUL(&(norbmat[1]),&nstate,&inormal,&(ovmat_eigv[1]),&nstate,&itransp,
             &(ov_scr[1]),&nstate,&nstate,&nstate,&nstate,
             &alpha,&beta);
  GEN_MATMUL(&(ovmat_eigv[1]),&nstate,&inormal,&(ov_scr[1]),&nstate,&inormal,
             &(norbmat[1]),&nstate,&nstate,&nstate,&nstate,
             &alpha,&beta);

  alpha = 1.0;  beta = 0.0;
  GEN_MATMUL(&(norbmati[1]),&nstate,&inormal,&(ovmat_eigv[1]),&nstate,&itransp,
             &(ov_scr[1]),&nstate,&nstate,&nstate,&nstate,
             &alpha,&beta);
  GEN_MATMUL(&(ovmat_eigv[1]),&nstate,&inormal,&(ov_scr[1]),&nstate,&inormal,
             &(norbmati[1]),&nstate,&nstate,&nstate,&nstate,
             &alpha,&beta);

/*========================================================================*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Rotate a CP vector (nstate x ncoef) into a different basis               */
/*==========================================================================*/

void cp_rotate_vector_dvr(double *creal,int icoef_form,
                          double *trans_mat,int *ioff,double *c1_temp,
                          CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*           Local variables                      */

  int ntot,i;
  int icoef_form_tmp;
  double alpha,beta;

  int nstate       = cp_comm_state_pkg->nstate;
  int ncoef        = cp_comm_state_pkg->nstate_ncoef_proc;
  int np_states    = cp_comm_state_pkg->num_proc;
  int myid_state   = cp_comm_state_pkg->myid;

/*========================================================================*/
/* I) Checks */
  if(np_states>1 && icoef_form !=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed form\n");
    printf("on state processor %d in cp_rotate_vector \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/*========================================================================*/
/* II) Rotation */

  alpha = 1.0;  beta = 0.0;
  ntot = ncoef*nstate;
  icoef_form_tmp = icoef_form;

  for(i=1;i<=ntot;i++){c1_temp[i]=creal[i];}

  cp_rotation_prim_dvr(c1_temp,icoef_form_tmp,creal,icoef_form,trans_mat,
                       alpha,beta,ioff,cp_comm_state_pkg);

/*========================================================================*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/* Rotation primitive: C2 = alpha*MxC1 + beta*C2                            */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_rotation_prim_dvr(double *c1,int i1coef_form,
                          double *c2,int i2coef_form,
                          double *matrix,double alpha_in,double beta_in,
                          int *ioff_st,CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int ioff,joff,is,js,iii;
  int ind;

  int np_states             = cp_comm_state_pkg->num_proc;
  int myid_state            = cp_comm_state_pkg->myid;
  int nstate                = cp_comm_state_pkg->nstate;
  int nstate_max            = cp_comm_state_pkg->nstate_max;
  int ncoef                 = cp_comm_state_pkg->ncoef;
  int itransp = 0;
  int inormal = 1;
  double beta = beta_in;
  double alpha=alpha_in;

/*========================================================================*/
/* I)  Serial */

  if(np_states == 1){

    GEN_MATMUL(&(c1[1]),&ncoef,&inormal,&(matrix[1]),&nstate,&inormal,
               &(c2[1]),&ncoef,&ncoef,&nstate,&nstate,
               &alpha,&beta);

  }/*endif*/

/*========================================================================*/
/* II) Parallel */

  if(np_states > 1){
 
    if((i1coef_form+i2coef_form) !=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_rotation_prim \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/

    cp_rotate_prim_par_dvr(c1,i1coef_form,c2,i2coef_form,matrix,
                           alpha_in,beta_in,cp_comm_state_pkg);

  }/* endif */

/*========================================================================*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/* Parallel Rotation primitive: C2 = alpha*MxC1 + beta*C2                   */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_rotate_prim_par_dvr(double *c1,int i1coef_form,
                            double *c2,int i2coef_form,
                            double *matrix,double alpha,double beta,
                            CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  double wgt_use,sum,alpha_loc,beta_loc;
  int ioff,joff,is,js,ig,ioff_c,i;
  int ind,jnd,iproc,iii,nread,isum;
  int itransp = 0;
  int inormal = 1;

  int ncoef_tot;
  int nstate                = cp_comm_state_pkg->nstate;
  int nstate_max            = cp_comm_state_pkg->nstate_max;
  int ncoef                 = cp_comm_state_pkg->ncoef;
  int nstate_proc           = cp_comm_state_pkg->nstate_proc;
  int nstate_proc_max       = cp_comm_state_pkg->nstate_proc_max;
  int nstate_ncoef_proc_max = cp_comm_state_pkg->nstate_ncoef_proc_max;
  int nstate_ncoef_proc     = cp_comm_state_pkg->nstate_ncoef_proc;
  int num_proc              = cp_comm_state_pkg->num_proc;
  int myid                  = cp_comm_state_pkg->myid;
  MPI_Comm comm  = cp_comm_state_pkg->comm;
  MPI_Comm world = cp_comm_state_pkg->world;

/*========================================================================*/
/* I) Check the forms */

  if((i1coef_form+i2coef_form) !=2){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed form\n");
    printf("on state processor %d in cp_rotate_prim_par \n",myid);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/*========================================================================*/
/* 0) Initialize matrix                                                   */

  ncoef_tot = nstate_ncoef_proc*nstate;
  for(is=1;is <= ncoef_tot;is++){c2[is] *= beta;}

/*========================================================================*/
/* II) Matrix multiply                                                    */

#ifdef HAND_ROT
  ioff = 0;
  for(is=1;is<=nstate;is++){
    joff = 0;
    for(js=1;js<=nstate;js++){
      ind  = (is-1)*nstate + js;
      for(ig=1;ig<=nstate_ncoef_proc ;ig++){
        c2[(ioff+ig)] += alpha*matrix[ind]*c1[(joff+ig)];
      }
      joff += nstate_ncoef_proc;
    }
    ioff += nstate_ncoef_proc;
  }
#else
  alpha_loc = alpha;
  beta_loc  = beta;
  GEN_MATMUL(&(c1[1]),&nstate_ncoef_proc,&inormal,&(matrix[1]),
             &nstate,&inormal,&(c2[1]),&nstate_ncoef_proc,
             &nstate_ncoef_proc,&nstate,&nstate,&alpha_loc,&beta_loc);

#endif

/*========================================================================*/
}/* end routine */
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rotate_occ_shuffle_dvr(double *creal,int icoef_form,
                            double *ovmat,double *ovmat_scr,
                            double *occ,double *rocc_sum,
                            int *ioff,double scale_fact,
                            CP_COMM_STATE_PKG *cp_comm_state_pkg)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int i,idiag,iii,j,iocc,ind1;
  int occ_temp;
  double alpha,beta;
  int np_states = cp_comm_state_pkg->num_proc;
  int myid      = cp_comm_state_pkg->myid;
  int nstate    = cp_comm_state_pkg->nstate;

/*========================================================================*/
/* I) Check the forms */

  if(np_states>1){
    if(icoef_form  !=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in rotate_occ_shuffle \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
  }/*endif*/

/*========================================================================*/
/* II) Calculate overlap matrix of nonorthogonal orbitals                 */

  alpha = scale_fact;
  beta = 1.0;
  cp_ovlap_mat_same_dvr(creal,icoef_form,ovmat,ovmat_scr,alpha,beta,ioff,
                        cp_comm_state_pkg);

/*========================================================================*/
/* III) Rejigger the occupation numbers                                   */

  for(idiag = 1; idiag <= nstate; idiag++){
    ind1 = (idiag -1 )*nstate + idiag;
    occ_temp    = NINT(ovmat[ind1]);
    occ[idiag]  = occ_temp;
  }/*endfor*/

/*  recalculate rocc_sum_up */
  iocc=0;
  for(i=1;i<= nstate ;i++){
    for(j=1;j<= nstate ;j++){
      iocc++;
      rocc_sum[iocc] = 1.0/(occ[i]+occ[j]);
    }/*endfor i*/
  }/* endfor j*/

/*========================================================================*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/* Normalize a CP vector */
/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_normalize_dvr(double *creal,int icoef_form,
                      double *occ,int *ioff,
                      double *norm_save,double *norm_tmp,
                      CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int iii,is,icf,i1,i2;
  double norm,scale,rsnorm;
  int ncoef      = cp_comm_state_pkg->nstate_ncoef_proc;
  int nstate     = cp_comm_state_pkg->nstate;
  int myid_state = cp_comm_state_pkg->myid;
  int np_states  = cp_comm_state_pkg->num_proc;
  MPI_Comm comm   = cp_comm_state_pkg->comm;

/*========================================================================*/
/* 0) Checks                                                              */

  if(icoef_form !=1&&np_states>1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed form\n");
    printf("on state processor %d in cp_normalize \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/*========================================================================*/
/* I) Calculate the norms                                                 */


  for(is=1;is<=nstate;is++){

    norm = 0.0;
    i1 = 1 + ioff[is];
    i2 = ncoef + ioff[is];
    for(icf=i1;icf<=i2; icf++){
      norm += creal[icf]*creal[icf];
    }
    norm_save[is] = norm;
  }/* endfor is*/

  if(np_states>1){
    for(is=1;is <= nstate;is++){norm_tmp[is] = norm_save[is];}
    Allreduce(&(norm_tmp[1]),&(norm_save[1]),nstate,MPI_DOUBLE,MPI_SUM,0,comm);
  }/*endif*/

/*========================================================================*/
/* I) Normalize each state                                                */

  for(is=1;is<=nstate;is++){
    rsnorm = sqrt(occ[is]/norm_save[is]);
    i1 = 1 + ioff[is];
    i2 = ncoef + ioff[is];
    for(icf=i1;icf<=i2; icf++){
      creal[icf] *= rsnorm;
    }/*endfor*/
  }/* endfor is*/

/*========================================================================*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* construct and add the KS contribution to the forces                      */
/*==========================================================================*/

void cp_add_ksmat_force_dvr(double *creal,int icoef_form,int icoef_orth,
                            double *fcreal,int ifcoef_form,int ifcoef_orth,
                            double *ksmat,double *ks_scr,
                            int *ioff,int cp_lsda,int cp_min,
                            double *occ,double scale_factor,
                            CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int nstate2;
  int iii,i,is,js,icf,jcf;
  int ind,ind1,ind2;
  int one=1;
  int koff;
  double alpha,beta;
  double factor = scale_factor/2.0;

  double tot1,tot2;

/*             Local pointer declarations                                */
   int np_states    = cp_comm_state_pkg->num_proc;
   int nstate       = cp_comm_state_pkg->nstate;
   int myid_state   = cp_comm_state_pkg->myid;
   int ncoef        = cp_comm_state_pkg->nstate_ncoef_proc;


/*========================================================================*/
/* I) Check the form of coefs and forces                                  */

  if(icoef_orth!=1 || ifcoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs and coef forces must be in orthogonal form    \n");
    printf("on state processor %d in cp_add_ksmat_force   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1){
   if(icoef_form!=1 || ifcoef_form!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs and coef forces must be in transposed form \n");
    printf("on state processor %d in cp_add_ksmat_force   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
  }/*endif*/

/*========================================================================*/
/* II) Compute KS matrix from the coef/orbital overlap  KS=F*C            */

  nstate2 = nstate*nstate;
  for(is=1;is<=nstate2; is++){ksmat[is] = 0.0;}

  alpha = 1.0;
  beta  = 1.0;
  cp_ovlap_mat_diff_dvr(creal,icoef_form,fcreal,ifcoef_form,
                        ksmat,ks_scr,alpha,beta,ioff,cp_comm_state_pkg);

/*========================================================================*/
/* III) If minimization: scale the matrix                                 */


  if(cp_min==1){
    alpha = 2.0;
    koff = 0;
    for(is=1;is<=nstate; is++){
      for(js=1;js<=nstate; js++){
        ksmat[(koff+js)] *= (alpha/occ[is]);
      }/* endfor */
      koff += nstate;
    }/* endfor */
  }/* endif */

/*========================================================================*/
/* IV) Add Kohn-Sham matrix contribution to force F -= KS*c               */

  alpha = -1.0;  beta = 1.0;
  if(cp_lsda == 1 && cp_min == 0){alpha = -2.0;}
  if(cp_lsda == 0 && nstate == 1 && occ[1] == 1.0 && cp_min == 0)  {alpha = -2.0;}

   alpha *= factor;
  cp_rotation_prim_dvr(creal,icoef_form,fcreal,ifcoef_form,ksmat,
                       alpha,beta,ioff,cp_comm_state_pkg);

/*========================================================================*/
  }/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Rotate back from a basis in which the overlap matrix is diagonal matrix  */
/* of occupation numbers to the original basis                              */
/*==========================================================================*/

void cp_rotate_gen_nonortho_dvr(double *creal,int icoef_form,int *icoef_orth,
                                double *norbmat,int *ioff,double *c1_temp,
                                CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*         Local Variables           */

  int myid_state = cp_comm_state_pkg->myid;
  int np_states  = cp_comm_state_pkg->num_proc;

/*========================================================================*/
/*  0) Checks                                                             */

  if((*icoef_orth)!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs must be in orthogonal form \n");
    printf("on state processor %d in cp_rotate_gen_nonortho   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1){
   if(icoef_form!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed form \n");
    printf("on state processor %d in cp_rotate_gen_nonortho  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
  }/*endif*/

/*========================================================================*/
/*  I) Rotate using the transformation matrix and flip the orth flag     */

  cp_rotate_vector_dvr(creal,icoef_form,norbmat,ioff,c1_temp,cp_comm_state_pkg);

  *icoef_orth = 0;

/*========================================================================*/
  }/*end routine*/
/*========================================================================*/


