/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: cp_con_utils.c                               */
/*                                                                          */
/* Contains utilities for wave function constraints and norb                */
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

#define DEBUG_OFF

/*==========================================================================*/
/* Fwd Transpose controller: nstate_proc x ncoef ->  nstate x ncoef_proc    */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_transpose_fwd_dvr(double *cre,int *icoef_form,double *c1_temp,
                          CP_COMM_STATE_PKG *cp_comm_state_pkg,
                          int nkf1,int nkf2,int nkf3)

/*======================================================================*/
/*                Begin Routine */
     {/*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */
#include "../typ_defs/typ_mask.h"

 int iii;
 int nstate                = cp_comm_state_pkg->nstate;
 int nstate_max            = cp_comm_state_pkg->nstate_max;
 int ncoef                 = cp_comm_state_pkg->ncoef;
 int nstate_proc           = cp_comm_state_pkg->nstate_proc;
 int nstate_proc_max       = cp_comm_state_pkg->nstate_proc_max;
 int nstate_proc_min       = cp_comm_state_pkg->nstate_proc_min;
 int nstate_ncoef_proc_max = cp_comm_state_pkg->nstate_ncoef_proc_max;
 /* int nstate_ncoef_proc_min = cp_comm_state_pkg->nstate_ncoef_proc_min;*/
 int nstate_ncoef_proc     = cp_comm_state_pkg->nstate_ncoef_proc;
 int num_proc              = cp_comm_state_pkg->num_proc;
 int myid                  = cp_comm_state_pkg->myid;
 MPI_Comm comm  = cp_comm_state_pkg->comm;
 int icoef_form_tmp;

/*========================================================================*/
/* I) Check the form                                                      */

  if(num_proc==1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("There is no need to transpose the vectors on 1 processor\n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/
  if((*icoef_form==1)){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in normal form\n");
    printf("on state processor %d in cp_transpose_fwd \n",myid);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/
/*========================================================================*/
/* II) Transpose the vector and flip the flag */

  icoef_form_tmp = 0;

  cp_transpose_fwd_prim_dvr(cre,&icoef_form_tmp,c1_temp,nstate,nstate_max,ncoef,
                            nstate_proc,nstate_proc_max,nstate_proc_min,
                            nstate_ncoef_proc_max,nstate_ncoef_proc,
                            nkf1,nkf2,nkf3,num_proc,myid,comm);
  (*icoef_form) = 1;

/*========================================================================*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/* Bck Transpose controller: nstate_proc x ncoef ->  nstate x ncoef_proc    */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_transpose_bck_dvr(double *cre,int *icoef_form,double *c1_temp,
                          CP_COMM_STATE_PKG *cp_comm_state_pkg,
                          int nkf1,int nkf2,int nkf3)

/*======================================================================*/
/*                Begin Routine */
     {/*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */
#include "../typ_defs/typ_mask.h"

 int nstate                = cp_comm_state_pkg->nstate;
 int nstate_max            = cp_comm_state_pkg->nstate_max;
 int ncoef                 = cp_comm_state_pkg->ncoef;
 int nstate_proc           = cp_comm_state_pkg->nstate_proc;
 int nstate_proc_max       = cp_comm_state_pkg->nstate_proc_max;
 int nstate_proc_min       = cp_comm_state_pkg->nstate_proc_min;
 int nstate_ncoef_proc_max = cp_comm_state_pkg->nstate_ncoef_proc_max;
 /*int nstate_ncoef_proc_min = cp_comm_state_pkg->nstate_ncoef_proc_min;*/
 int nstate_ncoef_proc     = cp_comm_state_pkg->nstate_ncoef_proc;
 int num_proc              = cp_comm_state_pkg->num_proc;
 int myid                  = cp_comm_state_pkg->myid;


 MPI_Comm comm  = cp_comm_state_pkg->comm;
 int icoef_form_tmp;

/*========================================================================*/
/* I) Check the form                                                      */

  if(num_proc==1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("There is no need to transpose the vectors on 1 processor\n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/
  if((*icoef_form==0)){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in normal form\n");
    printf("on state processor %d in cp_transpose_bck \n",myid);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/*========================================================================*/
/* II) Transpose the vector and flip the flag                             */

  icoef_form_tmp = 1;

  cp_transpose_bck_prim_dvr(cre,&icoef_form_tmp,c1_temp,nstate,nstate_max,ncoef,
                            nstate_proc,nstate_proc_max,nstate_proc_min,
                            nstate_ncoef_proc_max, nstate_ncoef_proc,
                            nkf1,nkf2,nkf3,num_proc,myid,comm);

  (*icoef_form) = 0;

/*========================================================================*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/* Perform the forward transpose: nstate_proc x ncoef-> nstate x ncoef_proc */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_transpose_fwd_prim_dvr(double *c,int *icoef_form,double *c_temp,
                               int nstate,int nstate_max,int ncoef,
                               int nstate_proc, int nstate_proc_max,
                               int nstate_proc_min,
                               int nstate_ncoef_proc_max,
                               int nstate_ncoef_proc,
                               int nkf1,int nkf2,int nkf3,
                               int num_proc,int myid,MPI_Comm comm)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */

#include "../typ_defs/typ_mask.h"
  int iremc,irems,irem,ig,is,ioff;
  int i,j,icount;
  int ioff_c,ioff_s;
  int iproc,iii;
  int sendcounts,recvcounts;
  int *ioffv,*inum;
  int irm,idv;

/*========================================================================*/
/*              Incoming variable declarations                            */
/*                                                                        */
/* nstate            = Total number of states                             */
/* nstate_proc       = number of states on this processor                 */
/* nstate_proc_max   = maximum number of states on any processor          */
/* ncoef             = Total number of coefficients in a state            */
/* nstate_ncoef_proc = Number of coefficients in a state on this processor*/
/*                     in the transposed data                             */
/* nstate_ncoef_proc_max = Maximum number of coefficients in a state on   */
/*                          any processesor in the transposed data        */
/* nstate_max        = nstate_proc_max*num_proc                           */
/* c                 = nstate_proc x ncoef array of coefficients          */
/* c_temp            = transposed data stored as nstate x ncoef_proc_max  */
/* ct_temp           = scratch space to help make transposed data         */
/* nscr_size         = size of scratch nstate_ncoef_proc_max*nstate_max   */

/*========================================================================*/
/* 0) Check the forms */

  if((*icoef_form) !=0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in normal form\n");
    printf("on state processor %d in cp_transpose_fwd_prim \n",myid);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  *icoef_form = 1;

/*========================================================================*/
/* I) Rearrange the coefficient data in c                                 */

  inum  = (int *) cmalloc(num_proc*sizeof(int))-1;
  ioffv = (int *) cmalloc(num_proc*sizeof(int))-1;

  iremc = ((nkf2*nkf3) % num_proc);
  irems = (nstate%num_proc);

  idv = (nkf2*nkf3)/num_proc;
  irm = (nkf2*nkf3)%num_proc;

  for(i=0; i < num_proc; i++){
    inum[i+1] = ( i < irm ? (idv+1)*nkf1 : idv*nkf1 ); 
  }

  if((iremc != 0) || (irems != 0) ){

    for(i=1;i<=nstate_ncoef_proc_max*nstate_max;i++){c_temp[i]=0.0;}

    for(is=1;is<=nstate_proc;is++){
      ioff_c = (is-1)*ncoef;
      icount = 0;
      for(iproc=1;iproc<=num_proc;iproc++){
        for(ig=1;ig<=inum[iproc];ig++){

          ioff = (iproc-1)*nstate_ncoef_proc_max*nstate_proc_max
                  + (is-1)*nstate_ncoef_proc_max;
          icount++;
          c_temp[(ioff+ig)] = c[ioff_c+icount];
        }/*endfor*/
      }/*endfor*/
    }/*endfor*/

    for(i=1;i<=nstate_ncoef_proc_max*nstate_max;i++){c[i]=c_temp[i];}

  }/*endif*/

/*========================================================================*/
/* II) Transpose the coef data                                            */

  sendcounts = nstate_ncoef_proc_max*nstate_proc_max;
  recvcounts = nstate_ncoef_proc_max*nstate_proc_max;

  Alltoall(&c[1],sendcounts,MPI_DOUBLE,&c_temp[1],recvcounts,MPI_DOUBLE,comm);


  irem  = MAX(iremc,irems);

  for(i=1; i <= irems; i++){ ioffv[i] = nstate_proc_max; }
  for(i=irems+1; i <= num_proc; i++){ ioffv[i] = nstate_proc_min; }

/*========================================================================*/
/* V) Internal rearrangement of coeff data                                */

  irems = (nstate%num_proc);
  irem  = MAX(iremc,irems);

  if((irems != 0) || (iremc != 0)){

    for(i=1; i <= irems; i++){ ioffv[i] = nstate_proc_max; }
    for(i=irems+1; i <= num_proc; i++){ ioffv[i] = nstate_proc_min; }
    for(i=1;i<=nstate_ncoef_proc_max*nstate_max;i++){c[i]=0.0;}

    icount = 0;
    for(iproc=1; iproc <= num_proc; iproc++){
      for(is=1; is <= ioffv[iproc]; is++){
        icount++;
        ioff_s = (icount-1)*nstate_ncoef_proc;
        ioff_c = (iproc-1)*nstate_ncoef_proc_max*nstate_proc_max
                 + (is-1)*nstate_ncoef_proc_max;
        for(ig=1; ig <= nstate_ncoef_proc; ig++){
          c[(ioff_s+ig)] = c_temp[(ioff_c+ig)];
        }
      }
    }

  }/*endif: remainder */

  cfree(&inum[1]);
  cfree(&ioffv[1]);

/*========================================================================*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/* Perform the back transpose: nstate x ncoef_proc -> nstate_proc x ncoef   */
/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

void cp_transpose_bck_prim_dvr(double *c,int *icoef_form,double *c_temp,
                               int nstate,int nstate_max,int ncoef,
                               int nstate_proc,int nstate_proc_max,
                               int nstate_proc_min,
                               int nstate_ncoef_proc_max,
                               int nstate_ncoef_proc,
                               int nkf1,int nkf2,int nkf3,
                               int num_proc,int myid,MPI_Comm comm)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */

#include "../typ_defs/typ_mask.h"

  int nfull,i,j;
  int idv,irm;
  int irems,iremc,irem,ig,is,ioff;
  int ioff_c,iproc,iii;
  int sendcounts,recvcounts;
  int icount;
  int *ioffv,*inum;

/*========================================================================*/
/*              Incoming variable declarations                            */
/*                                                                        */
/* nstate            = Total number of states                             */
/* nstate_proc       = number of states on this processor                 */
/* nstate_proc_max   = maximum number of states on any processor          */
/* ncoef             = Total number of coefficients in a state            */
/* nstate_ncoef_proc = Number of coefficients in a state on this processor*/
/*                     in the transposed data                             */
/* nstate_ncoef_proc_max = Maximum number of coefficients in a state on   */
/*                          any processesor in the transposed data        */
/* nstate_max        = nstate_proc_max*num_proc                           */
/* c                 = nstate_proc x ncoef array of coefficients          */
/* c_temp            = transposed data stored as nstate x ncoef_proc_max  */
/* ct_temp           = scratch space to help make transposed data         */
/* nscr_size         = size of scratch nstate_ncoef_proc_max*nstate_max   */
/*========================================================================*/
/* 0) Check the forms */

  if((*icoef_form) !=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed form\n");
    printf("on state processor %d in cp_transpose_bck_prim \n",myid);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  *icoef_form = 0;

/*========================================================================*/
/*========================================================================*/
/* I) Internal rearrangement of coeff data                                */

  irems = (nstate % num_proc);
  iremc = ((nkf2*nkf3)% num_proc);
  irem  = MAX(irems,iremc);

  ioffv = (int *) cmalloc(num_proc*sizeof(int))-1;
  inum  = (int *) cmalloc(num_proc*sizeof(int))-1;

  if((irems != 0)||(iremc != 0)){

    ioffv[1] = 0;
    for(i=1; i<= irems; i++){
      ioffv[i+1] = ioffv[i] + nstate_proc_max*nstate_ncoef_proc;
      inum[i]    = nstate_proc_max*nstate_ncoef_proc;
    }
    for(i=irems+1; i < num_proc; i++){
      ioffv[i+1] = ioffv[i] + nstate_proc_min*nstate_ncoef_proc;
      inum[i]    = nstate_proc_min*nstate_ncoef_proc;
    }
    inum[num_proc] = nstate_proc_min*nstate_ncoef_proc;

/* 1) copy data into temp array */

    for(i=1;i <= (num_proc*nstate_proc_max*nstate_ncoef_proc_max);i++){ 
      c_temp[i] = 0.00;
    }

    for(i=1; i<= num_proc; i++){
      ioff = (i-1)*nstate_proc_max*nstate_ncoef_proc_max;
      for(j=1; j <= inum[i]; j++){
        c_temp[(ioff+j)] = c[(ioffv[i]+j)];
      }
    }

/* 2) copy back into c     */

    nfull = nstate_ncoef_proc_max*nstate_max;
    for(ig=1;ig<=nfull;ig++){c[ig] = c_temp[ig];}

  }/*endif: remainder */

/*======================================================================*/
/* II) Send the transformed position data                               */

  sendcounts = nstate_ncoef_proc_max*nstate_proc_max;
  recvcounts = nstate_ncoef_proc_max*nstate_proc_max;

  Alltoall(&c[1],sendcounts,MPI_DOUBLE,&c_temp[1],recvcounts,
                MPI_DOUBLE,comm);

/*=======================================================================*/
/* III) Extract the transformed position data                            */

  idv   = (nkf2*nkf3)/num_proc;
  irm   = (nkf2*nkf3)%num_proc;

  for(i = 0; i < num_proc; i++){
    inum[i+1] = ( i < irm ? (idv+1)*nkf1 : idv*nkf1);
  }

  for(i=1;i<=nstate_ncoef_proc_max*nstate_max;i++){c[i]=0.0;}

  for(is=1; is <= nstate_proc; is++){
    icount = 0;
    ioff_c = (is-1)*ncoef;
    for(iproc=1;iproc<=num_proc;iproc++){
      ioff = (is-1)*inum[iproc]
           + (iproc-1)*nstate_proc_max*nstate_ncoef_proc_max;
      for(i=1;i<= inum[iproc];i++){
        icount++;
        c[ioff_c+icount] = c_temp[ioff+i];
      }/*endfor*/
    }/*endfor*/
  }

  cfree(&ioffv[1]);
  cfree(&inum[1]);

/*========================================================================*/
}/*end routine*/
/*========================================================================*/

