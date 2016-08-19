#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "mpi.h"


#define DEBUG
#define Alltoall MPI_Alltoall
#define Barrier MPI_Barrier

void cp_state_gvec_trans_fwd(double *,double *,double *,
                             int ,int ,int ,int ,int ,
                             int ,int ,int ,int ,MPI_Comm );


void cp_state_gvec_trans_bck(double *,double *,double *,
                             int ,int ,int ,int ,int ,
                             int ,int ,int ,int ,MPI_Comm );



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

int main (int argc, char *argv[])

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */
  
  double frac,whole;
  double *c;
  double *c_temp,*ct_temp;
  int nmall;
  int nstate_proc,nstate_proc_max,nstate_proc_st;
  int nstate_ncoef_proc,nstate_ncoef_proc_max;
  int myid,num_proc;
  int iproc,iii;
  int ioff_c;
  int irem,idiv;
  int ig,is,ioff;
  MPI_Comm world;
/*---------------------*/
  int ncoef       = 7;
  int nstate      = 17;
  int nstate_max;

/*========================================================================*/
/* 0) Initialize MPI                                                      */

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  world = MPI_COMM_WORLD;

/*========================================================================*/
/* 0) Calculate some constants                                            */

/*---------------------------*/
/* state per processor stuff */
  idiv             = nstate/num_proc;
  irem             = (nstate % num_proc);
  nstate_proc_max  = (irem > 0 ? idiv+1 : idiv);
  nstate_proc      = (myid < irem ? idiv+1 : idiv);
  if(myid <= irem) {
    nstate_proc_st = myid*(idiv+1)+1;
  } else {
    nstate_proc_st = irem*(idiv+1) + (myid-irem)*idiv+1;
  }/*endif*/
  nstate_max      = nstate_proc_max*num_proc;

/*--------------------------*/
/* coef per processor stuff */
  idiv                  = ncoef/num_proc;
  irem                  = (ncoef % num_proc);
  nstate_ncoef_proc_max = (irem > 0 ? idiv+1 : idiv);
  nstate_ncoef_proc     = (myid < irem ? idiv+1 : idiv);

/*========================================================================*/
/* I) Malloc the memory                                                   */

  nmall = nstate_ncoef_proc_max*nstate_max;
  c          = (double *)malloc(nmall*sizeof(double))-1;
  c_temp     = (double *)malloc(nmall*sizeof(double))-1;
  ct_temp    = (double *)malloc(nmall*sizeof(double))-1;

/*========================================================================*/
/* II) Fill the pos vector structure with model data                      */

  ioff=0;
  for(is=1;is<=nstate_proc;is++){
    for(ig=1;ig<=ncoef;ig++){
      frac = (double)ig;
      frac /= 10;
      whole = (double)(nstate_proc_st+is-1);
      c[(ig+ioff)] = frac + whole;
    }/*endfor*/
   ioff += ncoef;
  }/*endfor*/

#ifdef DEBUG
  if(myid==0){
    printf("This is x\n");
  }/*endif*/
  for(iproc=1;iproc<=num_proc;iproc++){
    if(myid==iproc-1){
     for(is=1;is<=nstate_proc;is++){
       ioff_c = (is-1)*ncoef;
       for(ig=1;ig<=ncoef;ig++){
        printf("%d %d %d %g\n",iproc,is,ig,c[(ig+ioff_c)]);
       }/*endfor*/
     }/*endfor*/
    }/*endif*/
    scanf("%d",&iii);
    Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif

/*========================================================================*/
/* III) Transpose the data */


  cp_state_gvec_trans_fwd(c,c_temp,ct_temp,nstate,nstate_max,ncoef,
                          nstate_proc,nstate_proc_max,
                          nstate_ncoef_proc_max,nstate_ncoef_proc,
                          num_proc,myid,world);

  cp_state_gvec_trans_bck(c,c_temp,ct_temp,nstate,nstate_max,ncoef,
                          nstate_proc,nstate_proc_max,
                          nstate_ncoef_proc_max,nstate_ncoef_proc,
                          num_proc,myid,world);

/*--------------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_state_gvec_trans_fwd(double *c,double *c_temp,double *ct_temp,
                             int nstate,int nstate_max,int ncoef,
                             int nstate_proc,int nstate_proc_max,
                             int nstate_ncoef_proc_max,int nstate_ncoef_proc,
                             int num_proc,int myid,MPI_Comm comm)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */
  
  double *c_temp_pt,*ct_temp_pt;
  int nfull,nblock;
  int irem,ig,is,ioff;
  int i,ioff_temp;
  int ioff_c;
  int joff;
  int iproc,itemp,iii;
  int proc_rem;
  int nstate_ncoef_proc_now;
  int sendcounts,recvcounts;

/*========================================================================*/
/*              Incoming variable declarations                            */

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
/* I) Extract the coefficient data from c                                 */

  for(i=1;i<=nstate_ncoef_proc_max*nstate_max;i++){c_temp[i]=0.0;}
  proc_rem = ncoef % num_proc;

  for(is=1;is<=nstate_proc;is++){
    ioff = 0;
    ioff_c = (is-1)*ncoef;
    for(iproc=1;iproc<=num_proc;iproc++){
      ioff_temp = (is-1)*nstate_ncoef_proc_max
                + (iproc-1)*(nstate_proc_max*nstate_ncoef_proc_max);
      nstate_ncoef_proc_now = nstate_ncoef_proc_max;
      if((iproc > proc_rem) && (proc_rem >0)) nstate_ncoef_proc_now--;
      for(ig=1;ig<=nstate_ncoef_proc_now;ig++){
        itemp = ig+ioff_temp;
        i     = ig+ioff;
        c_temp[itemp] = c[(i+ioff_c)];
      }/*endfor*/   
      ioff += nstate_ncoef_proc_now;
    }/*endfor*/   
  }/*endfor*/   

#ifdef DEBUG
  if(myid==0){
    printf("This is c_temp\n");
  }/*endif*/
  for(iproc=1;iproc<=num_proc;iproc++){
    if(myid==iproc-1){
     for(ig=1;ig<=nstate_ncoef_proc_max*nstate_max;ig++){
      printf("%d %d %g\n",iproc,ig,c_temp[ig]);
     }/*endfor*/
    }/*endif*/
    scanf("%d",&iii);
    Barrier(comm);
  }/*endfor*/
#endif

/*========================================================================*/
/* II) Transpose the coef data                                            */

  sendcounts = nstate_ncoef_proc_max*nstate_proc_max;
  recvcounts = nstate_ncoef_proc_max*nstate_proc_max;

  ct_temp_pt = ct_temp+1;
  c_temp_pt  = c_temp+1;
  Alltoall(c_temp_pt,sendcounts,MPI_DOUBLE,ct_temp_pt,recvcounts,
                MPI_DOUBLE,comm);

#ifdef DEBUG
  if(myid==0){
    printf("This is ct_temp\n");
  }/*endif*/
  for(is=1;is<=num_proc;is++){
    if(myid==is-1){
      for(ig=1;ig<=nstate_ncoef_proc_max*nstate_max;ig++){
         printf("%d %d %g\n",ig,is,ct_temp[ig]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&ig);
    Barrier(comm);
  }/*endfor*/
#endif
  

/*========================================================================*/
/* V) Internal rearrangement of coeff data                                */

  irem = (nstate % num_proc);

  if(irem != 0){

/*------------------------------------------------------------------------*/
/* A) copy full blocks                                                    */
    nfull = nstate_proc_max*nstate_ncoef_proc_max*irem;
    for(ig=1;ig<=nfull;ig++){c_temp[ig] = ct_temp[ig];}
/*------------------------------------------------------------------------*/
/* B) copy partial blocks                                                 */
    nblock = (nstate_proc_max-1)*nstate_ncoef_proc_max;
    ioff = nfull; joff=nfull;
    for(iproc=irem+1;iproc<=num_proc;iproc++){
     for(ig=1;ig<=nblock;ig++){c_temp[(ig+ioff)] = ct_temp[(ig+joff)];}
     ioff += nblock;joff += (nblock+nstate_ncoef_proc_max);
    }/*endfor*/

  }else{

/*------------------------------------------------------------------------*/
/* A) copy all blocks                                                     */
    nfull = nstate_ncoef_proc_max*nstate;
    for(ig=1;ig<=nfull;ig++){c_temp[ig] = ct_temp[ig];}

  }/*endif: remainder */
    
#ifdef DEBUG
  if(myid==0){
    printf("This is c_temp\n");
  }/*endif*/
  for(is=1;is<=num_proc;is++){
    if(myid==is-1){
      for(ig=1;ig<=nstate_ncoef_proc_max*nstate;ig++){
         printf("%d %d %g\n",ig,is,c_temp[ig]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&ig);
    Barrier(comm);
  }/*endfor*/
#endif
  

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/

void cp_state_gvec_trans_bck(double *c,double *c_temp,double *ct_temp,
                             int nstate,int nstate_max,int ncoef,
                             int nstate_proc,int nstate_proc_max,
                             int nstate_ncoef_proc_max,int nstate_ncoef_proc,
                             int num_proc,int myid,MPI_Comm comm)


/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */
  
  double *c_temp_pt,*ct_temp_pt;
  int nfull,nblock;
  int irem,ig,is,ioff;
  int i,ioff_temp;
  int ioff_c;
  int joff;
  int iproc,itemp,iii;
  int proc_rem;
  int nstate_ncoef_proc_now;
  int sendcounts,recvcounts;

/*========================================================================*/
/*              Incoming variable declarations                            */

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
/* I) Internal rearrangement of coeff data                                */

  irem = (nstate % num_proc);

  if(irem != 0){

/*------------------------------------------------------------------------*/
/* A) copy full blocks                                                    */
    nfull = nstate_proc_max*nstate_ncoef_proc_max*irem;
    for(ig=1;ig<=nfull;ig++){ct_temp[ig] = c_temp[ig];}
/*------------------------------------------------------------------------*/
/* B) copy partial blocks                                                 */
    nblock = (nstate_proc_max-1)*nstate_ncoef_proc_max;
    ioff = nfull; joff=nfull;
    for(iproc=irem+1;iproc<=num_proc;iproc++){
     for(ig=1;ig<=nblock;ig++){ct_temp[(ig+joff)] = c_temp[(ig+ioff)];}
     ioff += nblock;joff += (nblock+nstate_ncoef_proc_max);
    }/*endfor*/

  }else{

/*------------------------------------------------------------------------*/
/* A) copy all blocks                                                     */

    nfull = nstate_ncoef_proc_max*nstate;
    for(ig=1;ig<=nfull;ig++){ct_temp[ig] = c_temp[ig];}

  }/*endif: remainder */
    
#ifdef DEBUG
  if(myid==0){
    printf("This is ct_temp\n");
  }/*endif*/
  for(is=1;is<=num_proc;is++){
    if(myid==is-1){
      for(ig=1;ig<=nstate_ncoef_proc_max*nstate;ig++){
         printf("%d %d %g\n",ig,is,ct_temp[ig]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&ig);
    Barrier(comm);
  }/*endfor*/
#endif
  
/*========================================================================*/
/* II) Send the transformed position data                               */

  sendcounts = nstate_ncoef_proc_max*nstate_proc_max;
  recvcounts = nstate_ncoef_proc_max*nstate_proc_max;

  ct_temp_pt = ct_temp+1;
  c_temp_pt = c_temp+1;
  Alltoall(ct_temp_pt,sendcounts,MPI_DOUBLE,c_temp_pt,recvcounts,
                MPI_DOUBLE,comm);

#ifdef DEBUG
  if(myid==0){
    printf("This is c_temp\n");
  }/*endif*/
  for(is=1;is<=num_proc;is++){
    if(myid==is-1){
      for(ig=1;ig<=nstate_ncoef_proc_max*nstate_max;ig++){
         printf("%d %d %g\n",ig,is,c_temp[ig]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&ig);
    Barrier(comm);
  }/*endfor*/
#endif


/*========================================================================*/
/* III) Extract the transformed position data                               */

  for(i=1;i<=nstate_ncoef_proc_max*nstate_max;i++){c[i]=0.0;}
  proc_rem = ncoef % num_proc;

  for(is=1;is<=nstate_proc;is++){
    ioff = 0;
    ioff_c = (is-1)*ncoef;
    for(iproc=1;iproc<=num_proc;iproc++){
      ioff_temp = (is-1)*nstate_ncoef_proc_max
                + (iproc-1)*(nstate_proc_max*nstate_ncoef_proc_max);
      nstate_ncoef_proc_now = nstate_ncoef_proc_max;
      if((iproc>proc_rem)&&(proc_rem>0)) nstate_ncoef_proc_now--;
      for(ig=1;ig<=nstate_ncoef_proc_now;ig++){
        itemp = ig+ioff_temp;
        i     = ig+ioff;
        c[(i+ioff_c)] = c_temp[itemp];
      }/*endfor*/   
      ioff += nstate_ncoef_proc_now;
    }/*endfor*/   
  }/*endfor*/   

#ifdef DEBUG
  if(myid==0){
    printf("This is c\n");
  }/*endif*/
  for(iproc=1;iproc<=num_proc;iproc++){
    if(myid==iproc-1){
     for(is=1;is<=nstate_proc;is++){
       ioff_c = (is-1)*ncoef;
       for(ig=1;ig<=ncoef;ig++){
        printf("%d %d %d %g\n",iproc,is,ig,c[(ig+ioff_c)]);
       }/*endfor*/
     }/*endfor*/
    }/*endif*/
    scanf("%d",&iii);
    Barrier(comm);
  }/*endfor*/
#endif
/*========================================================================*/
   }/*end routine*/
/*========================================================================*/


















