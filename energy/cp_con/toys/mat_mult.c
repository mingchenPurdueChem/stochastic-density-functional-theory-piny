#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define MULTIPLY_OFF
#define MAXTIME 2147.483648
#define MINTIME (-MAXTIME)

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC  1000000
#endif

#define DEBUG_OFF 

void cputime(double *);


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

int main(int argc, char *argv[]) 

/*========================================================================*/
/* Begin Routine                                                          */
  {/* begin routine */
/*========================================================================*/
/*          Local Variables         */
  int ind,icomm,iproc,num_proc,myid,ntot;
  int nsend_a,nsend_b,nrecv_b,recv_id;
  int nm,nl_proc_a,nl_proc_max,nl_tot,nl_tot2;
  int nl_proc_a_st,nl_proc_a_end;
  int nl_proc_b_st_now,nl_proc_b_end_now;
  int send_id;
  int *nl_proc_b;
  int *nl_proc_b_st;
  int *nl_proc_b_end;
  int joff,jl,knd,jnd;
  int ncomm;
  int idiv,irem,ioff,il,im;
  int done,upper;
  double t0,t1,dt,dt_tmp;
  double frac,whole;
  double *a,*c;
  double *b_old,*b_new,*b;
  double c_tmp,*c_red;
  MPI_Request request;
  MPI_Request request1;
  MPI_Status stat;

/*========================================================================*/
/* 0) Initialize MPI                                                      */

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

/*========================================================================*/
/* I) Malloc memory                                                       */

   nl_tot = 128;
   nm = 50000;
   nl_proc_max = nl_tot/num_proc + 1;

   nl_tot2 = nl_tot*nl_tot;
   ncomm = ((num_proc % 2) == 0 ? num_proc/2 : (num_proc-1)/2);

   a             = (double *) malloc(nm*nl_proc_max*sizeof(double))-1;
   b             = (double *) malloc(nm*nl_proc_max*sizeof(double))-1;
   b_new         = (double *) malloc(nm*nl_proc_max*sizeof(double))-1;
   b_old         = (double *) malloc(nm*nl_proc_max*sizeof(double))-1;
   c             = (double *) malloc(nl_tot2*sizeof(double))-1;
   c_red         = (double *) malloc(nl_tot2*sizeof(double))-1;
   nl_proc_b     = (int *) malloc(ncomm*sizeof(int))-1;
   nl_proc_b_st  = (int *) malloc(ncomm*sizeof(int))-1;
   nl_proc_b_end = (int *) malloc(ncomm*sizeof(int))-1;


/*========================================================================*/
/* II) Find out amount of stuff on each proc and indices in state hierarchy */

   idiv = nl_tot/num_proc;
   irem = nl_tot % num_proc;
   nl_proc_a = (myid < irem ? idiv+1 : idiv);
   if(myid <= irem) {
    nl_proc_a_st = myid*(idiv+1)+1;
   } else {
    nl_proc_a_st = irem*(idiv+1) + (myid-irem)*idiv+1;
   }/*endif*/
   nl_proc_a_end = nl_proc_a_st + nl_proc_a -1;

/*========================================================================*/
/* III) Make the send/recv list (maps)                                    */

   send_id = myid+1;
   if(send_id > num_proc-1) send_id -= num_proc;

   idiv = nl_tot/num_proc;
   irem = nl_tot % num_proc;
   for(icomm=1;icomm<=ncomm;icomm++){
     recv_id   = myid-icomm; 
     if(recv_id < 0) recv_id += num_proc;
     nl_proc_b[icomm] = (recv_id < irem ? idiv+1 : idiv);
     if(recv_id <= irem) {
      nl_proc_b_st[icomm] = recv_id*(idiv+1)+1;
     } else {
      nl_proc_b_st[icomm] = irem*(idiv+1) + (recv_id-irem)*idiv+1;
     }/*endif*/
     nl_proc_b_end[icomm] = nl_proc_b_st[icomm] + nl_proc_b[icomm] -1;
   }/*endfor*/

#ifdef DEBUG
  if(myid==0){printf("Send/Recv information by shift\n");}
  for(iproc=0;iproc<num_proc;iproc++){
   if(myid==iproc) {
     for(im=1;im<=ncomm;im++){
       printf("send_id = %d nl_b[%d] = %d nl_b_st[%d] = %d myid = %d\n",
          send_id,im,nl_proc_b[im],im,nl_proc_b_st[im],myid);
     }/*endfor*/
   }/*endif*/
   scanf("%d",&il);
   MPI_Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif

/*========================================================================*/
/* II.5) Create model data                                                */

  ioff=0;
  for(il=1;il<=nl_proc_a;il++){
    for(im=1;im<=nm;im++){
      frac = (double)im;
      frac /= 10;
      whole = (double)(nl_proc_a_st+il-1);
      a[(im+ioff)] = frac + whole;
    }/*endfor*/
    ioff += nm;
  }/*endfor*/

#ifdef DEBUG
  if(myid==0){printf("Model data State.gpt\n");}
  for(iproc=0;iproc<num_proc;iproc++){
   if(myid==iproc) {
     for(im=1;im<=nm*nl_proc_a;im++){
       printf("a[%d] = %g   myid = %d\n",im,a[im],myid);
     }/*endfor*/
   }/*endif*/
   scanf("%d",&il);
   MPI_Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif

/*========================================================================*/
/* Zero the matrix */

 cputime(&t0);
 
  nl_tot2 = nl_tot*nl_tot;
  for(il=1;il<=nl_tot2;il++){c[il] = 0.0;}

/*========================================================================*/
/* Here begins the routine                                                */
/*========================================================================*/
/* I) Initialization */

  nsend_a   = nm*nl_proc_a;
  nrecv_b   = nsend_a;
  nl_proc_b_st_now = nl_proc_a_st;
  nl_proc_b_end_now = nl_proc_a_end;
  for(im=1;im<=nsend_a;im++){b_old[im]=a[im];}

/*========================================================================*/
/* II) Calculation loop */

  for(icomm=1;icomm<=(ncomm+1);icomm++){

    /*----------------*/
    /* copy bold to b */

    for(im=1;im<=nrecv_b;im++){b[im]=b_old[im];}
    nsend_b = nrecv_b;
    if(icomm<=ncomm){
     nrecv_b = nl_proc_b[icomm]*nm;
     if((myid % 2) == 0){
      MPI_Isend(&(b_old[1]),nsend_b,MPI_DOUBLE,send_id,send_id,
                MPI_COMM_WORLD,&request);
      MPI_Irecv(&(b_new[1]),nrecv_b,MPI_DOUBLE, 
                MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&request1);
     } else {
      MPI_Irecv(&(b_new[1]),nrecv_b,MPI_DOUBLE, 
                MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&request);
      MPI_Isend(&(b_old[1]),nsend_b,MPI_DOUBLE,send_id,send_id,
                MPI_COMM_WORLD,&request1);
     }/* endif myid */

#ifdef DEBUG
  MPI_Barrier(MPI_COMM_WORLD);
  if(myid==0){printf("Checking the %d send\n",icomm);}
  ntot = nm*nl_proc_b[icomm];
  for(iproc=0;iproc<num_proc;iproc++){
   if(myid==iproc) {
     for(im=1;im<=ntot;im++){
       printf("b_new[%d] = %g   myid = %d\n",im,b_new[im],myid);
     }/*endfor*/
   }/*endif*/
   scanf("%d",&il);
   MPI_Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif

    }/*endif*/
    MPI_Barrier(MPI_COMM_WORLD);

 done=0;
 if(icomm != ncomm+1 || myid < num_proc/2){
  ioff=0;
  for(il=nl_proc_a_st;il<=nl_proc_a_end;il++){
   joff=0;
   upper = nl_proc_b_end_now;
   if(nl_proc_a_end == nl_proc_b_end_now){upper = il;}
   for(jl=nl_proc_b_st_now;jl<=upper;jl++){
    knd = jl + (il-1)*nl_tot;
    for(im=1;im<=nm;im++){
     ind = im + ioff;
     jnd = im + joff;
     c[knd] += a[ind]*b[jnd];
    }/* endfor im */
    joff += nm;
   }/* endfor jl */
   ioff += nm;
  }/* endfor il */
 done=1;
 }/* endif */

#ifdef DEBUG
  MPI_Barrier(MPI_COMM_WORLD);
  if(myid==0){printf("Checking the %d multiply\n",icomm);}
  for(iproc=0;iproc<num_proc;iproc++){
   if(myid==iproc) {
       printf("myid = %d did %d nl_a_st %d nl_a_end %d nl_b_st %d nl_b_end %d\n",
        myid,done,nl_proc_a_st,nl_proc_a_end,nl_proc_b_st_now,nl_proc_b_end_now);
   }/*endif*/
   scanf("%d",&il);
   MPI_Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif


    /*-------------------*/
    /* copy bnew to bold */
    if(icomm<=ncomm){
      for(im=1;im<=nrecv_b;im++){b_old[im]=b_new[im];}
      nl_proc_b_st_now = nl_proc_b_st[icomm];
      nl_proc_b_end_now = nl_proc_b_end[icomm];
    }/*endif*/
  }/*endfor*/


/*========================================================================*/
/* V) Collect, symmetrize and reduce matrix */

 joff=0;
 for(il=1;il<=nl_tot;il++){
  ind  = (il-1)*nl_tot + il;
  c_red[ind] = c[ind];
  for(jl=il+1;jl<=nl_tot;jl++){
    ind  = (il-1)*nl_tot + jl;
    jnd =  (jl-1)*nl_tot + il;
    c_tmp = c[ind]+c[jnd];
    c_red[ind] = c_tmp;
    c_red[jnd] = c_tmp;
   }
 }

 MPI_Allreduce(&(c_red[1]),&(c[1]),nl_tot2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 

#ifdef DEBUG
 MPI_Barrier(MPI_COMM_WORLD);
 if(myid==0){
  for(il=1;il<=nl_tot2;il++){
   printf("c[%d] = %.12g\n",il,c[il]);
  }
 }
#endif

 cputime(&t1);
 printf("Time on proc %d  %.12g\n",myid,(t1-t0));
 dt_tmp = t1-t0;
 MPI_Reduce(&dt_tmp,&dt,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
 if(myid==0){
   printf("Maximum time = %.12g\n",dt);
 }

 MPI_Finalize();

/*========================================================================*/
}/* end main */
/*========================================================================*/



void cputime(double *time)
{
  static clock_t itime;
  static double to=0.,tn=0.;

  itime = clock();
  tn = (double)((double)itime/(double)CLOCKS_PER_SEC);

  if(to>tn){
    *time = (MAXTIME-to)+(tn-MINTIME);
  } else {
    *time = tn-to;
  }
  to = tn;
}

