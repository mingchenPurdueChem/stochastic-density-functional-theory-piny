#define HP_VECLIB
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "mpi.h"


#define MAXTIME 2147.483648
#define MINTIME (-MAXTIME)

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC  1000000
#endif

#define DEBUG_OFF 

void cputime(double *);



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void main (int argc, char *argv[])

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */

  
  double frac,whole;
  double *x;
  double *x_temp,*xt_temp,*x_temp_pt,*xt_temp_pt;
  double *c,*c_red;
  double t0,t1;
  double dt_tmp,dt;
  int js;
  int nfull,nblock;
  int irem,idiv;
  int ig,is,ioff,ioff_t,proc_rem;
  int i,ioff_temp;
  int nstate_max;
  int nstate_proc,nstate_ng_now,nstate_proc_max;
  int nstate_proc_st;
  int nstate_ng_proc;
  int ioff_x;
  int nstate2,nstate_ng_proc_use,joff,ind,jnd;
  int num_proc,iproc,itemp,iii;
  int sendcounts,recvcounts;
  int myid;
/*---------------------*/
  int ng_tot      = 50000;
  int nstate      = 128;

/*========================================================================*/
/* 0) Initialize MPI                                                      */

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

/*========================================================================*/
/* 0) Calculate some constants                                            */

  idiv            = nstate/num_proc;
  irem            = (nstate % num_proc);
  nstate_proc_max = (irem > 0 ? idiv+1 : idiv);
  nstate_max      = nstate_proc_max*num_proc;
  nstate_proc     = (myid < irem ? idiv+1 : idiv);
  nstate_ng_proc = ng_tot/num_proc;
  if( (ng_tot % num_proc) !=0){nstate_ng_proc++;}
  if(myid <= irem) {
    nstate_proc_st = myid*(idiv+1)+1;
   } else {
    nstate_proc_st = irem*(idiv+1) + (myid-irem)*idiv+1;
   }/*endif*/

 
/*========================================================================*/
/* I) Malloc the memory                                                   */

  x          = (double *)malloc((ng_tot*nstate_proc_max)*sizeof(double))-1;
  x_temp     = (double *)malloc((nstate_ng_proc*nstate_max)*sizeof(double))-1;
  xt_temp    = (double *)malloc((nstate_ng_proc*nstate_max)*sizeof(double))-1;
  c          = (double *)malloc((nstate*nstate)*sizeof(double))-1;
  c_red      = (double *)malloc((nstate*nstate)*sizeof(double))-1;

/*========================================================================*/
/* II) Fill the pos vector structure with model data                      */

  ioff=0;
  for(is=1;is<=nstate_proc;is++){
    for(ig=1;ig<=ng_tot;ig++){
      frac = (double)ig;
      frac /= 10;
      whole = (double)(nstate_proc_st+is-1);
      x[(ig+ioff)] = frac + whole;
    }/*endfor*/
   ioff += ng_tot;
  }/*endfor*/

#ifdef DEBUG
  if(myid==0){
    printf("This is x\n");
  }/*endif*/
  for(iproc=1;iproc<=num_proc;iproc++){
    if(myid==iproc-1){
     for(is=1;is<=nstate_proc;is++){
       ioff_x = (is-1)*ng_tot;
       for(ig=1;ig<=ng_tot;ig++){
        printf("%d %d %d %g\n",iproc,is,ig,x[(ig+ioff_x)]);
       }/*endfor*/
     }/*endfor*/
    }/*endif*/
    scanf("%d",&iii);
    MPI_Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif

/*========================================================================*/
/* III) Extract the position data from clatoms_pos                        */


 cputime(&t0);

 if(myid==0){
   printf("n_ng_proc %d n_max %d n_proc %d num_proc %d n_proc_max %d \n",
        nstate_ng_proc,nstate_max,nstate_proc,num_proc,nstate_proc_max);
   printf("ng_tot = %d \n",ng_tot);
 }

  for(i=1;i<=nstate_ng_proc*nstate_max;i++){x_temp[i]=0.0;}
  proc_rem = (ng_tot%num_proc);
  for(is=1;is<=nstate_proc;is++){
    ioff = 0;
    ioff_x = (is-1)*ng_tot;
    for(iproc=1;iproc<=num_proc;iproc++){
      ioff_temp = (is-1)*nstate_ng_proc
                + (iproc-1)*(nstate_proc_max*nstate_ng_proc);
      nstate_ng_now = nstate_ng_proc;
      if((iproc>proc_rem)&&(proc_rem>0)){nstate_ng_now = nstate_ng_proc-1;}
      for(ig=1;ig<=nstate_ng_now;ig++){
        itemp = ig+ioff_temp;
        i     = ig+ioff;
        x_temp[itemp] = x[(i+ioff_x)];
      }/*endfor*/   
      ioff += nstate_ng_now;
    }/*endfor*/   
  }/*endfor*/   

#ifdef DEBUG
  if(myid==0){
    printf("This is x_temp\n");
  }/*endif*/
  for(iproc=1;iproc<=num_proc;iproc++){
    if(myid==iproc-1){
     for(ig=1;ig<=nstate_ng_proc*nstate_max;ig++){
      printf("%d %d %g\n",iproc,ig,x_temp[ig]);
     }/*endfor*/
    }/*endif*/
    scanf("%d",&iii);
    MPI_Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif

/*========================================================================*/
/* IV) Send the position data                                             */

  sendcounts = nstate_ng_proc*nstate_proc_max;
  recvcounts = nstate_ng_proc*nstate_proc_max;

  xt_temp_pt = xt_temp+1;
  x_temp_pt = x_temp+1;
  MPI_Alltoall(x_temp_pt,sendcounts,MPI_DOUBLE,xt_temp_pt,recvcounts,
                MPI_DOUBLE,MPI_COMM_WORLD);

#ifdef DEBUG
  if(myid==0){
    printf("This is xt_temp\n");
  }/*endif*/
  for(is=1;is<=num_proc;is++){
    if(myid==is-1){
      for(ig=1;ig<=nstate_ng_proc*nstate_max;ig++){
         printf("%d %d %g\n",ig,is,xt_temp[ig]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&ig);
    MPI_Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif
  

/*========================================================================*/
/* V) Internal rearrangements                                             */

  if(irem>0){
   /* full blocks */
   nfull = nstate_proc_max*nstate_ng_proc*irem;
   for(ig=1;ig<=nfull;ig++){x_temp[ig] = xt_temp[ig];}
   /* partial blocks */
   nblock = (nstate_proc_max-1)*nstate_ng_proc;
   ioff = nfull; joff=nfull;
   for(iproc=irem+1;iproc<=num_proc;iproc++){
    for(ig=1;ig<=nblock;ig++){
     x_temp[(ig+ioff)] = xt_temp[(ig+joff)];
    }/*endfor*/
    ioff += nblock;joff += (nblock+nstate_ng_proc);
   }/*endfor*/
  }else{
   /* Only full blocks */
   nfull = nstate_ng_proc*nstate;
   for(ig=1;ig<=nfull;ig++){
    x_temp[ig] = xt_temp[ig];
   }/*endfor*/
  }/*endif*/
    
#ifdef DEBUG
  if(myid==0){
    printf("This is xt_temp\n");
  }/*endif*/
  for(is=1;is<=num_proc;is++){
    if(myid==is-1){
      for(ig=1;ig<=nstate_ng_proc*nstate_max;ig++){
         printf("%d %d %g\n",ig,is,x_temp[ig]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&ig);
    MPI_Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif
  

/*========================================================================*/
/* V) Perform the mini multiplies                                         */


  nstate2 = nstate*nstate;
  for(is=1;is<=nstate2;is++){c_red[is] = 0.0;c[is] = 0.0;}

  nstate_ng_proc_use = ng_tot/num_proc;
  irem               = ng_tot % num_proc;
  if(myid<irem){nstate_ng_proc_use++;}

  ioff = 0;
  for(is=1;is<=nstate;is++){
   joff = 0;
   for(js=1;js<=is;js++){
    ind  = (is-1)*nstate + js;
    for(ig=1;ig<=nstate_ng_proc_use;ig++){
      c_red[ind] += (x_temp[(joff+ig)]*x_temp[(ioff+ig)]);
    }/*endfor*/
    joff += nstate_ng_proc;
    jnd = (js-1)*nstate + is;
    c_red[jnd] = c_red[ind];
   }/*endfor*/
   ioff += nstate_ng_proc;
  }/*endfor*/


/*========================================================================*/
/* VI) Reduce the matrix                                                  */

 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Allreduce(&(c_red[1]),&(c[1]),nstate2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 

#ifdef DEBUG
 MPI_Barrier(MPI_COMM_WORLD);
 for(is=1;is<=num_proc;is++){
  if(myid==is-1){
   for(is=1;is<=nstate2;is++){
    printf("c[%d] = %.12g\n",is,c[is]);
   }/*endfor*/
  }/*endif*/
  scanf("%d",&ig);
  MPI_Barrier(MPI_COMM_WORLD);
 }/*endfor*/
#endif

 cputime(&t1);
 printf("Time on proc %d  %.12g\n",myid,(t1-t0));
 dt_tmp = t1-t0;
 MPI_Reduce(&dt_tmp,&dt,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
 if(myid==0){
   printf("Maximum time = %.12g\n",dt);
 }


/*========================================================================*/
/* VI) Finalize                                                           */

  MPI_Finalize();

/*========================================================================*/
   }/*end routine*/
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

