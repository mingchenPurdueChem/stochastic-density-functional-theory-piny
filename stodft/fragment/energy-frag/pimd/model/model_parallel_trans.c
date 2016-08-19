#define HP_VECLIB
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "mpi.h"

#define DEBUG 

typedef struct clatoms_pos{
  double *x;
} CLATOMS_POS;


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void main (int argc, char *argv[])

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */

  
  CLATOMS_POS *clatoms_pos;
  double frac,whole;
  double *clatoms_pos_x;
  double *x_temp,*xt_temp,*x_temp_pt,*xt_temp_pt;
  int iatm,ip,ioff,ioff_t,proc_rem;
  int i,ioff_temp;
/*  int natm_tot = 12;*/
/*  int pi_beads = 8;*/
  int natm_tot      = 7;
  int pi_beads      = 16;
  int pi_atm_proc,pi_atm_now;
  int pi_beads_proc;
  int *sendcounts,*senddspls,*recvcounts,*recvdspls;
  int num_proc,iproc,itemp,iii;
  int myid = 1;

/*========================================================================*/
/* 0) Initialize MPI                                                      */

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

/*========================================================================*/
/* 0) Calculate some constants                                            */

  pi_beads_proc = pi_beads/num_proc;
  if((pi_beads % num_proc) !=0){
   printf("Number of beads/number of processors not an integer\n");
   MPI_Finalize();
   exit(1);
  }/*endif*/
  pi_atm_proc = natm_tot/num_proc;
  if( (natm_tot % num_proc) !=0){pi_atm_proc++;}

/*========================================================================*/
/* I) Malloc the memory                                                   */

  clatoms_pos  = (CLATOMS_POS *)malloc(pi_beads_proc*sizeof(CLATOMS_POS))-1;
  for(ip=1;ip<=pi_beads_proc;ip++){
    clatoms_pos[ip].x  = (double *)malloc(natm_tot*sizeof(double))-1;
  }/*endfor*/

  x_temp     = (double *)malloc((pi_atm_proc*pi_beads)*sizeof(double))-1;
  xt_temp    = (double *)malloc((pi_atm_proc*pi_beads)*sizeof(double))-1;
  
  sendcounts     = (int *)malloc(num_proc*sizeof(int))-1;
  senddspls      = (int *)malloc(num_proc*sizeof(int))-1;
  recvcounts     = (int *)malloc(num_proc*sizeof(int))-1;
  recvdspls      = (int *)malloc(num_proc*sizeof(int))-1;

/*========================================================================*/
/* II) Fill the pos vector structure with model data                      */

  for(ip=1;ip<=pi_beads_proc;ip++){
    clatoms_pos_x = clatoms_pos[ip].x;
    for(iatm=1;iatm<=natm_tot;iatm++){
      frac = (double)iatm;
      frac /= 100;
      whole = (double)(myid*pi_beads_proc+ip);
      clatoms_pos_x[(iatm)] = frac + whole;
    }/*endfor*/
  }/*endfor*/

/*========================================================================*/
/* III) Extract the position data from clatoms_pos                        */

  for(i=1;i<=pi_atm_proc*pi_beads;i++){x_temp[i]=0.0;}
  proc_rem = (natm_tot%num_proc);
  for(ip=1;ip<=pi_beads_proc;ip++){
    ioff = 0;
    clatoms_pos_x = clatoms_pos[ip].x;
    for(iproc=1;iproc<=num_proc;iproc++){
      ioff_temp = (ip-1)*pi_atm_proc+(iproc-1)*(pi_beads_proc*pi_atm_proc);
      pi_atm_now = pi_atm_proc;
      if((iproc>proc_rem)&&(proc_rem>0)){pi_atm_now = pi_atm_proc-1;}
      for(iatm=1;iatm<=pi_atm_now;iatm++){
        itemp = iatm+ioff_temp;
        i     = iatm+ioff;
        x_temp[itemp] = clatoms_pos_x[i];
      }/*endfor*/   
      ioff += pi_atm_now;
    }/*endfor*/   
  }/*endfor*/   

#ifdef DEBUG
  if(myid==0){
    printf("This is x_temp\n");
  }/*endif*/
  for(iproc=1;iproc<=num_proc;iproc++){
    if(myid==iproc-1){
     for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
      printf("%d %d %g\n",iproc,iatm,x_temp[iatm]);
     }/*endfor*/
    }/*endif*/
    scanf("%d",&iii);
    MPI_Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif

/*========================================================================*/
/* IV) Send the position data                                             */

  for(ip=1;ip<=num_proc;ip++){
    sendcounts[ip] = pi_atm_proc*pi_beads_proc;
    senddspls[ip]  = (ip-1)*pi_atm_proc*pi_beads_proc;
    recvcounts[ip] = pi_atm_proc*pi_beads_proc;
    recvdspls[ip]  = (ip-1)*pi_atm_proc*pi_beads_proc;
  }/*endfor*/

  xt_temp_pt = xt_temp+1;
  x_temp_pt = x_temp+1;
  MPI_Alltoallv(x_temp_pt,&(sendcounts[1]),&(senddspls[1]),
                MPI_DOUBLE,xt_temp_pt,&(recvcounts[1]),
                &(recvdspls[1]), MPI_DOUBLE,MPI_COMM_WORLD);

#ifdef DEBUG
  if(myid==0){
    printf("This is xt_temp\n");
  }/*endif*/
  for(ip=1;ip<=num_proc;ip++){
    if(myid==ip-1){
      for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
         printf("%d %d %g\n",iatm,ip,xt_temp[iatm]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&iatm);
    MPI_Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif

/*========================================================================*/
/* V) Rearrange the position data                                         */

  for(ip=1;ip<=pi_beads;ip++){
    ioff_t = (ip-1)*pi_atm_proc;
    for(iatm=1;iatm<=pi_atm_proc;iatm++){
      ioff = (iatm-1)*pi_beads+ip;
      i    = iatm+ioff_t;
      x_temp[ioff] = xt_temp[i];
    }/*endfor*/
  }/*endfor*/

#ifdef DEBUG
  if(myid==0){
    printf("This is x_temp\n");
  }/*endif*/
  for(ip=1;ip<=num_proc;ip++){
    if(myid==ip-1){
      for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
         printf("%d %d %g\n",iatm,ip,x_temp[iatm]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&iatm);
    MPI_Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif

/*========================================================================*/
/* VI) Transform the position data                                         */

/*========================================================================*/
/* VII) Rearrange the transformed position data                             */


  for(ip=1;ip<=pi_beads;ip++){
    ioff_t = (ip-1)*pi_atm_proc;
    for(iatm=1;iatm<=pi_atm_proc;iatm++){
      ioff = (iatm-1)*pi_beads+ip;
      i    = iatm+ioff_t;
      xt_temp[i] = x_temp[ioff];
    }/*endfor*/
  }/*endfor*/

#ifdef DEBUG
  if(myid==0){
    printf("This is xt_temp\n");
  }/*endif*/
  for(ip=1;ip<=num_proc;ip++){
    if(myid==ip-1){
      for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
         printf("%d %d %g\n",iatm,ip,xt_temp[iatm]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&iatm);
    MPI_Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif

/*========================================================================*/
/* VIII) Send the transformed position data                               */

  for(ip=1;ip<=num_proc;ip++){
    sendcounts[ip] = pi_atm_proc*pi_beads_proc;
    senddspls[ip]  = (ip-1)*pi_atm_proc*pi_beads_proc;
    recvcounts[ip] = pi_atm_proc*pi_beads_proc;
    recvdspls[ip]  = (ip-1)*pi_atm_proc*pi_beads_proc;
  }/*endfor*/

  xt_temp_pt = xt_temp+1;
  x_temp_pt = x_temp+1;
  MPI_Alltoallv(xt_temp_pt,&(sendcounts[1]),&(senddspls[1]),
                MPI_DOUBLE,x_temp_pt,&(recvcounts[1]),
                &(recvdspls[1]), MPI_DOUBLE,MPI_COMM_WORLD);

#ifdef DEBUG
  if(myid==0){
    printf("This is x_temp\n");
  }/*endif*/
  for(ip=1;ip<=num_proc;ip++){
    if(myid==ip-1){
      for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
         printf("%d %d %g\n",iatm,ip,x_temp[iatm]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&iatm);
    MPI_Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif


/*========================================================================*/
/* IX) Extract the transformed position data                               */

  proc_rem = (natm_tot%num_proc);
  for(ip=1;ip<=pi_beads_proc;ip++){
    ioff = 0;
    clatoms_pos_x = clatoms_pos[ip].x;
    for(iproc=1;iproc<=num_proc;iproc++){
      ioff_temp = (ip-1)*pi_atm_proc+(iproc-1)*(pi_beads_proc*pi_atm_proc);
      pi_atm_now = pi_atm_proc;
      if((iproc>proc_rem)&&(proc_rem>0)){pi_atm_now = pi_atm_proc-1;}
      for(iatm=1;iatm<=pi_atm_now;iatm++){
        itemp = iatm+ioff_temp;
        i     = iatm+ioff;
        clatoms_pos_x[i] = x_temp[itemp];
      }/*endfor*/   
      ioff += pi_atm_now;
    }/*endfor*/   
  }/*endfor*/   

#ifdef DEBUG
  if(myid==0){
    printf("This is clatoms_pos_x\n");
  }/*endif*/
  for(iproc=1;iproc<=num_proc;iproc++){
    if(myid==iproc-1){
     for(ip=1;ip<=pi_beads_proc;ip++){
       clatoms_pos_x = clatoms_pos[ip].x;
       for(iatm=1;iatm<=natm_tot;iatm++){
        printf("%d %d %d %g\n",iproc,ip,iatm,clatoms_pos_x[iatm]);
       }/*endfor*/
     }/*endfor*/
    }/*endif*/
    scanf("%d",&iii);
    MPI_Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif


/*========================================================================*/
/* X) Finalize MPI                                                         */

  MPI_Finalize();

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/



