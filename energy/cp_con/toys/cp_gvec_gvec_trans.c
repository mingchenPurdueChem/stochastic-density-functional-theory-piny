#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "mpi.h"


#define DEBUG
#define Alltoall MPI_Alltoall
#define Barrier MPI_Barrier

void countkvec3d(int *,double ,int *,double *);

void setkvec3d_sm(int ,double ,int *,double *,
                  int *, int *, int *, int *, int *,int *);


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
  double hmat[10],hmati[10];
  int icoef_start;
  int *kastore,*kbstore,*kcstore;
  int *istart_fft;
  int *ibreak1,*ibreak2;
  int nfft_c;
  int nktot;
  int ncoef_proc,ncoef_proc_max;
  int myid,num_proc;
  int iproc,iii;
  int ig,ioff;
  int kmaxv[4];
  int ncoef;
  MPI_Comm world;
/*---------------------*/
/*========================================================================*/
/* 0) Initialize MPI                                                      */

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  world = MPI_COMM_WORLD;

/*========================================================================*/
/* 0) Calculate some constants                                            */

   for(i=1;i<=9;i++){hmat[i]=0.0;hmati[i]=0.0;}  

   hmat[1] = hmat[5] = hmat[9] = 10.0;
   hmati[1] = hmati[5] = hmati[9] = 1.0/hmat[1];
   ecut = 1.5;
   kmaxv[1]=30;
   kmaxv[2]=30;
   kmaxv[3]=30;
   
   countkvec3d(&nktot,&nfft_c,ecut,kmaxv,hmati);
   ncoef = nktot+1;

/*========================================================================*/
/* I) Malloc the memory                                                   */

  kastore = (int *) malloc(ncoef*sizeof(int))-1;
  kbstore = (int *) malloc(ncoef*sizeof(int))-1;
  kcstore = (int *) malloc(ncoef*sizeof(int))-1;
  istart_fft = (int *) malloc(nfft_c*sizeof(int))-1;
  ibreak1 = (int *) malloc(ncoef*sizeof(int))-1;
  ibreak2 = (int *) malloc(ncoef*sizeof(int))-1;
  c       = (double *) malloc(ncoef*sizeof(double))-1;
  c_temp  = (double *) malloc(ncoef*sizeof(double))-1;
  ct_temp = (double *) malloc(ncoef*sizeof(double))-1;

/*========================================================================*/
/* Fill the kvector arrays                                                */

  setkvec3d_sm(nktot,ecut,kmaxv,hmati,
               kastore, kbstore, kcstore, ibreak1, ibreak2,istart_fft);

/*========================================================================*/
/* coefs per proc                                                         */

  ncoef_proc = ncoef/num_proc;
  irem = (ncoef % num_proc);
  ncoef_proc_max = ncoef_proc;
  if(irem != 0) ncoef_proc_max++;
  if(myid < irem) ncoef_proc++;
  
  if(myid <= irem) icoef_start = ncoef_proc_max*myid + 1;
  if(myid > irem) icoef_start = ncoef_proc_max*irem 
                              + (myid-irem)*ncoef_proc+1;
  icoef_end = icoef_start + ncoef_proc - 1;
 
/*========================================================================*/
/* II) Fill the pos vector structure with model data                      */

  ioff = icoef_start-1;
  for(ig=1;ig<=ncoef_proc;ig++){
      c[ig] = 100000*kastore[(ig+ioff)] + 100*kbstore[(ig+ioff)]
            + kcstore[(ig+ioff)];
  }/*endfor*/

#ifdef DEBUG
  if(myid==0){
    printf("These are the coefficients\n");
  }/*endif*/
  for(iproc=1;iproc<=num_proc;iproc++){
    if(myid==iproc-1){
    for(ig=1;ig<=ncoef_proc;ig++){
        printf("%d %d %g\n",iproc,ig,c[ig]);
     }/*endfor*/
    }/*endif*/
    scanf("%d",&iii);
    Barrier(MPI_COMM_WORLD);
  }/*endfor*/
#endif

/*--------------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* COUNT THE NUMBER OF K VECTORS ON LARGE GRID */
/*==========================================================================*/

void countkvec3d(int *nktot,int *nfft_c_ret,
                 double ecut,int *kmaxv,double *hmatik)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int iii,icount;
  int nfft_c=0;
  int i, kbmin, kcmin, kbmax, kcmax, kamax, ka, kb, kc;
  double xk, yk, zk;
  double aka, akb, akc;
  double tpi,try;

/*==========================================================================*/
/* Count the kvectors */

  tpi = 2.0*M_PI;
  icount = 0;
  kamax = kmaxv[1];

/*********************************/

   for (i = 1; i <= kamax; ++i) {
    aka = (double) i;
    xk = aka * hmatik[1] * tpi;
    yk = aka * hmatik[4] * tpi;
    zk = aka * hmatik[7] * tpi;
    try = (xk * xk + yk * yk + zk * zk) * .5;
    if (try > ecut) {
      break;
    }
  }

  kamax = i - 1;

/***********************************/

  for (ka = 0; ka <= kamax; ++ka) {
    aka = (double) ka;
    kbmin = -kmaxv[2];
    if (ka == 0) {
      kbmin = 0;
    }
    for (i = kbmin; i <= 0; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      try = (xk * xk + yk * yk + zk * zk) * .5;
      if (try <= ecut) {
	break;
      }
    }

/*********************************/

    kbmin = i;
    for (i = 1; i <= kmaxv[2]; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      try = (xk * xk + yk * yk + zk * zk) * .5;
      if (try > ecut) {
	break;
      }
    }
    
    kbmax = i - 1;
    for (kb = kbmin; kb <= kbmax; ++kb) {
      
      akb = (double) kb;
      kcmin = -kmaxv[3];
      if (ka == 0 && kb == 0) {
	kcmin = 1;
      }
      for (i = kcmin; i <= 0; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try <= ecut) {
	  break;
	}
      }
/*********************************/

      kcmin = i;
      for (i = 1; i <= kmaxv[3]; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try > ecut) {
	  break;
	}
      }

      kcmax = i - 1;
      akc = (double) kcmin;
      if(kcmin <= kcmax) nfft_c++;
      for (kc = kcmin; kc <= kcmax; ++kc) {
	++icount;
      }
    }
  }
  *nktot = icount;
/*--------------------------------------------------------------------------*/
} /* countkvec3d */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* THIS SUBROUTINE DETERMINES THE ALLOWED SPHERICALLY TRUNCATED */
/* HALF SPACE K (.I.E. G) VECTORS GIVEN A CUTOFF. IT IS USED BY */
/* CP MODULES. */
/* SETS UP THE K VECTORS FOR A GIVEN CUTOFF AND SHAPE */
/*==========================================================================*/

void setkvec3d_sm(int nktot,double ecut,int *kmax_cp,double *hmatik,
                  int *kastore, int *kbstore, int *kcstore, 
		  int *ibreak1, int *ibreak2,int *istart_fft)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int i1, i2, i3;
  int nfft_c=0;

  int i, kbmin, kcmin, kbmax, kcmax, kamax, ka, kb, kc;
  double xk, yk, zk;
  int icount;
  double aka, akb, akc;
  double tpi, try;

/*==========================================================================*/
  /* SETUP THE KVECTORS */
  
    for(i=1;i<=nktot;i++){
      ibreak1[i] = 0;
      ibreak2[i] = 0;
    }


  tpi = M_PI * 2.;
  icount = 0;
  i1 = kmax_cp[1];
  for (i = 1; i <= i1; ++i) {
    aka = (double) i;
    xk = aka * hmatik[1] * tpi;
    yk = aka * hmatik[4] * tpi;
    zk = aka * hmatik[7] * tpi;
    try = (xk * xk + yk * yk + zk * zk) * .5;
    if (try > ecut) {
      break;
    }
  }

  kamax = i - 1;
  i1 = kamax;
  for (ka = 0; ka <= i1; ++ka) {
    aka = (double) ka;
    kbmin = -kmax_cp[2];
    if (ka == 0) {
      kbmin = 0;
    }
    for (i = kbmin; i <= 0; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      try = (xk * xk + yk * yk + zk * zk) * .5;
      if (try <= ecut) {
	break;
      }
    }
    kbmin = i;
    i2 = kmax_cp[2];
    for (i = 1; i <= i2; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      try = (xk * xk + yk * yk + zk * zk) * .5;
      if (try > ecut) {
	break;
      }
    }

    kbmax = i - 1;
    i2 = kbmax;
    for (kb = kbmin; kb <= i2; ++kb) {
      akb = (double) kb;
      kcmin = -kmax_cp[3];
      if (ka == 0 && kb == 0) {
	kcmin = 1;
      }
      for (i = kcmin; i <= 0; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try <= ecut) {
	  break;
	}
      }
      
      kcmin = i;
      i3 = kmax_cp[3];
      for (i = 1; i <= i3; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try > ecut) {
	  break;
	}
      }
      
      kcmax = i - 1;
      i3 = kcmax;
      for (kc = kcmin; kc <= i3; ++kc) {
	++icount;
	aka = (double) ka;
	akb = (double) kb;
	akc = (double) kc;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	kastore[icount] = ka;
	kbstore[icount] = kb;
	kcstore[icount] = kc;
	if (kc == kcmin) {
	  ibreak1[icount] = 1;
          nfft_c++;
          istart_fft[nfft_c] = icount;
	}
	if (kc < kcmax) {
	  ibreak2[icount] = 1;
	}
      }
    }
  }
  if(nktot!=icount){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Mismatch number of small kvectors\n");
       printf("%d vs %d\n",icount,nktot);
       printf("Contact technical support\n");
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
  }

  ++icount;
  kastore[icount] = 0;
  kbstore[icount] = 0;
  kcstore[icount] = 0;
/*--------------------------------------------------------------------------*/
} /* setkvec3d_sm */
/*==========================================================================*/

















