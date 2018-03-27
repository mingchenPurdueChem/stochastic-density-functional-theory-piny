/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: math_matrix                                  */
/*                                                                          */
/* Linear Algebra, timing ...                                               */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#define CLOCKS_PER_SEC_C  1000000
#define MAXTIME 2147.48



/*==========================================================================*/
/* Uniform random numbers */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

double ran_essl(double *qseed)

/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
  int n=1,ierr=0;
  double x;
/*==========================================================================*/

  DURAND(qseed,&n,&x,&ierr);
  if(ierr==1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Error in random number generator: Durand\n");
    if(n    < 0           ){printf("Parameter n=%d < 0\n"             ,n);}
    if(*qseed< 1.0        ){printf("Parameter qseed=%g < 1.0\n"       ,*qseed);}
    if(*qseed>2147483646.0){printf("Parameter qseed=%g > 2147483646\n",*qseed);}
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/  
  return x;

/*-------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/


/*===============================================================*/
/* Gaussian Random numbers */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void gaussran(int nran, int *iseed, int *iseed2, double *qseed, double gauss[])

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
   int i,iii,loop;
   double twopi,rad2,al,r,phi,arg;
/*========================================================================*/
/* I) Constants */
    twopi = 2.0*M_PI;
    rad2 = sqrt(2.0);
    loop = nran/2;

/*========================================================================*/
/* II) Make nran (or nran-1 if nran odd) Gaussian random numbers          */
    for(i=1;i<=loop;i++){
/*------------------------------------------------------------------------*/
/* A) uniform random numbers in r and phi */
       r   = ran_essl(qseed);
       r   = MAX(r,1e-30);
       r   = MIN(r,1.0);
       phi = ran_essl(qseed);
/*------------------------------------------------------------------------*/
/* B) Gaussify in x and y*/
       al  = sqrt(-log(r))*rad2;
       arg = twopi*phi;
       gauss[2*i-1] = al*cos(arg);
       gauss[2*i]   = al*sin(arg);
     }/*endfor*/
/*========================================================================*/
/* III) Make one more if nran is odd */
    if((nran % 2)!=0){
       r   = ran_essl(qseed);
       r   = MAX(r,1e-30);
       r   = MIN(r,1.0);
       phi = ran_essl(qseed);
       arg = twopi*phi;
       al  = sqrt(-log(r))*rad2;
       gauss[nran] = al*cos(arg);
     }/*endif*/

/*------------------------------------------------------------------------*/
}/*end routine*/
/*========================================================================*/

/*===============================================================*/
/* Gaussian Random numbers */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void gaussran2(int nran, int *iseed, int *iseed2, double *qseed, double gauss[])

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
   int i,iii,loop;
   double twopi,rad2,al,r,phi,arg;
/*========================================================================*/
/* I) Constants */
    twopi = 2.0*M_PI;
    rad2 = sqrt(2.0);
    loop = nran/2;

/*========================================================================*/
/* II) Make nran (or nran-1 if nran odd) Gaussian random numbers          */
    for(i=0;i<loop;i++){
/*------------------------------------------------------------------------*/
/* A) uniform random numbers in r and phi */
       r   = ran_essl(qseed);
       r   = MAX(r,1e-30);
       r   = MIN(r,1.0);
       phi = ran_essl(qseed);
/*------------------------------------------------------------------------*/
/* B) Gaussify in x and y*/
       al  = sqrt(-log(r))*rad2;
       arg = twopi*phi;
       gauss[2*i] = al*cos(arg);
       gauss[2*i+1]   = al*sin(arg);
     }/*endfor*/
/*========================================================================*/
/* III) Make one more if nran is odd */
    if((nran % 2)!=0){
       r   = ran_essl(qseed);
       r   = MAX(r,1e-30);
       r   = MIN(r,1.0);
       phi = ran_essl(qseed);
       arg = twopi*phi;
       al  = sqrt(-log(r))*rad2;
       gauss[nran-1] = al*cos(arg);
     }/*endif*/

/*------------------------------------------------------------------------*/
}/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* subroutine to time processes */
/*==========================================================================*/

void cputime(double *time)

/*==========================================================================*/
{
  long int itime;
  static double to=0.,tn=0.;

  itime = clock();
  tn = (double)((double)itime/(double)CLOCKS_PER_SEC_C);
  *time = tn;
  if(tn >= 0 && to >= 0){*time=tn;}
  if(tn < 0  && to >= 0){*time=MAXTIME*2.0+tn;}
  if(tn >= 0 && to <  0){*time=tn+MAXTIME;}
  if(tn <  0 && to <  0){*time=MAXTIME+tn;}

  to = tn;
}
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void  par_cpu_vomit(double cpu,MPI_Comm comm,int nproc,int myid,
                         char name[])

/*==========================================================================*/
{/*begin routine */
#include "../typ_defs/typ_mask.h"
 int i;
    for(i=0;i<nproc;i++){
      if(myid == i){
        printf("%s %d %.12g\n",name,myid,cpu);
      }
      Barrier(comm);
    }
 }/*end routine*/
/*==========================================================================*/



/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/

void matmul_2(double *a1, double *a2, double *a3, int n)
{
  int i, j, k;

  for (i = 1; i <= n; ++i) {
    for (j = 1; j <= n; ++j) {
      a3[(j + (i-1)*n)] = 0.;
    }
  }
  for (i = 1; i <= n; ++i) {
    for (j = 1; j <= n; ++j) {
      for (k = 1; k <= n; ++k) {
	a3[(j + (i-1)*n)] += a1[(j + (k-1)*n)] * a2[(k + (i-1)*n)];
      }
    }
  }
} /* matmul_2 */
/*===============================================================*/


/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
/* Square of a general nxn matrix */
/*===============================================================*/
void matmul_2s(double *a1, double *a3, int n)
{
  int i, j, k;

  /* Function Body */
  for (i = 1; i <= n; ++i) {
    for (j = 1; j <= n; ++j) {
      a3[(j + (i-1)*n)] = 0.;
    }
  }
  for (i = 1; i <= n; ++i) {
    for (j = 1; j <= n; ++j) {
      for (k = 1; k <= n; ++k) {
	a3[(j + (i-1)*n)] += a1[(j + (k-1)*n)] * a1[(k + (i-1)*n)];
      }
    }
  }
} /* matmul_2s */
/*===============================================================*/


/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
/* Matrix product for 3x3 matrices */
/*===============================================================*/
void matmul_3(double *a1, double *a2)
{
  int i, j, k;
  double a3[9]	/* was [3][3] */;

  /* Function Body */
  for (i = 1; i <= 3; ++i) {
    for (j = 1; j <= 3; ++j) {
      a3[(j + (i-1)*3)] = 0.;
    }
  }
  for (i = 1; i <= 3; ++i) {
    for (j = 1; j <= 3; ++j) {
      for (k = 1; k <= 3; ++k) {
	a3[(j + (i-1)*3)] += a1[(j + (k-1)*3)] * a2[(k + (i-1)*3)];
      }
    }
  }
  for (i = 1; i <= 3; ++i) {
    for (j = 1; j <= 3; ++j) {
      a1[(j + (i-1)*3)] = a3[(j + (i-1)*3)];
    }
  }
} /* matmul */
/*===============================================================*/


/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
/* Matrix product of transpose of A1 x transpose A2  for nxn matrices*/
/*===============================================================*/
void matmul_tt(double *a1, double *a2, double *a3, int n)
{
  /* Local variables */
  int i, j, k;
  
  /* Function Body */
  for (i = 1; i <= n; ++i) {
    for (j = 1; j <= n; ++j) {
      a3[(j + (i-1)*n)] = 0.;
    }
  }
  for (i = 1; i <= n; ++i) {
    for (j = 1; j <= n; ++j) {
      for (k = 1; k <= n; ++k) {
	a3[(j + (i-1)*n)] += a1[(k + (i-1)*n)] * a2[(j + (k-1)*n)];
      }
    }
  }
} /* matmul_tt */
/*===============================================================*/


/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
/* Matrix product of transpose of A1 x A2  for nxn matrices*/
/*===============================================================*/
void matmul_t(double *a1, double *a2, double *a3, int n)
{
  /* Local variables */
  int i, j, k;

  /* Function Body */
  for (i = 1; i <= n; ++i) {
    for (j = 1; j <= n; ++j) {
      a3[(j + (i-1) * n)] = 0.;
    }
  }
  for (i = 1; i <= n; ++i) {
    for (j = 1; j <= n; ++j) {
      for (k = 1; k <= n; ++k) {
	a3[(j + (i-1)*n)] += a1[(k + (i-1)*n)] * a2[(k + (j-1)*n)];
      }
    }
  }
} /* matmul_t */
/*===============================================================*/


/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
/* Matrix product of A1 x transpose A2 for nxn matrices */
/*===============================================================*/
void matmul_t2(double *a1, double *a2, double *a3, int n)
{
  /* Local variables */
  int i, j, k;


  /* Function Body */
  for (i = 1; i <= n; ++i) {
    for (j = 1; j <= n; ++j) {
      a3[(j + (i-1) * n)] = 0.;
    }
  }
  for (i = 1; i <= n; ++i) {
    for (j = 1; j <= n; ++j) {
      for (k = 1; k <= n; ++k) {
	a3[(j + (i-1)*n)] += a1[(i + (k-1)*n)] * a2[(j + (k-1)*n)];
      }
    }
  }
} /* matmul_t2 */
/*===============================================================*/


/*===============================================================*/
/* Diagonalize a 3x3 */
/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
void diag33(double *vmat, double *veig, double *veigv, double *fv1, 
            double *fv2)

{/*begin routine */
   int matz=1,ndiag=3,ndiagm=3;
   int ierr;
   RS(&ndiagm,&ndiag,&(vmat[1]),&(veig[1]),&matz,
       &(veigv[1]),&(fv1[1]),&(fv2[1]),&ierr);

}/* end routine */
/*===============================================================*/



/*===============================================================*/
/*  Inverse of a 3x3 */
/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/

void gethinv(double *hmat, double *hmati, double *deth, int iperd)

/*===============================================================*/
   {/*begin routine */
/*===============================================================*/
  double vol;
  int i;
/*===============================================================*/
/* gets inverse, hmati, of the iperd dimensional matrix hmat */
/* (stored as a 3 x 3) */

  *deth = 0.0;
  for(i=1;i<=9;i++){hmati[i]=0.0;}

/*===============================================================*/
/* Perd=3 */

  if (iperd == 3) {
    vol = (hmat[1] * (hmat[5] * hmat[9] - hmat[8] * hmat[6]) + 
	   hmat[4] * (hmat[8] * hmat[3] - hmat[2] * hmat[9]) + 
	   hmat[7] * (hmat[2] * hmat[6] - hmat[5] * hmat[3]));
    *deth = vol;
    hmati[1] = (hmat[5] * hmat[9] - hmat[8] * hmat[6]) / vol;
    hmati[5] = (hmat[1] * hmat[9] - hmat[7] * hmat[3]) / vol;
    hmati[9] = (hmat[1] * hmat[5] - hmat[4] * hmat[2]) / vol;
    hmati[4] = (hmat[7] * hmat[6] - hmat[4] * hmat[9]) / vol;
    hmati[2] = (hmat[3] * hmat[8] - hmat[2] * hmat[9]) / vol;
    hmati[7] = (hmat[4] * hmat[8] - hmat[7] * hmat[5]) / vol;
    hmati[3] = (hmat[2] * hmat[6] - hmat[3] * hmat[5]) / vol;
    hmati[8] = (hmat[7] * hmat[2] - hmat[8] * hmat[1]) / vol;
    hmati[6] = (hmat[3] * hmat[4] - hmat[6] * hmat[1]) / vol;
  }/*endif*/

/*===============================================================*/
/* Perd=2 */

  if (iperd == 2) {
    vol = hmat[1] * hmat[5] - hmat[4] * hmat[2];
    hmati[1] = hmat[5] / vol;
    hmati[5] = hmat[1] / vol;
    hmati[4] = -hmat[4] / vol;
    hmati[2] = -hmat[2] / vol;
    hmati[9] = 1. / hmat[9];
    *deth = vol * hmat[9];
  }/*endif*/

/*===============================================================*/
/* Perd=1,0,cluster_ewald */

  if(iperd <=1 || iperd==4) {
   *deth = hmat[1]*hmat[5]*hmat[9];
   hmati[1] = 1.0/hmat[1];
   hmati[5] = 1.0/hmat[5];
   hmati[9] = 1.0/hmat[9];
  }/*endif*/

/*===============================================================*/
/* Errors */

  if((*deth)==0.0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");    
    printf("The present volume is zero.                 \n");
    printf("If this is not an error in your input data, \n");
    printf("contact technical support                  \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  } /*endif*/

/*---------------------------------------------------------------*/
  } /* gethinv */
/*===============================================================*/



/*===============================================================*/
/*  Determinent of a 3x3 */
/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
double getdeth(double *hmat)
{
  
  double getdeth;
  /* gets det of hmat */
  
  getdeth = (hmat[1] * (hmat[5] * hmat[9] - hmat[8] * hmat[6]) + 
             hmat[4] * (hmat[8] * hmat[3] - hmat[2] * hmat[9]) + 
             hmat[7] * (hmat[2] * hmat[6] - hmat[5] * hmat[3]));
  return getdeth;

} /* getdeth */
/*===============================================================*/


/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
double ddot1(int n,double *a,int astep,double *b,int bstep)
{
  int i,j,k;
  double ddot1;

  ddot1 = 0.;
  for(i=1,j=1,k=1; i<=n ;i++,j+=astep,k+=bstep){
    ddot1 += a[j]*b[k];
  }

  return ddot1;
}
/*===============================================================*/

/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
double ddotBlasWrapper(int n,double *x,int indx,double *y,int indy){
  return ddot_(&n,x,&indx,y,&indy);
}

/*===============================================================*/

/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
void daxpyBlasWrapper(int n,double a,double *x,int indx,double *y,int indy){
  daxpy_(&n,&a,x,&indx,y,&indy);
}

/*===============================================================*/

/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
double dsum1(int n,double *a,int astep)
{
  int i,j;
  double dsum1;

  dsum1 = 0.;
  for(i=1,j=1; i<=n ;i++,j+=astep){
    dsum1 += a[j];
  }

  return dsum1;
}
/*===============================================================*/

/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
double dot(double *dot_1, double *dot_2)
{

    return dot_1[0]*dot_2[0]+dot_1[1]*dot_2[1]+dot_1[2]*dot_2[2];
}
/*==========================================================================*/

/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void cross_product(double *a, double *b, double *result)
{
    result[0] = a[1]*b[2]-a[2]*b[1];
    result[1] = a[2]*b[0]-a[0]*b[2];
    result[2] = a[0]*b[1]-a[1]*b[0];
}
/*==========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
double normalize3d(double *a)
{
  int i;
  double norm = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  a[0] /= norm;
  a[1] /= norm;
  a[2] /= norm;
  return norm;
}
/*==========================================================================*/


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
double dsymvWrapper(char uplo,int n, double alpha,double *A,int lda,double *x,int incx,
		    double beta,double *y,int incy)
{
  dsymv_(&uplo,&n,&alpha,A,&lda,x,&incx,&beta,y,&incy);

}

/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
double dgemvWrapper(char trans,int m,int n, double alpha,double *A,int lda,double *x,int incx,
                    double beta,double *y,int incy)
{
  dgemv_(&trans,&m,&n,&alpha,A,&lda,x,&incx,&beta,y,&incy);

}


/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
double gerf(double x)
{
/*===============================================================*/
/*  Local variables */

  double p=0.3614;
  double e1 = 0.2041422096422003, e2 = 0.1997535956961481;
  double e3 = 0.2213176596405576, e4 = 0.03360430734640255;
  double e5 = 0.4732592578721755, e6 =-0.509078520069735;
  double e7 = 0.6772631491947646, e8 =-0.369912979092217;
  double e9 = 0.06965131976970335;
  double eee,tt,gerf;

/*===============================================================*/
/* Calculate the error function */

   eee    = exp(-x*x);
   tt     = 1.0/(1.0+p*x);
   gerf   = 1.0 - ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                               +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;

/*===============================================================*/
   return gerf;
/*===============================================================*/
}/* end function */
/*===============================================================*/



/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
double gerfc(double x)
{
/*===============================================================*/
/*  Local variables */

  double p=0.3614;
  double e1 = 0.2041422096422003, e2 = 0.1997535956961481;
  double e3 = 0.2213176596405576, e4 = 0.03360430734640255;
  double e5 = 0.4732592578721755, e6 =-0.509078520069735;
  double e7 = 0.6772631491947646, e8 =-0.369912979092217;
  double e9 = 0.06965131976970335;
  double eee,tt,gerf,gerfc;

/*===============================================================*/
/* Calculate the error function */

   eee    = exp(-x*x);
   tt     = 1.0/(1.0+p*x);
   gerf   = 1.0 - ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                               +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
   gerfc = 1.0 - gerf;

/*===============================================================*/
   return gerfc;
/*===============================================================*/
}/* end function */
/*===============================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

double surf_corr(double x)
{
  return 2.0*sqrt(M_PI)*(x*sqrt(M_PI)*gerfc(x) - exp(-x*x));
}

/*===============================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

double dsurf_corr(double x)
{
  return 2.0*M_PI*gerfc(x);
}

/*===============================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

double d2surf_corr(double x)
{
  return -4.0*sqrt(M_PI)*exp(-x*x);
}

/*===============================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

int get_ind(int *index,int *nbin,int ndim){
  int i_dim;
  int index_c = 0;
  for(i_dim=0;i_dim<ndim;i_dim++){
     index_c += index[i_dim]*nbin[i_dim];
  }
  return index_c;
}
/*===============================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void get_ind_rev(int *index,int *n_bin,int i,int ndim){
  int div,res;
  int i_dim;
  res = i;
  for(i_dim=0;i_dim<ndim;i_dim++){
     index[i_dim] = (int)(res/n_bin[i_dim]);
     res = res%n_bin[i_dim];
  }
}
/*===============================================================*/

/***************************************************************************/
/* We will put some type of Fermi functions here.                          */
/***************************************************************************/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double fermiExpReal(double x,double mu,double beta){
  double xre = (double)x;
  return 1.0/(1.0+exp(beta*(x-mu)));

}
/*===============================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double complex fermiExpComplex(double complex x,double mu,double beta){
  return 1.0/(1.0+cexp(beta*(x-mu)));

}
/*===============================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double fermiErfcReal(double x,double mu,double beta){
  return 0.5*erfc(beta*(x-mu));
}
/*===============================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double gaussianReal(double x,double mu,double beta){
  return 0.3989422804014327*beta*exp(-beta*beta*(x-mu)*(x-mu));
}
/*===============================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double j0(double x){
  return sin(x)/x;
}
/*===============================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double j1(double x){
  return sin(x)/(x*x)-cos(x)/x;
}
/*===============================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double j2(double x){
  return (3.0/(x*x)-1.0)*sin(x)/x-3.0*cos(x)/(x*x);
}
/*===============================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double dj0(double x){
  return cos(x)/x-sin(x)/(x*x);
}
/*===============================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double dj1(double x){
  double x2 = x*x;
  double x3 = x2*x;
  return 2.0*cos(x)/x2+(1.0/x-2.0/x3)*sin(x);
}
/*===============================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double dj2(double x){
  double x2 = x*x;
  double x3 = x2*x;
  double x4 = x2*x2;
  return (-9.0/x4+4.0/x2)*sin(x)+(9.0/x3-1.0/x)*cos(x);
}
/*===============================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double dsysvWrapper(double *A,double *b,int n){
  char uplo = 'U';
  int nrhs = 1;
  int lda = n;
  int ldb = n;
  int i;
  int lwork = 64*n;
  int *ipiv = (int*)cmalloc(n*sizeof(int));
  double *work = (double*)cmalloc(lwork*sizeof(double));
  int info;

  double *ACpy = (double*)calloc(n*n,sizeof(double));
  //double *ACpy = (double*)cmalloc(n*n*sizeof(double));
  memcpy(ACpy,A,n*n*sizeof(double));
  //for(i=0;i<n*n;i++)ACpy[i] = A[i];
  dsysv_(&uplo,&n,&nrhs,ACpy,&lda,ipiv,b,&ldb,work,&lwork,&info);

  if(info>0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Error in solving the linear equations!\n");
    printf("The %i'th diag element is singular!\n",info);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(0);
  } 
  if(info<0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Error in solving the linear equations!\n");
    printf("The %i'th parameter had an illegal value!\n",info);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(0);
  } 
  
  free(ipiv);
  free(work);
  free(ACpy);
}
/*===============================================================*/

