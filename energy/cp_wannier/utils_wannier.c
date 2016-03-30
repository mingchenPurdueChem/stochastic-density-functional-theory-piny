/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: utils_wannier                                */
/*                                                                          */
/* This subprogram performs on the fly Wannier cp dynamics                  */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_wannier_cpcon_local.h"
#include "../proto_defs/proto_wannier_cpcon_entry.h"

#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0
#define its_max 200


static int minarg1, minarg2;
#define IMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1 < minarg2) ? (minarg1) : (minarg2))

static int maxarg1, maxarg2;
#define IMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1 > maxarg2) ? (maxarg1) : (maxarg2))


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comp_miller_weight(GENERAL_DATA *general_data, double *W)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/
/*            Local variable declarations                                    */
/*===========================================================================*/

#include "../typ_defs/typ_mask.h"

  int i,j,k;
  double *hmat_cp = general_data->cell.hmat_cp;
  double hvec1[3],hvec2[3],hvec3[3];

/*---------------------------------------------------------------------------*/
/* Assign cell vector                                                        */

  hvec1[0] = BOHR*hmat_cp[1];
  hvec1[1] = BOHR*hmat_cp[4];
  hvec1[2] = BOHR*hmat_cp[7];

  hvec2[0] = BOHR*hmat_cp[2];
  hvec2[1] = BOHR*hmat_cp[5];
  hvec2[2] = BOHR*hmat_cp[8];

  hvec3[0] = BOHR*hmat_cp[3];
  hvec3[1] = BOHR*hmat_cp[6];
  hvec3[2] = BOHR*hmat_cp[9];

/*----------------------------------------------------------------------------*/
/* Compute the Weight for I=1,,6                                              */
/*----------------------------------------------------------------------------*/
/* g_{1,2} = h(1)*h(2), g_{2,3} = h(2)*h(3), etc.. where h(i) is the ith      */
/* column of cell matrix                                                      */

  for (i=1;i<=6;i++) {
    W[i]=0.0;
  }

  for (j=0;j<3;j++){
    W[1]+=(hvec1[j]*(hvec1[j]-hvec2[j]-hvec3[j]));
    W[2]+=(hvec2[j]*(hvec2[j]-hvec1[j]-hvec3[j]));
    W[3]+=(hvec3[j]*(hvec3[j]-hvec1[j]-hvec2[j]));
    W[4]+=(hvec1[j]*hvec2[j]);
    W[5]+=(hvec1[j]*hvec3[j]);
    W[6]+=(hvec2[j]*hvec3[j]);
  }

/*============================================================================ */
}/* end routine */
/*============================================================================ */

void comp_tensor_z_dvr(GENERAL_DATA *general_data, CP *cp,
                            double ***Z_real, double ***Z_imag)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/

#include "../typ_defs/typ_mask.h"

  int myid           = cp->communicate.myid_state;
  int num_proc       = cp->communicate.np_states;
  MPI_Comm comm      = cp->communicate.comm_states;

  int nstate_up=cp->cpcoeffs_info.nstate_up;
  int nstate_dn=cp->cpcoeffs_info.nstate_dn;
  int nstate=nstate_up;

  double *hmat             = general_data->cell.hmat;
  double *hmati            = general_data->cell.hmati;
  double *hmat_cp          = general_data->cell.hmat_cp;
  double vol               = general_data->cell.vol_cp;

  int nkf1                 = cp->cp_para_fft_pkg3d_sm.nkf1;
  int nkf2                 = cp->cp_para_fft_pkg3d_sm.nkf2;
  int nkf3                 = cp->cp_para_fft_pkg3d_sm.nkf3;
  int skc_fft_ka_proc      = cp->cp_para_fft_pkg3d_sm.skc_fft_ka_proc;
  int ekc_fft_ka_proc      = cp->cp_para_fft_pkg3d_sm.ekc_fft_ka_proc;
  int skb_fft_ka_proc      = cp->cp_para_fft_pkg3d_sm.skb_fft_ka_proc;
  int ekb_fft_ka_proc      = cp->cp_para_fft_pkg3d_sm.ekb_fft_ka_proc;
  int nfft_proc            = cp->cp_para_fft_pkg3d_sm.nfft_proc;
  int nfft2                = cp->cp_para_fft_pkg3d_sm.nfft/2;

  double *dvrc             = cp->cpcoeffs_pos_dvr[1].dvrc_up;
  double dfact             = vol/((double)nfft2);

  /* This is a waste of memory, but scratch arrays are already allocated */
  double *dvrc_real        = cp->cpscr.cpscr_wave.cre_up;
  double *dvrc_imag        = cp->cpscr.cpscr_wave.cim_up;
  double *dvrc_real_tmp    = cp->cpscr.cpscr_wave.zfft;
  double *dvrc_imag_tmp    = cp->cpscr.cpscr_wave.zfft_tmp;

/*        Local Variable declarations             */

  int i,j,k,ka,kb,kc,kb_str,kb_end,is,lda,ldb,ldc,nrow,ncol,nvec,ngrid,index;
  double sa,sb,sc,da,db,dc,aka,akb,akc;
  double grid_x,grid_y,grid_z;
  double dxx,dyy,dzz,exp_gr_r,exp_gr_i,gr;
  double pi,tpi;
  double G_I[3];
  int g_hat_I[6][3];
  int itransp = 0;
  int inormal = 1;
  double beta=1.0;
  double alpha=1.0;

/*============================================================================*/
/* 0) constants for the grid set up                                           */

  pi = M_PI;
  tpi = 2.0*pi;

  dxx = hmat[1]/(2.0*(double) nkf1);
  dyy = hmat[5]/(2.0*(double) nkf2);
  dzz = hmat[9]/(2.0*(double) nkf3);

  da = 1.0/((double) nkf1);
  db = 1.0/((double) nkf2);
  dc = 1.0/((double) nkf3);

/*=============================================================================*/
/* I) Miller index                                                           */

  g_hat_I[0][0] = 1;
  g_hat_I[0][1] = 0;
  g_hat_I[0][2] = 0;
  g_hat_I[1][0] = 0;
  g_hat_I[1][1] = 1;
  g_hat_I[1][2] = 0;
  g_hat_I[2][0] = 0;
  g_hat_I[2][1] = 0;
  g_hat_I[2][2] = 1;
  g_hat_I[3][0] = 1;
  g_hat_I[3][1] = 1;
  g_hat_I[3][2] = 0;
  g_hat_I[4][0] = 1;
  g_hat_I[4][1] = 0;
  g_hat_I[4][2] = 1;
  g_hat_I[5][0] = 0;
  g_hat_I[5][1] = 1;
  g_hat_I[5][2] = 1;

  ngrid = (ekc_fft_ka_proc - skc_fft_ka_proc - 1)*nkf1*nkf2
           + (nkf2 - skb_fft_ka_proc + 1)*nkf1
           + (ekb_fft_ka_proc)*nkf1;

/*===========================================================================*/
/* II) Compute Z-tensor for each I             */

   for (i=0;i<3;i++){  /* There are total six miller indices
                          but assume orthorombic cell*/
     for (j=0;j<3;j++){
       G_I[j]=0.0;
     }

     for (j=1;j<=ngrid*nstate;j++){
       dvrc_real_tmp[j]=0.0;
       dvrc_imag_tmp[j]=0.0;
       dvrc_real[j]=0.0;
       dvrc_imag[j]=0.0;
     }

     aka = (double)(g_hat_I[i][0]);
     akb = (double)(g_hat_I[i][1]);
     akc = (double)(g_hat_I[i][2]);

     G_I[0] = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
     G_I[1] = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
     G_I[2] = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;

     /* compute \psi*e(igr) */
               
     index=0;
     for(kc=skc_fft_ka_proc;kc<=ekc_fft_ka_proc;kc++){
       kb_str = (kc==skc_fft_ka_proc ? skb_fft_ka_proc : 1);
       kb_end = (kc==ekc_fft_ka_proc ? ekb_fft_ka_proc : nkf2);
       sc = dc*((double)(kc-1)) - 0.5;

       for(kb=kb_str;kb<=kb_end;kb++){
         sb = db*((double)(kb-1)) - 0.5;
         for(ka=1;ka<=nkf1;ka++){
           sa = da*((double)(ka-1)) - 0.5;

           grid_x  = sa*hmat[1]+sb*hmat[4]+sc*hmat[7];
           grid_y  = sa*hmat[2]+sb*hmat[5]+sc*hmat[8];
           grid_z  = sa*hmat[3]+sb*hmat[6]+sc*hmat[9];

           grid_x += dxx;
           grid_y += dyy;
           grid_z += dzz;

           gr = grid_x*G_I[0]+grid_y*G_I[1]+grid_z*G_I[2];
           exp_gr_r = cos(gr);
           exp_gr_i = sin(gr);

           index += 1;
           for(is=1;is<=nstate;is++){
             dvrc_real[(is-1)*ngrid+index]=dvrc[(is-1)*ngrid+index]*exp_gr_r;
             dvrc_imag[(is-1)*ngrid+index]=dvrc[(is-1)*ngrid+index]*exp_gr_i;
           }
         }/*endfor ka*/
       }/*endfor kb*/
     }/*endfor kc*/

    /*compute z_tensor for ith miller index */

    lda=ngrid;
    ldb=ngrid;
    ldc=nstate;
    nrow = nstate;
    ncol = nstate;
    nvec = ngrid;

    /* GEN_MATMUL(A,LDA,T/N,B,LDB,T/N,C,LDC,M,N,K,alpha,beta);
       A=MxK, B=KxN, C=MxN, C=alpha*A(T/N)*B(T/N)+beta*C) */

    GEN_MATMUL(&(dvrc[1]),&lda,&itransp,&(dvrc_real[1]),&ldb,&inormal,
               &(dvrc_real_tmp[1]),&ldc,&nrow,&ncol,&nvec,&alpha,&beta);

    GEN_MATMUL(&(dvrc[1]),&lda,&itransp,&(dvrc_imag[1]),&ldb,&inormal,
               &(dvrc_imag_tmp[1]),&ldc,&nrow,&ncol,&nvec,&alpha,&beta);

    if(num_proc > 1){

      Barrier(comm);
      Allreduce(&(dvrc_real_tmp[1]),&(dvrc_real[1]),(nstate*nstate),MPI_DOUBLE,
                MPI_SUM,0,comm);

      Barrier(comm);
      Allreduce(&(dvrc_imag_tmp[1]),&(dvrc_imag[1]),(nstate*nstate),MPI_DOUBLE,
                MPI_SUM,0,comm);

      for(k=1;k<=nstate;k++){
        for(j=1;j<=nstate;j++){
          Z_real[i+1][j][k]=dvrc_real[(k-1)*nstate+j]*dfact/2.0;
          Z_imag[i+1][j][k]=dvrc_imag[(k-1)*nstate+j]*dfact/2.0;
        }
      }
    }else{
      for(k=1;k<=nstate;k++){
        for(j=1;j<=nstate;j++){
          Z_real[i+1][j][k]=dvrc_real_tmp[(k-1)*nstate+j]*dfact/2.0;
          Z_imag[i+1][j][k]=dvrc_imag_tmp[(k-1)*nstate+j]*dfact/2.0;
        }
      }
    }/*endif num_proc*/

  }/*end for i*/

/*===========================================================================*/
} /*end routine */
/*============================================================================*/

void comp_omega(GENERAL_DATA *general_data,CP *cp, double *X,
                double ***Z_real, double ***Z_imag, double *W,double *fp, int NDIM)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/
/*            Local variable declarations                                    */
/*===========================================================================*/

#include "../typ_defs/typ_mask.h"

  int nstate_up=cp->cpcoeffs_info.nstate_up;
  int nstate_dn=cp->cpcoeffs_info.nstate_dn;
  int nstate=nstate_up;
  int choice_diag = cp->cp_wannier.wan_diag_typ;
  int choice_func = cp->cp_wannier.wan_func_typ;


/* declartion of local variable */

  int I=6;
  int i,j,k;
  int flag,incx,incy,mm,nn,kk,lda,ldb,ldc;

  double norm,prod1,prod2,pi,fpisq,Omega_int,fp_now;
  double alpha,beta;
  char transa,transb;

  double **A        = cp->cpscr.cpscr_wannier.A;
  double **R_real   = cp->cpscr.cpscr_wannier.R_real;
  double **R_imag   = cp->cpscr.cpscr_wannier.R_imag;
  double **U_real   = cp->cpscr.cpscr_wannier.U_real;
  double **U_imag   = cp->cpscr.cpscr_wannier.U_imag;
  double **U_tmp1   = cp->cpscr.cpscr_wannier.U_tmp1;
  double **U_tmp2   = cp->cpscr.cpscr_wannier.U_tmp2;
  double *D         = cp->cpscr.cpscr_wannier.D;


/* ============================================================================ */
/*  0) Constant */

  pi= M_PI;
  fpisq = 4.0*pi*pi;

/*===============================================================================*/
/* I) Reassign X(1D) to A(2D)                                                   */

  flag=0;
  for(i=1;i<=nstate;i++){
    for(j=1;j<=nstate;j++){
      A[i][j]=X[(i-1)*nstate+j];
      if(A[i][j]!=0.0) flag=2;  /*Is this OK? comparison of real number!!*/
    }
  }

  if(flag==2){
    for(i=1;i<=nstate;i++){
       for (j=1;j<=nstate;j++){
         if(fabs(A[i][j])>=0.000000001) flag=1;
       }
    }
    flag=1; /*flag=1 (A is big) anyway just for now. */
  }

/*===============================================================================*/
/* II) diagonalize A matrix                                                      */

  if(flag==0){ /* A=0 */
    for(i=1;i<=nstate;i++){
      D[i]=0.0;
      for(j=1;j<=nstate;j++){
        R_imag[i][j]=0.0;
        U_imag[i][j]=0.0;
        if(i==j){
          R_real[i][j]=1.0;
          U_real[i][j]=1.0;
        }else{
          R_real[i][j]=0.0;
          U_real[i][j]=0.0;
        }
      }/*endfor j*/
    }/*endfor i*/
  }else if(flag==1){ /* A is large: direct diagonalization*/

    diagonalize_asymm_mat(A,nstate,D,U_real,U_imag,R_real,R_imag,choice_diag);

  }else{ /* A is small: Taylor expansion*/

    transa='N'; transb='N'; mm=nstate; nn=nstate; kk=nstate; lda=nstate;
    ldb=nstate; ldc=nstate; alpha=1.0; beta=0.0;

    /* first calculate the R_real=A*A, R_imag= 1/6 * A + 1/2 * I, U_real=A+I */

    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&A[1][1],&lda,&A[1][1],
           &ldb,&beta,&R_real[1][1],&ldc);

    for(i=1;i<=nstate;i++){
      for(j=1;j<=nstate;j++){
        if(i==j){
          R_imag[i][j]=0.5+A[i][j]/6.0;
          U_real[i][j]=1.0+A[i][j];
        }else{
          R_imag[i][j]=A[i][j]/6.0;
          U_real[i][j]=A[i][j];
        }
      }/*end for j*/
    }/*end for i*/

    /* Now U_real=U_real+R_real*R_imag = I+A+1/2*A*A+1/6*A*A*A : simple Taylor exp */
    beta=1.0;
    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&R_imag[1][1],&lda,&R_real[1][1],
           &ldb,&beta,&U_real[1][1],&ldc);

    /* Verify that U_real is indeed a unitary matrix */
    beta=0.0; transb='T';
    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&U_real[1][1],&lda,&U_real[1][1],
           &ldb,&beta,&R_real[1][1],&ldc);
    for (i=1;i<=nstate;i++) {
      R_real[i][i]-=1.0;
    }
    incx=incy=1;
    norm=sqrt(DDOT(&NDIM,&R_real[1][1],&incx,&R_real[1][1],&incy));
    printf("Using Taylor expansion the norm in omega is %g\n",norm);

  } /*end if*/

/*============================================================================*/
/* III) Calculate Omega                                                       */

  fp_now=0.0;

  for(i=1;i<=I;i++){
    if(W[i]!=0){
      /* first calculate the matrix Z*U. Do that by multiplying U*Z using a
         fortran routine. Next, calculate the inner product between the columns
         of (U*Z) and U   */

      transa='N'; transb='N'; mm=nstate; nn=kk=nstate; lda=nstate;
      ldb=nstate; ldc=nstate; alpha=1.0; beta=0.0; incx=nstate; incy=nstate;

      DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&U_real[1][1],&lda,&Z_real[i][1][1],
             &ldb,&beta,&U_tmp1[1][1],&ldc);

      DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&U_real[1][1],&lda,&Z_imag[i][1][1],
             &ldb,&beta,&U_tmp2[1][1],&ldc);

      Omega_int=0.0;
      for(j=1;j<=nstate;j++){
        prod1=DDOT(&lda,&U_real[1][j],&incx,&U_tmp1[1][j],&incy);
        prod2=DDOT(&lda,&U_real[1][j],&incx,&U_tmp2[1][j],&incy);

        switch(choice_func){
          case 1:
            Omega_int+=log((prod1*prod1)+(prod2*prod2));
            /* Omega_int+=-log((prod1*prod1)+(prod2*prod2));*/
          break;

          case 2:
            Omega_int+=(1.0-prod1*prod1-prod2*prod2);
          break;

          case 3:
            Omega_int+=sqrt((prod1*prod1)+(prod2*prod2));
            /*Omega_int+=1-sqrt((prod1*prod1)+(prod2*prod2));*/
          break;

          default:
            printf("@@@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@\n");
            printf(" Wrong choice for the functional should have been rejected\n");
            printf(" previously\n");
            printf("@@@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@\n");
            fflush(stdout);
            Finalize(); exit(1);
          break;
        }/*end switch*/
      }/*end for*/

    fp_now += (W[i]/fpisq*Omega_int);
    }/*end if*/
  }/*end for i*/

  *fp=fp_now;


/* ========================================================================== */
}/*end routine                                                               */
/* ========================================================================== */

void comp_domega(GENERAL_DATA *general_data,CP *cp, double *X,
                double *df, double ***Z_real, double ***Z_imag, double *W,int NDIM)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/
/*            Local variable declarations                                    */
/*===========================================================================*/

#include "../typ_defs/typ_mask.h"

  int nstate_up=cp->cpcoeffs_info.nstate_up;
  int nstate_dn=cp->cpcoeffs_info.nstate_dn;
  int nstate=nstate_up;
  int choice_diag = cp->cp_wannier.wan_diag_typ;
  int choice_func = cp->cp_wannier.wan_func_typ;

/*local variables*/

  int i,j,k,l,m,n,s,t,mm,nn,kk;
  int ii,jj,iii,ioff,ioff2,istarti,istartj,istartn,count,flag;
  int lda,ldb,ldc,incx,incy;
  int I=6;
  static int first_time=1;
  char transa,transb;
  double alpha,beta;
  double prod1,prod2,tmp1_real,tmp1_imag,tmp2_real,tmp2_imag,tmp,tmp1,tmp2;
  double MaxImagPart;

  double pi,fpisq;

  double **A        = cp->cpscr.cpscr_wannier.A;
  double **R_real   = cp->cpscr.cpscr_wannier.R_real;
  double **R_imag   = cp->cpscr.cpscr_wannier.R_imag;
  double **U_real   = cp->cpscr.cpscr_wannier.U_real;
  double **U_imag   = cp->cpscr.cpscr_wannier.U_imag;
  double **Tmp1_real   = cp->cpscr.cpscr_wannier.U_tmp1;
  double **Tmp1_imag   = cp->cpscr.cpscr_wannier.U_tmp2;
  double **Tmp2_real   = cp->cpscr.cpscr_wannier.U_tmp3;
  double **Tmp2_imag   = cp->cpscr.cpscr_wannier.U_tmp4;
  double *D         = cp->cpscr.cpscr_wannier.D;
  double **Bt_real   = cp->cpscr.cpscr_wannier.Bt_real;
  double **Bt_imag   = cp->cpscr.cpscr_wannier.Bt_imag;
  double **M_real   = cp->cpscr.cpscr_wannier.M_real;
  double *real   = cp->cpscr.cpscr_wannier.real;
  double *imag   = cp->cpscr.cpscr_wannier.imag;
  double *norm   = cp->cpscr.cpscr_wannier.norm;

/*============================================================================*/
/* I) reassign matrix A (1D->2D)             */

  pi=M_PI;
  fpisq=4.0*pi*pi;

  for(i=1;i<=nstate;i++){
    for (j=1;j<=nstate;j++){
      M_real[i][j]=0.0;
      A[i][j]=X[(i-1)*nstate+j];
    }
  }

/*==========================================================================*/
/* II) diagonalize matrix A                                                 */
/*     The matrix R will contain the eigenvectors of A as its columns.      */
/*     The vector D will contain the eigenvalues of A.                      */
/*     If A=0, then R=1, D=0.                                               */
/*     The matrix A is multiplied by sqrt(-1) =>hermitian matrix            */
/*     Also compute the matrix U=exp(A)=R * exp(D) * R and the matrix Bt    */
/*==========================================================================*/

  flag=0;  /*A=0*/
  for(i=1;i<=nstate;i++){
    for (j=1;j<=nstate;j++){
      if(A[i][j]!=0.0) flag=1; /*A=/=0*/
    }
  }

  if(flag==0){ /* A=0 */
    for(i=1;i<=nstate;i++){
      D[i]=0.0;
      for(j=1;j<=nstate;j++){
        R_imag[i][j]=0.0;
        U_imag[i][j]=0.0;
        if(i==j){
          R_real[i][j]=1.0;
          U_real[i][j]=1.0;
        }else{
          R_real[i][j]=0.0;
          U_real[i][j]=0.0;
        }
      }/*endfor j*/
    }/*endfor i*/
  }else{ /*flag==1 A is large: direct diagonalization*/

    diagonalize_asymm_mat(A,nstate,D,U_real,U_imag,R_real,R_imag,choice_diag);

  } /*no Taylor expansion option this time*/

/*===========================================================================*/
/* III) Compute matrix Bt                                                    */

  if(flag==0){
    for(i=1;i<=nstate;i++){
      for(j=i;j<=nstate;j++){
        Bt_imag[i][j]=0.0;
        Bt_real[i][j]=1.0;
        Bt_imag[j][i]=Bt_imag[i][j];
        Bt_real[j][i]=Bt_real[i][j];
      } /* end for j */
    } /* end for i */
  }else{ /* flag==1 */
    /*Note that the actual eigenvalues of A seems to be lambda_j=-i*D[j]*/
    for (i=1;i<=nstate;i++){
      for (j=1;j<=nstate;j++){
        if (D[i]==D[j]){
          Bt_real[j][i]=cos(D[i]);
          Bt_imag[j][i]=-sin(D[i]);
        }else{ /* if D[i]!=D[j] */
          Bt_real[j][i]=(sin(D[i])-sin(D[j]))/(D[i]-D[j]);
          Bt_imag[j][i]=(cos(D[i])-cos(D[j]))/(D[i]-D[j]);
        }
      } /*endfor j*/
    } /*end for i */
  } /* end if flag==0 */

/*=============================================================================*/
/* IV) Compute matrix M                                                        */

  for (m=1;m<=I;m++) {
    if (W[m]!=0.0){
      /* First, compute U*Z */
      /* Remember U is real at Gamma point. It should be the same for DVR.
         Also note that M_{st} is proportional to (UZ)_{st} */

      transa='N'; transb='N'; mm=nstate; nn=nstate; kk=nstate; lda=nstate;
      ldb=nstate; ldc=nstate; alpha=1.0; beta=0.0; incx=nstate; incy=nstate;

      DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&U_real[1][1],&lda,
             &Z_real[m][1][1],&ldb,&beta,&Tmp1_real[1][1],&ldc);
      DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&U_real[1][1],&lda,
             &Z_imag[m][1][1],&ldb,&beta,&Tmp1_imag[1][1],&ldc);

      for(j=1;j<=nstate;j++){
        /* compute Z_I,n */
        real[j]=DDOT(&lda,&U_real[1][j],&incx,&Tmp1_real[1][j],&incy);
        imag[j]=DDOT(&lda,&U_real[1][j],&incx,&Tmp1_imag[1][j],&incy);

        /* compute f'(|Z_I,n|^2) */
        switch(choice_func){
          case 1:
            norm[j]=1.0/(real[j]*real[j]+imag[j]*imag[j]);
          break;
          case 2:
            norm[j]=1.0;
          break;
          case 3:
            norm[j]=0.5/sqrt((real[j]*real[j]+imag[j]*imag[j]));
          break;
          default:
            printf("@@@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@\n");
            printf(" Mistake in omega.c. Wrong choice for the functional\n");
            printf(" should have been rejected previously\n");
            printf("@@@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@\n");
            fflush(stdout);
            Finalize(); exit(1);
          break;
        }/*end switch*/
      }/*end for j*/
      for(j=1;j<=nstate;j++){
        real[j]*=(4.0*norm[j]);
        imag[j]*=(4.0*norm[j]);
      }
      /*sum over n is not really a sum due to the delta_tn*/
      for (j=1;j<=nstate;j++){
        DSCAL(&lda,&real[j],&Tmp1_real[1][j],&incx);
        DSCAL(&lda,&imag[j],&Tmp1_imag[1][j],&incx);
      }

      /* Update: add M_st for this I to M_real. Note that M matrix is real*/
      lda=nstate*nstate; alpha=W[m]/fpisq; incx=1; incy=1;
      DAXPY(&lda,&alpha,&Tmp1_real[1][1],&incx,&M_real[1][1],&incy);
      DAXPY(&lda,&alpha,&Tmp1_imag[1][1],&incx,&M_real[1][1],&incy);

    }/*end if W*/
  }/*end m*/

/*===========================================================================*/
/* V) Compute the derivative: dOmega/dA_ij                                   */
/*  i.e. Calculate R{(R(*)M(t)R),B(t)}R(*)                                   */

/* 1. First compute R(*)M(t)=(R_r(t)-iR_i(t))M(t)(R_r+iR_i) */

  /* R_r(t)M(t)-iR_i(t)M(t) = Tmp1_real+iTmp1_imag */
  transa='T'; transb='T'; mm=nn=kk=lda=ldb=ldc=nstate; alpha=1.0; beta=0.0;
  DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&M_real[1][1],&lda,&R_real[1][1],
         &ldb,&beta,&Tmp1_real[1][1],&ldc);
  transa='T'; transb='T'; mm=nn=kk=lda=ldb=ldc=nstate; alpha=-1.0; beta=0.0;
  DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&M_real[1][1],&lda,&R_imag[1][1],
         &ldb,&beta,&Tmp1_imag[1][1],&ldc);

  /* R(*)M(t)= (Tmp1_r+iTmp1_i)(R_r+iR_i)
             = (Tmp1_r R_r - Tmp1_i R_i)+i(Tmp1_r R_i + Tmp1_i R_r)
             =  Tmp2_r + iTmp2_i */

  transa='N'; transb='N'; mm=nn=kk=lda=ldb=ldc=nstate; alpha=1.0; beta=0.0;
  DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&R_real[1][1],&lda,&Tmp1_real[1][1],
         &ldb,&beta,&Tmp2_real[1][1],&ldc);
  transa='N'; transb='N'; mm=nn=kk=lda=ldb=ldc=nstate; alpha=-1.0; beta=1.0;
  DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&R_imag[1][1],&lda,&Tmp1_imag[1][1],
         &ldb,&beta,&Tmp2_real[1][1],&ldc);
  transa='N'; transb='N'; mm=nn=kk=lda=ldb=ldc=nstate; alpha=1.0; beta=0.0;
  DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&R_real[1][1],&lda,&Tmp1_imag[1][1],
         &ldb,&beta,&Tmp2_imag[1][1],&ldc);
  transa='N'; transb='N'; mm=nn=kk=lda=ldb=ldc=nstate; alpha=1.0; beta=1.0;
  DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&R_imag[1][1],&lda,&Tmp1_real[1][1],
         &ldb,&beta,&Tmp2_imag[1][1],&ldc);

/* 2. Component-wise multiplication { } */

  component_wise_matrix_multiply(nstate,1,1,Tmp2_real,Tmp2_imag,Bt_real,Bt_imag,
                                 Tmp1_real,Tmp1_imag);

/* 3. R Tmp_1 R(*) */

  /* R x Tmp_1 = (R_r+iR_i)(Tmp1_r + iTmp1_i)
               = [R_r Tmp1_r - R_i Tmp1_i]+i[R_r Tmp1_i+R_i Tmp1_r]
               = Tmp2_r + i Tmp2_i */

  transa='N'; transb='N'; mm=nn=kk=lda=ldb=ldc=nstate; alpha=1.0; beta=0.0;
  DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&Tmp1_real[1][1],&lda,&R_real[1][1],
         &ldb,&beta,&Tmp2_real[1][1],&ldc);
  transa='N'; transb='N'; mm=nn=kk=lda=ldb=ldc=nstate; alpha=-1.0; beta=1.0;
  DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&Tmp1_imag[1][1],&lda,&R_imag[1][1],
         &ldb,&beta,&Tmp2_real[1][1],&ldc);

  transa='N'; transb='N'; mm=nn=kk=lda=ldb=ldc=nstate; alpha=1.0; beta=0.0;
  DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&Tmp1_real[1][1],&lda,&R_imag[1][1],
         &ldb,&beta,&Tmp2_imag[1][1],&ldc);
  transa='N'; transb='N'; mm=nn=kk=lda=ldb=ldc=nstate; alpha=1.0; beta=1.0;
  DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&Tmp1_imag[1][1],&lda,&R_real[1][1],
         &ldb,&beta,&Tmp2_imag[1][1],&ldc);

  /* Tmp2 x R(*) = (Tmp2_r+iTmp2_i)(R_r(t) -iR_i(t))
                 = [Tmp2_rR_r(t)+Tmp2_iR_i(t)] +i[Tmp2_iR_r(t)-Tmp2_rR_i(t)] */

  transa='T'; transb='N'; mm=nn=kk=lda=ldb=ldc=nstate; alpha=1.0; beta=0.0;
  DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&R_real[1][1],&lda,&Tmp2_real[1][1],
         &ldb,&beta,&Tmp1_real[1][1],&ldc);
  transa='T'; transb='N'; mm=nn=kk=lda=ldb=ldc=nstate; alpha=1.0; beta=1.0;
  DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&R_imag[1][1],&lda,&Tmp2_imag[1][1],
         &ldb,&beta,&Tmp1_real[1][1],&ldc);
  transa='T'; transb='N'; mm=nn=kk=lda=ldb=ldc=nstate; alpha=1.0; beta=0.0;
  DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&R_real[1][1],&lda,&Tmp2_imag[1][1],
         &ldb,&beta,&Tmp1_imag[1][1],&ldc);
  transa='T'; transb='N'; mm=nn=kk=lda=ldb=ldc=nstate; alpha=-1.0; beta=1.0;
  DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&R_imag[1][1],&lda,&Tmp2_real[1][1],
         &ldb,&beta,&Tmp1_imag[1][1],&ldc);

/* 4. Finish up the calculation */

  MaxImagPart=0.0;
  for (i=1;i<=nstate;i++){
    for (j=1;j<=nstate;j++){
      df[(i-1)*nstate+j]=-(Tmp1_real[j][i]-Tmp1_real[i][j]);
      if (MaxImagPart < fabs(Tmp1_imag[j][i]-Tmp1_imag[i][j])){
         MaxImagPart=fabs(Tmp1_imag[j][i]-Tmp1_imag[i][j]);
      }
    }
  }
  if (MaxImagPart>=0.0000001){
    printf("Warning: The imaginary part of dOmega is larger than 0.0000001 \n");
    printf("MaxImagPart=%- 15.10e\n",MaxImagPart);
    fflush(NULL);
  }

  first_time=0;

/* ========================================================================== */
}/*end routine                                                               */
/* ========================================================================== */

void comp_domega_zero(int nstate_up, int nstate_dn, int choice_func, double *df,
                      double ***Z_real, double ***Z_imag, double *W)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/
/*            Local variable declarations                                    */
/*===========================================================================*/

#include "../typ_defs/typ_mask.h"

  int nstate=nstate_up;

  /*local variables*/
  int i,j,k;
  int I=6;
  double pi, fpisq, dOmega;

/*============================================================================*/
/* I) Initialization                                                          */

  pi= M_PI;
  fpisq = 4.0*pi*pi;

  for (i=1;i<=nstate;i++){
    for (j=1;j<=nstate;j++){
      df[(i-1)*nstate+j]=0.0;
    }
  }
/*============================================================================*/
/* II) Compute dOmega in the limit of A=0                                     */
/* See Eq. (B1) of JCP paper                                                  */

  for (k=1;k<=I;k++) {
    if (W[k]!=0.0) {
    for (i=1;i<=nstate;i++) {
      for (j=1;j<=nstate;j++) {
        switch(choice_func){
          case 1:
            dOmega=(Z_real[k][i][j]*Z_real[k][j][j]+Z_imag[k][i][j]*Z_imag[k][j][j])
                  /(Z_real[k][j][j]*Z_real[k][j][j]+Z_imag[k][j][j]*Z_imag[k][j][j])
                 - (Z_real[k][i][j]*Z_real[k][i][i]+Z_imag[k][i][j]*Z_imag[k][i][i])
                  /(Z_real[k][i][i]*Z_real[k][i][i]+Z_imag[k][i][i]*Z_imag[k][i][i]);
          break;

          case 2:
            if (i==j) {
              dOmega=0.0;
            }else{
              dOmega=(Z_real[k][i][j]*Z_real[k][j][j]+Z_imag[k][i][j]*Z_imag[k][j][j])
                    -(Z_real[k][i][j]*Z_real[k][i][i]+Z_imag[k][i][j]*Z_imag[k][i][i]);
            }
          break;

          case 3:
            dOmega=(Z_real[k][i][j]*Z_real[k][j][j]+Z_imag[k][i][j]*Z_imag[k][j][j])
          /(2.0*sqrt(Z_real[k][j][j]*Z_real[k][j][j]+Z_imag[k][j][j]*Z_imag[k][j][j]))
                 - (Z_real[k][i][j]*Z_real[k][i][i]+Z_imag[k][i][j]*Z_imag[k][i][i])
          /(2.0*sqrt(Z_real[k][i][i]*Z_real[k][i][i]+Z_imag[k][i][i]*Z_imag[k][i][i]));
          break;

          default:
            printf("@@@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@\n");
            printf(" Wrong choice for the functional should have been rejected\n");
            printf(" previously\n");
            printf("@@@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@\n");
            fflush(stdout);
            Finalize(); exit(1);
          break;
        }

        df[(i-1)*nstate+j]+=(-4.0*W[k]/fpisq*dOmega);

      }
    }
    } /*end if*/
  }

/* ========================================================================== */
}/*end routine                                                               */
/* ========================================================================== */

void comp_ddomega_zero(int nstate_up, int nstate_dn, int choice_func, double *ddf,
                      double ***Z_real, double ***Z_imag, double *W)


/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/
/*            Local variable declarations                                    */
/*===========================================================================*/


#include "../typ_defs/typ_mask.h"

  int nstate=nstate_up;
  int i,j,k;
  int I=6;

  double pi, fpisq, ddOmega;
  double real_ii, real_ij, real_ji, real_jj;
  double imag_ii, imag_ij, imag_ji, imag_jj;

/*============================================================================*/
/* I) Initialization                                                          */

  pi= M_PI;
  fpisq = 4.0*pi*pi;

  for (i=1;i<=nstate;i++){
    for (j=1;j<=nstate;j++){
      ddf[(i-1)*nstate+j]=0.0;
    }
  }

/*============================================================================*/
/* II) Compute ddOmega in the limit of A=0                                    */
/* The expression in Berghold' paper is WRONG!                                */

  for (k=1;k<=I;k++) {
    if (W[k]!=0.0) {
    for (i=1;i<=nstate;i++) {
      for (j=1;j<=nstate;j++) {
        real_ii=Z_real[k][i][i];
        real_ij=Z_real[k][i][j];
        real_ji=Z_real[k][j][i];
        real_jj=Z_real[k][j][j];
        imag_ii=Z_imag[k][i][i];
        imag_ij=Z_imag[k][i][j];
        imag_ji=Z_imag[k][j][i];
        imag_jj=Z_imag[k][j][j];
        switch(choice_func){
          case 1:
            if (i==j) {
              ddOmega=0.0;
            }else{
              ddOmega=-((real_ii*real_jj+imag_ii*imag_jj)
                         +2.0*(real_ij*real_ij+imag_ij*imag_ij))
                      *(1.0/(real_ii*real_ii+imag_ii*imag_ii)
                       +1.0/(real_jj*real_jj+imag_jj*imag_jj)) + 2.0
                      +4.0*(real_ij*real_jj+imag_ij*imag_jj)
                          *(real_ij*real_jj+imag_ij*imag_jj)
                          /pow((real_jj*real_jj+imag_jj*imag_jj),3.0)
                      -4.0*(real_ij*real_ii+imag_ij*imag_ii)
                          *(real_ij*real_ii+imag_ij*imag_ii)
                          /pow((real_ii*real_ii+imag_ii*imag_ii),3.0);
            }
          break;

         case 2:
           if (i==j) {
             ddOmega=0.0;
           }else{
             ddOmega= (real_ii*real_ii+imag_ii*imag_ii+
                        real_jj*real_jj+imag_jj*imag_jj)
                 -2.0*(real_ii*real_jj+imag_ii*imag_jj)
                 -4.0*(real_ij*real_ij+imag_ij*imag_ij);
           }
         break;

         case 3:
           if (i==j) {
             ddOmega=0.0;
           }else{
             ddOmega=-( (real_ii*real_jj+imag_ii*imag_jj)
                       +2.0*(real_ij*real_ij+imag_ij*imag_ij))
                     *(0.5/sqrt(real_ii*real_ii+imag_ii*imag_ii) +
                       0.5/sqrt(real_jj*real_jj+imag_jj*imag_jj))
                     +0.5*sqrt(real_ii*real_ii+imag_ii*imag_ii)
                     +0.5*sqrt(real_jj*real_jj+imag_jj*imag_jj)
                     +0.5*(real_ij*real_jj+imag_ij*imag_jj)
                         *(real_ij*real_jj+imag_ij*imag_jj)
                         /pow((real_jj*real_jj+imag_jj*imag_jj),2.0)
                     -0.5*(real_ij*real_ii+imag_ij*imag_ii)
                         *(real_ij*real_ii+imag_ij*imag_ii)
                         /pow((real_ii*real_ii+imag_ii*imag_ii),2.0);
            }
          break;

          default:
            printf("@@@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@\n");
            printf(" Wrong choice for the functional should have been rejected\n");
            printf(" previously\n");
            printf("@@@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@\n");
            fflush(stdout);
            Finalize(); exit(1);
          break;
        }
        ddf[(i-1)*nstate+j]+=(4.0*W[k]/fpisq*ddOmega);
      }
    }
    }/*end if*/
  }

/* ========================================================================== */
}/*end routine                                                               */
/* ========================================================================== */

void diagonalize_asymm_mat(double **A,int nstate,double *D, double **U_real,
                           double **U_imag, double **R_real, double **R_imag,
                           int choice_diag)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/
/*            Local variable declarations                                    */
/*===========================================================================*/

#include "../typ_defs/typ_mask.h"

/* declaration of local variables */

  int i,j,k,ii,jj,kk;
  int INFO,N,LDZ,LDA,count1,count,flag,IL,IU,M;
  double tmp,VL,VU,ABSTOL;  /* for safety ABSTOL=DLAMCH(CMACH) with CMACH='S'; */

  int *ISUPPZ;
  double *DD,*EE;

  zomplex *TAU,**iTimesA,**Z;

  char UP;          /* zhetrd uses only the upper part of iTimesA in diagonalization */
  char JOBZ,COMPZ;  /* zstedc and zsteqr will calculate BOTH the eigenvalues and
                      eigenvectors of iTimesA */
  char RANGE;       /* all eigenvalues and all eigenvectors */
  char CMACH;       /* CMACH='S' can be used such that DLAMCH(CMACH) returns a
                      "safe minimum, i.e: 1/sfmin will not overflow base */
  int NB;           /* optimal blocksize in ZHETRD ; not optimized ! */

  int LWORK1=nstate*nstate;/* optimal value for LWORK is returned by ZHETRD in WORK[0]*/
  int LWORK2=2*nstate*nstate;  /* optimal LWORK in divide and conquer and zsteqr ? */
  int LWORK3= nstate*nstate+64*nstate; /*  optimal LWORK for the RRR alg ? */
  int LRWORK2=5*nstate*nstate+24*nstate+1; /* optimal LRWORK for divide and conquer ? */
  int LRWORK3= 200*nstate;   /* optimal LRWORK for the RRR algorithm ? */
  int LIWORK2=5*nstate*nstate+6*nstate+6;  /* optimal LIWORK for divide and conquer */
  int LIWORK3=100*nstate;       /* optimal LIWORK for the RRR algorithm */

  zomplex *WORK1,*WORK2,*WORK3;
  double *RWORK2,*RWORK1,*RWORK3;
  int *IWORK2,*IWORK3;

  int NDIM = nstate*nstate;
  static int first_time=1;
  static double *old,*oldA;

/*===============================================================================*/
/* 0) first time memory allocation                                               */

  ISUPPZ = (int *) malloc((size_t)(2*nstate*sizeof(int)));
  DD = (double *) malloc((size_t)(nstate*sizeof(double)));
  EE = (double *) malloc((size_t)((nstate-1)*sizeof(double)));

  iTimesA=cmall_zomp2(0,nstate-1,0,nstate-1);
  Z = cmall_zomp2(0,nstate-1,0,nstate-1);
  TAU = (zomplex *) malloc((size_t)((nstate-1)*sizeof(zomplex)));

  WORK1 = (zomplex *) malloc((size_t)(LWORK1*sizeof(zomplex)));
  WORK2 = (zomplex *) malloc((size_t)(LWORK2*sizeof(zomplex)));
  WORK3 = (zomplex *) malloc((size_t)(LWORK3*sizeof(zomplex)));
  
  RWORK1 = (double *) malloc((size_t)((2*nstate-2)*sizeof(double)));
  RWORK2 = (double *) malloc((size_t)(LRWORK2*sizeof(double)));
  RWORK3 = (double *) malloc((size_t)(LRWORK3*sizeof(double)));

  IWORK2 = (int *) malloc((size_t)(LIWORK2*sizeof(int))); 
  IWORK3 = (int *) malloc((size_t)(LIWORK3*sizeof(int)));

  if (first_time==1){
    oldA=(double *)cmalloc(NDIM*sizeof(double)); /* No -1 here !!! */
    old=(double *)cmalloc(5*NDIM*sizeof(double)); /* so it starts from 0 */
  }

  count=0;
  count1=0;
  flag=1;

  if(first_time==0){
    for(ii=0;ii<nstate;ii++){
      for(jj=0;jj<nstate;jj++){
        if(A[jj+1][ii+1]!= oldA[count1]) flag=0;
        count1 += 1;
      }
    }
  }

/*=============================================================================*/
/* I) NOT the first time AND  A==oldA  */

  if((flag==1) && (first_time==0)){ /* A==oldA */
    for (i=1;i<=nstate;i++){
      D[i]=old[count];
      count +=1;
    }
    for (i=1;i<=nstate;i++) {
      for (j=1;j<=nstate;j++){
        R_real[i][j]=old[count];
        count += 1;
        R_imag[i][j]=old[count];
        count += 1;
      }
    }
    for (i=1;i<=nstate;i++){
      for (j=1;j<=nstate;j++){
        U_real[i][j]=old[count];
        count += 1;
        U_imag[i][j]=old[count];
        count += 1;
      }
    }
    return;
  }

/*==========================================================================*/
/* II) ALL other cases   */

  count1 = 0;

  for(ii=0;ii<nstate;ii++){
    for(jj=0;jj<nstate;jj++){
      iTimesA[ii][jj].re=0.0;
      /* the matrix iTimesA is i times A, i=sqrt(-1).*/
      oldA[count1]=iTimesA[ii][jj].im=A[jj+1][ii+1];
      /*why transposed? Eigenvectors are saved in transposed form. See below*/

      count1 += 1;

    }
  }

  count=0;
  UP='U'; N=nstate; LDA=LDZ=nstate; CMACH='S';
  ABSTOL=0.00000001;
  JOBZ='V'; RANGE='A'; COMPZ='V';

  switch(choice_diag){
    case 1 :

      ZHETRD(&UP,&N,&iTimesA[0][0],&LDA,&DD[0],&EE[0],&TAU[0],&WORK1[0],&LWORK1,&INFO);
      if (INFO!=0) {
        printf("Error: the %d argument in zhetrd had an illegal value\n",-INFO);
      }
      ZUNGTR(&UP,&N,&iTimesA[0][0],&LDA,&TAU[0],&WORK1[0],&LWORK1,&INFO);
      if (INFO!=0){
        printf("Error: the %d argument in zungtr had an illegal value\n",-INFO);
        fflush(NULL);
        exit(1);
      }
      /* The classical QR algorithm ; usually the slowest of the 3 choices,
         but the most robust */
      ZSTEQR(&COMPZ,&N,&DD[0],&EE[0],&iTimesA[0][0],&LDA,&RWORK1[0],&INFO);
    break;

    case 2 :
      /* the divide and conquer algorithm */
      ZHEEVD(&JOBZ,&UP,&N,&iTimesA[0][0],&LDA,&DD[0],&WORK2[0],&LWORK2,&RWORK2[0],
              &LRWORK2,&IWORK2[0],&LIWORK2,&INFO);
    break;

    case 3 :
     /* the relatively robust method -- should be the fastest,
        but sometimes produces NAN's */
      ZHEEVR(&JOBZ,&RANGE,&UP,&N,&iTimesA[0][0],&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,&DD[0],
            &Z[0][0],&LDZ,&ISUPPZ[0],&WORK3[0],&LWORK3,&RWORK3[0],&LRWORK3,&IWORK3[0],
            &LIWORK3,&INFO);
    break;

    default :
      printf("@@@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("Error in the choice for the diagonalization algorithm !!!\n");
      printf("The legal choices are :                                  \n");
      printf(" 1. The QR method --- used in Lapack 1.0                 \n");
      printf(" 2. The divide and conquer method --- used in Lapack 2.0 \n");
      printf(" 3. The relatively robust representation --- Lapack 3.0  \n");
      printf("@@@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@@@@\n");

      fflush(stdout);
      Finalize();
      exit(1);
    break;
  }

/*--------------------------------------------------------------------------*/
/* III) Error message from LAPACK routine  */

  if (INFO!=0){
    printf("@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf(" Error in asymetric matrix diagonalization    \n");
    printf(" LAPACK routine: zsteqr_ (QR), zheevd_ (DQ), zheevr_ (R) \n");
    if (INFO<0){
      printf("Error: the %d argument in the diagonalization function\n",-INFO);
      printf("       had an illegal value\n");
    }else{
      printf("The diagonalization algorithm failed\n");
    }
    printf("@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(NULL);
    exit(1);
  }

/*==============================================================================*/
/* IV) Update old[]                                                             */
/* Note that DD contains the eigenvalues and iTimesA contains eigenvectors now  */

  for(i=1;i<=nstate;i++){
    old[count]=DD[i-1];
    D[i]=DD[i-1];
    count += 1; /* count was not reset for this for loop*/
  }

  if (choice_diag!=3){
    for (i=1;i<=nstate;i++){
      for (j=1;j<=nstate;j++){
        old[count]=R_real[i][j]=iTimesA[j-1][i-1].re;
        count +=1;
        old[count]=R_imag[i][j]=iTimesA[j-1][i-1].im;
        count +=1;
      }
    }
  }else{
    for (i=1;i<=nstate;i++){
      for (j=1;j<=nstate;j++){
        old[count]=R_real[i][j]=Z[j-1][i-1].re;
        count +=1;
        old[count]=R_imag[i][j]=Z[j-1][i-1].im;
        count +=1;
      }
    }
  }/*end if*/

/*========================================================================*/
/* V) Compute U=exp(A)=R*exp(-i*D)*R                                      */

  for (i=1;i<=nstate;i++){
    for (j=1;j<=nstate;j++){
      U_real[i][j]=U_imag[i][j]=0.0;
      for (k=1;k<=nstate;k++){
        U_real[i][j]+= cos(D[k])*(R_real[i][k]*R_real[j][k]
                                 +R_imag[i][k]*R_imag[j][k])
                      +sin(D[k])*(R_real[j][k]*R_imag[i][k]
                                 -R_imag[j][k]*R_real[i][k]);
        U_imag[i][j]+= cos(D[k])*(R_real[j][k]*R_imag[i][k]
                                 -R_imag[j][k]*R_real[i][k])
                      -sin(D[k])*(R_real[i][k]*R_real[j][k]
                                 +R_imag[i][k]*R_imag[j][k]);
      } /* end for k */
      old[count]=U_real[i][j];
      count +=1;
      old[count]=U_imag[i][j];
      count +=1;
    } /* end for j */
  }/*endfor i*/

 first_time=0;

/* ========================================================================== */
}/*end routine                                                               */
/* ========================================================================== */

void  component_wise_matrix_multiply(int nstate,int flag1,int flag2,
           double **A1_real, double **A1_imag, double **A2_real, double **A2_imag,
           double **A3_real, double **A3_imag)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/
/*            Local variable declarations                                    */
/*===========================================================================*/
  int i,j;

  for (i=1;i<=nstate;i++){
    for (j=1;j<=nstate;j++){
      if (flag1==1 && flag2==1){
        A3_real[i][j]=A1_real[i][j]*A2_real[i][j]-A1_imag[i][j]*A2_imag[i][j];
        A3_imag[i][j]=A1_real[i][j]*A2_imag[i][j]+A1_imag[i][j]*A2_real[i][j];
      }else if (flag1==1 && flag2==2){
        A3_real[i][j]=A1_real[i][j]*A2_real[j][i]+A1_imag[i][j]*A2_imag[j][i];
        A3_imag[i][j]=-A1_real[i][j]*A2_imag[j][i]+A1_imag[i][j]*A2_real[j][i];
      }else if (flag1==2 && flag2==1){
        A3_real[i][j]=A1_real[j][i]*A2_real[i][j]+A1_imag[j][i]*A2_imag[i][j];
        A3_imag[i][j]=A1_real[j][i]*A2_imag[i][j]-A1_imag[j][i]*A2_real[i][j];
      }else if (flag1==2 && flag2==2){
        A3_real[i][j]=A1_real[j][i]*A2_real[j][i]-A1_imag[j][i]*A2_imag[j][i];
        A3_imag[i][j]=-A1_real[j][i]*A2_imag[j][i]-A1_imag[j][i]*A2_real[j][i];
      }else {
        printf("@@@@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        printf(" Error in component-wise multiplication. The legal option for\n");
        printf(" second and third integers when multiplying matrices with    \n");
        printf(" complex numbers are 1 or 2:                                 \n");
        printf("     1 represent normal matrix,      \n");
        printf("     2 represents its complex conjugate\n");
        printf("@@@@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        exit(1);
      }
    }/*endfor j*/
  }/*endfor i*/

/* ========================================================================== */
}/*end routine                                                               */
/* ========================================================================== */


void comp_tensor_z(GENERAL_DATA *general_data, CP *cp, double *W,
                            double ***Z_real, double ***Z_imag)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/

#include "../typ_defs/typ_mask.h"

  int cp_lsda        = cp->cpopts.cp_lsda;
  int myid           = cp->cp_comm_state_pkg_up.myid;
  int num_proc       = cp->cp_comm_state_pkg_up.num_proc;
  MPI_Comm world     = cp->cp_comm_state_pkg_up.world;
  MPI_Comm comm      = cp->cp_comm_state_pkg_up.comm;

  int nstate_up=cp->cpcoeffs_info.nstate_up;                  
  int nstate_dn=cp->cpcoeffs_info.nstate_dn;                 
  int nstate=nstate_up;                                   
  int ncoef=cp->cpcoeffs_info.ncoef;

  int nk1 = cp->cpewald.kmax_cp[1];     
  int nk2 = cp->cpewald.kmax_cp[2];   
  int nk3 = cp->cpewald.kmax_cp[3];     

  int nstate_max            = cp->cp_comm_state_pkg_up.nstate_max;
  int nstate_proc           = cp->cp_comm_state_pkg_up.nstate_proc;
  int nstate_proc_max       = cp->cp_comm_state_pkg_up.nstate_proc_max;
  int nstate_ncoef_proc_max = cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
  int nstate_ncoef_proc     = cp->cp_comm_state_pkg_up.nstate_ncoef_proc;
  int icoef_start           = cp->cp_comm_state_pkg_up.icoef_start;

  int istate_up_st          = cp->cpcoeffs_info.istate_up_st;
  int istate_up_end         = cp->cpcoeffs_info.istate_up_end;
  int istate_up_st_recv, istate_up_end_recv;

  static int ifirst_time = 1;
  static int max_buffer;
  static int ***triplet;
  static char *buffer;
  static int *count_proc;
  static double ***creal_imag_f_b;

  int i,j,k,l,m,n,g1,g2,g3;
  int I=6;
  int count,pos,position,s1,s2,isend,iproc_send, iproc_recv,rounds;

  double *creal1,*creal2,*cimag1,*cimag2;

/*=============================================================================*/
/* 0) Check   */

  if(cp_lsda==1){
    if(myid==0){
      printf("@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf(" Wannier dynamics has not been implemented for \n");
      printf(" spin-polarized DFT (LSDA) \n");
      printf("@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@@\n");
    }
    fflush(stdout);
    Finalize();
    exit(1);
  }

/*==============================================================================*/
/* I) First time allocation  */

  if (ifirst_time==1){

    triplet= cmall_itens3(-nk1-2,nk1+2,-nk2-2,nk2+2,-nk3-2,nk3+2);

    /* max number of buffer space neeeded to transmit data */
    max_buffer=(nstate/num_proc+2)*(nstate/num_proc+2)*(num_proc/2+3)
               *(sizeof(int)+6*sizeof(double));
    if ((W[4]>0.0) || (W[5]>0.0) || (W[6]>0.0)) {
      max_buffer*=2;
    }
    buffer=(char *)cmalloc(max_buffer*sizeof(char));

    /* count[i]=how many elements are received from processor i */
    if (num_proc>1) {
      count_proc = cmall_ivector(1,num_proc-1);
      creal_imag_f_b=cmall_tens3(1,4,1,3,1,ncoef*nstate_proc_max);
    }else{
      creal_imag_f_b=cmall_tens3(1,4,1,3,1,ncoef*nstate);
    }
  }/*endif*/

/*===========================================================================*/
/* II) Initialization   */

  for (k=1;k<=I;k++) {
    for (i=0;i<nstate;i++) {
      for (j=0;j<nstate;j++) {
        Z_real[k][i+1][j+1]=Z_imag[k][i+1][j+1]=0.0;
      }
    }
  }

  creal1=cp->cpcoeffs_pos[1].cre_up;
  creal2=cp->cpcoeffs_pos[1].cre_up;
  cimag1=cp->cpcoeffs_pos[1].cim_up;
  cimag2=cp->cpcoeffs_pos[1].cim_up;

  for (i=-nk1-2;i<=nk1+2;i++) {
    for (j=-nk2-2;j<=nk2+2;j++) {
      for (k=-nk3-2;k<=nk3+2;k++) {
        triplet[i][j][k]=0;
      }
    }
  }

  for (i=1;i<=ncoef;i++) {
    g1=cp->cpewald.kastr_sm[i];
    g2=cp->cpewald.kbstr_sm[i];
    g3=cp->cpewald.kcstr_sm[i];
    triplet[g1][g2][g3]=i;
  }
/*==========================================================================*/
/* III) Rearrange coefficient along x-direction */


  order_coef (creal1,cimag1,creal_imag_f_b[1],creal_imag_f_b[2],creal_imag_f_b[3],
              creal_imag_f_b[4], cp,triplet, ncoef,nstate_proc,nstate_proc_max);


/*==========================================================================*/
/* IV) calculate tensor elements with my coefficients */


  istate_up_st_recv = istate_up_st;
  istate_up_end_recv = istate_up_end;

  comp_tensor_z_prim(Z_real,Z_imag,istate_up_st,istate_up_end,istate_up_st_recv,
                     istate_up_end_recv,ncoef,nstate,creal1,cimag1,creal_imag_f_b[1],
                     creal_imag_f_b[2],creal_imag_f_b[3],creal_imag_f_b[4]);


  if(num_proc>1){

/*============================================================================*/
/* V) calculate tensor elements with coefficients on other PE */

    rounds=(num_proc%2==0)?(num_proc+1)/2:(num_proc-1)/2;


    for (isend=1;isend<=rounds;isend++) {

      /* send and receive coefficient */

      iproc_send=(myid+1)%num_proc;
      iproc_recv=(num_proc+myid-1)%num_proc;

      Barrier(comm);
      Sendrecv_replace(&istate_up_st_recv,1,MPI_INT,iproc_send,myid,iproc_recv,
                       iproc_recv,comm);
      Sendrecv_replace(&istate_up_end_recv,1,MPI_INT,iproc_send,myid,iproc_recv,
                       iproc_recv,comm);

      for(k=1;k<=4;k++){
        for(j=1;j<=3;j++){
          Barrier(comm);
          Sendrecv_replace(&creal_imag_f_b[k][j][1],nstate_proc_max*ncoef,MPI_DOUBLE,
                           iproc_send,myid,iproc_recv,iproc_recv,comm);
        }
      }

      /* compute z-element using received coefficients */


      comp_tensor_z_prim(Z_real,Z_imag,istate_up_st,istate_up_end,istate_up_st_recv,
                           istate_up_end_recv,ncoef,nstate,creal1,cimag1,creal_imag_f_b[1],
                           creal_imag_f_b[2],creal_imag_f_b[3],creal_imag_f_b[4]);
    }
 
/*==================================================================================*/
/* VI) Send the non-zero elements of Z_real and Z_imag on processor 0 */

  
    count=position=0;

    /* data packing */

    if (myid!=0) {
      for (i=istate_up_st-1;i<istate_up_st-1+nstate_proc;i++){
        for (j=0;j<istate_up_st-1;j++) {
          if ((Z_real[1][i+1][j+1]!=0.0) || (Z_imag[1][i+1][j+1]!=0.0) ||
              (Z_real[2][i+1][j+1]!=0.0) || (Z_imag[2][i+1][j+1]!=0.0) ||
              (Z_real[3][i+1][j+1]!=0.0) || (Z_imag[3][i+1][j+1]!=0.0)) {
            count++;
            pos=IMIN(i,j)*nstate+IMAX(i,j);
            Pack(&pos,1,MPI_INT,buffer,max_buffer,&position,comm);
            for (k=1;k<=I;k++) if (W[k]!=0.0) {
              Pack(&Z_real[k][i+1][j+1],1,MPI_DOUBLE,buffer,max_buffer,&position,comm);
              Pack(&Z_imag[k][i+1][j+1],1,MPI_DOUBLE,buffer,max_buffer,&position,comm);
            } /* end for k */
          } /* end if Z_real, .. */
        } /* end for j */

        for (j=i;j<istate_up_st-1+nstate_proc;j++) {
          if ((Z_real[1][i+1][j+1]!=0.0) || (Z_imag[1][i+1][j+1]!=0.0) ||
              (Z_real[2][i+1][j+1]!=0.0) || (Z_imag[2][i+1][j+1]!=0.0) ||
              (Z_real[3][i+1][j+1]!=0.0) || (Z_imag[3][i+1][j+1]!=0.0)) {
            count++;
            pos=IMIN(i,j)*nstate+IMAX(i,j);
            Pack(&pos,1,MPI_INT,buffer,max_buffer,&position,comm);
            for (k=1;k<=I;k++) if (W[k]!=0.0) {
              Pack(&Z_real[k][i+1][j+1],1,MPI_DOUBLE,buffer,max_buffer,&position,comm);
              Pack(&Z_imag[k][i+1][j+1],1,MPI_DOUBLE,buffer,max_buffer,&position,comm);
            } /* end for k */
          } /* end if Z_real, .. */
        } /* end for j */

        for (j=istate_up_st-1+nstate_proc;j<nstate;j++) {
          if ((Z_real[1][i+1][j+1]!=0.0) || (Z_imag[1][i+1][j+1]!=0.0) ||
              (Z_real[2][i+1][j+1]!=0.0) || (Z_imag[2][i+1][j+1]!=0.0) ||
              (Z_real[3][i+1][j+1]!=0.0) || (Z_imag[3][i+1][j+1]!=0.0)) {
            count++;
            pos=IMIN(i,j)*nstate+IMAX(i,j);
            Pack(&pos,1,MPI_INT,buffer,max_buffer,&position,comm);
            for (k=1;k<=I;k++) {
              if (W[k]!=0.0) {
                Pack(&Z_real[k][i+1][j+1],1,MPI_DOUBLE,buffer,max_buffer,&position,comm);
                Pack(&Z_imag[k][i+1][j+1],1,MPI_DOUBLE,buffer,max_buffer,&position,comm);
              }
            } /* end for k */
          } /* end if Z_real, .. */
        } /* end for j */
      } /* end for i*/
    } /* myid!=0 */

    /* broadcasting */

    if (myid!=0) {
      Send(&count,1,MPI_INT,0,myid,comm);
    }else {
      for (i=1;i<=num_proc-1;i++) {
        Recv(&count_proc[i],1,MPI_INT,i,i,comm);
      }
    }

    /* send the packed elements and unpacking*/

    if (myid!=0) {
      Send(buffer,position,MPI_PACKED,0,myid,comm);
    }else {
      for (i=1;i<=num_proc-1;i++){
        position=0;
        Recv(buffer,max_buffer,MPI_PACKED,i,i,comm);
        for (j=1;j<=count_proc[i];j++){
          Unpack(buffer,max_buffer,&position,&pos,1,MPI_INT,comm);
          m=pos/nstate;
          n=pos%nstate;
          for (k=1;k<=I;k++) {
            if (W[k]!=0.0) {
              Unpack(buffer,max_buffer,&position,&Z_real[k][n+1][m+1],1,MPI_DOUBLE,comm);
              Unpack(buffer,max_buffer,&position,&Z_imag[k][n+1][m+1],1,MPI_DOUBLE,comm);
              Z_real[k][m+1][n+1]=Z_real[k][n+1][m+1];
              Z_imag[k][m+1][n+1]=Z_imag[k][n+1][m+1];
            }
          }
        } /* end for j */
      } /* end for i */
    }/*endif myid*/

  }/*endif parallel */

  ifirst_time=0;

/* ========================================================================== */
}/*end routine                                                               */
/* ========================================================================== */


void comp_tensor_z_prim(double ***Z_real,double ***Z_imag,int s1_start,
                     int s1_end, int s2_start,int s2_end, int ncoef, int nstate,
                     double *crealF, double *cimagF,double **crealf, 
                     double **cimagf, double **creals, double **cimags)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/
/*            Local variable declarations                                    */
/*===========================================================================*/

#include "../typ_defs/typ_mask.h"

  char transa,transb;
  int mm,nn,kk,lda,ldb,ldc,i,j,k;
  double alpha,beta;
  int direction;
  double *creal1,*cimag1,*creal2,*cimag2,*creal3,*cimag3;

  for (direction=1;direction<=3;direction++){
    creal1=crealF; 
    cimag1=cimagF;
    creal2=crealf[direction]; 
    cimag2=cimagf[direction];
    creal3=creals[direction]; 
    cimag3=cimags[direction];

    /*----------------------------------------------------------------*/
    /*                                                           T    */
    /* compute Z[direction]=(creal1-i cimag1) * (creal2+i cimag2)     */
    /*                                                                */
    /* Since Fortran and C have different convention, it should look  */
    /* like                                                           */
    /*                  T                                             */  
    /* (creal2+i cimag2) * (creal1-i cimag1)                          */
    /*----------------------------------------------------------------*/

    transa='T'; transb='N'; mm=s2_end-s2_start+1; nn=s1_end-s1_start+1; kk=ncoef; 
    lda=ncoef; ldb=ncoef; ldc=nstate;

    alpha=0.5; beta=0.0;
    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&creal2[1],&lda,&creal1[1],&ldb,&beta,
          &Z_real[direction][s1_start][s2_start],&ldc);
    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&cimag2[1],&lda,&creal1[1],&ldb,&beta,
          &Z_imag[direction][s1_start][s2_start],&ldc);
  
    alpha=0.5; beta=1.0;
    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&cimag2[1],&lda,&cimag1[1],&ldb,&beta,
          &Z_real[direction][s1_start][s2_start],&ldc);

    alpha=-0.5; beta=1.0;
    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&creal2[1],&lda,&cimag1[1],&ldb,&beta,
          &Z_imag[direction][s1_start][s2_start],&ldc);

    /*----------------------------------------------------------------*/
    /*                                                          T     */
    /* compute Z[direction]+=(creal1+i cimag1)*(creal3-i cimag3)      */
    /*                            T                                   */
    /* Looks like (creal3-icimag3)  *(creal1+i cimag1)                */
    /*----------------------------------------------------------------*/

    transa='T'; transb='N'; mm=s2_end-s2_start+1; nn=s1_end-s1_start+1; kk=ncoef; 
    lda=ncoef; ldb=ncoef; ldc=nstate;

    alpha=0.5; beta=1.0;
    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&creal3[1],&lda,&creal1[1],&ldb,&beta,
          &Z_real[direction][s1_start][s2_start],&ldc);
    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&creal3[1],&lda,&cimag1[1],&ldb,&beta,
          &Z_imag[direction][s1_start][s2_start],&ldc);
    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&cimag3[1],&lda,&cimag1[1],&ldb,&beta,
          &Z_real[direction][s1_start][s2_start],&ldc);

    alpha=-0.5;beta=1.0;
    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&cimag3[1],&lda,&creal1[1],&ldb,&beta,
          &Z_imag[direction][s1_start][s2_start],&ldc);
  }


  for (direction=1;direction<=3;direction++){
    for (i=s1_start;i<=s1_end;i++){
      for (j=s2_start;j<=s2_end;j++){
        if (Z_real[direction][j][i]==0.0 && Z_imag[direction][j][i]==0.0){
          Z_real[direction][j][i]=Z_real[direction][i][j]; 
          Z_imag[direction][j][i]=Z_imag[direction][i][j];
        }
      }
    }
  }

/* ========================================================================== */
}/*end routine                                                                */
/* ========================================================================== */


void order_coef(double *creal1,double *cimag1,double **creal2,double **cimag2, 
                double **creal3, double **cimag3, CP *cp,int ***triplet,
                int ncoef,int nstate1,int nstate2)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/
/*            Local variable declarations                                    */
/*===========================================================================*/

#include "../typ_defs/typ_mask.h"

  int i,j,g1,g2,g3,J1,J2,J3,direction;
  double *p1_real,*p1_imag,*p2_real,*p2_imag,*p3_real,*p3_imag;

  if (nstate1>nstate2){
    printf("Error in order_coef: nstate1 is greater than nstate2\n");
    Finalize(); 
    exit(1);
  }
 
  for (direction=1;direction<=3;direction++){
    p1_real=creal1;
    p2_real=creal2[direction];
    p3_real=creal3[direction];
    p1_imag=cimag1;
    p2_imag=cimag2[direction];
    p3_imag=cimag3[direction];

    for (j=1;j<=nstate1;j++){

      for (i=1;i<=ncoef;i++){ 
        g1=cp->cpewald.kastr_sm[i]; 
        g2=cp->cpewald.kbstr_sm[i]; 
        g3=cp->cpewald.kcstr_sm[i];

        switch(direction){
          case 1:
            J1=triplet[g1+1][g2][g3]; 
            J2=triplet[g1-1][g2][g3]; 
            J3=triplet[1-g1][-g2][-g3];
          break;
          case 2:
            J1=triplet[g1][g2+1][g3]; 
            J2=triplet[g1][g2-1][g3]; 
            J3=triplet[-g1][1-g2][-g3];
          break;
          case 3:
            J1=triplet[g1][g2][g3+1]; 
            J2=triplet[g1][g2][g3-1]; 
            J3=triplet[-g1][-g2][1-g3];
          break;
          default:
            printf("Error in order_coef: illegal direction!\n");
            Finalize(); exit(1);
        }

        if (J1){
          p2_real[i]=p1_real[J1];
          p2_imag[i]=p1_imag[J1];
        }else{
          p2_real[i]=p2_imag[i]=0.0;
        }

        if (J2){
          p3_real[i]=p1_real[J2];
          p3_imag[i]=p1_imag[J2];
        }else if (J3){
          p3_real[i]=p1_real[J3];
          p3_imag[i]=-p1_imag[J3];
        }else{
          p3_real[i]=p3_imag[i]=0.0;
        }
      }/*endfor ncoef*/

      p3_real[ncoef]=p3_imag[ncoef]=0.0;

      p1_real+=ncoef;
      p1_imag+=ncoef;
      p2_real+=ncoef;
      p2_imag+=ncoef;
      p3_real+=ncoef;
      p3_imag+=ncoef;

    }/*endfor state*/
  }/*endfor direction*/

/* ========================================================================== */
}/*end routine                                                                */
/* ========================================================================== */

