/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: analysis_cp                                  */
/*                                                                          */
/* This subprogram performs on the fly analysis of CP data                  */
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
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_wannier_cpcon_local.h"
#include "../proto_defs/proto_wannier_cpcon_entry.h"

#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0
#define its_max 200

#define DEBUG_WANNIER_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calcul_wannier(GENERAL_DATA *general_data,CP *cp, double ***Z_real,
                    double ***Z_imag, double *W, double **wan_cent,
                    int dip_calc_flag)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/

#include "../typ_defs/typ_mask.h"

  int myid         = cp->communicate.myid_state;
  int num_proc     = cp->communicate.np_states;
  MPI_Comm comm    = cp->communicate.comm_states;

  int itime = general_data->timeinfo.itime;

  int icoef_form_up    = cp->cpcoeffs_pos[1].icoef_form_up;
  int icoef_form_dn    = cp->cpcoeffs_pos[1].icoef_form_dn;
  int cp_wan_min_on    = cp->cpopts.cp_wan_min_opt;
  int cp_wan_on        = cp->cpopts.cp_wan_opt;
  int *ioff_st         = cp->cpcoeffs_info.ioff_upt;
  int cp_lsda          = cp->cpopts.cp_lsda;
  int nstate_up        = cp->cpcoeffs_info.nstate_up;
  int nstate_dn        = cp->cpcoeffs_info.nstate_dn;
  int ncoef_up         = cp->cp_comm_state_pkg_up.nstate_ncoef_proc;
  int ncoef_dn         = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc;
  int nstate=nstate_up;
  int tot_beads=cp->cpcoeffs_info.pi_beads;

  double alpha=1.0, beta=0.0;

  int NDIM=nstate*nstate;
  double *elementA = cp->electronic_properties.elementA;
  double *U_final = cp->cpscr.cpscr_wannier.U_final;

  int i,j,k;

/*=======================================================================*/
/* 0) Parallel checks                                         */

  if(num_proc>1){
    if(icoef_form_up!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up Coef forces are not in normal form \n");
      printf("on state processor %d in calcul_wannier \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
      if(icoef_form_dn!=0){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Dn Coef forces are not in normal form \n");
        printf("on state processor %d in calcul_wannier \n",myid);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
    }/*endif*/
  }/*endif*/

  if(tot_beads != 1){
    if(myid==0){
      printf("@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("  No Path-Integral MD yet for Wannier/Diople calculation\n");
      printf("@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    }
    fflush(stdout);
    Finalize();
    exit(1);
  }

/*=============================================================================*/
/* I) Initialization                                                              */

  for (j=1;j<=nstate;j++){
    for (k=1;k<=3;k++){
      wan_cent[j][k]=0.0;
    }
  }
  for (i=1;i<=NDIM;i++){
    elementA[i]=0.0;
  }

/*==============================================================================*/
/* II) Compute z_I,ij tensor elements                                           */

  if(dip_calc_flag==0 ){
    comp_tensor_z(general_data,cp,W,Z_real,Z_imag);
  }


/*==============================================================================*/
/* III) Wannier minimization for proc=0                                            */

  if(num_proc > 1){
    Barrier(comm);
  }

  if(myid ==0){

    if(cp_wan_min_on==1 ){
      /*Full minimization: LDFP Minimization*/
      lbfg_min(general_data,cp,elementA,NDIM,Z_real,Z_imag,W);
    }else if(cp_wan_on==1){
      /*on-the_fly wannier: NEWTON Minimization*/
      newton(general_data,cp,elementA,NDIM,Z_real,Z_imag,W);
    }

    comp_wannier_center(general_data,cp,elementA,Z_real,Z_imag,
                         W,wan_cent,U_final);
  }

/*===========================================================================*/
/* IV) Broadcast the optimized A and Wannier centers                         */

  if(num_proc > 1){

    Barrier(comm);

    Bcast(&elementA[1],NDIM,MPI_DOUBLE,0,comm);
    Bcast(&U_final[1],NDIM,MPI_DOUBLE,0,comm);
    Bcast(&wan_cent[1][1],3*nstate,MPI_DOUBLE,0,comm);
  }

/*============================================================================*/
/* V) Rotate orbital, velocities and force                                    */

  if(num_proc>1){
    control_coef_transpose_fwd(cp,1);
  }

  cp_rotate_all(cp->cpcoeffs_pos[1].cre_up, cp->cpcoeffs_pos[1].cim_up,
                cp->cpcoeffs_pos[1].icoef_form_up,
                cp->cpcoeffs_pos[1].vcre_up, cp->cpcoeffs_pos[1].vcim_up,
                cp->cpcoeffs_pos[1].ivcoef_form_up,
                cp->cpcoeffs_pos[1].fcre_up, cp->cpcoeffs_pos[1].fcim_up,
                cp->cpcoeffs_pos[1].ifcoef_form_up,
                U_final, cp->cpcoeffs_info.ioff_upt,
                cp->cpscr.cpscr_wave.cre_up,cp->cpscr.cpscr_wave.cim_up,
                &(cp->cp_comm_state_pkg_up));

/*===========================================================================*/
} /*end routine */
/*============================================================================*/

void calcul_initial_wannier(GENERAL_DATA *general_data,CP *cp)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/
/*            Local variable declarations                                    */
/*===========================================================================*/
#include "../typ_defs/typ_mask.h"

  int i,j,k;

  int myid         = cp->communicate.myid_state;
  int num_proc     = cp->communicate.np_states;
  MPI_Comm comm    = cp->communicate.comm_states;

  int itime = general_data->timeinfo.itime;

  int icoef_form_up    = cp->cpcoeffs_pos[1].icoef_form_up;
  int icoef_form_dn    = cp->cpcoeffs_pos[1].icoef_form_dn;
  int cp_wan_min_on    = cp->cpopts.cp_wan_min_opt;
  int cp_wan_on        = cp->cpopts.cp_wan_opt;
  int *ioff_st         = cp->cpcoeffs_info.ioff_upt;
  int cp_lsda          = cp->cpopts.cp_lsda;
  int nstate_up        = cp->cpcoeffs_info.nstate_up;
  int nstate_dn        = cp->cpcoeffs_info.nstate_dn;
  int ncoef_up         = cp->cp_comm_state_pkg_up.nstate_ncoef_proc;
  int ncoef_dn         = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc;
  int nstate=nstate_up;
  int tot_beads=cp->cpcoeffs_info.pi_beads;
  int iwrite_init_wcent = cp->cpopts.iwrite_init_wcent;
  int iwrite_init_worb  = cp->cpopts.iwrite_init_worb;

  double alpha=1.0, beta=0.0;

  int NDIM=nstate*nstate;
  double *elementA = cp->electronic_properties.elementA;
  double *U_final = cp->cpscr.cpscr_wannier.U_final;

  int I=6;
  double *W         = cp->electronic_properties.weight;
  double ***Z_real  = cp->electronic_properties.Z_real;
  double ***Z_imag  = cp->electronic_properties.Z_imag;
  double **wan_cent = cp->electronic_properties.wannier_cent;;

/*=======================================================================*/
/* 0) Parallel checks                                         */

  if(num_proc>1){
    if(icoef_form_up!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up Coef forces are not in trasposed form \n");
      printf("on state processor %d in calcul_wannier_init \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
      if(icoef_form_dn!=1){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Dn Coef forces are not in transposed form \n");
        printf("on state processor %d in calcul_wannier_init \n",myid);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
    }/*endif*/
  }/*endif*/

  if(tot_beads != 1){
    if(myid==0){
      printf("@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("  No Path-Integral MD yet for Wannier/Diople calculation\n");
      printf("@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    }
    fflush(stdout);
    Finalize();
    exit(1);
  }

/*===============================================================================*/
/* I) Untraspose into normal form  */

  if(num_proc > 1){
    control_coef_transpose_bck(cp,1);
  }

/*=============================================================================*/
/* II) Initialization                                                              */

  for (j=1;j<=nstate;j++){
    for (k=1;k<=3;k++){
      wan_cent[j][k]=0.0;
    }
  }
  for (i=1;i<=NDIM;i++){
    elementA[i]=0.0;
  }

/*==============================================================================*/
/* III) Compute Miller indices and z_I,ij tensor elements */

  comp_miller_weight(general_data,W);

  comp_tensor_z(general_data,cp,W,Z_real,Z_imag);


/*==============================================================================*/
/* IV) Wannier minimization for proc=0                                            */

  if(num_proc > 1){
    Barrier(comm);
  }

  if(myid ==0){

    /*Full minimization: LDFP Minimization*/
    lbfg_min(general_data,cp,elementA,NDIM,Z_real,Z_imag,W);

    comp_wannier_center(general_data,cp,elementA,Z_real,Z_imag,
                         W,wan_cent,U_final);

  }

/*===========================================================================*/
/* V) Broadcast the optimized A and Wannier centers                         */

  if(num_proc > 1){

    Barrier(comm);

    Bcast(&elementA[1],NDIM,MPI_DOUBLE,0,comm);
    Bcast(&U_final[1],NDIM,MPI_DOUBLE,0,comm);
    Bcast(&wan_cent[1][1],3*nstate,MPI_DOUBLE,0,comm);
  }

/*============================================================================*/
/* VI) Rotate orbital, velocities and force                                    */

  if(num_proc>1){
    control_coef_transpose_fwd(cp,1);
  }

  cp_rotate_all(cp->cpcoeffs_pos[1].cre_up, cp->cpcoeffs_pos[1].cim_up,
                cp->cpcoeffs_pos[1].icoef_form_up,
                cp->cpcoeffs_pos[1].vcre_up, cp->cpcoeffs_pos[1].vcim_up,
                cp->cpcoeffs_pos[1].ivcoef_form_up,
                cp->cpcoeffs_pos[1].fcre_up, cp->cpcoeffs_pos[1].fcim_up,
                cp->cpcoeffs_pos[1].ifcoef_form_up,
                U_final, cp->cpcoeffs_info.ioff_upt,
                cp->cpscr.cpscr_wave.cre_up,cp->cpscr.cpscr_wave.cim_up,
                &(cp->cp_comm_state_pkg_up));

/*============================================================================*/
/* V) print out initial Wannier center and orbital */

  if(iwrite_init_wcent==1){
    if(myid==0){
      printf("\n");
      printf("     WFC X(A)          WFC Y(A)          WFC Z(A)     State\n\n");
      for (i=1;i<=nstate;i++){
        printf("  %- 18.5e%- 18.5e%- 18.5e%-8d\n", 
                 wan_cent[i][1], wan_cent[i][2], wan_cent[i][3],i);
      }
    }
  }
  printf("\n");

  if(iwrite_init_worb==1){
    write_wannier_orb(general_data, cp);
  }

/*===========================================================================*/
} /*end routine */
/*============================================================================*/


void lbfg_min(GENERAL_DATA *general_data, CP *cp, double *A,
               int NDIM, double ***Z_real, double ***Z_imag, double *W)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/

#include "../typ_defs/typ_mask.h"

  double tol=general_data->minopts.tol_wan_coef;
  int nstate_up=cp->cpcoeffs_info.nstate_up;
  int nstate=nstate_up;

/* local variable declaration */

  int i,j,k;
  int its,iflag,flag,info,diagco;
  int nn=NDIM,mm=5; /*number of corrections ?!*/
  double fp,xtol=0.000000000000001;
  double *g    = cp->cpscr.cpscr_wannier.g;
  double *diag = cp->cpscr.cpscr_wannier.diag;
  double *scr  = cp->cpscr.cpscr_wannier.scr;
  int  *iprint = cp->cpscr.cpscr_wannier.iprint;

  iprint[1]=1;

/*============================================================================*/
/* I) compute Omega and dOmega                                               */

  comp_omega(general_data,cp,A,Z_real,Z_imag,W,&fp,NDIM);

  comp_domega(general_data,cp,A,g,Z_real,Z_imag,W,NDIM);

/*============================================================================*/
/* II) Minimize the functional (omega)                                       */

  its=0;
  iflag=0;
  flag=0;
  diagco=0;
  while (its<its_max){

    LBFGS(&nn,&mm,&A[1],&fp,&g[1],&diagco,&diag[1],&iprint[1],&tol,&xtol,
           &scr[1],&iflag,&info);

    switch(iflag){
      case -3 :
        printf("Error in lbfg_min calling lbfgs: nn or mm are not positive!\n");
        fflush(stdout);
        Finalize();
        exit(0);
      break;

      case -2 :
        printf("Error: An element of the approximate inverse matrix diag\n");
        printf("in lbfg_min calling lbfgs is not positive!\n");
        fflush(stdout);
        Finalize();
        exit(0);
      break;

      case -1 :
        printf("The line search MCSRCH in lbfgs failed. INFO=%d\n",info);
        fflush(stdout);
        Finalize();
        exit(1);
      break;

      case 2 :
        printf("You must provide the approximate diagonal for the inverse hessian\n");
        fflush(stdout);
        Finalize();
        exit(0);
      break;

      case 1 :
        comp_omega(general_data,cp,A,Z_real,Z_imag,W,&fp,NDIM);
        comp_domega(general_data,cp,A,g,Z_real,Z_imag,W,NDIM);
      break;

      case 0 :
        its=its_max;
        flag=1;
      break;

      default :
        printf("Erroneous iflag returned by lbfgs in lbfg_min\n");
        fflush(stdout);
        Finalize();
        exit(0);
      break;
    } /* end switch */

    its++;
  } /*end while*/

/*===============================================================================*/
/* III) Maximum iteration is reached without convergence   */

  if (its>its_max && flag==0){
    printf("Error: lbfgs did not converge after 200 iterations\n");
    fflush(stdout);
    Finalize();
    exit(0);
  }

/*===========================================================================*/
} /*end routine */
/*============================================================================*/

void newton(GENERAL_DATA *general_data, CP *cp, double *A,
               int NDIM, double ***Z_real, double ***Z_imag, double *W)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/

#include "../typ_defs/typ_mask.h"

  double tol=general_data->minopts.tol_wan_coef;
  int nstate_up=cp->cpcoeffs_info.nstate_up;
  int nstate_dn=cp->cpcoeffs_info.nstate_dn;
  int nstate=nstate_up;
  int choice_func = cp->cp_wannier.wan_func_typ;

/* local variable declaration */

  int i,j,k;
  int its,iflag,flag,info,diagco;
  int nn=NDIM,mm=5; /*number of corrections ?!*/
  double *g  = cp->cpscr.cpscr_wannier.g;
  double *gg = cp->cpscr.cpscr_wannier.gg;
  double fp,fpp,fpm,sum;
  double xtol=0.000000000000001;
  double delta= 0.00001;

/*============================================================================*/
/* I) compute Omega and check the initial value of the functional            */

  comp_omega(general_data,cp,A,Z_real,Z_imag,W,&fp,NDIM);
  printf("Before Newton Update: Omega_i = =%- 15.10e\n",fp);

/*============================================================================*/
/* II) Compute domega_zero and ddomega_zero                                  */

  comp_domega_zero(nstate_up,nstate_dn,choice_func,g,Z_real,Z_imag,W);

  sum=0.0;
  for(i=1;i<=NDIM;i++){
    sum+= g[i]*g[i];
  }
  sum=sqrt(sum);
  printf("Before Newton Update: Norm(dOmega)= %- 15.10e\n",sum);

  comp_ddomega_zero(nstate_up,nstate_dn,choice_func,gg,Z_real,Z_imag,W);

#ifdef DEBUG_WANNIER

  for (i=1;i<=nstate;i++){
    for (j=1;j<=nstate;j++){
      A[(i-1)*nstate+j]+=delta;
      if (i!=j) A[(j-1)*nstate+i]-=delta;
      comp_omega(general_data,cp,A,Z_real,Z_imag,W,&fpp,NDIM);
      A[(i-1)*nstate+j]-=(2.0*delta);
      if (i!=j) A[(j-1)*nstate+i]+=(2.0*delta);
      comp_omega(general_data,cp,A,Z_real,Z_imag,W,&fpm,NDIM);
      A[(i-1)*nstate+j]+=delta;
      if (i!=j) A[(j-1)*nstate+i]-=delta;
      printf("fd_omega_zero[%d][%d]: Numerical= %- 15.10e, Analytic= %- 15.10e\n",
              i,j,(fpp-fpm)/(2.0*delta),g[(i-1)*nstate+j]);
      printf("sd_omega_zero[%d][%d]: Numerical= %- 15.10e, Analytic= %- 15.10e\n",
              i,j,(fpp+fpm-2.0*fp)/(delta*delta),gg[(i-1)*nstate+j]);

    }
  }  
#endif

/*============================================================================*/
/* IV) Update A matrix                                                        */

  for (i=1;i<=nstate;i++){
    for (j=1;j<=nstate;j++){
      if (i==j) {
        A[(i-1)*nstate+j]=0.0;
      }else{
        A[(i-1)*nstate+j]=(-g[(i-1)*nstate+j]/gg[(i-1)*nstate+j]);
      }
    }
  }

  comp_omega(general_data,cp,A,Z_real,Z_imag,W,&fp,NDIM);
  printf("After Newton Update: Omega_i = =%- 15.10e\n",fp);

  comp_domega(general_data,cp,A,g,Z_real,Z_imag,W,NDIM);

  sum=0.0;
  for(i=1;i<=NDIM;i++){
    sum+= g[i]*g[i];
  }
  sum=sqrt(sum);
  printf("After Newton Update: Norm(dOmega)= %- 15.10e\n",sum);

/*===========================================================================*/
} /*end routine */
/*============================================================================*/

void  comp_wannier_center(GENERAL_DATA *general_data,CP *cp,
                          double *x, double ***Z_real, double ***Z_imag,
                          double *W, double **wan_cent, double *U_final)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/

#include "../typ_defs/typ_mask.h"

  int nstate_up=cp->cpcoeffs_info.nstate_up;
  int nstate_dn=cp->cpcoeffs_info.nstate_dn;
  int nstate=nstate_up;
  int I=3;      /* maximum number of non-zero Miller indices: Changed from I=6!! */
                /* assumme cubic*/
  int choice_diag = cp->cp_wannier.wan_diag_typ;

  int NDIM=nstate*nstate;

  int i,j,kk,l,mm,nn,r,s,t,m,n,k,flag,lda,ldb,ldc,incx,incy;
  double alpha,beta,prod,distance,norm,pi;
  char transa,transb;

  double **A        = cp->cpscr.cpscr_wannier.A;
  double **R_real   = cp->cpscr.cpscr_wannier.R_real;
  double **R_imag   = cp->cpscr.cpscr_wannier.R_imag;
  double **U_real   = cp->cpscr.cpscr_wannier.U_real;
  double **U_imag   = cp->cpscr.cpscr_wannier.U_imag;
  double **U_tmp1   = cp->cpscr.cpscr_wannier.U_tmp1;
  double **phi      = cp->cpscr.cpscr_wannier.phi;
  double **HMatrix  = cp->cpscr.cpscr_wannier.HMatrix;
  double *real      = cp->cpscr.cpscr_wannier.real;
  double *imag      = cp->cpscr.cpscr_wannier.imag;
  double *D         = cp->cpscr.cpscr_wannier.D;

/*==================================================================================*/
/* I) Initialization */

  pi=M_PI;

  for(i=1;i<=nstate;i++){
    for (j=1;j<=nstate;j++){
      A[i][j]=x[(i-1)*nstate+j];
    }
  }

/*=================================================================================*/
/* II) Diagonalize the matrix A   */

  flag=0;
  for(i=1;i<=nstate;i++){
    for (j=1;j<=nstate;j++){
      if(A[i][j]!=0.0) flag=2;
    }
  }
  if (flag==2){
    for(i=1;i<=nstate;i++){
      for (j=1;j<=nstate;j++){
        if(fabs(A[i][j])>=0.1) flag=1; /* 0.1 is much bigger than before*/
      }
    }
    /* comment out the following line if Taylor expansion is OK */
    flag=1; /* always direct diagonalization if A =/= 0 */
  }

  if (flag==0){
    for(i=1;i<=nstate;i++){
      D[i]=0.0;
      for (j=1;j<=nstate;j++){
        R_imag[i][j]=0.0;
        U_imag[i][j]=0.0;
        if (i==j){
          R_real[i][j]=1.0;
          U_real[i][j]=1.0;
        }else{
          R_real[i][j]=0.0;
          U_real[i][j]=0.0;
        }
      } /* end for j */
    } /* endfor i; */
  }else if (flag==1){ /* too large elements to use Taylor summation */

    diagonalize_asymm_mat(A,nstate,D,U_real,U_imag,R_real,R_imag,choice_diag);

  }else{
    /* Taylor summation of order 4 via 2 matrix multiplication
       A_1=A; A_2=A_1*A_1; F_1=1/24 * A_2 + 1/6 * A_1 + I
       F_2=A_2*F_1+A_1+I;       F_2 is the desired result */

    transa='N'; transb='N'; mm=nstate; nn=nstate; kk=nstate; lda=nstate;
    ldb=nstate; ldc=nstate; alpha=1.0; beta=0.0;

    /* step 1: calculate R_real=A*A,
                         R_imag=1/24 * R_real + 1/6 * A_1 + 1/2 * I, U_real=A+I */
    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&R_real[1][1],&lda,&A[1][1],&ldb,
           &beta,&A[1][1],&ldc);
    for (i=1;i<=nstate;i++){
      for (j=1;j<=nstate;j++){
        if (i==j){
          R_imag[i][j]=0.5+A[i][j]/6.0;
          U_real[i][j]=1.0+A[i][j];
        }else{
          R_imag[i][j]=A[i][j]/6.0;
          U_real[i][j]=A[i][j];
        }
      }
    }

    /* step 2: calculate U_real=U_real+R_real*R_imag */
    beta=1.0;
    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&U_real[1][1],&lda,&R_imag[1][1],
           &ldb,&beta,&R_real[1][1],&ldc);

    /* Verify that U_real is indeed a unitary matrix */

    beta=0.0; transb='T';
    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&R_real[1][1],&lda,&U_real[1][1],
           &ldb,&beta,&U_real[1][1],&ldc);
    for (i=1;i<=nstate;i++) R_real[i][i]-=1.0;
    incx=incy=1;
    norm=sqrt(DDOT(&NDIM,&R_real[1][1],&incx,&R_real[1][1],&incy));
    /* only pe=0 is supposed to print out */
    printf("Using Taylor expansion the norm in computeWannierCenters is %g\n",norm);
  }/* endif flag*/

/*==================================================================================*/
/* III) Calculate Matrix Phi */

  transa='N'; transb='N'; mm=nstate; nn=nstate; kk=nstate; lda=nstate;
  ldb=nstate; ldc=nstate; alpha=1.0; beta=0.0; incx=nstate; incy=nstate;

  for (i=1;i<=3;i++){ /* I guess only rectangular cell is allowed for wan center */

    /* step 1. Calculate matrix Z*U */
    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&U_real[1][1],&lda,&Z_real[i][1][1],
           &ldb,&beta,&U_tmp1[1][1],&ldc);
    for (j=1;j<=nstate;j++){
      real[j]=DDOT(&lda,&U_real[1][j],&incx,&U_tmp1[1][j],&incy);
    }

    DGEMM(&transa,&transb,&mm,&nn,&kk,&alpha,&U_real[1][1],&lda,&Z_imag[i][1][1],
           &ldb,&beta,&U_tmp1[1][1],&ldc);

    for (j=1;j<=nstate;j++){
      imag[j]=DDOT(&lda,&U_real[1][j],&incx,&U_tmp1[1][j],&incy);
    }

    /* step 2. compute Phi matrix */
    for (j=1;j<=nstate;j++){
      phi[j][i]= atan2(imag[j],real[j]); /* sign is changed from - to + */
    }
  } /*end for i*/

/*==================================================================================*/
/* IV) Compute Wannier Centers */

  HMatrix[1][1]=BOHR*general_data->cell.hmat_cp[1];
  HMatrix[2][1]=BOHR*general_data->cell.hmat_cp[4];
  HMatrix[3][1]=BOHR*general_data->cell.hmat_cp[7];


  HMatrix[1][2]=BOHR*general_data->cell.hmat_cp[2];
  HMatrix[2][2]=BOHR*general_data->cell.hmat_cp[5];
  HMatrix[3][2]=BOHR*general_data->cell.hmat_cp[8];


  HMatrix[1][3]=BOHR*general_data->cell.hmat_cp[3];
  HMatrix[2][3]=BOHR*general_data->cell.hmat_cp[6];
  HMatrix[3][3]=BOHR*general_data->cell.hmat_cp[9];

  for (n=1;n<=nstate;n++){
    distance=0.0;
    for(i=1;i<=3;i++){
      wan_cent[n][i]=0.0;
      for (k=1;k<=3;k++){
        wan_cent[n][i]+=(HMatrix[i][k]*phi[n][k]/(2.0*pi)); /* in AA*/
      }
    }
  }

/*---------------------------------------------------------------------------*/
/* Return final U matrix                                                     */

  n=0;
  for(i=1;i<=nstate;i++){
    for(j=1;j<=nstate;j++){
      n += 1;
      U_final[n]=U_real[j][i];
    }
  }

/*===========================================================================*/
} /*end routine */
/*============================================================================*/

void write_wannier_orb(GENERAL_DATA *general_data,CP *cp)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/

#include "../typ_defs/typ_mask.h"

  int myid         = cp->communicate.myid_state;
  int num_proc     = cp->communicate.np_states;
  MPI_Comm comm    = cp->communicate.comm_states;

  int iwrite_st       = cp->cpopts.iwrite_init_state;
  int cp_ngrid_skip   = cp->cpopts.cp_ngrid_skip;

  int ncoef           = cp->cpcoeffs_info.nstate_ncoef_proc_max_up;

  int nfft              = cp->cp_para_fft_pkg3d_lg.nfft;
  int nfft_proc         = cp->cp_para_fft_pkg3d_lg.nfft_proc;
  int nfft2             = nfft/2;
  int nfft2_proc        = nfft_proc/2;
  int skc_fft_ka_proc   = cp->cp_para_fft_pkg3d_lg.skc_fft_ka_proc;
  int ekc_fft_ka_proc   = cp->cp_para_fft_pkg3d_lg.ekc_fft_ka_proc;
  int skb_fft_ka_proc   = cp->cp_para_fft_pkg3d_lg.skb_fft_ka_proc;
  int ekb_fft_ka_proc   = cp->cp_para_fft_pkg3d_lg.ekb_fft_ka_proc;

  double *creal         = cp->cpcoeffs_pos[1].cre_up;
  double *cimag         = cp->cpcoeffs_pos[1].cim_up;

  double *zfft          = cp->cpscr.cpscr_wave.zfft;
  double *zfft_tmp      = cp->cpscr.cpscr_wave.zfft_tmp;
  double *cp_wan_orb    = cp->cpscr.cpscr_wave.zfft_tmp; /* not a bug*/
  int *kmax_cp          = cp->cpewald.kmax_cp;

  double *hmat          = general_data->cell.hmat;

  FILE *fp_wan;
  int nkf1,nkf2,nkf3;
  int kb_str,kb_end;
  int i,j,k, iproc, ka,kb,kc,kap,kbp,kcp,ioff,skc_use;
  int skc_fft_ka_proc0,ekc_fft_ka_proc0,skb_fft_ka_proc0,ekb_fft_ka_proc0;
  double sa,sb,sc,da,db,dc,x,y,z,scale_fact;


/*=============================================================================*/
/* I) Transform the specified orbital into real-space */

  ioff=  (iwrite_st-1)*ncoef;

  sngl_pack_coef(&creal[ioff],&cimag[ioff],zfft,&(cp->cp_para_fft_pkg3d_sm));

  para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,&(cp->cp_para_fft_pkg3d_sm));

  for(i=1;i<=nfft2_proc;i++){
    cp_wan_orb[i] = zfft[2*i-1];
  }

/*==========================================================================*/
/* II) Check commensurability                                                   */

  nkf1 = 4*(kmax_cp[1]+1);
  nkf2 = 4*(kmax_cp[2]+1);
  nkf3 = 4*(kmax_cp[3]+1);

  da = 1.0/((double) nkf1);
  db = 1.0/((double) nkf2);
  dc = 1.0/((double) nkf3);

  scale_fact=sqrt((double)(nkf1*nkf2*nkf3));
  scale_fact = 1.0/scale_fact;

  if(((nkf1 % cp_ngrid_skip) != 0) || ((nkf2 % cp_ngrid_skip) != 0) ||
    ((nkf3 % cp_ngrid_skip) != 0)){
    if(myid==0){
      printf("@@@@@@@@@@@@@@@@@@@@-ERROR-@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("The grid skip value, cp_ngrid_skip, must be commensurate\n");
      printf("with the number of points along each direction\n");
      printf("nx = %d, ny = %d, nz = %d, grid skip = %d\n",nkf1,nkf2,nkf3,
               cp_ngrid_skip);
      printf("@@@@@@@@@@@@@@@@@@@@-ERROR-@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    }/* endif myid */
    if(num_proc>1){Barrier(comm);}
    Finalize();
    exit(1);
  }/* endif */

/*==========================================================================*/
/* III) Open Wannier orbital file    */

  if(myid==0){
    fp_wan = cfopen(general_data->filenames.worbname,"w");
  }

/*==========================================================================*/
/* IV) Write orbital ( Inspired from Elf writing routine) */

  if(num_proc >1){

    for(iproc=0;iproc<num_proc;iproc++){
      if(iproc==0){
        if(myid == 0){

          fprintf(fp_wan, "%d  %d  %d\n",
                  nkf1/cp_ngrid_skip,nkf2/cp_ngrid_skip,nkf3/cp_ngrid_skip);

          /* Assumes that skc_fft_ka_proc=1 for the 0th processor */
          i=1;
          for(kc=skc_fft_ka_proc;kc<=ekc_fft_ka_proc;kc+=cp_ngrid_skip){
            kb_str = (kc==skc_fft_ka_proc ? skb_fft_ka_proc : 1);
            kb_end = (kc==ekc_fft_ka_proc ? ekb_fft_ka_proc : nkf2);
            for(kb=kb_str;kb<=kb_end;kb+=cp_ngrid_skip){
              for(ka=1;ka<=nkf1;ka+=cp_ngrid_skip){

                sa = da*((double)(ka-1));
                sb = db*((double)(kb-1));
                sc = dc*((double)(kc-1));
                x  = sa*hmat[1]+sb*hmat[4]+sc*hmat[7];
                y  = sa*hmat[2]+sb*hmat[5]+sc*hmat[8];
                z  = sa*hmat[3]+sb*hmat[6]+sc*hmat[9];

                fprintf(fp_wan,"%.8g  %.8g  %.8g  %.10g\n",z,y,x,cp_wan_orb[i]*scale_fact);
                i += cp_ngrid_skip;
              }/* endfor */
              i += (cp_ngrid_skip-1)*nkf1;
            }/* endfor */
            i += (cp_ngrid_skip-1)*nkf1*nkf2;
          }/* endfor */
          skc_fft_ka_proc0 = skc_fft_ka_proc;
          ekc_fft_ka_proc0 = ekc_fft_ka_proc;
          skb_fft_ka_proc0 = skb_fft_ka_proc;
          ekb_fft_ka_proc0 = ekb_fft_ka_proc;
        }/* endif myid == 0 */
      }/* endif iproc == 0 */

      if(num_proc>1){Barrier(comm);}

      if(iproc != 0){
        if(myid == iproc){
          Ssend(&cp_wan_orb[1],nfft2_proc,MPI_DOUBLE,0,0,comm);
          Ssend(&skc_fft_ka_proc,1,MPI_INT,0,1,comm);
          Ssend(&ekc_fft_ka_proc,1,MPI_INT,0,2,comm);
          Ssend(&skb_fft_ka_proc,1,MPI_INT,0,3,comm);
          Ssend(&ekb_fft_ka_proc,1,MPI_INT,0,4,comm);
        } /* endif myid == iproc */
        if(myid == 0){
          Recv(&cp_wan_orb[1],nfft2_proc,MPI_DOUBLE,MPI_ANY_SOURCE,0,comm);
          Recv(&skc_fft_ka_proc,1,MPI_INT,MPI_ANY_SOURCE,1,comm);
          Recv(&ekc_fft_ka_proc,1,MPI_INT,MPI_ANY_SOURCE,2,comm);
          Recv(&skb_fft_ka_proc,1,MPI_INT,MPI_ANY_SOURCE,3,comm);
          Recv(&ekb_fft_ka_proc,1,MPI_INT,MPI_ANY_SOURCE,4,comm);
          skc_use = (kb>nkf2 ? kc : kc-cp_ngrid_skip);
          i=(i%cp_ngrid_skip)+1;
          for(kc=skc_use;kc<=ekc_fft_ka_proc;kc+=cp_ngrid_skip){
            kb_str = (kc==skc_use ? kb : 1);
            kb_str = (kb>nkf2 ? kb-nkf2 : kb);
            kb_end = (kc==ekc_fft_ka_proc ? ekb_fft_ka_proc : nkf2);
            for(kb=kb_str;kb<=kb_end;kb+=cp_ngrid_skip){
              for(ka=1;ka<=nkf1;ka+=cp_ngrid_skip){

                sa = da*((double)(ka-1)) ;
                sb = db*((double)(kb-1)) ;
                sc = dc*((double)(kc-1)) ;
                x  = sa*hmat[1]+sb*hmat[4]+sc*hmat[7];
                y  = sa*hmat[2]+sb*hmat[5]+sc*hmat[8];
                z  = sa*hmat[3]+sb*hmat[6]+sc*hmat[9];

                fprintf(fp_wan,"%.8g  %.8g  %.8g  %.10g\n",z,y,x,cp_wan_orb[i]*scale_fact);

                i += cp_ngrid_skip;
              }/* endfor */
              i += (cp_ngrid_skip-1)*nkf1;
            }/* endfor */
            i += (cp_ngrid_skip-1)*nkf1*nkf2;
          }/* endfor */
        }/* endif myid == 0 */
      }/* endif iproc != 0 */

      if(num_proc>1){Barrier(comm);}

    }/* endfor iproc */

  }else{ /*serial*/

    fprintf(fp_wan, "%d  %d  %d\n",
                   nkf3/cp_ngrid_skip,nkf2/cp_ngrid_skip,nkf1/cp_ngrid_skip);

    for(kc=1;kc<=nkf3;kc+=cp_ngrid_skip){
      for(kb=1;kb<=nkf2;kb+=cp_ngrid_skip){
        for(ka=1;ka<=nkf1;ka+=cp_ngrid_skip){
          i = (ka-1) + (kb-1)*nkf1 + (kc-1)*nkf1*nkf2 + 1;

          sa = da*((double)(ka-1)) ;
          sb = db*((double)(kb-1)) ;
          sc = dc*((double)(kc-1)) ;
          x  = sa*hmat[1]+sb*hmat[4]+sc*hmat[7];
          y  = sa*hmat[2]+sb*hmat[5]+sc*hmat[8];
          z  = sa*hmat[3]+sb*hmat[6]+sc*hmat[9];

          fprintf(fp_wan,"%.8g  %.8g  %.8g  %.10g\n",z,y,x,cp_wan_orb[i]*scale_fact);

        }/* endfor */
     }/* endfor */
    }/* endfor */

  }/*end if parallel */

/*==========================================================================*/
/* V) close the file and restore the fft parameters for root*/

  if(myid==0){

    fclose(fp_wan);

    if(num_proc > 1){
      cp->cp_para_fft_pkg3d_lg.skc_fft_ka_proc = skc_fft_ka_proc0;
      cp->cp_para_fft_pkg3d_lg.ekc_fft_ka_proc = ekc_fft_ka_proc0;
      cp->cp_para_fft_pkg3d_lg.skb_fft_ka_proc = skb_fft_ka_proc0;
      cp->cp_para_fft_pkg3d_lg.ekb_fft_ka_proc = ekb_fft_ka_proc0;
    }
  }


/*===========================================================================*/
} /*end routine */
/*============================================================================*/


