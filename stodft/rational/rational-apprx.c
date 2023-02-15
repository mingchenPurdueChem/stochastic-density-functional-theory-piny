/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: rational.c                                     */
/*                                                                          */
/* This routine initialize the rational approximation part.                 */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_stodft_local.h"

#include "../proto_defs/proto_rational_elliptic.h"
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void test(double m){
  //m = 0.999701;
  printf("Complete_Elliptic_Integral_First_Kind \n");
  printf("integral: %lg \n", Complete_Elliptic_Integral_First_Kind('m', m));
}
/*===============================================================*/


double complex fun_test(double complex x, double dmu, double beta, double epsilon) {
  double complex arg_raw, fout;
  double arg_raw_re, arg_raw_im, exp_re;
  arg_raw = beta*(x + dmu);
  arg_raw_re = creal(arg_raw);
  if(arg_raw_re < -200.0) arg_raw_re = -200.0;
  if(arg_raw_re > 1.0e10) arg_raw_re =  1.0e10;
  arg_raw_im = cimag(arg_raw);
  exp_re = exp(-arg_raw_re);
  fout = (1.0 + epsilon) * exp_re / ((1.0 + epsilon)*exp_re + cos(arg_raw_im) + sin(arg_raw_im) * I);
  return fout;
}

/*========================================================================================*/
void solve_shifted_eqn_cocg_mu( CP *cp, CLASS *class, GENERAL_DATA *general_data, 
                            double complex *z, double complex *x, double *rhs, int nz, int ndim, int id) { 
/*========================================================================================*/
/* =------------------------------------------------------------------------------------= */
#include "../typ_defs/typ_mask.h"

  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[1]);
  CLATOMS_POS*  clatoms_pos     = &(class->clatoms_pos[1]);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  COMMUNICATE *communicate      = &(cp->communicate);

  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;

  int nfft          = cp_para_fft_pkg3d_lg->nfft;
  int nfft2         = nfft/2;

  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  double *wfUpRe1 = stodftCoefPos->wfUpRe1;
  double *wfUpIm1 = stodftCoefPos->wfUpIm1;
  double *wfDnRe1 = stodftCoefPos->wfDnRe1;
  double *wfDnIm1 = stodftCoefPos->wfDnIm1;
  int spinFlag;

  int iState,iCoeff, iOff, iCoeffStart,index1,index2;

/* =------------------------------------------------------------------------------------= */
int  itermax, iter_old;
double threshold;
double complex z_seed;

double complex *v12, *v2, *r_l;

int status[3];

int i, j, iter, jiter, dim, zdim;

double scale1, scale2;

scale1 = -0.5;
scale2 = -1.0;

zdim = nz;
dim = nfft2; // ndim;
itermax =  9000;
threshold = 1.0E-5;
spinFlag = 0;


v12 = (double complex*)malloc((dim)*sizeof(double complex));
v2 = (double complex*)malloc((dim)*sizeof(double complex));
r_l = (double complex*)malloc((dim)*sizeof(double complex));
/* =------------------- */

printf("--- tessssssssssssst solve_shifted_eqn start ---- \n");

genNoiseOrbitalRealRational(cp,cpcoeffs_pos, v2, id); 


//double sum;
//sum = 0.0;
//
//for (i=0; i < nfft2; i++){
//  sum += creal(v2[i])*creal(v2[i]);
//}
//printf("--- sum here = %lg \n", sum);

for (i=0; i < nfft2; i++){
  rhs[i] = creal(v2[i]);
}


printf("#####  CG Iteration 1 #####\n");

printf("zdim %i dim %i ndim %i \n", zdim, dim, ndim);

printf("#####  CG Iteration 0 #####\n");

komega_cocg_init(&ndim, &dim, &nz, x, z, &itermax, &threshold, NULL);

printf("#####  CG Iteration  #####\n");

for (iter = 0; iter < itermax; iter++){

  for (i = 0; i < dim; i++) {
    r_l[i] = v2[i];
  }

  calcCoefForceWrapSCFReal(class,general_data,cp,cpcoeffs_pos,clatoms_pos,v2,v12,spinFlag);


  komega_cocg_update(v12, v2, x, r_l, status);

  printf(" DEBUG : %i %i %i %i %lg \n", iter, status[0], status[1], status[2], creal(v12[0]));

  if(status[0] < 0) break;

}

switch(status[1]) {
  case (0) :
    printf("  Converged in iteration %d \n", abs(status[0]));
    break;
  case (1) :
    printf("  Not Converged in iteration %d \n", abs(status[0]));
    break;
  case (2) :
    printf("  Alpha becomes infinity %d \n", abs(status[0]));
    break;
  case (3) :
    printf("  Pi_seed becomes zero %d \n", abs(status[0]));
    break;
  case (4) :
    printf("  Residual & Shadow residual are orthogonal %d \n", abs(status[0]));
    break;
}

komega_cocg_finalize();


printf("--- tessssssssssssst solve_shifted_eqn end ---- \n");
}
/*==========================================================================*/
/*========================================================================================*/
/*========================================================================================*/
void solve_shifted_eqn_cocg( CP *cp, CLASS *class, GENERAL_DATA *general_data, 
                            double complex *z, double complex *x, int nz, int ndim, int id) { 
/*========================================================================================*/
/* =------------------------------------------------------------------------------------= */
#include "../typ_defs/typ_mask.h"

  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[1]);
  CLATOMS_POS*  clatoms_pos     = &(class->clatoms_pos[1]);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  COMMUNICATE *communicate      = &(cp->communicate);

  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;

  int nfft          = cp_para_fft_pkg3d_lg->nfft;
  int nfft2         = nfft/2;

  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  double *wfUpRe1 = stodftCoefPos->wfUpRe1;
  double *wfUpIm1 = stodftCoefPos->wfUpIm1;
  double *wfDnRe1 = stodftCoefPos->wfDnRe1;
  double *wfDnIm1 = stodftCoefPos->wfDnIm1;
  int spinFlag;

  int iState,iCoeff, iOff, iCoeffStart,index1,index2;

/* =------------------------------------------------------------------------------------= */
int  itermax, iter_old;
double threshold;
double complex z_seed;

double complex *v12, *v2, *r_l;

int status[3];

int i, j, iter, jiter, dim, zdim;

double scale1, scale2;

scale1 = -0.5;
scale2 = -1.0;

zdim = nz;
dim = nfft2; // ndim;
itermax =  9000;
threshold = 1.0E-5;
spinFlag = 0;


v12 = (double complex*)malloc((dim)*sizeof(double complex));
v2 = (double complex*)malloc((dim)*sizeof(double complex));
r_l = (double complex*)malloc((dim)*sizeof(double complex));
/* =------------------- */

printf("--- tessssssssssssst solve_shifted_eqn start ---- \n");

genNoiseOrbitalRealRational(cp,cpcoeffs_pos, v2, id); 

printf("#####  CG Iteration 1 #####\n");

printf("zdim %i dim %i ndim %i \n", zdim, dim, ndim);

printf("#####  CG Iteration 0 #####\n");

komega_cocg_init(&ndim, &dim, &nz, x, z, &itermax, &threshold, NULL);

printf("#####  CG Iteration  #####\n");

for (iter = 0; iter < itermax; iter++){

  for (i = 0; i < dim; i++) {
    r_l[i] = v2[i];
  }

  calcCoefForceWrapSCFReal(class,general_data,cp,cpcoeffs_pos,clatoms_pos,v2,v12,spinFlag);


  komega_cocg_update(v12, v2, x, r_l, status);

  printf(" DEBUG : %i %i %i %i %lg \n", iter, status[0], status[1], status[2], creal(v12[0]));

  if(status[0] < 0) break;

}

switch(status[1]) {
  case (0) :
    printf("  Converged in iteration %d \n", abs(status[0]));
    break;
  case (1) :
    printf("  Not Converged in iteration %d \n", abs(status[0]));
    break;
  case (2) :
    printf("  Alpha becomes infinity %d \n", abs(status[0]));
    break;
  case (3) :
    printf("  Pi_seed becomes zero %d \n", abs(status[0]));
    break;
  case (4) :
    printf("  Residual & Shadow residual are orthogonal %d \n", abs(status[0]));
    break;
}

komega_cocg_finalize();


printf("--- tessssssssssssst solve_shifted_eqn end ---- \n");
}
/*==========================================================================*/
/*==========================================================================*/



/*===============================================================*/
void solve_shifted_eqn( CP *cp, CLASS *class, GENERAL_DATA *general_data, double complex *z, double complex *x, int nz, int ndim, int sign) { 
//                      int *ndim, int *nz, double complex *ham, double complex *rhs, double complex *z, double complex *x  ) {

/* =------------------- */
#include "../typ_defs/typ_mask.h"

  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[1]);
  CLATOMS_POS*  clatoms_pos     = &(class->clatoms_pos[1]);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  COMMUNICATE *communicate      = &(cp->communicate);


  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;

  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  double *wfUpRe1 = stodftCoefPos->wfUpRe1;
  double *wfUpIm1 = stodftCoefPos->wfUpIm1;
  double *wfDnRe1 = stodftCoefPos->wfDnRe1;
  double *wfDnIm1 = stodftCoefPos->wfDnIm1;


  int iState,iCoeff, iOff, iCoeffStart,index1,index2;

/* =------------------- */
int  itermax, iter_old;
double threshold;
double complex z_seed;

double complex *v12, *v2, *r_l;
double complex *v14, *v4;

int status[3];

int i, j, iter, jiter, dim, zdim;

double scale1, scale2;

scale1 = -0.5;
scale2 = -1.0;
//scale1 = -0.25;
//scale2 = -0.50;


// ndim = 2*(numCoeff - 1) + 1
zdim = nz;
dim = ndim;
itermax = 990000;
threshold = 1.0E-5;

v12 = (double complex*)malloc((dim)*sizeof(double complex));
v2 = (double complex*)malloc((dim)*sizeof(double complex));
v14 = (double complex*)malloc((dim)*sizeof(double complex));
v4 = (double complex*)malloc((dim)*sizeof(double complex));
r_l = (double complex*)malloc((dim)*sizeof(double complex));
/* =------------------- */

//z[0] = 0.0 + 0.0*I;
//for (i = 0; i < nz; i++ ) {
//  z[i] = creal(z[i]) + 0.0*I;
//  printf("Z VALUES are %lg, %lg\n", creal(z[i]), cimag(z[i]));
//}


  printf("--- tessssssssssssst solve_shifted_eqn start ---- \n");

genNoiseOrbitalReal(cp,cpcoeffs_pos); 


printf("#####  CG Iteration 1 #####\n");

if(sign == 1){

  for (i = 1; i < numCoeff; i++ ) {
    v2[i-1] = cre_up[i] + cim_up[i] * I;
    v2[i-1 + numCoeff] = cre_up[i] - cim_up[i]*I;
    //
    v4[i-1] = cre_up[i] - cim_up[i] * I;
    v4[i-1 + numCoeff] = cre_up[i] + cim_up[i]*I;
  }
  v2[numCoeff-1] = cre_up[numCoeff] + cim_up[numCoeff] * I;
  //
  v4[numCoeff-1] = cre_up[numCoeff] - cim_up[numCoeff] * I;

}
else {

  for (i = 1; i < numCoeff; i++ ) {
    v2[i-1] = conj( cre_up[i] + cim_up[i] * I);
    v2[i-1 + numCoeff] = conj(cre_up[i] - cim_up[i]*I);
    //
    v4[i-1] = conj(cre_up[i] - cim_up[i] * I);
    v4[i-1 + numCoeff] = conj(cre_up[i] + cim_up[i]*I);
  }
  v2[numCoeff-1] = conj(cre_up[numCoeff] + cim_up[numCoeff] * I);
  //
  v4[numCoeff-1] = conj(cre_up[numCoeff] - cim_up[numCoeff] * I);

}

//print before for checking
//    for (i = 1; i < numCoeff; i++) {
//      printf("before %i %lg %lg %lg %lg \n", i, creal(v2[i-1]), cimag(v2[i-1]), creal(v2[i-1 + numCoeff]), cimag(v2[i-1 + numCoeff]));
//    }

printf("#####  CG Iteration 0 #####\n");

komega_bicg_init(&ndim, &dim, &nz, x, z, &itermax, &threshold, NULL);

printf("#####  CG Iteration  #####\n");

for (iter = 0; iter < itermax; iter++){

  for (i = 0; i < dim; i++) {
    r_l[i] = v2[i];
  }
/********************************/
    for (i = 1; i < numCoeff; i++) {
      cre_up[i]  = 0.5* creal(v2[i-1] + v2[i-1+numCoeff] );
      cim_up[i+numCoeff]  = 0.5* creal(-v2[i-1] + v2[i-1+numCoeff] ); 
      cre_up[i+numCoeff]  = 0.5* cimag(v2[i-1] + v2[i-1+numCoeff] ); 
      cim_up[i]  = 0.5* cimag(v2[i-1] - v2[i-1+numCoeff] );
    }
    cre_up[numCoeff]  = creal(v2[numCoeff-1]);
    cim_up[numCoeff]  = 0.0;
    cre_up[2*numCoeff]  = cimag(v2[numCoeff-1]);
    cim_up[2*numCoeff]  = 0.0;
/********************************/
    calcCoefForceWrapSCF(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

    for (i = 1; i < numCoeff; i++) {
      fcre_up[i] *= scale1;
      fcim_up[i] *= scale1;
    }
    fcre_up[numCoeff] *= scale2;
    fcim_up[numCoeff] = 0.0;
    for (i = 1; i < numCoeff; i++) {
      fcre_up[numCoeff+i] *= scale1;
      fcim_up[numCoeff+i] *= scale1;
    }
    fcre_up[2*numCoeff] *= scale2;
    fcim_up[2*numCoeff] = 0.0;
    
    for (i = 1; i < numCoeff; i++) {
      v12[i-1] = fcre_up[i]-fcim_up[i+numCoeff] + (fcim_up[i] + fcre_up[i+numCoeff])* I;
      v12[i-1 + numCoeff] = fcre_up[i] + fcim_up[i+numCoeff] + (-fcim_up[i] + fcre_up[i+numCoeff])*I;
    }
    v12[numCoeff-1] = fcre_up[numCoeff] + fcre_up[2*numCoeff] * I;

/********************************/
    for (i = 1; i < numCoeff; i++) {
      cre_up[i]  = 0.5* creal(v4[i-1] + v4[i-1+numCoeff] );
      cim_up[i+numCoeff]  = 0.5* creal(-v4[i-1] + v4[i-1+numCoeff] ); 
      cre_up[i+numCoeff]  = 0.5* cimag(v4[i-1] + v4[i-1+numCoeff] ); 
      cim_up[i]  = 0.5* cimag(v4[i-1] - v4[i-1+numCoeff] );
    }
    cre_up[numCoeff]  = creal(v4[numCoeff-1]);
    cim_up[numCoeff]  = 0.0;
    cre_up[2*numCoeff]  = cimag(v4[numCoeff-1]);
    cim_up[2*numCoeff]  = 0.0;

    /*
    for (i = 1; i < numCoeff; i++) {
      cre_up[i]  = creal(v4[i-1]);
      cim_up[i]  = cimag(v4[i-1]);
    }
    cre_up[numCoeff]  = creal(v4[numCoeff-1]);
    cim_up[numCoeff]  = cimag(v4[numCoeff-1]);
    */
/********************************/
    calcCoefForceWrapSCF(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

    for (i = 1; i < numCoeff; i++) {
      fcre_up[i] *= scale1;
      fcim_up[i] *= scale1;
    }
    fcre_up[numCoeff] *= scale2;
    fcim_up[numCoeff] = 0.0;
    for (i = 1; i < numCoeff; i++) {
      fcre_up[numCoeff+i] *= scale1;
      fcim_up[numCoeff+i] *= scale1;
    }
    fcre_up[2*numCoeff] *= scale2;
    fcim_up[2*numCoeff] = 0.0;
    
    for (i = 1; i < numCoeff; i++) {
      v14[i-1] = fcre_up[i]-fcim_up[i+numCoeff] + (fcim_up[i] + fcre_up[i+numCoeff])* I;
      v14[i-1 + numCoeff] = fcre_up[i] + fcim_up[i+numCoeff] + (-fcim_up[i] + fcre_up[i+numCoeff])*I;
    }
    v14[numCoeff-1] = fcre_up[numCoeff] + fcre_up[2*numCoeff] * I;


    /*
    for (i = 1; i < numCoeff; i++) {
      v14[i-1] = scale1*fcre_up[i] + scale1*fcim_up[i] * I;
      v14[i-1 + numCoeff] = scale1*fcre_up[i] - scale1*fcim_up[i]*I;
    }
    v14[numCoeff-1] = scale2*fcre_up[numCoeff] + scale2*fcim_up[numCoeff] * I;
    */

/********************************/
  komega_bicg_update(v12, v2, v14, v4, x, r_l, status);
//printf("v2 %i real %lg imag %lg , real %lg imag %lg \n", iter, creal(v2[1]), cimag(v2[1]), creal(v2[100]), cimag(v2[100]));

/*for (i = 1; i < numCoeff; i++ ) {
  printf("v2 %i real %lg imag %lg , real %lg imag %lg \n", iter, creal(v2[i-1]), cimag(v2[i-1]), creal(v2[i-1+numCoeff]), cimag(v2[i-1+numCoeff]));
  printf("v4 %i real %lg imag %lg , real %lg imag %lg \n", iter, creal(v4[i-1]), cimag(v4[i-1]), creal(v4[i-1+numCoeff]), cimag(v4[i-1+numCoeff]));
}*/

  printf(" DEBUG : %i %i %i %i %lg \n", iter, status[0], status[1], status[2], creal(v12[0]));

  if(status[0] < 0) break;

}

switch(status[1]) {
  case (0) :
    printf("  Converged in iteration %d \n", abs(status[0]));
    break;
  case (1) :
    printf("  Not Converged in iteration %d \n", abs(status[0]));
    break;
  case (2) :
    printf("  Alpha becomes infinity %d \n", abs(status[0]));
    break;
  case (3) :
    printf("  Pi_seed becomes zero %d \n", abs(status[0]));
    break;
  case (4) :
    printf("  Residual & Shadow residual are orthogonal %d \n", abs(status[0]));
    break;
}

komega_bicg_finalize();

/*
double complex **x_tmp;

x_tmp = (double complex**)malloc((2*nz)*sizeof(double complex*));

for (int i=0; i< nz; i++){
  x_tmp[i] = (double complex*)malloc((numCoeff)*sizeof(double complex));
}

// memcpy(x_tmp, x, ndim*nz*sizeof(double complex));

for  (int j =0; j < ndim*nz; j++){
  x_tmp[j] = x[j];
}

int id;
for (int j =0; j < nz; j++){
  //for (int i = 0; i < numCoeff; i++){
   // sum_p_p = sum_p_p + fun_p[j] * x_p_p[j*ndim + i];
    id = j*ndim;
    for (i = 1; i < numCoeff; i++) {
      //cre_up[i]  = 0.5* creal(v2[i-1] + v2[i-1+numCoeff] );
      //cim_up[i+numCoeff]  = 0.5* creal(-v2[i-1] + v2[i-1+numCoeff] ); 
      //cre_up[i+numCoeff]  = 0.5* cimag(v2[i-1] + v2[i-1+numCoeff] ); 
      //cim_up[i]  = 0.5* cimag(v2[i-1] - v2[i-1+numCoeff] );
      x[i-1 + id] = 0.5* creal(x_tmp[i-1 + id] + x_tmp[i-1+numCoeff + id] ) + 0.5* cimag(x_tmp[i-1 + id] - x_tmp[i-1+numCoeff + id] )*I;
      x[i-1 + numCoeff + id] = 0.5* cimag(x_tmp[i-1 + id] + x_tmp[i-1+numCoeff + id] ) + 0.5* creal(-x_tmp[i-1 + id] + x_tmp[i-1+numCoeff + id] )*I;
    }
    x[numCoeff-1 + id] = creal(x_tmp[numCoeff-1 + id]) + 0.0*I;
    x[2*numCoeff-1 + id] = cimag(x_tmp[numCoeff-1 +id]) + 0.0*I;
    //cre_up[numCoeff]  = creal(v2[numCoeff-1]);
    //cim_up[numCoeff]  = 0.0;
    //cre_up[2*numCoeff]  = cimag(v2[numCoeff-1]);
    //cim_up[2*numCoeff]  = 0.0;

  //}
}
*/


/*
//cross check the solution
    for (i = 1; i < numCoeff; i++) {
      cre_up[i]  = 0.5* creal(x[i-1] + x[i-1+numCoeff] );
      cim_up[i+numCoeff]  = 0.5* creal(-x[i-1] + x[i-1+numCoeff] ); 
      cre_up[i+numCoeff]  = 0.5* cimag(x[i-1] + x[i-1+numCoeff] ); 
      cim_up[i]  = 0.5* cimag(x[i-1] - x[i-1+numCoeff] );
    }
    cre_up[numCoeff]  = creal(x[numCoeff-1]);
    cim_up[numCoeff]  = 0.0;
    cre_up[2*numCoeff]  = cimag(x[numCoeff-1]);
    cim_up[2*numCoeff]  = 0.0;

    calcCoefForceWrapSCF(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

    for (i = 1; i < numCoeff; i++) {
      fcre_up[i] *= scale1;
      fcim_up[i] *= scale1;
    }
    fcre_up[numCoeff] *= scale2;
    fcim_up[numCoeff] = 0.0;
    for (i = 1; i < numCoeff; i++) {
      fcre_up[numCoeff+i] *= scale1;
      fcim_up[numCoeff+i] *= scale1;
    }
    fcre_up[2*numCoeff] *= scale2;
    fcim_up[2*numCoeff] = 0.0;
    
    for (i = 1; i < numCoeff; i++) {
      v12[i-1] = fcre_up[i]-fcim_up[i+numCoeff] + (fcim_up[i] + fcre_up[i+numCoeff])* I;
      v12[i-1 + numCoeff] = fcre_up[i] + fcim_up[i+numCoeff] + (-fcim_up[i] + fcre_up[i+numCoeff])*I;
    }
    v12[numCoeff-1] = fcre_up[numCoeff] + fcre_up[2*numCoeff] * I;

genNoiseOrbitalReal(cp,cpcoeffs_pos); 

if(sign == 1){

  for (i = 1; i < numCoeff; i++ ) {
    v2[i-1] = cre_up[i] + cim_up[i] * I;
    v2[i-1 + numCoeff] = cre_up[i] - cim_up[i]*I;
    //
    v4[i-1] = cre_up[i] - cim_up[i] * I;
    v4[i-1 + numCoeff] = cre_up[i] + cim_up[i]*I;
  }
  v2[numCoeff-1] = cre_up[numCoeff] + cim_up[numCoeff] * I;
  //
  v4[numCoeff-1] = cre_up[numCoeff] - cim_up[numCoeff] * I;

}
else {

  for (i = 1; i < numCoeff; i++ ) {
    v2[i-1] = conj( cre_up[i] + cim_up[i] * I);
    v2[i-1 + numCoeff] = conj(cre_up[i] - cim_up[i]*I);
    //
    v4[i-1] = conj(cre_up[i] - cim_up[i] * I);
    v4[i-1 + numCoeff] = conj(cre_up[i] + cim_up[i]*I);
  }
  v2[numCoeff-1] = conj(cre_up[numCoeff] + cim_up[numCoeff] * I);
  //
  v4[numCoeff-1] = conj(cre_up[numCoeff] - cim_up[numCoeff] * I);

}

//print before for checking
    for (i = 1; i < numCoeff; i++) {
      printf("before %i %lg %lg %lg %lg \n", i, creal(v2[i-1]), cimag(v2[i-1]), creal(v2[i-1 + numCoeff]), cimag(v2[i-1 + numCoeff]));
      //printf("after %i %lg %lg %lg %lg \n", i, creal(v12[i-1]), cimag(v12[i-1]), creal(v12[i-1 + numCoeff]), cimag(v12[i-1 + numCoeff]));
      printf("after %i %lg %lg %lg %lg \n", i, creal(v12[i-1]-z[0]*x[i-1]), cimag(v12[i-1]-z[0]*x[i-1]), creal(v12[i-1 + numCoeff]-z[0]*x[i-1+numCoeff]), cimag(v12[i-1 + numCoeff]-z[0]*x[i-1+numCoeff]));
    }
    //v12[numCoeff-1] = fcre_up[numCoeff] + fcre_up[2*numCoeff] * I;
*/

/*
if(communicate->myid_state  == 0){
FILE *fp;
fp = fopen("out-sto-orb.dat", "w");

  for(iState=0; iState<numStateUpProc; iState++){
    iOff = iState*numCoeff;
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      fprintf(fp, "%.8f %.8f \n",  cre_up[iOff + iCoeff], cim_up[iOff + iCoeff]);
    }
  }

}
*/
/*==========================================================================*/
/* 1) Calculate the H(norm)|phi> */
/*
  calcCoefForceWrapSCF(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

      memcpy(&cre_up[1],&fcre_up[1],numCoeffUpTotal*sizeof(double));
      memcpy(&cim_up[1],&fcim_up[1],numCoeffUpTotal*sizeof(double));

  calcCoefForceWrapSCF(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

      memcpy(&cre_up[1],&fcre_up[1],numCoeffUpTotal*sizeof(double));
      memcpy(&cim_up[1],&fcim_up[1],numCoeffUpTotal*sizeof(double));

if(communicate->myid_state  == 0){
FILE *fp1;
fp1 = fopen("out-sto-orb1.dat", "w");

  for(iState=0; iState<numStateUpProc; iState++){
    iOff = iState*numCoeff;
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      fprintf(fp1, "%.8f %.8f \n",  cre_up[iOff + iCoeff], cim_up[iOff + iCoeff]);
    }
  }

}
*/
  printf("--- tessssssssssssst solve_shifted_eqn end ---- \n");
/*
int  itermax, iter_old;
double threshold;
double complex z_seed;

double complex *v12, *v2, *r_l;

int status[3];

int i, j, iter, jiter, dim, zdim;

zdim = *nz;
dim = *ndim;
itermax = 5000;
threshold = 1.0E-5;

v12 = (double complex*)malloc((dim)*sizeof(double complex));
v2 = (double complex*)malloc((dim)*sizeof(double complex));
r_l = (double complex*)malloc((dim)*sizeof(double complex));


for (i = 0; i < dim; i++ ) {
  v2[i] = rhs[i];
}


komega_cocg_init(ndim, &dim, nz, x, z, &itermax, &threshold, NULL);


printf("#####  CG Iteration starts  #####\n");
for (iter = 0; iter < itermax; iter++){

  for (i = 0; i < dim; i++) {
    r_l[i] = v2[i];
  }

  for (i = 0; i < dim; i++){
    v12[i] = 0.0 + 0.0*I;
    for (j = 0; j < dim; j++) {
      v12[i] = v12[i] + ham[i + dim * j] * v2[j];
    }
  }

  komega_cocg_update(v12, v2, x, r_l, status);

  printf(" DEBUG : %i %i %i %i %lg \n", iter, status[0], status[1], status[2], creal(v12[0]));

  if(status[0] < 0) break;

} // for i
printf("#####  CG Iteration finished  #####\n");

switch(status[1]) {
  case (0) :
    printf("  Converged in iteration %d \n", abs(status[0]));
    break;
  case (1) :
    printf("  Not Converged in iteration %d \n", abs(status[0]));
    break;
  case (2) :
    printf("  Alpha becomes infinity %d \n", abs(status[0]));
    break;
  case (3) :
    printf("  Pi_seed becomes zero %d \n", abs(status[0]));
    break;
  case (4) :
    printf("  Residual & Shadow residual are orthogonal %d \n", abs(status[0]));
    break;
}

komega_cocg_finalize();
*/
}
/*===============================================================*/
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void filterRational(CP *cp,CLASS *class,GENERAL_DATA *general_data,
			  int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  #include "../typ_defs/typ_mask.h"

  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CLATOMS_POS*  clatoms_pos     = &(class->clatoms_pos[ip_now]);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  COMMUNICATE *communicate      = &(cp->communicate);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm = &(cp->cp_sclr_fft_pkg3d_sm);

  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

  NEWTONINFO *newtonInfo = stodftInfo->newtonInfo;

  int expanType	     = stodftInfo->expanType;
  int polynormLength = stodftInfo->polynormLength;
  int numChemPot     = stodftInfo->numChemPot;
  int smearOpt       = stodftInfo->smearOpt;
  int filterDiagFlag = stodftInfo->filterDiagFlag;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int cpLsda	     = cpopts->cp_lsda;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int myidState       = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int numThreads = cp_sclr_fft_pkg3d_sm->numThreads;
  int pseudoRealFlag = cp->pseudo.pseudoReal.pseudoRealFlag;
  int imu,iCoeff,iPoly,indexStart,iState;
  int iOff2;
  int startIndex;
  int storeChebyMomentsFlag = stodftInfo->storeChebyMomentsFlag;
  MPI_Comm comm_states   =    communicate->comm_states;


  double energyDiff  = stodftInfo->energyDiff;
  double energyMin   = stodftInfo->energyMin;
  double energyMax   = stodftInfo->energyMax;
  double energyMean  = stodftInfo->energyMean;
  double scale       = newtonInfo->scale;
  double polyCoeff;
  double timeProc,timeTot;
  double *sampPoint = (double*)newtonInfo->sampPoint;
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;
  double *entropyUpRe = stodftCoefPos->entropyUpRe;
  double *entropyUpIm = stodftCoefPos->entropyUpIm;
  double *entropyDnRe = stodftCoefPos->entropyDnRe;
  double *entropyDnIm = stodftCoefPos->entropyDnIm;
  double *entropyExpanCoeff = stodftCoefPos->entropyExpanCoeff;
  double *chebyMomentsUp = stodftCoefPos->chebyMomentsUp;
  double *chebyMomentsDn = stodftCoefPos->chebyMomentsDn;
  double *wfUpRe1,*wfUpIm1;
  double *wfDnRe1,*wfDnIm1;
  double *wfUpRe0,*wfUpIm0;
  double *wfDnRe0,*wfDnIm0;


  double *expanCoeff = (double*)stodftCoefPos->expanCoeff;
  double timeStart,timeEnd; 
  double timeStart2,timeEnd2;
  double timeStart3,timeEnd3;
  double deltaTime = 0.0;
  double deltaTime2 = 0.0;
  double dot;


  int nfft          = cp_para_fft_pkg3d_lg->nfft;
  int nfft2         = nfft/2;
  

/*****************************************************/
// define variables, move to a different place later !!
//
  //FERMIFUNC fermiFunction;
  //fermiFunction = &fermiExpComplex;
  double K, K_prim;
  double mA, MA, ktmp, k, kinv, m, dmu, epsilon;
  int ntgrid;
  //double dx;
  //int ngrid;
  double energymax;
  double complex sum;
/*****************************************************/
 
//if(myidState==0){ 
  printf("from RAAAAAAAAAAAAAAAATIONALLLLLLL \n");
  printf("Energy Max= %.16lg.\n", stodftInfo->energyMax);
  printf("Energy Min= %.16lg.\n", stodftInfo->energyMin);
  //stodftInfo->beta = 6000.0;
  printf("Beta= %.16lg.\n", stodftInfo->beta);
  //stodftInfo->chemPotTrue = 0.3177870; // 0.22182526; // 0.23124281;//0.243;
  printf("Correct Chemical Potential = %.16lg.\n", stodftInfo->chemPotTrue);


  energymax = stodftInfo->energyMax -  stodftInfo->chemPotTrue;
  mA = (M_PI*M_PI)/(stodftInfo->beta*stodftInfo->beta);
  MA = energymax*energymax + mA;
  //MA = stodftInfo->energyMax*stodftInfo->energyMax + mA;
  ktmp = sqrt(MA/mA);
  k = (ktmp - 1.0)/(ktmp + 1.0);
  kinv = 1.0/k;
  m = k*k;
  dmu = 0.0005;
  epsilon = 0.001;

  ntgrid=100; // 100; // 0;
//  ngrid = 10000;
  
//  dx = (stodftInfo->energyMax - stodftInfo->energyMin)/ngrid;

  printf(" mA = %lg , MA = %lg \n", mA, MA);
  printf(" k = %lg , kinv =  %lg , m =  %lg \n", k, kinv, m);

  if (m > 0.9){
    K = Complete_Elliptic_Integral_First_Kind('m', m);
  }
  else {
    K = Complete_Elliptic_Integral_First_Kind('m', m);
  }

  K_prim = Complete_Elliptic_Integral_First_Kind('m', 1 - m);

  printf("m= %lg, K= %lg, K_prim = %lg \n", m, K, K_prim);

  printf("from RAAAAAAAAAAAAAAAATIONALLLLLLL \n");


  printf("from RAAAAAAAAAAAAAAAATIONALLLLLLL \n");
/*****************************************************/
double complex *tgrid;

tgrid = (double complex*)malloc((ntgrid)*sizeof(double complex));

double *tgrid_re, *tgrid_im;
double *sn_re, *cn_re, *dn_re;
double *sn_tmp, *cn_tmp, *dn_tmp;
double *cn_im, *dn_im;

double complex *sn_im;
double complex *sn_c, *cn_c, *dn_c;


tgrid_re = (double*)malloc((ntgrid)*sizeof(double));
tgrid_im = (double*)malloc((ntgrid)*sizeof(double));
sn_re = (double*)malloc((ntgrid)*sizeof(double));
cn_re = (double*)malloc((ntgrid)*sizeof(double));
dn_re = (double*)malloc((ntgrid)*sizeof(double));
sn_tmp = (double*)malloc((ntgrid)*sizeof(double));
cn_tmp = (double*)malloc((ntgrid)*sizeof(double));
dn_tmp = (double*)malloc((ntgrid)*sizeof(double));
cn_im = (double*)malloc((ntgrid)*sizeof(double));
dn_im = (double*)malloc((ntgrid)*sizeof(double));

sn_im = (double complex*)malloc((ntgrid)*sizeof(double complex));
sn_c = (double complex*)malloc((ntgrid)*sizeof(double complex));
cn_c = (double complex*)malloc((ntgrid)*sizeof(double complex));
dn_c = (double complex*)malloc((ntgrid)*sizeof(double complex));

double complex *z;
double *z_re, *z_im;

z = (double complex*)malloc((ntgrid)*sizeof(double complex));
z_re = (double*)malloc((ntgrid)*sizeof(double));
z_im = (double*)malloc((ntgrid)*sizeof(double));

double complex *ksi_p, *ksi_m;
//double *ksi_p_re, *ksi_p_im;
//double *ksi_m_re, *ksi_m_im;

ksi_p = (double complex*)malloc((ntgrid)*sizeof(double complex));
ksi_m = (double complex*)malloc((ntgrid)*sizeof(double complex));
//ksi_p_re = (double*)malloc((ntgrid)*sizeof(double));
//ksi_p_im = (double*)malloc((ntgrid)*sizeof(double));
//ksi_m_re = (double*)malloc((ntgrid)*sizeof(double));
//ksi_m_im = (double*)malloc((ntgrid)*sizeof(double));

double complex *fun_p, *fun_m;
fun_p = (double complex*)malloc((ntgrid)*sizeof(double complex));
fun_m = (double complex*)malloc((ntgrid)*sizeof(double complex));

/*****************************************************/
for (int i = 0; i < ntgrid; i++){
  tgrid[i] = -K + 2.0*(i + 0.5)*K/ntgrid + K_prim*0.5 * I ; // 0.0 + (double) i; 
}

for (int i = 0; i < ntgrid; i++){
  tgrid_re[i] = creal(tgrid[i]);
  tgrid_im[i] = cimag(tgrid[i]);
}

for (int i = 0; i < ntgrid; i++){
  Jacobi_sn_cn_dn(tgrid_re[i], 'm', m, &sn_re[i], &cn_re[i], &dn_re[i] );
  Jacobi_sn_cn_dn(tgrid_im[i], 'm', 1-m, &sn_tmp[i], &cn_tmp[i], &dn_tmp[i] );
  sn_im[i] = sn_tmp[i]/cn_tmp[i] * I;
  cn_im[i] = 1 / cn_tmp[i];
  dn_im[i] = dn_tmp[i] / cn_tmp[i];
}

for (int i = 0; i < ntgrid; i++){
  sn_c[i] = (sn_re[i]*sn_re[i] - sn_im[i]*sn_im[i]) / (sn_re[i]*cn_im[i]*dn_im[i] - sn_im[i]*cn_re[i]*dn_re[i]);
  cn_c[i] = (sn_re[i]*cn_re[i]*dn_im[i] - sn_im[i]*cn_im[i]*dn_re[i])/(sn_re[i]*cn_im[i]*dn_im[i] - sn_im[i]*cn_re[i]*dn_re[i]);
  dn_c[i] = (sn_re[i]*cn_im[i]*dn_re[i] - sn_im[i]*cn_re[i]*dn_im[i])/(sn_re[i]*cn_im[i]*dn_im[i] - sn_im[i]*cn_re[i]*dn_re[i]);
}

for (int i = 0; i < ntgrid; i++){
  z[i] = sqrt(mA*MA) * (kinv + sn_c[i]) / (kinv - sn_c[i]) ;
  z_re[i] = creal(z[i]);
  z_im[i] = cimag(z[i]);
}

for (int i = 0; i < ntgrid; i++){
  ksi_p[i] = csqrt(z[i] - mA);
  ksi_m[i] = - csqrt(z[i] - mA);
  //ksi_p_re[i] = creal(ksi_p[i]);
  //ksi_p_im[i] = cimag(ksi_p[i]);
  //ksi_m_re[i] = creal(ksi_m[i]);
  //ksi_m_im[i] = cimag(ksi_m[i]);
}

for (int i = 0; i < ntgrid; i++){
    
  //fun_p[i] = fermiFunction(ksi_p[i], stodftInfo->chemPotTrue, stodftInfo->beta);
  //fun_m[i] = fermiFunction(ksi_m[i], stodftInfo->chemPotTrue, stodftInfo->beta);
  fun_p[i] = fun_test(ksi_p[i], dmu, stodftInfo->beta, epsilon);
  fun_m[i] = fun_test(ksi_m[i], dmu, stodftInfo->beta, epsilon);
  fun_p[i] = fun_p[i]*cn_c[i]*dn_c[i] / ((kinv - sn_c[i])*(kinv - sn_c[i])) / ksi_p[i];
  fun_m[i] = fun_m[i]*cn_c[i]*dn_c[i] / ((kinv - sn_c[i])*(kinv - sn_c[i])) / ksi_m[i];
}

/******************************************************************************/
double complex *zseed, *x;
int ndim, nl, nz;

ndim = nfft2; // 2*(numCoeff-1) + 1; //TODO change here
nl =ndim;
//nz = ntgrid;
nz = 2*ntgrid;


zseed = (double complex*)malloc((nz)*sizeof(double complex));

x = (double complex*)malloc((nl*nz)*sizeof(double complex));
/******************************************************************************/


for (int j =0; j < ntgrid; j++){
  zseed[j] =  ksi_p[j]; 
  zseed[j + ntgrid] =  ksi_m[j];
 
  //z_p[j] =  ksi_p[j] + stodftInfo->chemPotTrue; //TODO check
  //z_m[j] =  ksi_m[j] + stodftInfo->chemPotTrue;
}


//Barrier(comm_states);
//printf("@@@@@@@@@@@@@@@@@@@@_forced_stop__@@@@@@@@@@@@@@@@@@@@\n");
//fflush(stdout);
//exit(1);

//for (int i = 0; i < nz; i++ ) {
//  printf("BEFORE Z VALUES are %lg, %lg\n", creal(zseed[i]), cimag(zseed[i]));
//}


double *frhs;

frhs = (double*)malloc((ndim)*sizeof(double));

double preRat = -2 * K *sqrt(mA*MA) / M_PI / ntgrid * kinv;

printf("before preRat = %lg %lg %lg %lg %lg %lg %i \n", preRat, K, mA, MA, M_PI, kinv, ntgrid);
//printf("preRat %lg\n",preRat);
//printf("before preRat %lg %lg %lg \n",preRat, creal(zseed[1]), cimag(zseed[1]));


printf("numStateUpProc %i %i %i %i \n", numStateUpProc, myidState, numProcStates, numThreads);


/******************************************************************************/
printf("Start Calculating ChemicalPotential with Rational Approximation\n");
/******************************************************************************/

calcChemPotRational(cp, class, general_data, ip_now);

Barrier(comm_states);

/******************************************************************************/
printf("Done Calculating ChemicalPotential with Rational Approximation\n");
/******************************************************************************/



/******************************************************************************/
printf("Start Filtering ChemicalPotential with Rational Approximation\n");
/******************************************************************************/
for(iState=0;iState<numStateUpProc;iState++){

  printf("myidState %i %i %i \n", myidState, numStateUpProc, iState);
  solve_shifted_eqn_cocg( cp, class, general_data, zseed, x, nz, ndim, iState);
  
  for (int i = 0; i<ndim;i++){
    sum = 0.0 + 0.0 *I ;
    for (int j =0; j < ntgrid; j++){
      sum = sum + fun_p[j] * x[j*ndim + i] + fun_m[j] * x[(j+ntgrid)*ndim + i];
    }
    frhs[i] = preRat*cimag(sum);
  }
  
  printf("==== final results === \n");
  
  for (int j =1; j <= numCoeff; j++){
    printf("actuai %i %i  %.8f  %.8f \n", myidState, j,  stoWfUpRe[0][j + (iState*numCoeff)],  stoWfUpIm[0][j + (iState*numCoeff)] );
  }
  
  rhsReal(class, general_data, cp, cpcoeffs_pos, clatoms_pos, frhs, iState);

}
/******************************************************************************/
/******************************************************************************/
printf("Done Filtering ChemicalPotential with Rational Approximation\n");
/******************************************************************************/

Barrier(comm_states);

printf("@@@@@@@@@@@@@@@@@@@@_forced_stop__@@@@@@@@@@@@@@@@@@@@\n");
fflush(stdout);
exit(1);

/*****************************************************/
/*****************************************************/

  // performance
  stodftInfo->cputime1 = 0.0;
  stodftInfo->cputime2 = 0.0;
  stodftInfo->cputime3 = 0.0;
  stodftInfo->cputime4 = 0.0;
  stodftInfo->cputime5 = 0.0;
  stodftInfo->cputime6 = 0.0;
  stodftInfo->cputime7 = 0.0;
  cp_sclr_fft_pkg3d_sm->cputime = 0.0;
  for(imu=0;imu<100;imu++)stodftInfo->cputime_new[imu] = 0.0;

  omp_set_num_threads(numThreads);


/*==========================================================================*/
/* 0) Allocate temp memories */

  stodftCoefPos->wfUpRe1 = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  stodftCoefPos->wfUpIm1 = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  wfUpRe1 = stodftCoefPos->wfUpRe1;
  wfUpIm1 = stodftCoefPos->wfUpIm1;
  if(cpLsda==1&&numStateDnProc!=0){
    stodftCoefPos->wfDnRe1 = (double*)cmalloc(numCoeffDnTotal*sizeof(double))-1;
    stodftCoefPos->wfDnIm1 = (double*)cmalloc(numCoeffDnTotal*sizeof(double))-1;
    wfDnRe1 = stodftCoefPos->wfDnRe1;
    wfDnIm1 = stodftCoefPos->wfDnIm1;
  }

  
/*==========================================================================*/
/* 1) Copy the initial stochastic orbital */

  //debug 
  //storeChebyMomentsFlag = 1;

  if(storeChebyMomentsFlag==1){
    printf("%p\n",stodftCoefPos->chebyMomentsUp);
    stodftCoefPos->chebyMomentsUp = (double*)cmalloc((polynormLength+1)*sizeof(double));
    //stodftCoefPos->chebyMomentsUp = (double*)crealloc(stodftCoefPos->chebyMomentsUp,
    //                                               (polynormLength+1)*sizeof(double));
    chebyMomentsUp = stodftCoefPos->chebyMomentsUp;
    if(cpLsda==1&&numStateDnProc!=0){
      stodftCoefPos->chebyMomentsDn = (double*)crealloc(stodftCoefPos->chebyMomentsDn,
                                                   (polynormLength+1)*sizeof(double));
      chebyMomentsDn = stodftCoefPos->chebyMomentsDn;
    }
    //stodftCoefPos->chebyMomentsUp = (double*)cmalloc((polynormLength+1)*sizeof(double));
    for(iPoly=0;iPoly<polynormLength;iPoly++)chebyMomentsUp[iPoly] = 0.0;
    if(cpLsda==1&&numStateDnProc!=0){
      for(iPoly=0;iPoly<polynormLength;iPoly++)chebyMomentsDn[iPoly] = 0.0;
    }
    stodftCoefPos->wfUpRe0 = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
    stodftCoefPos->wfUpIm0 = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
    wfUpRe0 = stodftCoefPos->wfUpRe0;
    wfUpIm0 = stodftCoefPos->wfUpIm0;
    if(cpLsda==1&&numStateDnProc!=0){
      stodftCoefPos->wfDnRe1 = (double*)cmalloc(numCoeffDnTotal*sizeof(double))-1;
      stodftCoefPos->wfDnIm1 = (double*)cmalloc(numCoeffDnTotal*sizeof(double))-1;
      wfDnRe0 = stodftCoefPos->wfDnRe0;
      wfDnIm0 = stodftCoefPos->wfDnIm0;
    }
    memcpy(&wfUpRe0[1],&cre_up[1],numCoeffUpTotal*sizeof(double));
    memcpy(&wfUpIm0[1],&cim_up[1],numCoeffUpTotal*sizeof(double));
    if(cpLsda==1&&numStateDnProc!=0){
      memcpy(&wfDnRe0[1],&cre_dn[1],numCoeffDnTotal*sizeof(double));
      memcpy(&wfDnIm0[1],&cim_dn[1],numCoeffDnTotal*sizeof(double));
    }
    for(iState=0;iState<numStateUpProc;iState++){
      iOff2 = iState*numCoeff;
      dot = 2.0*ddotBlasWrapperThreads(numCoeff,&wfUpRe0[iOff2+1],1,
            &wfUpRe0[iOff2+1],1,numThreads)+
            2.0*ddotBlasWrapperThreads(numCoeff,&wfUpIm0[iOff2+1],1,
            &wfUpIm0[iOff2+1],1,numThreads)-
            wfUpRe0[iOff2+numCoeff]*wfUpRe0[iOff2+numCoeff];
      chebyMomentsUp[0] += dot;
    }
    if(cpLsda==1&&numStateDnProc!=0){
      for(iState=0;iState<numStateDnProc;iState++){
        iOff2 = iState*numCoeff;
        dot = 2.0*ddotBlasWrapperThreads(numCoeff,&cre_dn[iOff2+1],1,
              &cre_dn[iOff2+1],1,numThreads)+
              2.0*ddotBlasWrapperThreads(numCoeff,&cim_dn[iOff2+1],1,
              &cim_dn[iOff2+1],1,numThreads)-
              cre_dn[iOff2+numCoeff]*cre_dn[iOff2+numCoeff];
        chebyMomentsDn[0] += dot;
      }
    }
  }
 
  //printf("1111111 cheby 0 %lg\n",chebyMomentsUp[0]);

  for(imu=0;imu<numChemPot;imu++){
    #pragma omp parallel for private(iCoeff)
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      stoWfUpRe[imu][iCoeff] = expanCoeff[imu]*cre_up[iCoeff];
      stoWfUpIm[imu][iCoeff] = expanCoeff[imu]*cim_up[iCoeff];
    }//endfor iCoeff
    if(cpLsda==1&&numStateDnProc!=0){
      #pragma omp parallel for private(iCoeff)
      for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
	stoWfDnRe[imu][iCoeff] = expanCoeff[imu]*cre_dn[iCoeff];
	stoWfDnIm[imu][iCoeff] = expanCoeff[imu]*cim_dn[iCoeff];
      }//endfor iCoeff      
    }//endif 
  }//endfor imu

  if(smearOpt>0&&filterDiagFlag==0){
    #pragma omp parallel for private(iCoeff)
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      entropyUpRe[iCoeff] = entropyExpanCoeff[0]*cre_up[iCoeff];
      entropyUpIm[iCoeff] = entropyExpanCoeff[0]*cim_up[iCoeff];
      //if(isnormal(entropyUpRe[iCoeff])==0||isnormal(entropyUpIm[iCoeff])==0){
      //  printf("bbbbbb %i %lg %lg %lg\n",iCoeff,entropyUpRe[iCoeff],entropyUpIm[iCoeff],entropyExpanCoeff[imu]);
      //}

    }//endfor iCoeff
    if(cpLsda==1&&numStateDnProc!=0){
      #pragma omp parallel for private(iCoeff)
      for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
        entropyDnRe[iCoeff] = entropyExpanCoeff[0]*cre_dn[iCoeff];
        entropyDnIm[iCoeff] = entropyExpanCoeff[0]*cim_dn[iCoeff];
      }//endfor iCoeff
    }//endif cpLsda
  }//endif smearOpt

  memcpy(&wfUpRe1[1],&cre_up[1],numCoeffUpTotal*sizeof(double));
  memcpy(&wfUpIm1[1],&cim_up[1],numCoeffUpTotal*sizeof(double));
  if(cpLsda==1&&numStateDnProc!=0){
    memcpy(&wfDnRe1[1],&cre_dn[1],numCoeffDnTotal*sizeof(double));
    memcpy(&wfDnIm1[1],&cim_dn[1],numCoeffDnTotal*sizeof(double));
  }   

/*==========================================================================*/
/* 1) Loop over all polynomial terms (iPoly=0<=>polynomial order 1) */

  timeStart = omp_get_wtime();

  for(iPoly=1;iPoly<polynormLength;iPoly++){
    if(iPoly%1000==0&&myidState==0){
      printf("%lg%% ",iPoly*100.0/polynormLength);
      fflush(stdout);
    }
    timeStart3 = omp_get_wtime();  
    normHCheby(cp,class,general_data,
                 cpcoeffs_pos,clatoms_pos,iPoly,1);
    timeEnd3 = omp_get_wtime();
    deltaTime2 += timeEnd3-timeStart3;
    timeStart2 = omp_get_wtime();

    if(storeChebyMomentsFlag==1){
      for(iState=0;iState<numStateUpProc;iState++){
        iOff2 = iState*numCoeff;
        dot = 2.0*ddotBlasWrapperThreads(numCoeff,&wfUpRe0[iOff2+1],1,
              &cre_up[iOff2+1],1,numThreads)+
              2.0*ddotBlasWrapperThreads(numCoeff,&wfUpIm0[iOff2+1],1,
              &cim_up[iOff2+1],1,numThreads)-
              wfUpRe0[iOff2+numCoeff]*cre_up[iOff2+numCoeff];
        chebyMomentsUp[iPoly] += dot;
      }
      //printf("1111111 cheby %i %lg\n",iPoly,chebyMomentsUp[iPoly]);
      if(cpLsda==1&&numStateDnProc!=0){
        for(iState=0;iState<numStateDnProc;iState++){
          iOff2 = iState*numCoeff;
          dot = 2.0*ddotBlasWrapperThreads(numCoeff,&wfDnRe0[iOff2+1],1,
                &cre_dn[iOff2+1],1,numThreads)+
                2.0*ddotBlasWrapperThreads(numCoeff,&wfDnIm0[iOff2+1],1,
                &cim_dn[iOff2+1],1,numThreads)-
                wfDnRe0[iOff2+numCoeff]*cre_dn[iOff2+numCoeff];
          chebyMomentsDn[iPoly] += dot;
        }
      }
    }

    for(imu=0;imu<numChemPot;imu++){
      polyCoeff = expanCoeff[iPoly*numChemPot+imu];
      //omp_set_num_threads(numThreads);
      #pragma omp parallel for private(iCoeff)
      for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
	stoWfUpRe[imu][iCoeff] += polyCoeff*cre_up[iCoeff];	
        stoWfUpIm[imu][iCoeff] += polyCoeff*cim_up[iCoeff];                       

      }//endfor iCoeff
      if(cpLsda==1&&numStateDnProc!=0){
	#pragma omp parallel for private(iCoeff)
        for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
	  stoWfDnRe[imu][iCoeff] += polyCoeff*cre_dn[iCoeff];                     
	  stoWfDnIm[imu][iCoeff] += polyCoeff*cim_dn[iCoeff];
        }//endfor iCoeff        
      }//endif 
    }//endfor imu
    if(smearOpt>0&&filterDiagFlag==0){ 
      polyCoeff = entropyExpanCoeff[iPoly];
      #pragma omp parallel for private(iCoeff)
      for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
        entropyUpRe[iCoeff] += polyCoeff*cre_up[iCoeff];
        entropyUpIm[iCoeff] += polyCoeff*cim_up[iCoeff];
        //if(isnormal(entropyUpRe[iCoeff])==0||isnormal(entropyUpIm[iCoeff])==0){
        //  printf("bbbbbb %i %lg %lg %lg\n",iCoeff,entropyUpRe[iCoeff],entropyUpIm[iCoeff],polyCoeff);
        //}
      }
      if(cpLsda==1&&numStateDnProc!=0){
        #pragma omp parallel for private(iCoeff)
        for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
          entropyDnRe[iCoeff] += polyCoeff*cre_dn[iCoeff];
          entropyDnIm[iCoeff] += polyCoeff*cim_dn[iCoeff];         
        }
      }
    }
    timeEnd2 = omp_get_wtime();
    deltaTime += timeEnd2-timeStart2;
  }//endfor iPoly
  timeEnd = omp_get_wtime();
  if(myidState==0)printf("\n");
  timeProc = timeEnd-timeStart;
  if(numProcStates>1)Reduce(&timeProc,&timeTot,1,MPI_DOUBLE,MPI_SUM,0,comm_states);
  else timeTot = timeProc;
  /*
  if(myidState==0){
    if(pseudoRealFlag==0){
      printf("Average Filter time is %lg\n",timeTot/numProcStates);
      printf("0th process filter time is %lg\n",timeProc);
      printf("non-local pp time %.8lg\n",stodftInfo->cputime1);
      printf("calc force(fft) time %.8lg\n",stodftInfo->cputime2);
      printf("FFTW3D time %.8lg\n",cp_sclr_fft_pkg3d_sm->cputime);
      printf("non-local pp matrix %.8lg\n",stodftInfo->cputime3);
      printf("non-local pp energy %.8lg\n",stodftInfo->cputime4);
      printf("non-local pp coef force %.8lg\n",stodftInfo->cputime5);
    }
    else{
      printf("Average Filter time is %lg\n",timeTot/numProcStates);
      //printf("0th process filter time is %lg\n",timeProc);
      printf("Nlpp part 1 %.8lg\n",stodftInfo->cputime0);
      printf("Nlpp part 2 %.8lg\n",stodftInfo->cputime1);
      printf("Apply KS pot %.8lg\n",stodftInfo->cputime2);
      printf("Pack fft %.8lg\n",stodftInfo->cputime3);
      printf("Unpack fft %.8lg\n",stodftInfo->cputime4);
      printf("Kinetic %.8lg\n",stodftInfo->cputime5);
      printf("FFTW3D time %.8lg\n",cp_sclr_fft_pkg3d_sm->cputime);
      printf("FFTW3D to r pre time %.8lg\n",cp_sclr_fft_pkg3d_sm->cputime1);
      printf("FFTW3D to r post time %.8lg\n",cp_sclr_fft_pkg3d_sm->cputime2);
      printf("FFTW3D to g pre time %.8lg\n",cp_sclr_fft_pkg3d_sm->cputime3);
      printf("FFTW3D to g post time %.8lg\n",cp_sclr_fft_pkg3d_sm->cputime4);
      printf("Scale H|phi> time %.8lg\n",stodftInfo->cputime7);
      printf("Accumulate P(H)|phi> time %.8lg\n",deltaTime);
    }
  }
  */
  if(numProcStates>1)Barrier(comm_states);

  if(myidState==0){
    printf("Process ID %i filter time %.8lg total NormH time %.8lg Accumulate P(H)|phi> time %.8lg\n",myidState,timeProc,deltaTime2,deltaTime);
    printf("Process ID %i Nlpp %.8lg Apply-KS-pot %.8lg Pack-fft %.8lg Unpack-fft %.8lg kinetic %.8lg scale-H|phi> %.8lg\n",myidState,stodftInfo->cputime_new[0],stodftInfo->cputime3,stodftInfo->cputime2,stodftInfo->cputime4,stodftInfo->cputime5,stodftInfo->cputime7);
    printf("Process ID Nlpp-part1 %.8lg Nlpp-part2 %.8lg\n",stodftInfo->cputime0,stodftInfo->cputime1);
    printf("Process ID %i FFTW3D time %.8lg FFTW3D-to-r-pre %.8lg FFTW3D-to-r-post %.8lg FFTW3D-to-g-pre %.8lg FFTW3D-to-g-post %.8lg\n",myidState,cp_sclr_fft_pkg3d_sm->cputime,cp_sclr_fft_pkg3d_sm->cputime1,cp_sclr_fft_pkg3d_sm->cputime2,cp_sclr_fft_pkg3d_sm->cputime3,cp_sclr_fft_pkg3d_sm->cputime4);
    printf("Process ID %i cpy-wf %.8lg %.8lg add force %.8lg %.8lg prepare-loc %.8lg zero %.8lg\n",myidState,stodftInfo->cputime_new[1],stodftInfo->cputime_new[3],stodftInfo->cputime_new[2],stodftInfo->cputime_new[4],stodftInfo->cputime_new[5],stodftInfo->cputime_new[6]);
   fflush(stdout);
  }

  if(numProcStates>1)Barrier(comm_states);

  //debug
  /*
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    printf("stowftesttest %.16lg %.16lg\n",stoWfUpRe[0][iCoeff],stoWfUpIm[0][iCoeff]);
  }
  */

  free(&wfUpRe1[1]);
  free(&wfUpIm1[1]);
  if(cpLsda==1&&numStateDnProc!=0){
    free(&wfDnRe1[1]);
    free(&wfDnIm1[1]);    
  }
  if(storeChebyMomentsFlag==1){
    free(&wfUpRe0[1]);
    free(&wfUpIm0[1]);
    if(cpLsda==1&&numStateDnProc!=0){
      free(&wfDnRe0[1]);
      free(&wfDnIm0[1]);    
    }
  }

/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcChemPotRational(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                          int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This routine first generate all chebyshev momentum. The Chebyshev	 */
/* is half of the newton polynoimial length (so make the newton one even */
/* ). Then calculate all chebyshev moments. Then do optimization	 */
/* Attention: after calculating all chebyshev moments, check whether the */
/* initial chem pots DO give you the > and < # of electrons		 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CLATOMS_POS*  clatoms_pos     = &(class->clatoms_pos[ip_now]);
  COMMUNICATE *communicate      = &(cp->communicate);

  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;

  int iPoly,iState,iChem;
  int iScf = stodftInfo->iScf;
  int polynormLength = stodftInfo->polynormLength;
  int numChebyMoments = (polynormLength%2==0)?(polynormLength/2+1):((polynormLength+1)/2);
  int numFFTGridMutpl = 32;
  int numChebyGridInit = polynormLength*numFFTGridMutpl;
  int numChebyGrid;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int cpLsda         = cpopts->cp_lsda;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int numProcStates   = communicate->np_states;
  int myidState	      = communicate->myid_state;
  int numStateStoUp = stodftInfo->numStateStoUp;
  int numStateStoDn = stodftInfo->numStateStoDn;
  int printChebyMoment = stodftInfo->printChebyMoment;
  MPI_Comm comm_states   =    communicate->comm_states;

  double chemPotDiff = 1000.0;
  double numElecTrue = stodftInfo->numElecTrue;
  double numElecTol = 1.0e-11*numElecTrue;
  double chemPotMin,chemPotMax;
  double chemPotInit = stodftInfo->chemPotInit;
  double gapInit = stodftInfo->gapInit;
  double chemPotNew,chemPotOld;
  double numElecMin,numElecMax;
  double numElecNew,numElecOld;
  double Smin = chebyshevInfo->Smin;
  double Smax = chebyshevInfo->Smax;
  double energyDiff = stodftInfo->energyDiff;
  
  double *chebyCoeffs = (double*)cmalloc(polynormLength*sizeof(double));
  double *chebyMomentsTemp;
  double *chemPot = stodftCoefPos->chemPot;
  double *chebyMomentsUp,*chebyMomentsDn;

  fftw_complex *chebyCoeffsFFT,*funValGridFFT;


  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  int nfft          = cp_para_fft_pkg3d_lg->nfft;
  int nfft2         = nfft/2;


  int ntgrid;
  ntgrid=100;
  double complex *zseed;
  zseed = (double complex*)malloc((2*ntgrid)*sizeof(double complex));
  double preRat;

double complex *fun_p, *fun_m;
fun_p = (double complex*)malloc((ntgrid)*sizeof(double complex));
fun_m = (double complex*)malloc((ntgrid)*sizeof(double complex));


init_zseed(cp, zseed, ntgrid, &preRat, fun_p, fun_m);
printf("preRat %lg %lg %lg \n", preRat, creal(zseed[1]), cimag(zseed[1]));


/******************************************************************************/
printf("Start ChemicalPotential with Rational Approximation\n");
/******************************************************************************/
double complex *x;
x = (double complex*)malloc((nfft2*2*ntgrid)*sizeof(double complex));
double complex sum;

double *frhs;
double *rhs;

double dsum;
double NElecTot;
NElecTot = 0.0;
dsum = 0.0;

frhs = (double*)malloc((nfft2)*sizeof(double));
rhs = (double*)malloc((nfft2)*sizeof(double));

printf(" rhoRealGridNum %i %i \n", stodftInfo->rhoRealGridNum, stodftInfo->rhoRealGridTot);

for(iState=0;iState<numStateUpProc;iState++){

  printf("myidState %i %i %i %i %i \n", myidState, numStateUpProc, iState, numProcStates, numStateStoUp);
  solve_shifted_eqn_cocg_mu( cp, class, general_data, zseed, x, rhs, 2*ntgrid, nfft2, iState);

  for (int i = 0; i<nfft2;i++){
    sum = 0.0 + 0.0 *I ;
    for (int j =0; j < ntgrid; j++){
      sum = sum + fun_p[j] * x[j*nfft2 + i] + fun_m[j] * x[(j+ntgrid)*nfft2 + i];
    }
    frhs[i] = preRat*cimag(sum);
    dsum = dsum + 2.0*frhs[i]*rhs[i];
  }

  //printf("==== final results === \n");

 // for (int j =1; j <= numCoeff; j++){
 //   printf("actuai %i %i  %.8f  %.8f \n", myidState, j,  stoWfUpRe[0][j + (iState*numCoeff)],  stoWfUpIm[0][j + (iState*numCoeff)] );
 // }

  //rhsReal(class, general_data, cp, cpcoeffs_pos, clatoms_pos, frhs, iState);

}

printf("b==== final results === %lg  \n", dsum*(1.0/stodftInfo->rhoRealGridTot));
dsum = (dsum/numStateStoUp)*(1.0/stodftInfo->rhoRealGridTot);
printf("a==== final results === %lg  \n", dsum);

if(numProcStates>1)Reduce(&dsum,&NElecTot,1,MPI_DOUBLE,MPI_SUM,0,comm_states);

printf("==== final results === %lg %lg \n", dsum, NElecTot);

/******************************************************************************/
/******************************************************************************/
printf("Done ChemicalPotential with Rational Approximation\n");
/******************************************************************************/



printf("@@@@@@@@@@@@@@@@@@@@_forced_stop__@@@@@@@@@@@@@@@@@@@@\n");
fflush(stdout);
exit(1);

 
/*==========================================================================*/
/* I) Allocate memories */
   
  chebyshevInfo->scale = (Smax-Smin)/energyDiff;

  stodftInfo->numChebyMoments = numChebyMoments;
  stodftCoefPos->chebyMomentsUp = (double*)cmalloc((polynormLength+1)*sizeof(double));
  chebyMomentsUp = stodftCoefPos->chebyMomentsUp;
  if(cpLsda==1&&numStateDnProc!=0){
    stodftCoefPos->chebyMomentsDn = (double*)cmalloc((polynormLength+1)*sizeof(double));
    chebyMomentsDn = stodftCoefPos->chebyMomentsDn;
  }
  if(myidState==0){
    numChebyGrid = roundFFT(numChebyGridInit);
    stodftInfo->numChebyGrid = numChebyGrid;
    stodftCoefPos->chebyCoeffsFFT = fftw_malloc(numChebyGrid*sizeof(fftw_complex));
    stodftCoefPos->funValGridFFT = fftw_malloc(numChebyGrid*sizeof(fftw_complex));
    chebyCoeffsFFT = stodftCoefPos->chebyCoeffsFFT;
    funValGridFFT = stodftCoefPos->funValGridFFT;
    stodftInfo->fftwPlanForward = fftw_plan_dft_1d(numChebyGrid,funValGridFFT,chebyCoeffsFFT,
                                    FFTW_FORWARD,FFTW_MEASURE);
  }

  //debug
  /*
  stodftInfo->numElecStoWf = (double*)cmalloc((polynormLength+1)*numStateUpProc*sizeof(double));
  double *numElecStoWf = stodftInfo->numElecStoWf;
  for(iPoly=0;iPoly<=polynormLength;iPoly++){
    for(iState=0;iState<numStateUpProc;iState++){
      numElecStoWf[iPoly*numStateUpProc+iState] = 0.0;
    }
  }
  */

/*==========================================================================*/
/* II) Calculate Chebyshev Moments */

  if(myidState==0)printf("Start Calculating Chebyshev Moments\n");
#ifdef FAST_FILTER   
  calcChebyMomentsFake(cp,class,general_data,ip_now);
#else
  calcChebyMoments(cp,class,general_data,ip_now);
#endif
  if(numProcStates>1){
    if(myidState==0){
      chebyMomentsTemp = (double*)cmalloc((polynormLength+1)*sizeof(double));
      for(iPoly=0;iPoly<=polynormLength;iPoly++)chebyMomentsTemp[iPoly] = 0.0;
    }
    Reduce(chebyMomentsUp,chebyMomentsTemp,polynormLength+1,
	    MPI_DOUBLE,MPI_SUM,0,comm_states);
    if(myidState==0){
      memcpy(chebyMomentsUp,chebyMomentsTemp,(polynormLength+1)*sizeof(double));
    }
    if(cpLsda==1&&numStateDnProc!=0){
      if(myidState==0){
	for(iPoly=0;iPoly<=polynormLength;iPoly++)chebyMomentsTemp[iPoly] = 0.0;
      }
      Reduce(chebyMomentsDn,chebyMomentsTemp,polynormLength+1,
             MPI_DOUBLE,MPI_SUM,0,comm_states);
      if(myidState==0){
	memcpy(chebyMomentsDn,chebyMomentsTemp,(polynormLength+1)*sizeof(double));
      }//endif
    }//endif
  }//endif
  if(myidState==0)printf("Finish Calculating Chebyshev Moments\n");
  if(numProcStates>1)Barrier(comm_states);


/*==========================================================================*/
/* II) Solve N(mu)=Ne */

  if(myidState==0){
    printf("Start Calculating Chemical Potential\n");
    for(iPoly=0;iPoly<=polynormLength;iPoly++){
      //printf("chebyMomentsUp %i %lg\n",iPoly,chebyMomentsUp[iPoly]);
      chebyMomentsUp[iPoly] /= (double)numStateStoUp;
      //printf("chebyMomentsUp %i %lg\n",iPoly,chebyMomentsUp[iPoly]);
    }
    //printf("chebyMomentsUp %lg %lg %lg\n",chebyMomentsUp[0],chebyMomentsUp[1],chebyMomentsUp[2]);
    if(cpLsda==1&&numStateDnProc!=0){
      for(iPoly=0;iPoly<=polynormLength;iPoly++){
	chebyMomentsDn[iPoly] /= (double)numStateStoDn;
      }      
    }
    
    if(stodftInfo->printChebyMoment==1){
      FILE *filecheby = fopen("cheby-moment-output","w");
      for(iPoly=0;iPoly<=polynormLength;iPoly++){
        fprintf(filecheby,"%.16lg\n",chebyMomentsUp[iPoly]);
      }
      fclose(filecheby);
    }
    //exit(0);
    
    /*
    FILE *filecheby = fopen("cheby-moment","r");
    for(iPoly=0;iPoly<=polynormLength;iPoly++){
      fscanf(filecheby,"%lg",&chebyMomentsUp[iPoly]);
    }
    fclose(filecheby);
    */
    
    //debug
    /*
    double numChemTest = 1000;
    double deltChemPot = gapInit/numChemTest;
    int iChem;
    for(iChem=0;iChem<numChemTest;iChem++){
      chemPotNew = chemPotMin+(iChem+0.5)*deltChemPot;
      numElecNew = calcNumElecCheby(cp,chemPotNew,chebyCoeffs);
      printf("Scaleeee %.16lg %.16lg\n",chemPotNew,numElecNew);
    }
    //exit(0);
    */
    
    printf("2222222 Ne test %.16lg\n",calcNumElecCheby(cp,0.2472244378617802,chebyCoeffs));
    chemPotMin = chemPotInit-gapInit*0.5;
    chemPotMax = chemPotInit+gapInit*0.5;
    numElecMin = calcNumElecCheby(cp,chemPotMin,chebyCoeffs);
    numElecMax = calcNumElecCheby(cp,chemPotMax,chebyCoeffs);
    printf("numEmin %.16lg numEmax %.16lg\n",numElecMin,numElecMax);
    while(numElecMax<numElecTrue){
      printf("numEmin %lg numEmax %lg\n",numElecMin,numElecMax);
      chemPotMin = chemPotMax;
      chemPotMax += gapInit;
      numElecMin = calcNumElecCheby(cp,chemPotMin,chebyCoeffs);
      numElecMax = calcNumElecCheby(cp,chemPotMax,chebyCoeffs);
    }
    while(numElecMin>numElecTrue){
      printf("numEmin %lg numEmax %lg\n",numElecMin,numElecMax);
      chemPotMax = chemPotMin;
      chemPotMin -= gapInit;
      numElecMin = calcNumElecCheby(cp,chemPotMin,chebyCoeffs);
      numElecMax = calcNumElecCheby(cp,chemPotMax,chebyCoeffs);
    }
    printf("numEmin %lg numEmax %lg\n",numElecMin,numElecMax);
    chemPotNew = (numElecTrue-numElecMin)*(chemPotMax-chemPotMin)/(numElecMax-numElecMin)+
		  chemPotMin;  
    numElecNew = calcNumElecCheby(cp,chemPotNew,chebyCoeffs); 
    while(fabs(numElecNew-numElecTrue)>numElecTol){
      if(numElecNew>numElecTrue){
	chemPotMax = chemPotNew;
	numElecMax = numElecNew;
      }
      if(numElecNew<numElecTrue){
	chemPotMin = chemPotNew;
	numElecMin = numElecNew;
      }
      /*
      chemPotNew = (numElecTrue-numElecMin)*(chemPotMax-chemPotMin)/(numElecMax-numElecMin)+
		    chemPotMin;
      */
      chemPotNew = 0.5*(chemPotMin+chemPotMax);
      numElecNew = calcNumElecCheby(cp,chemPotNew,chebyCoeffs);
      //printf("chemPotNew %lg numElecNew %lg\n",chemPotNew,numElecNew);
    }//endwhile
    printf("Finish Calculating Chemical Potential\n");
    printf("The correct chemical potential is %.16lg Ne %.16lg DNe %.16lg\n",chemPotNew,numElecNew,
            fabs(numElecNew-numElecTrue));


    //DEBUG
    /*
    chemPotNew = 0.2223443349358951;

    printf("Finish Calculating Chemical Potential\n");
    printf("The correct chemical potential is %.16lg Ne %.16lg DNe %.16lg\n",chemPotNew,numElecNew,
	    fabs(numElecNew-numElecTrue));
    */
    //chemPotMin = chemPotInit-gapInit*0.5;
    //chemPotMax = chemPotInit+gapInit*0.5;
    //test DOS
    /*
    double chemPot1 = chemPotNew-0.0734;
    double chemPot2 = chemPotNew+0.0734;
    int numChemPotTest = 100;
    double dmu = (chemPot2-chemPot1)/numChemPotTest;
    double *numElecMu = (double*)cmalloc(numChemPotTest*sizeof(double));
    double x,dos;
    double dmuInv = 0.5/dmu;
    for(iChem=0;iChem<numChemPotTest;iChem++){
      x = chemPot1+iChem*dmu;
      numElecMu[iChem] = calcNumElecCheby(cp,x,chebyCoeffs);
      printf("%i Nemuuuu %.10lg %.10lg\n",iScf,x,numElecMu[iChem]);
    }
    for(iChem=1;iChem<numChemPotTest-1;iChem++){
      x = chemPot1+iChem*dmu;
      dos = (numElecMu[iChem+1]-numElecMu[iChem-1])*dmuInv;
      printf("%i dossssss %.10lg %.10lg\n",iScf,x,dos);
    }
    free(&numElecMu[0]);
    */
  }
  if(numProcStates>1){
    Barrier(comm_states);
    Bcast(&chemPotNew,1,MPI_DOUBLE,0,comm_states);
  }
  
  chemPot[0] = chemPotNew;
  stodftInfo->chemPotTrue = chemPotNew; // another backup
  if(myidState==0){
    printf("Correct Chemical Potential is %.16lg\n",chemPotNew);
    fflush(stdout);
  }

  free(chebyCoeffs);
  free(stodftCoefPos->chebyMomentsUp);
  if(cpLsda==1&&numStateDnProc!=0)free(stodftCoefPos->chebyMomentsDn);
  if(numProcStates>1&&myidState==0)free(chebyMomentsTemp);
  //exit(0);

  if(myidState==0){
    fftw_destroy_plan(stodftInfo->fftwPlanForward);
    fftw_free(chebyCoeffsFFT);
    fftw_free(funValGridFFT);
  }
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*========================================================================================*/
/*========================================================================================*/
void init_zseed( CP *cp, double complex *zseed, int ntgrid, double *preRat, double complex *fun_p, double complex *fun_m ) {

#include "../typ_defs/typ_mask.h"


STODFTINFO *stodftInfo        = cp->stodftInfo;

/*         Local Variable declarations                                   */
double K, K_prim;
double mA, MA, ktmp, k, kinv, m;
double dmu, epsilon;
//int ntgrid;
double energymax;

printf("from init_zseed \n");
printf("Energy Max= %.16lg.\n", stodftInfo->energyMax);
printf("Energy Min= %.16lg.\n", stodftInfo->energyMin);
//stodftInfo->beta = 6000.0;
printf("Beta= %.16lg.\n", stodftInfo->beta);
stodftInfo->chemPotTrue = 0.26615; // 0.22182526; // 0.23124281;//0.243;
printf("Correct Chemical Potential = %.16lg.\n", stodftInfo->chemPotTrue);


  energymax = stodftInfo->energyMax -  stodftInfo->chemPotTrue;
  mA = (M_PI*M_PI)/(stodftInfo->beta*stodftInfo->beta);
  MA = energymax*energymax + mA;
//MA = stodftInfo->energyMax*stodftInfo->energyMax + mA;
  ktmp = sqrt(MA/mA);
  k = (ktmp - 1.0)/(ktmp + 1.0);
  kinv = 1.0/k;
  m = k*k;
  dmu = 0.0005;
  epsilon = 0.001;


//  ntgrid=100; // 100; // 0;

  printf(" mA = %lg , MA = %lg \n", mA, MA);
  printf(" k = %lg , kinv =  %lg , m =  %lg \n", k, kinv, m);

  if (m > 0.9){
    K = Complete_Elliptic_Integral_First_Kind('m', m);
  }
  else {
    K = Complete_Elliptic_Integral_First_Kind('m', m);
  }

  K_prim = Complete_Elliptic_Integral_First_Kind('m', 1 - m);

  printf("m= %lg, K= %lg, K_prim = %lg \n", m, K, K_prim);

  printf("from init_zseed %i \n", ntgrid);


/*****************************************************/
double complex *tgrid;

double *tgrid_re, *tgrid_im;
double *sn_re, *cn_re, *dn_re;
double *sn_tmp, *cn_tmp, *dn_tmp;
double *cn_im, *dn_im;

double complex *sn_im;
double complex *sn_c, *cn_c, *dn_c;

double complex *z;
double *z_re, *z_im;

double complex *ksi_p, *ksi_m;

tgrid = (double complex*)malloc((ntgrid)*sizeof(double complex));
tgrid_re = (double*)malloc((ntgrid)*sizeof(double));
tgrid_im = (double*)malloc((ntgrid)*sizeof(double));
sn_re = (double*)malloc((ntgrid)*sizeof(double));
cn_re = (double*)malloc((ntgrid)*sizeof(double));
dn_re = (double*)malloc((ntgrid)*sizeof(double));
sn_tmp = (double*)malloc((ntgrid)*sizeof(double));
cn_tmp = (double*)malloc((ntgrid)*sizeof(double));
dn_tmp = (double*)malloc((ntgrid)*sizeof(double));
cn_im = (double*)malloc((ntgrid)*sizeof(double));
dn_im = (double*)malloc((ntgrid)*sizeof(double));

sn_im = (double complex*)malloc((ntgrid)*sizeof(double complex));
sn_c = (double complex*)malloc((ntgrid)*sizeof(double complex));
cn_c = (double complex*)malloc((ntgrid)*sizeof(double complex));
dn_c = (double complex*)malloc((ntgrid)*sizeof(double complex));


z = (double complex*)malloc((ntgrid)*sizeof(double complex));
z_re = (double*)malloc((ntgrid)*sizeof(double));
z_im = (double*)malloc((ntgrid)*sizeof(double));


ksi_p = (double complex*)malloc((ntgrid)*sizeof(double complex));
ksi_m = (double complex*)malloc((ntgrid)*sizeof(double complex));

/*****************************************************/

for (int i = 0; i < ntgrid; i++){
  tgrid[i] = -K + 2.0*(i + 0.5)*K/ntgrid + K_prim*0.5 * I ; // 0.0 + (double) i; 
}

for (int i = 0; i < ntgrid; i++){
  tgrid_re[i] = creal(tgrid[i]);
  tgrid_im[i] = cimag(tgrid[i]);
}

for (int i = 0; i < ntgrid; i++){
  Jacobi_sn_cn_dn(tgrid_re[i], 'm', m, &sn_re[i], &cn_re[i], &dn_re[i] );
  Jacobi_sn_cn_dn(tgrid_im[i], 'm', 1-m, &sn_tmp[i], &cn_tmp[i], &dn_tmp[i] );
  sn_im[i] = sn_tmp[i]/cn_tmp[i] * I;
  cn_im[i] = 1 / cn_tmp[i];
  dn_im[i] = dn_tmp[i] / cn_tmp[i];
}

for (int i = 0; i < ntgrid; i++){
  sn_c[i] = (sn_re[i]*sn_re[i] - sn_im[i]*sn_im[i]) / (sn_re[i]*cn_im[i]*dn_im[i] - sn_im[i]*cn_re[i]*dn_re[i]);
  cn_c[i] = (sn_re[i]*cn_re[i]*dn_im[i] - sn_im[i]*cn_im[i]*dn_re[i])/(sn_re[i]*cn_im[i]*dn_im[i] - sn_im[i]*cn_re[i]*dn_re[i]);
  dn_c[i] = (sn_re[i]*cn_im[i]*dn_re[i] - sn_im[i]*cn_re[i]*dn_im[i])/(sn_re[i]*cn_im[i]*dn_im[i] - sn_im[i]*cn_re[i]*dn_re[i]);
}

for (int i = 0; i < ntgrid; i++){
  z[i] = sqrt(mA*MA) * (kinv + sn_c[i]) / (kinv - sn_c[i]) ;
  z_re[i] = creal(z[i]);
  z_im[i] = cimag(z[i]);
}

for (int i = 0; i < ntgrid; i++){
  ksi_p[i] = csqrt(z[i] - mA);
  ksi_m[i] = - csqrt(z[i] - mA);
}


for (int i = 0; i < ntgrid; i++){
  fun_p[i] = fun_test(ksi_p[i], dmu, stodftInfo->beta, epsilon)*fun_test(ksi_p[i], dmu, stodftInfo->beta, epsilon);
  fun_m[i] = fun_test(ksi_m[i], dmu, stodftInfo->beta, epsilon)*fun_test(ksi_m[i], dmu, stodftInfo->beta, epsilon);
  fun_p[i] = fun_p[i]*cn_c[i]*dn_c[i] / ((kinv - sn_c[i])*(kinv - sn_c[i])) / ksi_p[i];
  fun_m[i] = fun_m[i]*cn_c[i]*dn_c[i] / ((kinv - sn_c[i])*(kinv - sn_c[i])) / ksi_m[i];
}

for (int j =0; j < ntgrid; j++){
  zseed[j] =  ksi_p[j];
  zseed[j + ntgrid] =  ksi_m[j];
}

*preRat = -2 * K *sqrt(mA*MA) / M_PI / ntgrid * kinv;

//printf("preRat = %lg %lg %lg %lg %lg %lg %i \n", preRat, K, mA, MA, M_PI, kinv, ntgrid);

free(tgrid);
free(tgrid_re);
free(tgrid_im);
free(sn_re);
free(cn_re);
free(dn_re);
free(sn_tmp);
free(cn_tmp);
free(dn_tmp);
free(cn_im);
free(dn_im);
free(sn_im);
free(sn_c);
free(cn_c);
free(dn_c);
free(z);
free(z_re);
free(z_im);
free(ksi_p);
free(ksi_m);


}
/*========================================================================================*/
/*========================================================================================*/


