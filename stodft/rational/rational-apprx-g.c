/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: rational.c                                     */
/*                                                                          */
/* This routine initialize the rational approximation part in G space.      */
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


double complex fermi_fun_g(double complex x, double dmu, double beta, double epsilon) {
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
void solve_shifted_eqn_cocg_g( CP *cp, CLASS *class, GENERAL_DATA *general_data, KOMEGAINFO *komegaInfo, int id, int is_mu_calc) { 
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
  RATIONALINFO *rationalInfo    = stodftInfo->rationalInfo;
 // KOMEGAINFO *komegaInfo = rationalInfo->komegaInfo;
  //komegaInfo = (KOMEGAINFO*)malloc(sizeof(KOMEGAINFO));
  MPI_Comm comm_states   =    communicate->comm_states;


  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

  int numThreads = communicate->numThreads;

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

  int itermax = rationalInfo->itermax;
  int ntgrid = rationalInfo->ntgrid;
  double threshold = rationalInfo->threshold;

  double *creRev12 = rationalInfo->creRev12;
  double *cimRev12 = rationalInfo->cimRev12;
  double *creImv12 = rationalInfo->creImv12;
  double *cimImv12 = rationalInfo->cimImv12;

  double *creRev2 = rationalInfo->creRev2;
  double *cimRev2 = rationalInfo->cimRev2;
  double *creImv2 = rationalInfo->creImv2;
  double *cimImv2 = rationalInfo->cimImv2;

  double *creRex = rationalInfo->creRex;
  double *cimRex = rationalInfo->cimRex;
  double *creImx = rationalInfo->creImx;
  double *cimImx = rationalInfo->cimImx;

  double *creRer_l = rationalInfo->creRer_l;
  double *cimRer_l = rationalInfo->cimRer_l;
  double *creImr_l = rationalInfo->creImr_l;
  double *cimImr_l = rationalInfo->cimImr_l;

  double *crerhs = rationalInfo->crerhs;
  double *cimrhs = rationalInfo->cimrhs;

  //double complex *v12 = rationalInfo->v12;
  //double complex *v2 = rationalInfo->v2;
  //double complex *r_l = rationalInfo->r_l;
  double complex *zseed = rationalInfo->zseed;
  //double complex *x = rationalInfo->x;
  //double *rhs = rationalInfo->rhs;
/* =------------------------------------------------------------------------------------= */
  int status[3];
  double timeStart1, timeEnd1;
  double timeStart2, timeEnd2;
  double t1, t2, tot, tot1;
  int i, j, iter, jiter, ndim, nz;

  double sum1, sum2, sum3, sum4, sum;


  tot = 0.0;
  tot1 = 0.0;
  nz = 2*ntgrid;
  ndim = numCoeff; //nfft2; // ndim TODO;
  spinFlag = 0;

  int ncoef;
  ncoef = ndim;

  komegaInfo->normfact = nfft2; 
  //threshold = threshold*nfft2;
  //printf("nfft2 %i %lg \n", nfft2, threshold);
  printf("--- Starting solving shifted COCG eqn  ---- numThreads %i numCoeff %i \n", numThreads, numCoeff);

  timeStart1 = omp_get_wtime();

 
/* =------------------------------------------------------------------------------------= */
    #pragma omp parallel for private(i) //TODO CHECK index etc. 
    for (i = 0; i < ncoef; i++) {
      creRev2[i] = cre_up[i+1];
      cimRev2[i] = cim_up[i+1];
      creImv2[i] = 0.0; //cre_up[ncoef+i+1];
      cimImv2[i] = 0.0; //cim_up[ncoef+i+1];
    }
    //printf("creRev2 %lg  %lg \n", creRev2[0], creRev2[100]);
    //printf("cimRev2 %lg  %lg \n", creRev2[0], creRev2[100]);
    //printf("creImv2 %lg  %lg \n", creRev2[0], creRev2[100]);
    //printf("cimImv2 %lg  %lg \n", creRev2[0], creRev2[100]);
/* =------------------------------------------------------------------------------------= */
  
  // get stochastic orbital from cre_up and cim_up //TODO
  //genNoiseOrbitalRealRational(cp,cpcoeffs_pos, v2, id); 
  
  
//TODO
//  if(is_mu_calc == 1) {
//    #pragma omp parallel for private(i)
//    for (i=0; i < nfft2; i++){
//      rhs[i] = creal(v2[i]);
//    }
//  }
  
  //printf("zdim %i dim %i ndim %i \n", zdim, dim, ndim);
  
  timeEnd1 = omp_get_wtime(); 
  timeStart2 = omp_get_wtime();
  //komega_cocg_init(&ndim, &ndim, &nz, x, zseed, &itermax, &threshold, NULL);
  komega_COCG_init_g(komegaInfo, ndim, ndim, nz, creRex, cimRex, creImx, cimImx, zseed, itermax, threshold);
  for (iter = 0; iter < itermax; iter++){
 

/* =------------------------------------------------------------------------------------= */
    #pragma omp parallel for private(i) //TODO memcpy? 
    for (i = 0; i < ndim; i++) {
      //r_l[i] = v2[i];
      creRer_l[i] = creRev2[i];
      cimRer_l[i] = cimRev2[i];
      creImr_l[i] = creImv2[i];
      cimImr_l[i] = cimImv2[i];
    }

    

/* =------------------------------------------------------------------------------------= */

    t1 = omp_get_wtime(); 
    calcCoefForceWrapSCFReal_g(class,general_data,cp,cpcoeffs_pos,clatoms_pos, creRev2, cimRev2, creImv2, cimImv2,
                            creRev12, cimRev12, creImv12, cimImv12, spinFlag);
    t2 = omp_get_wtime(); 
    tot = tot + (t2-t1); 
  
    //printf("creRev12 %lg  %lg \n", creRev12[0], creRev12[100]);
    //printf("cimRev12 %lg  %lg \n", creRev12[0], creRev12[100]);
    //printf("creImv12 %lg  %lg \n", creRev12[0], creRev12[100]);
    //printf("cimImv12 %lg  %lg \n", creRev12[0], creRev12[100]);

    t1 = omp_get_wtime(); 
    //komega_cocg_update(v12, v2, x, r_l, status);
    komega_COCG_update_g(komegaInfo, creRev12, cimRev12, creImv12, cimImv12,
                       creRev2, cimRev2, creImv2, cimImv2,
                       creRex, cimRex, creImx, cimImx,
                       creRer_l, cimRer_l, creImr_l, cimImr_l,
                     // v12, v2, x, r_l, 
                      status, numThreads);
    t2 = omp_get_wtime(); 
    tot1 = tot1 + (t2-t1); 
 
  sum1 = 0.0;
  sum2 = 0.0;
  sum3 = 0.0;
  sum4 = 0.0;
  for(i=0;i<ncoef-1;i++){
    sum1 += creRex[i]*creRex[i]+cimRex[i]*cimRex[i];
    sum2 += creImx[i]*creImx[i]+cimImx[i]*cimImx[i];
    sum3 += creRex[i]*creImx[i]+cimRex[i]*cimImx[i];
  }
  sum1 = sum1*2.0+creRex[ncoef-1]*creRex[ncoef-1];
  sum2 = sum2*2.0+creImx[ncoef-1]*creImx[ncoef-1];
  sum3 = sum3*2.0+creRex[ncoef-1]*creImx[ncoef-1];
  sum = sum1-sum2+2.0*sum3*I;


   // printf("last term %lg %lg \n", creRex[ncoef-1], creImx[ncoef-1] ); 
   //printf("x dotp x %lg %lg \n", creal(sum), cimag(sum));
    //printf(" DEBUG : %i %i %i %i %lg \n", iter, status[0], status[1], status[2], creal(v12[0]));
  
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
  
  //komega_cocg_finalize();
  komega_COCG_finalize_g(komegaInfo);
  
  timeEnd2 = omp_get_wtime(); 
  printf("--- Finished solving shifted COCG eqn  ---- %lg %lg %lg %lg \n", timeEnd1-timeStart1, timeEnd2-timeStart2, tot, tot1);



}
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void filterRational_g(CP *cp,CLASS *class,GENERAL_DATA *general_data, KOMEGAINFO *komegaInfo,
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
  RATIONALINFO *rationalInfo = stodftInfo->rationalInfo;

  //KOMEGAINFO *komegaInfo = rationalInfo->komegaInfo; // = (KOMEGAINFO*)cmalloc(sizeof(KOMEGAINFO));
  //rationalInfo->komegaInfo = komegaInfo; 

  int expanType      = stodftInfo->expanType;
  int polynormLength = stodftInfo->polynormLength;
  int numChemPot     = stodftInfo->numChemPot;
  int smearOpt       = stodftInfo->smearOpt;
  int filterDiagFlag = stodftInfo->filterDiagFlag;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numStateStoUp = stodftInfo->numStateStoUp;
  int numStateStoDn = stodftInfo->numStateStoDn;
  int numCoeff       = cpcoeffs_info->ncoef;
  int cpLsda         = cpopts->cp_lsda;
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

  double *chemPot = stodftCoefPos->chemPot;


  double *expanCoeff = (double*)stodftCoefPos->expanCoeff;
  //double timeStart,timeEnd;
  //double timeStart2,timeEnd2;
  //double timeStart3,timeEnd3;
  double deltaTime = 0.0;
  double deltaTime2 = 0.0;
  double dot;


  int nfft          = cp_para_fft_pkg3d_lg->nfft;
  int nfft2         = nfft/2;
  double complex sum, sum_0, sum_p, sum_m;
  double dsum_0, dsum_p, dsum_m; 
  int ntgrid = rationalInfo->ntgrid;

  double complex *ksi_p = rationalInfo->ksi_p;
  double complex *ksi_m = rationalInfo->ksi_m;
  double complex *fun_p = rationalInfo->fun_p;
  double complex *fun_m = rationalInfo->fun_m;
  double complex *rat_fact_p = rationalInfo->rat_fact_p;
  double complex *rat_fact_m = rationalInfo->rat_fact_m;
  double dmu = rationalInfo->dmu;
  double epsilon = rationalInfo->epsilon; 
   
  double complex *fun_p_0 = rationalInfo->fun_p_0;
  double complex *fun_m_0 = rationalInfo->fun_m_0;
  double complex *fun_p_p = rationalInfo->fun_p_p;
  double complex *fun_m_p = rationalInfo->fun_m_p;
  double complex *fun_p_m = rationalInfo->fun_p_m;
  double complex *fun_m_m = rationalInfo->fun_m_m;

  double complex *x = rationalInfo->x;
  //double preRat = rationalInfo->preRat;
  double *frhs = rationalInfo->frhs;
  double *rhs = rationalInfo->rhs;
  double numElecMin,numElecMax, numElecTot;

  double dNdm, fmu;
  double maxmu = rationalInfo->maxmu;
  double large_dmu = rationalInfo->large_dmu;
  double small_dmu = rationalInfo->small_dmu;
  double numElecTrue = stodftInfo->numElecTrue;
  double chemPotNew;
  double numElecNew;
  double chemPotMin,chemPotMax;
  double numElecTol = 1.0e-8*numElecTrue;
  double delta_mu;

  double timeStart1, timeEnd1;
  double timeStart2, timeEnd2;
  double timeStart3, timeEnd3;
  double t1, t2, t3, t4, t5, t6, t7, t8, t9;

  double sum_reRe, sum_imRe, sum_reIm, sum_imIm, fact1_real, fact1_imag, fact2_real, fact2_imag;
  double dsum_re, dsum_im;

  double *creRex = rationalInfo->creRex;
  double *cimRex = rationalInfo->cimRex;
  double *creImx = rationalInfo->creImx;
  double *cimImx = rationalInfo->cimImx;

  numElecTot = 0.0;
  numElecMin = 0.0;
  numElecMax = 0.0;

  printf("start tessssssssssssst RA \n");

  printf("myidState = %i %i %lg %lg %i %lg %lg %lg %lg %i %lg \n", 
   myidState,rationalInfo->ntgrid, rationalInfo->dmu, rationalInfo->threshold,
   rationalInfo->itermax, rationalInfo->large_dmu,
   rationalInfo->small_dmu, rationalInfo->maxmu, rationalInfo->epsilon, expanType,
   rationalInfo->init_mu);

/******************************************************************************/
printf("Starting Filtering with Rational Approximation\n");
/******************************************************************************/

t1 = omp_get_wtime();

init_zseed_g(cp, 0);

dsum_0 = 0.0;
dsum_p = 0.0;
dsum_m = 0.0;
dsum_re = 0.0;
dsum_im = 0.0;

t2 = omp_get_wtime();

//printf(" dmu, stodftInfo->beta, epsilon %lg %lg %lg \n", dmu, stodftInfo->beta, epsilon );
//printf("==== FINAL RESULTS   === %i %lg %lg %lg %lg \n", myidState, creal(ksi_p[10]), creal(ksi_m[10]), cimag(ksi_p[10]), cimag(ksi_m[10]));

for (int i = 0; i < ntgrid; i++){
  fun_p[i] = fermi_fun_g(ksi_p[i], 0.0, stodftInfo->beta, epsilon);
  fun_m[i] = fermi_fun_g(ksi_m[i], 0.0, stodftInfo->beta, epsilon);

  fun_p_0[i] = fun_p[i]*fun_p[i];
  fun_m_0[i] = fun_m[i]*fun_m[i];
  //double check the sign 
  fun_p_p[i] = fermi_fun_g(ksi_p[i], -dmu, stodftInfo->beta, epsilon)*fermi_fun_g(ksi_p[i], -dmu, stodftInfo->beta, epsilon);
  fun_m_p[i] = fermi_fun_g(ksi_m[i], -dmu, stodftInfo->beta, epsilon)*fermi_fun_g(ksi_m[i], -dmu, stodftInfo->beta, epsilon);
  fun_p_m[i] = fermi_fun_g(ksi_p[i], dmu, stodftInfo->beta, epsilon)*fermi_fun_g(ksi_p[i], dmu, stodftInfo->beta, epsilon);
  fun_m_m[i] = fermi_fun_g(ksi_m[i], dmu, stodftInfo->beta, epsilon)*fermi_fun_g(ksi_m[i], dmu, stodftInfo->beta, epsilon);
}

t3 = omp_get_wtime();
/*
printf("==== final RESULTS   === %i %lg %lg %lg %lg \n", myidState, creal(fun_p[10]), creal(fun_m[10]), creal(fun_p_0[10]), creal(fun_m_0[10]));
printf("==== final RESULTS I   === %i %lg %lg %lg %lg \n", myidState, cimag(fun_p[10]), cimag(fun_m[10]), cimag(fun_p_0[10]), cimag(fun_m_0[10]));

printf("==== final results:   === %i %lg %lg %lg %lg \n", myidState, creal(fun_p_p[10]), creal(fun_m_p[10]), creal(fun_p_m[10]), creal(fun_m_m[10]));
printf("==== final results: I   === %i %lg %lg %lg %lg \n", myidState, cimag(fun_p_p[10]), cimag(fun_m_p[10]), cimag(fun_p_m[10]), cimag(fun_m_m[10]));


  printf("end tessssssssssssst RA \n");
  Barrier(comm_states);
  printf("@@@@@@@@@@@@@@@@@@@@_forced_stop__@@@@@@@@@@@@@@@@@@@@\n");
  fflush(stdout);
  exit(1);
*/

for(iState=0;iState<numStateUpProc;iState++){

  //printf("myidState %i %i %i \n", myidState, numStateUpProc, iState);
  solve_shifted_eqn_cocg_g( cp, class, general_data, komegaInfo, iState, 1);


  printf("last term %lg %lg \n", creRex[numCoeff-1], creImx[numCoeff-1] );
  //for (int i = 0; i<nfft2;i++){
  for (int i = 0; i<numCoeff;i++){
    //sum   = 0.0 + 0.0 *I ;
    sum_0 = 0.0 + 0.0 *I ;
    sum_p = 0.0 + 0.0 *I ;
    sum_m = 0.0 + 0.0 *I ;
    sum_reRe = 0.0;
    sum_imRe = 0.0;
    sum_reIm = 0.0;
    sum_imIm = 0.0;
    for (int j =0; j < ntgrid; j++){
      //sum = sum + fun_p[j] * rat_fact_p[j] * x[j*nfft2 + i] 
      //          + fun_m[j] * rat_fact_m[j] * x[(j+ntgrid)*nfft2 + i];

      fact1_real = creal(fun_p_0[j] * rat_fact_p[j]);
      fact1_imag = cimag(fun_p_0[j] * rat_fact_p[j]);
      fact2_real = creal(fun_m_0[j] * rat_fact_m[j]);
      fact2_imag = cimag(fun_m_0[j] * rat_fact_m[j]);

      sum_reRe += fact1_real*creRex[j*numCoeff + i] - fact1_imag*creImx[j*numCoeff + i] +
	      fact2_real*creRex[(j+ntgrid)*numCoeff + i] - fact2_imag*creImx[(j+ntgrid)*numCoeff + i]; 

      sum_imRe += fact1_real*cimRex[j*numCoeff + i] - fact1_imag*cimImx[j*numCoeff + i] +
	      fact2_real*cimRex[(j+ntgrid)*numCoeff + i] - fact2_imag*cimImx[(j+ntgrid)*numCoeff + i]; 

      sum_reIm += fact1_real*creImx[j*numCoeff + i] + fact1_imag*creRex[j*numCoeff + i] +
	      fact2_real*creImx[(j+ntgrid)*numCoeff + i] + fact2_imag*creRex[(j+ntgrid)*numCoeff + i]; 

      sum_imIm += fact1_real*cimImx[j*numCoeff + i] + fact1_imag*cimRex[j*numCoeff + i] +
	      fact2_real*cimImx[(j+ntgrid)*numCoeff + i] + fact2_imag*cimRex[(j+ntgrid)*numCoeff + i]; 

//      sum_0 = sum_0 + fun_p_0[j] * rat_fact_p[j] * x[j*nfft2 + i] 
//                    + fun_m_0[j] * rat_fact_m[j] * x[(j+ntgrid)*nfft2 + i];
//      sum_p = sum_p + fun_p_p[j] * rat_fact_p[j] * x[j*nfft2 + i] 
//                    + fun_m_p[j] * rat_fact_m[j] * x[(j+ntgrid)*nfft2 + i];
//      sum_m = sum_m + fun_p_m[j] * rat_fact_p[j] * x[j*nfft2 + i] 
//                   + fun_m_m[j] * rat_fact_m[j] * x[(j+ntgrid)*nfft2 + i];
    }
    //frhs[i] = rationalInfo->preRat*cimag(sum);
    //printf("i = %i, %lg %lg \n", i, 
    stoWfUpRe[0][i+1] = rationalInfo->preRat*sum_reIm;
    stoWfUpIm[0][i+1] = rationalInfo->preRat*sum_imIm;

    printf("i = %i, %lg %lg \n", i, rationalInfo->preRat*sum_reIm, rationalInfo->preRat*sum_imIm);
    // get stochastic orbital from cre_up and cim_up to rhs //TODO
    //dsum_0 += 2.0 * rationalInfo->preRat * cimag(sum_0) * rhs[i];
    if(i < numCoeff -1) {
        dsum_0 += 2.0 * (stoWfUpRe[0][i+1] * cre_up[i+1] + stoWfUpIm[0][i+1] * cim_up[i+1]);
        //dsum_0 += 2.0 * (stoWfUpRe[0][i+1] * stoWfUpRe[0][i+1] + stoWfUpIm[0][i+1] * stoWfUpIm[0][i+1]);
    }
    else{
        dsum_0 += stoWfUpRe[0][i+1] * cre_up[i+1];     
        //dsum_0 += stoWfUpRe[0][i+1] * stoWfUpRe[0][i+1];     
    }
    //dsum_re += 2.0 * stoWfUpRe[0][i+1] * cre_up[i+1];
    //dsum_im += 2.0 * stoWfUpIm[0][i+1] * cim_up[i+1];

    //dsum_p += 2.0 * rationalInfo->preRat * cimag(sum_p) * rhs[i];
    //dsum_m += 2.0 * rationalInfo->preRat * cimag(sum_m) * rhs[i];
  }

  //rhsReal(class, general_data, cp, cpcoeffs_pos, clatoms_pos, frhs, iState);
}


t4 = omp_get_wtime();

dsum_0 = 2.0*(dsum_0/numStateStoUp);
//dsum_0 = (dsum_0/numStateStoUp)*(1.0/stodftInfo->rhoRealGridTot);
//dsum_p = (dsum_p/numStateStoUp)*(1.0/stodftInfo->rhoRealGridTot);
//dsum_m = (dsum_m/numStateStoUp)*(1.0/stodftInfo->rhoRealGridTot);

printf("==== final results: dsum  === %i %lg %lg %lg %lg \n", myidState, dsum_0, dsum_p, dsum_m, rationalInfo->preRat);

if(numProcStates>1)Reduce(&dsum_0,&numElecTot,1,MPI_DOUBLE,MPI_SUM,0,comm_states);
if(numProcStates>1)Reduce(&dsum_p,&numElecMax,1,MPI_DOUBLE,MPI_SUM,0,comm_states);
if(numProcStates>1)Reduce(&dsum_m,&numElecMin,1,MPI_DOUBLE,MPI_SUM,0,comm_states);

if(myidState == 0)printf("==== final results: numElec  === %.10lg %.10lg %.10lg\n", numElecTot, numElecMax, numElecMin);



   if(numProcStates>1){
     Barrier(comm_states);
     Bcast(&numElecTot,1,MPI_DOUBLE,0,comm_states);
     Bcast(&numElecMax,1,MPI_DOUBLE,0,comm_states);
     Bcast(&numElecMin,1,MPI_DOUBLE,0,comm_states);
   }


  printf("end checking \n");
  Barrier(comm_states);
  printf("@@@@@@@@@@@@@@@@@@@@_forced_stop__@@@@@@@@@@@@@@@@@@@@\n");
  fflush(stdout);
  exit(1);
/*=====================================================================================*/


t5 = omp_get_wtime();

if ((numElecMax > numElecTrue) && (numElecTrue > numElecMin)){

    //chemPotNew = ((stodftInfo->chemPotTrue + dmu ) + (stodftInfo->chemPotTrue - dmu) ) / 2.0;
    chemPotNew = stodftInfo->chemPotTrue;

    numElecNew = numElecTot; 
    //numElecNew = calcNumberElecRational_g(cp, chemPotNew);

    chemPotMin = stodftInfo->chemPotTrue - dmu;
    chemPotMax = stodftInfo->chemPotTrue + dmu;

    //printf("test MU dmu %lg \n", calcNumberElecRational_g(cp, dmu));
    //printf("test MU -dmu %lg \n", calcNumberElecRational_g(cp, -dmu));

    while(fabs(numElecNew-numElecTrue)>numElecTol){
      if(numElecNew>numElecTrue){
        chemPotMax = chemPotNew;
        numElecMax = numElecNew;
      }
      if(numElecNew<numElecTrue){
        chemPotMin = chemPotNew;
        numElecMin = numElecNew;
      }
      chemPotNew = 0.5*(chemPotMin+chemPotMax);
      //numElecNew = calcNumElecCheby(cp,chemPotNew,chebyCoeffs);
      numElecNew = calcNumberElecRational_g(cp, stodftInfo->chemPotTrue - chemPotNew);

      if(myidState == 0)printf(" while loop %lg %lg %lg %lg %lg %lg %lg \n", chemPotMax, chemPotMin, chemPotNew, numElecMax, numElecMin, numElecNew, stodftInfo->chemPotTrue - chemPotNew );
    }





/*************************************************/
   if(numProcStates>1){
     Barrier(comm_states);
     Bcast(&chemPotNew,1,MPI_DOUBLE,0,comm_states);
   }

   delta_mu = stodftInfo->chemPotTrue - chemPotNew;

   chemPot[0] = chemPotNew;

   stodftInfo->chemPotTrue = chemPotNew;

   if(myidState==0){
     printf("Correct Chemical Potential is %.16lg\n",chemPotNew);
     fflush(stdout);
   }

/*************************************************/
//applyFilterRational_g(cp, class, general_data, delta_mu, ip_now);
/*************************************************/


/*************************************************/
/*************************************************/
// DELETE LATER
/*************************************************/
  init_zseed_g(cp, 0);
  
  for (int i = 0; i < ntgrid; i++){
    fun_p[i] = fermi_fun_g(ksi_p[i], 0.0, stodftInfo->beta, epsilon);
    fun_m[i] = fermi_fun_g(ksi_m[i], 0.0, stodftInfo->beta, epsilon);
  }
  
  for(iState=0;iState<numStateUpProc;iState++){
  
    //printf("myidState %i %i %i \n", myidState, numStateUpProc, iState);
    solve_shifted_eqn_cocg_g( cp, class, general_data, komegaInfo, iState, 0);
  
    for (int i = 0; i<nfft2;i++){
      sum   = 0.0 + 0.0 *I ;
      for (int j =0; j < ntgrid; j++){
        sum = sum + fun_p[j] * rat_fact_p[j] * x[j*nfft2 + i] 
                  + fun_m[j] * rat_fact_m[j] * x[(j+ntgrid)*nfft2 + i];
  
      }
      frhs[i] = rationalInfo->preRat*cimag(sum);
    }
  
    rhsReal(class, general_data, cp, cpcoeffs_pos, clatoms_pos, frhs, iState);
  }
/******************************************************************************/

printf("here num elec %.16lg \n", calcNumberElecRational_g(cp, 0.0));

/*************************************************/
// DELETE UPTO HERE
/*************************************************/
/*************************************************/

}
else {

/******************************************************************************/
/******************************************************************************/
   chemPotNew = stodftInfo->chemPotTrue;
   
   dNdm = (numElecMax - numElecMin)/(2.0*small_dmu);

   fmu = -(numElecTot - numElecTrue)/(numElecTrue*numElecTrue)*dNdm;

   if(myidState == 0)printf("dNdm fmu fmu*large_dmu %lg %lg %lg %lg \n", dNdm, fmu, 
             fmu*large_dmu, -(numElecTot - numElecTrue)/(numElecTrue*numElecTrue) );

   if (fabs(fmu*large_dmu) > maxmu ){
     if(fmu>0) chemPotNew = chemPotNew + maxmu;
     else chemPotNew = chemPotNew - maxmu;
   }
   else {
     chemPotNew = chemPotNew + (fmu*large_dmu);
   }

/******************************************************************************/
   if(numProcStates>1){
     Barrier(comm_states);
     Bcast(&chemPotNew,1,MPI_DOUBLE,0,comm_states);
   }

   delta_mu = stodftInfo->chemPotTrue - chemPotNew;

   chemPot[0] = chemPotNew;

   stodftInfo->chemPotTrue = chemPotNew;

   if(myidState==0){
     printf("Correct Chemical Potential is %.16lg\n",chemPotNew);
     fflush(stdout);
   }
/******************************************************************************/
  init_zseed_g(cp, 0);
  
  for (int i = 0; i < ntgrid; i++){
    fun_p[i] = fermi_fun_g(ksi_p[i], 0.0, stodftInfo->beta, epsilon);
    fun_m[i] = fermi_fun_g(ksi_m[i], 0.0, stodftInfo->beta, epsilon);
  }
  
  for(iState=0;iState<numStateUpProc;iState++){
  
    //printf("myidState %i %i %i \n", myidState, numStateUpProc, iState);
    solve_shifted_eqn_cocg_g( cp, class, general_data, komegaInfo, iState, 0);
  
    for (int i = 0; i<nfft2;i++){
      sum   = 0.0 + 0.0 *I ;
      for (int j =0; j < ntgrid; j++){
        sum = sum + fun_p[j] * rat_fact_p[j] * x[j*nfft2 + i] 
                  + fun_m[j] * rat_fact_m[j] * x[(j+ntgrid)*nfft2 + i];
  
      }
      frhs[i] = rationalInfo->preRat*cimag(sum);
    }
  
    rhsReal(class, general_data, cp, cpcoeffs_pos, clatoms_pos, frhs, iState);
  }
/******************************************************************************/

printf("here num elec %.16lg \n", calcNumberElecRational_g(cp, 0.0));

}

/*
   if(numProcStates>1){
     Barrier(comm_states);
     Bcast(&chemPotNew,1,MPI_DOUBLE,0,comm_states);
   }

   delta_mu = stodftInfo->chemPotTrue - chemPotNew;

   chemPot[0] = chemPotNew;

   stodftInfo->chemPotTrue = chemPotNew;

   if(myidState==0){
     printf("Correct Chemical Potential is %.16lg\n",chemPotNew);
     fflush(stdout);
   }


t6 = omp_get_wtime();

  printf("delta_mu is %lg \n", delta_mu);
*/
/******************************************************************************/


/******************************************************************************/
//printf("Finished Filtering with Rational Approximation times %lg %lg \n", timeEnd2 - timeStart2, timeEnd1 - timeStart1);
/******************************************************************************/
/******************************************************************************/

/*
if(fabs(delta_mu) >  small_dmu) {

  init_zseed_g(cp, 0);
  
  for (int i = 0; i < ntgrid; i++){
    fun_p[i] = fermi_fun_g(ksi_p[i], 0.0, stodftInfo->beta, epsilon);
    fun_m[i] = fermi_fun_g(ksi_m[i], 0.0, stodftInfo->beta, epsilon);
  }
  
  for(iState=0;iState<numStateUpProc;iState++){
  
    //printf("myidState %i %i %i \n", myidState, numStateUpProc, iState);
    solve_shifted_eqn_cocg_g( cp, class, general_data, komegaInfo, iState, 0);
  
    for (int i = 0; i<nfft2;i++){
      sum   = 0.0 + 0.0 *I ;
      for (int j =0; j < ntgrid; j++){
        sum = sum + fun_p[j] * rat_fact_p[j] * x[j*nfft2 + i] 
                  + fun_m[j] * rat_fact_m[j] * x[(j+ntgrid)*nfft2 + i];
  
      }
      frhs[i] = rationalInfo->preRat*cimag(sum);
    }
  
    rhsReal(class, general_data, cp, cpcoeffs_pos, clatoms_pos, frhs, iState);
  }
}
else {

applyFilterRational_g(cp, class, general_data, delta_mu, ip_now);

}

t7 = omp_get_wtime();
*/
/******************************************************************************/
printf("Finished Filtering with Rational Approximation times %lg %lg %lg %lg %lg %lg \n", t2-t1, t3-t2, t4-t3, t5-t4, t6-t5, t7-t6);
/******************************************************************************/
/*
  printf("end tessssssssssssst RA \n");
  Barrier(comm_states);
  printf("@@@@@@@@@@@@@@@@@@@@_forced_stop__@@@@@@@@@@@@@@@@@@@@\n");
  fflush(stdout);
  exit(1);
*/
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/

void applyFilterRational_g(CP *cp, CLASS *class, GENERAL_DATA *general_data, double delta_mu, int ip_now) {

STODFTINFO *stodftInfo        = cp->stodftInfo;
CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
COMMUNICATE *communicate      = &(cp->communicate);
CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
CLATOMS_POS *clatoms_pos      = &(class->clatoms_pos[ip_now]);
RATIONALINFO *rationalInfo    = stodftInfo->rationalInfo;
PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

int numStateUpProc = cpcoeffs_info->nstate_up_proc;
int numCoeff       = cpcoeffs_info->ncoef;

  int nfft          = cp_para_fft_pkg3d_lg->nfft;
  int nfft2         = nfft/2;

  int ntgrid = rationalInfo->ntgrid;
  double complex *fun_p = rationalInfo->fun_p;
  double complex *fun_m = rationalInfo->fun_m;
  double complex *ksi_p = rationalInfo->ksi_p;
  double complex *ksi_m = rationalInfo->ksi_m;
  double complex *rat_fact_p = rationalInfo->rat_fact_p;
  double complex *rat_fact_m = rationalInfo->rat_fact_m;
  //double complex *x = rationalInfo->x;
  double *frhs = rationalInfo->frhs;
  double epsilon = rationalInfo->epsilon;


  double *creRex = rationalInfo->creRex;
  double *cimRex = rationalInfo->cimRex;
  double *creImx = rationalInfo->creImx;
  double *cimImx = rationalInfo->cimImx;

//TODO CHECK
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;

  int  iState, i, j;
  double complex sum;

  double sum_cre, sum_cim;

printf(" from applying delta_mu = %lg \n", delta_mu);

/*==========================================================================*/
for (int i = 0; i < ntgrid; i++){
  fun_p[i] = fermi_fun_g(ksi_p[i], delta_mu, stodftInfo->beta, epsilon);
  fun_m[i] = fermi_fun_g(ksi_m[i], delta_mu, stodftInfo->beta, epsilon);
}

double *scale1_re = (double*)cmalloc(ntgrid*sizeof(double));
double *scale1_im = (double*)cmalloc(ntgrid*sizeof(double));
double *scale2_re = (double*)cmalloc(ntgrid*sizeof(double));
double *scale2_im = (double*)cmalloc(ntgrid*sizeof(double));
for(j=0;j<ntgrid;j++){
  scale1_re[j] = creal(fun_p[j] * rat_fact_p[j]);
  scale1_im[j] = cimag(fun_p[j] * rat_fact_p[j]);
  scale2_re[j] = creal(fun_m[j] * rat_fact_m[j]);
  scale2_im[j] = cimag(fun_m[j] * rat_fact_m[j]);
}

for(iState=0;iState<numStateUpProc;iState++){
  /*----------rational-g----------*/
  //for (int i = 0; i<nfft2;i++){
  //  sum   = 0.0 + 0.0 *I ;
  //  for (int j =0; j < ntgrid; j++){
  //    sum = sum + fun_p[j] * rat_fact_p[j] * x[j*nfft2 + i] 
  //              + fun_m[j] * rat_fact_m[j] * x[(j+ntgrid)*nfft2 + i];

  //  }
  //  frhs[i] = rationalInfo->preRat*cimag(sum);
  //}

  for(i=0;i<numCoeff;i++){
    sum_cre = 0.0;
    sum_cim = 0.0;
    for(j=0;j<ntgrid;j++){
       sum_cre += scale1_re[j]*creImx[j*numCoeff+i]+
                  scale1_im[j]*creRex[j*numCoeff+i]+
                  scale2_re[j]*creImx[(j+ntgrid)*numCoeff+i]+
                  scale2_im[j]*creRex[(j+ntgrid)*numCoeff+i];
       sum_cim += scale1_re[j]*cimImx[j*numCoeff+i]+
                  scale1_im[j]*cimRex[j*numCoeff+i]+
                  scale2_re[j]*cimImx[(j+ntgrid)*numCoeff+i]+
                  scale2_im[j]*cimRex[(j+ntgrid)*numCoeff+i];
    }
    cre_up[i+1] = sum_cre*rationalInfo->preRat;
    cim_up[i+1] = sum_cim*rationalInfo->preRat;
  }

  /*----------end rational-g----------*/

  //rhsReal(class, general_data, cp, cpcoeffs_pos, clatoms_pos, frhs, iState);
}
free(scale1_re);
free(scale1_im);
free(scale2_re);
free(scale2_im);


/*==========================================================================*/

}


/*==========================================================================*/
double calcNumberElecRational_g(CP *cp, double delta_mu) {

STODFTINFO *stodftInfo        = cp->stodftInfo;
CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
COMMUNICATE *communicate      = &(cp->communicate);
RATIONALINFO *rationalInfo = stodftInfo->rationalInfo;
PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

int numStateUpProc = cpcoeffs_info->nstate_up_proc;
int numStateStoUp = stodftInfo->numStateStoUp;
int numProcStates = communicate->np_states;
int myidState       = communicate->myid_state;
MPI_Comm comm_states   =    communicate->comm_states;

  int nfft          = cp_para_fft_pkg3d_lg->nfft;
  int nfft2         = nfft/2;

  int ntgrid = rationalInfo->ntgrid;
  double complex *ksi_p = rationalInfo->ksi_p;
  double complex *ksi_m = rationalInfo->ksi_m;
  double complex *fun_p_p = rationalInfo->fun_p_p;
  double complex *fun_m_p = rationalInfo->fun_m_p;
  double complex *rat_fact_p = rationalInfo->rat_fact_p;
  double complex *rat_fact_m = rationalInfo->rat_fact_m;
  double complex *x = rationalInfo->x;
  double *rhs = rationalInfo->rhs;
  double epsilon = rationalInfo->epsilon;


  double dsum_p;
  int  iState;
  double complex sum_p;
  double numElecMax;

printf("delta_mu = %lg \n", delta_mu);

dsum_p = 0.0;

for (int i = 0; i < ntgrid; i++){
  fun_p_p[i] = fermi_fun_g(ksi_p[i], delta_mu, stodftInfo->beta, epsilon)*fermi_fun_g(ksi_p[i], delta_mu, stodftInfo->beta, epsilon);
  fun_m_p[i] = fermi_fun_g(ksi_m[i], delta_mu, stodftInfo->beta, epsilon)*fermi_fun_g(ksi_m[i], delta_mu, stodftInfo->beta, epsilon);
}


for(iState=0;iState<numStateUpProc;iState++){

  for (int i = 0; i<nfft2;i++){
    sum_p = 0.0 + 0.0 *I ;
    for (int j =0; j < ntgrid; j++){
      sum_p = sum_p + fun_p_p[j] * rat_fact_p[j] * x[j*nfft2 + i] 
                    + fun_m_p[j] * rat_fact_m[j] * x[(j+ntgrid)*nfft2 + i];
    }
    dsum_p += 2.0 * rationalInfo->preRat * cimag(sum_p) * rhs[i];
  }

}

dsum_p = (dsum_p/numStateStoUp)*(1.0/stodftInfo->rhoRealGridTot);

if(numProcStates>1)Reduce(&dsum_p,&numElecMax,1,MPI_DOUBLE,MPI_SUM,0,comm_states);


   if(numProcStates>1){
     Barrier(comm_states);
     Bcast(&numElecMax,1,MPI_DOUBLE,0,comm_states);
   }

if(myidState == 0)printf("==== final results: numElec  === %.10lg \n", numElecMax);
return numElecMax;
}


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcChemPotRational_g(CP *cp,CLASS *class,GENERAL_DATA *general_data, KOMEGAINFO *komegaInfo,
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
  RATIONALINFO *rationalInfo    = stodftInfo->rationalInfo;

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

  int ntgrid = rationalInfo->ntgrid;
  double complex *fun_p = rationalInfo->fun_p;
  double complex *fun_m = rationalInfo->fun_m;
  double complex *rat_fact_p = rationalInfo->rat_fact_p;
  double complex *rat_fact_m = rationalInfo->rat_fact_m;
  double complex *x = rationalInfo->x;
  //double preRat = rationalInfo->preRat;
  double *frhs = rationalInfo->frhs;
  double *rhs = rationalInfo->rhs;
  double maxmu = rationalInfo->maxmu;
  double large_dmu = rationalInfo->large_dmu;
  double small_dmu = rationalInfo->small_dmu;

  double complex sum;
  double dsum;
  double NElecTot;
  double dNdm, fmu;

  NElecTot = 0.0;
/******************************************************************************/
printf("Start ChemicalPotential with Rational Approximation\n");
/******************************************************************************/

chemPotNew = stodftInfo->chemPotTrue;

/******************************************************************************/
init_zseed_g(cp, 1);

printf("myidStates %i %lg  \n", myidState, rationalInfo->preRat);
//for (int j =0; j < ntgrid; j++){
// printf("myidState %i %lg %lg \n", myidState,creal(rationalInfo->fun_m[j]), cimag(rationalInfo->fun_m[j]));
//}
/*
Barrier(comm_states);
printf("@@@@@@@@@@@@@@@@@@@@_forced_stop__@@@@@@@@@@@@@@@@@@@@\n");
fflush(stdout);
exit(1);
*/
dsum = 0.0;
for(iState=0;iState<numStateUpProc;iState++){

  //printf("myidState %i %i %i %i %i \n", myidState, numStateUpProc, iState, numProcStates, numStateStoUp);
  solve_shifted_eqn_cocg_g( cp, class, general_data, komegaInfo, iState, 1);

  for (int i = 0; i<nfft2;i++){
    sum = 0.0 + 0.0 *I ;
    for (int j =0; j < ntgrid; j++){
      sum = sum + fun_p[j] * x[j*nfft2 + i] + fun_m[j] * x[(j+ntgrid)*nfft2 + i];
    }
    frhs[i] = rationalInfo->preRat*cimag(sum);
    dsum = dsum + 2.0*frhs[i]*rhs[i];
  }

}

dsum = (dsum/numStateStoUp)*(1.0/stodftInfo->rhoRealGridTot);
if(numProcStates>1)Reduce(&dsum,&NElecTot,1,MPI_DOUBLE,MPI_SUM,0,comm_states);
printf("==== final results: NElecTot === %lg %lg \n", dsum, NElecTot);
/******************************************************************************/
/******************************************************************************/
stodftInfo->chemPotTrue = chemPotNew - small_dmu;
init_zseed_g(cp, 1);

dsum = 0.0;
for(iState=0;iState<numStateUpProc;iState++){

  //printf("myidState %i %i %i %i %i \n", myidState, numStateUpProc, iState, numProcStates, numStateStoUp);
  solve_shifted_eqn_cocg_g( cp, class, general_data, komegaInfo, iState, 1);

  for (int i = 0; i<nfft2;i++){
    sum = 0.0 + 0.0 *I ;
    for (int j =0; j < ntgrid; j++){
      sum = sum + fun_p[j] * x[j*nfft2 + i] + fun_m[j] * x[(j+ntgrid)*nfft2 + i];
    }
    frhs[i] = rationalInfo->preRat*cimag(sum);
    dsum = dsum + 2.0*frhs[i]*rhs[i];
  }

}

dsum = (dsum/numStateStoUp)*(1.0/stodftInfo->rhoRealGridTot);
if(numProcStates>1)Reduce(&dsum,&numElecMin,1,MPI_DOUBLE,MPI_SUM,0,comm_states);
printf("==== final results: numElecMin === %lg %lg \n", dsum, numElecMin);
/******************************************************************************/
/******************************************************************************/
stodftInfo->chemPotTrue = chemPotNew + small_dmu;
init_zseed_g(cp, 1);

dsum = 0.0;
for(iState=0;iState<numStateUpProc;iState++){

  //printf("myidState %i %i %i %i %i \n", myidState, numStateUpProc, iState, numProcStates, numStateStoUp);
  solve_shifted_eqn_cocg_g( cp, class, general_data, komegaInfo, iState, 1);

  for (int i = 0; i<nfft2;i++){
    sum = 0.0 + 0.0 *I ;
    for (int j =0; j < ntgrid; j++){
      sum = sum + fun_p[j] * x[j*nfft2 + i] + fun_m[j] * x[(j+ntgrid)*nfft2 + i];
    }
    frhs[i] = rationalInfo->preRat*cimag(sum);
    dsum = dsum + 2.0*frhs[i]*rhs[i];
  }

}

dsum = (dsum/numStateStoUp)*(1.0/stodftInfo->rhoRealGridTot);
if(numProcStates>1)Reduce(&dsum,&numElecMax,1,MPI_DOUBLE,MPI_SUM,0,comm_states);
printf("==== final results: numElecMax === %lg %lg \n", dsum, numElecMax);
/******************************************************************************/
stodftInfo->chemPotTrue = chemPotNew;

   dNdm = (numElecMax - numElecMin)/(2.0*small_dmu);
   //fmu = -(NElecTot - numElecTrue)*dNdm; 
   fmu = -(NElecTot - numElecTrue)/(numElecTrue*numElecTrue)*dNdm; 

   if (fabs(fmu*large_dmu) > maxmu ){
     if(fmu>0) chemPotNew = chemPotNew + maxmu;
     else chemPotNew = chemPotNew - maxmu;
   } 
   else {
     chemPotNew = chemPotNew + (fmu*large_dmu);
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
/******************************************************************************/
printf("Finished ChemicalPotential with Rational Approximation\n");
/******************************************************************************/
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/
/*========================================================================================*/
void init_zseed_g( CP *cp, int is_fun_sq) {

#include "../typ_defs/typ_mask.h"


STODFTINFO *stodftInfo        = cp->stodftInfo;
RATIONALINFO *rationalInfo    = stodftInfo->rationalInfo;

double K = rationalInfo->K;
double K_prim = rationalInfo->K_prim;
double mA = rationalInfo->mA;
double MA = rationalInfo->MA;
double k = rationalInfo->k;
double kinv = rationalInfo->kinv;
double m = rationalInfo->m;
double dmu = rationalInfo->dmu;
double epsilon = rationalInfo->epsilon;
int ntgrid = rationalInfo->ntgrid;
double energymax;
double ktmp;

double complex *tgrid = rationalInfo->tgrid;
double complex *sn_im = rationalInfo->sn_im;
double complex *sn_c = rationalInfo->sn_c;
double complex *cn_c = rationalInfo->cn_c;
double complex *dn_c = rationalInfo->dn_c;
double complex *z = rationalInfo->z;
double complex *ksi_p = rationalInfo->ksi_p;
double complex *ksi_m = rationalInfo->ksi_m;
double complex *fun_p = rationalInfo->fun_p;
double complex *fun_m = rationalInfo->fun_m;
double complex *rat_fact_p = rationalInfo->rat_fact_p;
double complex *rat_fact_m = rationalInfo->rat_fact_m;

double *tgrid_re  = rationalInfo->tgrid_re;
double *tgrid_im  = rationalInfo->tgrid_im;
double *sn_re  = rationalInfo->sn_re;
double *cn_re  = rationalInfo->cn_re;
double *dn_re  = rationalInfo->dn_re;
double *cn_im  = rationalInfo->cn_im;
double *dn_im  = rationalInfo->dn_im;
double *sn_tmp  = rationalInfo->sn_tmp;
double *cn_tmp  = rationalInfo->cn_tmp;
double *dn_tmp  = rationalInfo->dn_tmp;
double *z_re  = rationalInfo->z_re;
double *z_im  = rationalInfo->z_im;

double complex *zseed = rationalInfo->zseed;
//double preRat = rationalInfo->preRat; 

printf("Energy Max= %.16lg.\n", stodftInfo->energyMax);
printf("Energy Min= %.16lg.\n", stodftInfo->energyMin);
printf("Beta= %.16lg.\n", stodftInfo->beta);
printf("Correct Chemical Potential = %.16lg.\n", stodftInfo->chemPotTrue);


  energymax = stodftInfo->energyMax -  stodftInfo->chemPotTrue;
  mA = (M_PI*M_PI)/(stodftInfo->beta*stodftInfo->beta);
  MA = energymax*energymax + mA;
//MA = stodftInfo->energyMax*stodftInfo->energyMax + mA;
  ktmp = sqrt(MA/mA);
  k = (ktmp - 1.0)/(ktmp + 1.0);
  kinv = 1.0/k;
  m = k*k;

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

/*
for (int i = 0; i < ntgrid; i++){
  if(is_fun_sq == 1) {
    fun_p[i] = fermi_fun_g(ksi_p[i], dmu, stodftInfo->beta, epsilon)*fermi_fun_g(ksi_p[i], dmu, stodftInfo->beta, epsilon);
    fun_m[i] = fermi_fun_g(ksi_m[i], dmu, stodftInfo->beta, epsilon)*fermi_fun_g(ksi_m[i], dmu, stodftInfo->beta, epsilon);
  }else{
    fun_p[i] = fermi_fun_g(ksi_p[i], dmu, stodftInfo->beta, epsilon);
    fun_m[i] = fermi_fun_g(ksi_m[i], dmu, stodftInfo->beta, epsilon);
  }
  fun_p[i] = fun_p[i]*cn_c[i]*dn_c[i] / ((kinv - sn_c[i])*(kinv - sn_c[i])) / ksi_p[i];
  fun_m[i] = fun_m[i]*cn_c[i]*dn_c[i] / ((kinv - sn_c[i])*(kinv - sn_c[i])) / ksi_m[i];
}
*/

for (int j =0; j < ntgrid; j++){
  zseed[j] =  ksi_p[j];
  zseed[j + ntgrid] =  ksi_m[j];
}

for (int i = 0; i < ntgrid; i++){
  rat_fact_p[i] = cn_c[i]*dn_c[i] / ((kinv - sn_c[i])*(kinv - sn_c[i])) / ksi_p[i];
  rat_fact_m[i] = cn_c[i]*dn_c[i] / ((kinv - sn_c[i])*(kinv - sn_c[i])) / ksi_m[i];
}

rationalInfo->preRat = -2 * K *sqrt(mA*MA) / M_PI / ntgrid * kinv;
printf("myidStatesin %lg  \n", rationalInfo->preRat);
}
/*========================================================================================*/
/*========================================================================================*/
/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                               KOMEGA LIBRARY                             */
/*                             Solve shifted COCG                           */
/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
void komega_COCG_init_g(KOMEGAINFO *komegaInfo, int ndim0, int nl0, int nz0, 
                       double *creRex, double *cimRex, double *creImx, double *cimImx,
                       double complex *z0, int itermax0, double threshold0){

  int i;
  //KOMEGAINFO *komegaInfo;
  //komegaInfo = (KOMEGAINFO*)malloc(sizeof(KOMEGAINFO));
  
  komegaInfo->ndim = ndim0;
  komegaInfo->nl = nl0;
  komegaInfo->nz = nz0;
  komegaInfo->itermax = itermax0;
  komegaInfo->threshold = threshold0;
  komegaInfo->almost0 = 1.0E-50;

  //printf("%i %i %i %i %lg\n", komegaInfo->ndim, komegaInfo->nl, komegaInfo->nz, komegaInfo->itermax, komegaInfo->threshold);

  komegaInfo->z = (double complex*)malloc((komegaInfo->nz)*sizeof(double complex));  
  komegaInfo->pi = (double complex*)malloc((komegaInfo->nz)*sizeof(double complex));  
  komegaInfo->pi_old = (double complex*)malloc((komegaInfo->nz)*sizeof(double complex));  
  komegaInfo->lz_conv = (int*)malloc((komegaInfo->nz)*sizeof(int));  
/*----------------------------------------------------*/
//  komegaInfo->p = (double complex*)malloc((komegaInfo->nz*komegaInfo->nl)*sizeof(double complex));  
//  komegaInfo->v3 = (double complex*)malloc((komegaInfo->ndim)*sizeof(double complex));  

//  double complex *v3 = komegaInfo->v3;
//  double complex *p = komegaInfo->p;

  komegaInfo->creRep = (double*)malloc((komegaInfo->nz*komegaInfo->nl)*sizeof(double));
  komegaInfo->cimRep = (double*)malloc((komegaInfo->nz*komegaInfo->nl)*sizeof(double));
  komegaInfo->creImp = (double*)malloc((komegaInfo->nz*komegaInfo->nl)*sizeof(double));
  komegaInfo->cimImp = (double*)malloc((komegaInfo->nz*komegaInfo->nl)*sizeof(double));

  komegaInfo->creRev3 = (double*)malloc((komegaInfo->ndim)*sizeof(double));  
  komegaInfo->cimRev3 = (double*)malloc((komegaInfo->ndim)*sizeof(double));  
  komegaInfo->creImv3 = (double*)malloc((komegaInfo->ndim)*sizeof(double));  
  komegaInfo->cimImv3 = (double*)malloc((komegaInfo->ndim)*sizeof(double));  

  double *creRep = komegaInfo->creRep;
  double *cimRep = komegaInfo->cimRep;
  double *creImp = komegaInfo->creImp;
  double *cimImp = komegaInfo->cimImp;

  double *creRev3 = komegaInfo->creRev3;
  double *cimRev3 = komegaInfo->cimRev3;
  double *creImv3 = komegaInfo->creImv3;
  double *cimImv3 = komegaInfo->cimImv3;
/*----------------------------------------------------*/
  double complex *z = komegaInfo->z;
  double complex *pi = komegaInfo->pi;
  double complex *pi_old = komegaInfo->pi_old;
  int *lz_conv = komegaInfo->lz_conv; 
  
  #pragma omp parallel for private(i)
  for ( i = 0; i < komegaInfo->nz; i++){
    z[i] = z0[i];
  }

  //printf("inside %lg %lg \n", creal(z0[0]) , creal(z[0]));

/*----------------------------------------------------*/
  //TODO2023: g-space
//  #pragma omp parallel for private(i)
//  for ( i = 0; i < komegaInfo->ndim; i++){
//    v3[i] = 0.0 + 0.0*I;
//  }

//  #pragma omp parallel for private(i)
//  for ( i = 0; i < komegaInfo->nl*komegaInfo->nz; i++){
//    p[i] = 0.0 + 0.0*I;
//    x[i] = 0.0 + 0.0*I;
//  }


  #pragma omp parallel for private(i)
  for ( i = 0; i < komegaInfo->ndim; i++){
    creRev3[i] = 0.0;
    cimRev3[i] = 0.0;
    creImv3[i] = 0.0;
    cimImv3[i] = 0.0;
  }


  #pragma omp parallel for private(i)
  for ( i = 0; i < komegaInfo->nl*komegaInfo->nz; i++){
    creRep[i] = 0.0;
    cimRep[i] = 0.0;
    creImp[i] = 0.0;
    cimImp[i] = 0.0;
    creRex[i] = 0.0;
    cimRex[i] = 0.0;
    creImx[i] = 0.0;
    cimImx[i] = 0.0;
  }


/*----------------------------------------------------*/

  #pragma omp parallel for private(i)
  for ( i = 0; i < komegaInfo->nz; i++){
    pi[i] = 1.0 + 0.0*I;
    pi_old[i] = 1.0 + 0.0*I;
  }
 
  komegaInfo->rho = 1.0 + 0.0*I; 
  komegaInfo->alpha = 1.0 + 0.0*I; 
  komegaInfo->beta = 0.0 + 0.0*I; 

  komegaInfo->iz_seed = 0;
  komegaInfo->z_seed = z[komegaInfo->iz_seed];
  komegaInfo->iter = 0;

  #pragma omp parallel for private(i)
  for ( i = 0; i < komegaInfo->nz; i++){
    lz_conv[i] = 0;
  }
  


} /* End komega_COCG_init */
/*==========================================================================*/

/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
void komega_COCG_finalize_g(KOMEGAINFO *komegaInfo){

  // TODO2023: g-space
  free(komegaInfo->z);
  free(komegaInfo->pi);
  free(komegaInfo->lz_conv);
  free(komegaInfo->pi_old);

 // free(komegaInfo->p);
 // free(komegaInfo->v3);

  free(komegaInfo->creRep);
  free(komegaInfo->cimRep);
  free(komegaInfo->creImp);
  free(komegaInfo->cimImp);

  free(komegaInfo->creRev3);
  free(komegaInfo->cimRev3);
  free(komegaInfo->creImv3);
  free(komegaInfo->cimImv3);


} /* End komega_COCG_finalize */
/*==========================================================================*/

/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
void komega_COCG_update_g(KOMEGAINFO *komegaInfo, 
                       double *creRev12, double *cimRev12, double *creImv12, double *cimImv12,
                       double *creRev2, double *cimRev2, double *creImv2, double *cimImv2,
                       double *creRex, double *cimRex, double *creImx, double *cimImx,
                       double *creRer_l, double *cimRer_l, double *creImr_l, double *cimImr_l,
// double complex *v12, double complex *v2, double complex *x, double complex *r_l, 
                       int *status, int numThreads){

  double complex rho_old, alpha_denom;
  double cdotp_sq;
  double complex cdotp;
  double complex conts1, conts2;
  int i,j;
  double complex *pi = komegaInfo->pi; 
  int *lz_conv = komegaInfo->lz_conv; 
  double complex sum;

  double sum1, sum2, sum3, sum4;
  double zreal, zimag;
/*------------------------------------------------------*/
 // double complex *v3 = komegaInfo->v3; 

  double *creRev3 = komegaInfo->creRev3;
  double *cimRev3 = komegaInfo->cimRev3;
  double *creImv3 = komegaInfo->creImv3;
  double *cimImv3 = komegaInfo->cimImv3;

  double conts1_re, conts1_im, alpha_re, alpha_im, conts2_re, conts2_im;
/*------------------------------------------------------*/
  int numCoeff    = komegaInfo->ndim; //TODO CHECK

  double t1, t2, t3, t4, t5, t6;

  t1 = omp_get_wtime();
  //printf("here 0 inside COCG %i \n", numThreads);

  komegaInfo->iter = komegaInfo->iter + 1;

  for (i = 0; i < 3; i++){
    status[i] = 0;
  }


  rho_old = komegaInfo->rho;
  /*------- v2\dot v2-----------*/

  //sum = 0.0 +0.0*I;
  sum1 = 0.0;
  sum2 = 0.0;
  sum3 = 0.0;
  sum4 = 0.0;
  for(i=0;i<numCoeff-1;i++){
    sum1 += creRev2[i]*creRev2[i]+cimRev2[i]*cimRev2[i];
    sum2 += creImv2[i]*creImv2[i]+cimImv2[i]*cimImv2[i];
    sum3 += creRev2[i]*creImv2[i]+cimRev2[i]*cimImv2[i];
  }
  sum1 = sum1*2.0+creRev2[numCoeff-1]*creRev2[numCoeff-1];
  sum2 = sum2*2.0+creImv2[numCoeff-1]*creImv2[numCoeff-1];
  sum3 = sum3*2.0+creRev2[numCoeff-1]*creImv2[numCoeff-1];
  sum = sum1-sum2+2.0*sum3*I;

  //#pragma omp parallel for reduction(+:sum) private(i)
  //for (i = 0; i < komegaInfo->ndim; i++){
  //  sum += v2[i]*v2[i];
  //}
  /*-------end v2\dot v2-----------*/


  komegaInfo->rho = sum;
  //printf("here 1 %lg %lg \n", creal(komegaInfo->rho), cimag(komegaInfo->rho));

  if(komegaInfo->iter == 1){
    komegaInfo->beta =  0.0 +0.0*I;
  }
  else{
    komegaInfo->beta = komegaInfo->rho / rho_old;
  }

  //printf("here 2 %lg %lg \n", creal(komegaInfo->beta), cimag(komegaInfo->beta));
  //printf("here 2 \n");
  
  /*-------update v12--------------*/
  zreal = creal(komegaInfo->z_seed);
  zimag = cimag(komegaInfo->z_seed);
  for(i=0;i<numCoeff;i++){
    creRev12[i] = (zreal*creRev2[i]-zimag*creImv2[i])-creRev12[i];
    cimRev12[i] = (zreal*cimRev2[i]-zimag*cimImv2[i])-cimRev12[i];
    creImv12[i] = (zreal*creImv2[i]+zimag*creRev2[i])-creImv12[i];
    cimImv12[i] = (zreal*cimImv2[i]+zimag*cimRev2[i])-cimImv12[i];
  }


  //#pragma omp parallel for private(i)
  //for (i = 0; i < komegaInfo->ndim; i++){
  //  v12[i] = komegaInfo->z_seed * v2[i] - v12[i];
  //}
  /*-------end update v12--------------*/


  komegaInfo->alpha_old = komegaInfo->alpha;

  /*------- v2\dot v12-----------*/

  sum1 = 0.0;
  sum2 = 0.0;
  sum3 = 0.0;
  sum4 = 0.0;
  for(i=0;i<numCoeff-1;i++){
    sum1 += creRev2[i]*creRev12[i]+cimRev2[i]*cimRev12[i];
    sum2 += creImv2[i]*creImv12[i]+cimImv2[i]*cimImv12[i];
    sum3 += creRev2[i]*creImv12[i]+cimRev2[i]*cimImv12[i];
    sum4 += creRev12[i]*creImv2[i]+cimRev12[i]*cimImv2[i];
  }
  sum1 = sum1*2.0+creRev2[numCoeff-1]*creRev12[numCoeff-1];
  sum2 = sum2*2.0+creImv2[numCoeff-1]*creImv12[numCoeff-1];
  sum3 = sum3*2.0+creRev2[numCoeff-1]*creImv12[numCoeff-1];
  sum4 = sum4*2.0+creRev12[numCoeff-1]*creImv2[numCoeff-1];

  cdotp = sum1-sum2+(sum3+sum4)*I;
  //cdotp = 0.0 +0.0*I;
  //#pragma omp parallel for reduction(+:cdotp) private(i)
  //for (i = 0; i < komegaInfo->ndim; i++){
  //  cdotp += v2[i]*v12[i];
  //}

  /*-------end v2\dot v12-----------*/

  //printf("cdotp %lg %lg %lg %lg \n", sum1, sum2, sum3, sum4);
  //printf("here 2 %lg %lg \n", creal(cdotp), cimag(cdotp));

  //printf("here a %lg %lg \n", creal(komegaInfo->alpha), cimag(komegaInfo->alpha));
  //printf("here b %lg %lg \n", creal(komegaInfo->beta), cimag(komegaInfo->beta));
  //printf("here r %lg %lg \n", creal(komegaInfo->rho), cimag(komegaInfo->rho));

  alpha_denom = cdotp - (komegaInfo->beta * komegaInfo->rho / komegaInfo->alpha) ;
  //printf("here 3 %lg %lg \n", creal(alpha_denom), cimag(alpha_denom));
  //printf("alpha_denom %lg \n", cabs(alpha_denom));
  //printf("here 3 %lg %lg \n", creal(komegaInfo->rho), cimag(komegaInfo->rho));
  //printf("komegaInfo->rho %lg \n", cabs(komegaInfo->rho));

  if( cabs(alpha_denom) < komegaInfo->almost0){
    status[1] = 2;
  }
  else if( cabs(komegaInfo->rho) < komegaInfo->almost0 ){
    status[1] = 4;
  }

  komegaInfo->alpha = komegaInfo->rho / alpha_denom;

  //printf("here 4 %lg %lg \n ", creal(komegaInfo->alpha), cimag(komegaInfo->alpha));

  t2 = omp_get_wtime();
  /* call Shifted equation */

  komega_COCG_shiftedeqn_g(komegaInfo, creRer_l, cimRer_l, creImr_l, cimImr_l,
                           creRex, cimRex, creImx, cimImx,
                           numThreads);

  t3 = omp_get_wtime();
  /* Update residual */

  conts1 = (1.0 + komegaInfo->alpha * komegaInfo->beta / komegaInfo->alpha_old);
  conts2 = komegaInfo->alpha * komegaInfo->beta / komegaInfo->alpha_old;

  //printf("here 5 %lg %lg  %lg %lg \n ", creal(conts1), cimag(conts1), creal(conts2), cimag(conts2));
  /*-------update v12--------------*/
  conts1_re = creal(conts1);
  conts1_im = cimag(conts1);
  alpha_re = creal(komegaInfo->alpha);
  alpha_im = cimag(komegaInfo->alpha);
  conts2_re = creal(conts2);
  conts2_im = cimag(conts2);
  //printf("here debug sum %lg %lg %lg %lg\n",
  //       creRev12[numCoeff-1],creImv12[numCoeff-1],
  //       creRev2[numCoeff-1],creImv2[numCoeff-1]);
  //printf("alpha check %lg %lg\n",alpha_re,alpha_im);
  //printf("last coeff check %lg\n",-alpha_re*creImv12[numCoeff-1]-alpha_im*creRev12[numCoeff-1]);
  double tmp1, tmp2;
  for(i=0;i<numCoeff;i++){ //TODO CHECK
    tmp1 = creRev12[i];
    tmp2 = cimRev12[i];
    creRev12[i] = conts1_re*creRev2[i]-conts1_im*creImv2[i]-
                  alpha_re*creRev12[i]+alpha_im*creImv12[i]-
                  conts2_re*creRev3[i]+conts2_im*creImv3[i];

    cimRev12[i] = conts1_re*cimRev2[i]-conts1_im*cimImv2[i]-
                  alpha_re*cimRev12[i]+alpha_im*cimImv12[i]-
                  conts2_re*cimRev3[i]+conts2_im*cimImv3[i];

    creImv12[i] = conts1_re*creImv2[i]+conts1_im*creRev2[i]-
                  alpha_re*creImv12[i]-alpha_im*tmp1-
                  conts2_re*creImv3[i]-conts2_im*creRev3[i];
    //creImv12[i] = -alpha_re*creImv12[i]-alpha_im*creRev12[i];

    cimImv12[i] = conts1_re*cimImv2[i]+conts1_im*cimRev2[i]-
                  alpha_re*cimImv12[i]-alpha_im*tmp2-
                  conts2_re*cimImv3[i]-conts2_im*cimRev3[i];
  }

  //printf("debug %lg %lg %lg %lg\n", creRev12[numCoeff-1], cimRev12[numCoeff-1], creImv12[numCoeff-1], cimImv12[numCoeff-1]);

  /*------- v12\dot v12-----------*/

  //sum = 0.0 +0.0*I;
  sum1 = 0.0;
  sum2 = 0.0;
  sum3 = 0.0;
  sum4 = 0.0;
  for(i=0;i<numCoeff-1;i++){
    sum1 += creRev12[i]*creRev12[i]+cimRev12[i]*cimRev12[i];
    sum2 += creImv12[i]*creImv12[i]+cimImv12[i]*cimImv12[i];
    sum3 += creRev12[i]*creImv12[i]+cimRev12[i]*cimImv12[i];
  }
  sum1 = sum1*2.0+creRev12[numCoeff-1]*creRev12[numCoeff-1];
  sum2 = sum2*2.0+creImv12[numCoeff-1]*creImv12[numCoeff-1];
  sum3 = sum3*2.0+creRev12[numCoeff-1]*creImv12[numCoeff-1];
  sum = sum1-sum2+2.0*sum3*I;

 // printf("here debug %lg %lg %lg \n", creal(sum), cimag(sum), sum3);
  //printf("debug %lg %lg %lg %lg\n", creRev12[numCoeff-1], cimRev12[numCoeff-1], creImv12[numCoeff-1], cimImv12[numCoeff-1]);
  //#pragma omp parallel for reduction(+:sum) private(i)
  //for (i = 0; i < komegaInfo->ndim; i++){
  //  sum += v2[i]*v2[i];
  //}
  /*-------end v2\dot v2-----------*/



  //#pragma omp parallel for private(i)
  //for (i = 0; i < komegaInfo->ndim; i++){
  //  v12[i] = conts1*v2[i] - komegaInfo->alpha * v12[i] - conts2*v3[i];
  //}
  /*-------end update v12--------------*/

  
  /*-------copy v2 and v12---------*/
  //#pragma omp parallel for private(i)
  //for (i = 0; i < komegaInfo->ndim; i++){
  //  v3[i] = v2[i];
  //}
  //#pragma omp parallel for private(i)
  //for (i = 0; i < komegaInfo->ndim; i++){
  //  v2[i] = v12[i];
  //}
  memcpy(creRev3,creRev2,sizeof(double)*numCoeff);
  memcpy(cimRev3,cimRev2,sizeof(double)*numCoeff);
  memcpy(creImv3,creImv2,sizeof(double)*numCoeff);
  memcpy(cimImv3,cimImv2,sizeof(double)*numCoeff);

  memcpy(creRev2,creRev12,sizeof(double)*numCoeff);
  memcpy(cimRev2,cimRev12,sizeof(double)*numCoeff);
  memcpy(creImv2,creImv12,sizeof(double)*numCoeff);
  memcpy(cimImv2,cimImv12,sizeof(double)*numCoeff);


  /*-------end copy v2 and v12---------*/


  t4 = omp_get_wtime();

  /* Seed Switching  */
  komega_COCG_seed_switch_g(komegaInfo, creRev2, cimRev2, creImv2, cimImv2,status);

  t5 = omp_get_wtime();


  /*------- v2\dot v2-----------*/

  //sum = 0.0 +0.0*I;
  sum1 = 0.0;
  sum2 = 0.0;
  sum3 = 0.0;
  sum4 = 0.0;
  for(i=0;i<numCoeff-1;i++){
    sum1 += creRev2[i]*creRev2[i]+cimRev2[i]*cimRev2[i];
    sum2 += creImv2[i]*creImv2[i]+cimImv2[i]*cimImv2[i];
    sum3 += creRev2[i]*creImv2[i]+cimRev2[i]*cimImv2[i];
  }
  sum1 = sum1*2.0+creRev2[numCoeff-1]*creRev2[numCoeff-1];
  sum2 = sum2*2.0+creImv2[numCoeff-1]*creImv2[numCoeff-1];
  sum3 = sum3*2.0+creRev2[numCoeff-1]*creImv2[numCoeff-1];
  sum = sum1-sum2+2.0*sum3*I;


  //printf("v2 check %lg %lg \n", creal(sum), cimag(sum));

  /* Convergence check  */
  /*--------------c dot p-----------------*/
  //cdotp = 0.0 + 0.0*I;
  //#pragma omp parallel for reduction(+:cdotp) private(i) 
  //for (i = 0; i < komegaInfo->ndim; i++){
  // cdotp += v2[i]* (creal(v2[i]) - cimag(v2[i])*I);
  //}
  sum1 = 0.0;
  sum2 = 0.0;
  for(i=0;i<numCoeff-1;i++){
    sum1 += creRev2[i]*creRev2[i]+cimRev2[i]*cimRev2[i];
    sum2 += creImv2[i]*creImv2[i]+cimImv2[i]*cimImv2[i];
  }
  sum1 = 2.0*sum1+creRev2[numCoeff-1]*creRev2[numCoeff-1];
  sum2 = 2.0*sum2+creImv2[numCoeff-1]*creImv2[numCoeff-1];
  cdotp = sum1+sum2;
 
  cdotp = cdotp*creal(komegaInfo->normfact);
 
  /*--------------end c dot p-----------------*/
    //printf("cdotp %lg %lg \n", creal(cdotp), cimag(cdotp));
    //printf("sum %lg %lg \n", sum1, sum2);

    //v12[0] = sqrt(creal(cdotp)) + 0.0*I; //TODO
    //komegaInfo->resnorm = creal(v12[0]);
    cdotp_sq = sqrt(creal(cdotp));
    komegaInfo->resnorm = cdotp_sq;

    //printf("komegaInfo->normfact = %i \n", komegaInfo->normfact);
    //printf("v12[0] %lg %lg \n", creal(v12[0]), cimag(v12[0]));
   //printf("cdotp_sq %lg \n", cdotp_sq);
   //printf("cabs(cdotp_sq/pi[i]) %lg \n", cabs(cdotp_sq/pi[i]));
  
  //#pragma omp parallel for private(i)
  for (i = 0; i < komegaInfo->nz; i++){
    //if( cabs(v12[0]/pi[i]) < komegaInfo->threshold ) lz_conv[i] = 1;
    if( cabs(cdotp_sq/pi[i]) < komegaInfo->threshold ) lz_conv[i] = 1;
    //printf("cabs(v12[0]/pi[i]) %i %lg \n", i, cabs(cdotp_sq/pi[i]));
  }

  if( cdotp_sq < komegaInfo->threshold){
     /* Converged */
     status[0] = - komegaInfo->iter;
     status[1] = 0;
  }
  else if(komegaInfo->iter == komegaInfo->itermax){
     /* NOT Converged in itermax */
     status[0] = - komegaInfo->iter;
     status[1] = 1;
  }
  else if(status[1] == 2){
     /* alpha becomes infinite */
     status[0] = - komegaInfo->iter;
  }
  else if(status[1] == 3){
     /* pi_seed becomes zero */
     status[0] = - komegaInfo->iter;
  }
  else if(status[1] == 4){
     /* rho becomes zero */
     status[0] = - komegaInfo->iter;
  }
  else {
     /* Continue */
     status[0] = komegaInfo->iter;
     status[1] = 0;
  }

  t6 = omp_get_wtime(); 

  //printf("TIMES %lg %lg %lg %lg %lg \n", t2-t1, t3-t2, t4-t3, t5-t4, t6-t5);

}/* End komega_COCG_update */
/*==========================================================================*/

/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
void komega_COCG_shiftedeqn_g(KOMEGAINFO *komegaInfo,
                              double *creRer_l, double *cimRer_l, double *creImr_l, double *cimImr_l,
                              double *creRex, double *cimRex, double *creImx, double *cimImx,
                              int numThreads) {

int *lz_conv = komegaInfo->lz_conv; 
double complex *z = komegaInfo->z; 
//double complex *p = komegaInfo->p; 
double complex *pi = komegaInfo->pi; 
double complex *pi_old = komegaInfo->pi_old; 
double complex pi_new;
double complex cons;
double complex const1;
double complex const2;
double complex ptemp;
int iz, i;

  double *creRep = komegaInfo->creRep;
  double *cimRep = komegaInfo->cimRep;
  double *creImp = komegaInfo->creImp;
  double *cimImp = komegaInfo->cimImp;

double t1, t2, t3, t4, tsum, tsum2;
double const1_re, const1_im, const2_re, const2_im, cons_re, cons_im;
double ptemp_cre_re, ptemp_cre_im, ptemp_cim_re, ptemp_cim_im;

int numCoeff = komegaInfo->nl;

tsum = 0.0;
tsum2 = 0.0;

  //omp_set_num_threads(4);

  //printf(" threads from COCG_shiftedeqn %i \n", numThreads);

 // #pragma omp parallel
 //   {
 //       printf("Hello from process: %d\n", omp_get_thread_num());
 //   }

  for ( iz = 0; iz < komegaInfo->nz; iz++){
    if(lz_conv[iz] == 1) continue;
  

    pi_new = (1.0 + komegaInfo->alpha * (z[iz] - komegaInfo->z_seed)) * pi[iz]
            - komegaInfo->alpha * komegaInfo->beta / komegaInfo->alpha_old * (pi_old[iz] - pi[iz]);
  
    t1 = omp_get_wtime();
    //printf(" %i %lg %lg \n", iz, creal(pi_new), cimag(pi_new));
    const1 = (pi_old[iz] / pi[iz])*(pi_old[iz] / pi[iz]) * komegaInfo->beta;
    const2 = 1.0/pi[iz];
    cons = pi[iz]/ pi_new * komegaInfo->alpha;
    
    //printf("consts %i %lg %lg %lg %lg \n", iz, creal(const1), cimag(const1), creal(const2), cimag(const2));
    //printf("const2 %i %lg %lg \n", iz, creal(const2), cimag(const2));
    //printf("cons %i %lg %lg \n", iz, creal(cons), cimag(cons));
    /*------updating shifted eq-----------*/
    //#pragma omp parallel for private(i) 
    //for ( i = 0; i < komegaInfo->nl; i++){
    //  //p[i*komegaInfo->nz + iz] = r_l[i] / pi[iz] + const1 * p[i*komegaInfo->nz + iz];
    //  //p[i + iz*komegaInfo->nl] = const2 * r_l[i]  + const1 * p[i + iz*komegaInfo->nl];
    //  ptemp = const2 * r_l[i]  + const1 * p[i + iz*komegaInfo->nl];
    //  x[i + iz*komegaInfo->nl] = x[i + iz*komegaInfo->nl] + cons * ptemp;
    //  p[i + iz*komegaInfo->nl] = ptemp;
    //}
    const1_re = creal(const1);
    const1_im = cimag(const1);
    const2_re = creal(const2);
    const2_im = cimag(const2);
    cons_re = creal(cons);
    cons_im = cimag(cons);
    for(i=0;i<numCoeff;i++){
      ptemp_cre_re = const2_re*creRer_l[i]-const2_im*creImr_l[i]+
                     const1_re*creRep[i+iz*komegaInfo->nl]-const1_im*creImp[i+iz*komegaInfo->nl];

      ptemp_cre_im = const2_re*creImr_l[i]+const2_im*creRer_l[i]+
                     const1_re*creImp[i+iz*komegaInfo->nl]+const1_im*creRep[i+iz*komegaInfo->nl];

      ptemp_cim_re = const2_re*cimRer_l[i]-const2_im*cimImr_l[i]+
                     const1_re*cimRep[i+iz*komegaInfo->nl]-const1_im*cimImp[i+iz*komegaInfo->nl];

      ptemp_cim_im = const2_re*cimImr_l[i]+const2_im*cimRer_l[i]+
                     const1_re*cimImp[i+iz*komegaInfo->nl]+const1_im*cimRep[i+iz*komegaInfo->nl];

      creRex[i+iz*komegaInfo->nl] += cons_re*ptemp_cre_re-cons_im*ptemp_cre_im;
      creImx[i+iz*komegaInfo->nl] += cons_re*ptemp_cre_im+cons_im*ptemp_cre_re;
      cimRex[i+iz*komegaInfo->nl] += cons_re*ptemp_cim_re-cons_im*ptemp_cim_im;
      cimImx[i+iz*komegaInfo->nl] += cons_re*ptemp_cim_im+cons_im*ptemp_cim_re;
 
      creRep[i+iz*komegaInfo->nl] = ptemp_cre_re;
      creImp[i+iz*komegaInfo->nl] = ptemp_cre_im;
      cimRep[i+iz*komegaInfo->nl] = ptemp_cim_re;
      cimImp[i+iz*komegaInfo->nl] = ptemp_cim_im;

    }

    /*-----end updating shifted eq-----------*/
    
    t2 = omp_get_wtime();
    tsum = tsum + t2 - t1;
 
    t3 = omp_get_wtime();
    //cons = pi[iz]/ pi_new * komegaInfo->alpha;
    //#pragma omp parallel for private(i)
    //for ( i = 0; i < komegaInfo->nl; i++){
    //  //x[i*komegaInfo->nz + iz] = x[i*komegaInfo->nz + iz] + cons * p[i*komegaInfo->nz + iz];
    //  x[i + iz*komegaInfo->nl] = x[i + iz*komegaInfo->nl] + cons * p[i + iz*komegaInfo->nl];
    //}
    
  
    t4 = omp_get_wtime();
    tsum2 = tsum2 + t4 - t3;

    pi_old[iz] = pi[iz];
    pi[iz] = pi_new;
 
    //printf("pi_old[iz] %i %lg %lg  \n", iz, creal(pi_old[iz]), cimag(pi_old[iz])); 
  }

  //printf("times %lg %lg \n", tsum, tsum2);
  //printf(" %i %lg %lg \n", iz, creal(x[0]), cimag(x[0]));
  //printf(" %i %lg %lg \n", iz, creal(x[4999]), cimag(x[4999]));
  //printf(" %i %lg %lg \n", iz, creal(x[komegaInfo->nz*komegaInfo->nl-1]), cimag(x[komegaInfo->nz*komegaInfo->nl-1]));

} /* END komega_COCG_shiftedeqn */
/*==========================================================================*/


/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
void komega_COCG_seed_switch_g(KOMEGAINFO *komegaInfo,
                               double *creRev2, double *cimRev2, double *creImv2, double *cimImv2,
                               int *status) {

double complex *z = komegaInfo->z; 
//double complex *v3 = komegaInfo->v3; 
double complex *pi = komegaInfo->pi; 
double complex *pi_old = komegaInfo->pi_old; 
int *lz_conv = komegaInfo->lz_conv;
//int iz_seed = komegaInfo->iz_seed;
int i;
double complex scale;
double minimum;
int location;
double scale_re, scale_im;

double *creRev3 = komegaInfo->creRev3;
double *cimRev3 = komegaInfo->cimRev3;
double *creImv3 = komegaInfo->creImv3;
double *cimImv3 = komegaInfo->cimImv3;

int numCoeff = komegaInfo->ndim; //TODO CHECK

  //printf("%lg %lg %lg \n", cabs(pi[0]), creal(pi[0]), cimag(pi[0]) );

  //status[2] = MINLOC(ABS(pi(1:nz)), 1, .NOT. lz_conv(1:nz));

    minimum = cabs(pi[0]);
    location = 0; 
    for ( i = 0 ; i < komegaInfo->nz ; i++ ) 
    {
        if ( (cabs(pi[i]) < minimum) && (lz_conv[i] == 0) ) 
        {
           minimum = cabs(pi[i]);
           location = i;
        }
    }

    //printf("loc %i %lg \n",  location, minimum);

    status[2] = location;


  if ( cabs(pi[status[2]]) < komegaInfo->almost0 ) {
    status[1] = 3;
  }


  if(status[2] != komegaInfo->iz_seed) {

    komegaInfo->iz_seed = status[2];
    komegaInfo->z_seed = z[komegaInfo->iz_seed];

    komegaInfo->alpha = komegaInfo->alpha * pi_old[komegaInfo->iz_seed] / pi[komegaInfo->iz_seed];
    komegaInfo->rho = komegaInfo->rho / (pi_old[komegaInfo->iz_seed]*pi_old[komegaInfo->iz_seed]);

    scale = 1.0 / pi[komegaInfo->iz_seed];

    //printf("iz_seed  %i %i \n", iz_seed, komegaInfo->iz_seed);
    //printf("scale %lg %lg \n", creal(scale), cimag(scale));
    /*----------scale v2-----------------*/
    //#pragma omp parallel for private(i)
    //for (i = 0; i < komegaInfo->ndim; i++){
    //  v2[i] = scale*v2[i];
    //}
    double tmp1, tmp2;
    scale_re = creal(scale);
    scale_im = cimag(scale);
    for(i=0;i<numCoeff;i++){
      tmp1 = creRev2[i];
      tmp2 = cimRev2[i];
      creRev2[i] = scale_re*creRev2[i]-scale_im*creImv2[i];
      creImv2[i] = scale_re*creImv2[i]+scale_im*tmp1;
      cimRev2[i] = scale_re*cimRev2[i]-scale_im*cimImv2[i];
      cimImv2[i] = scale_re*cimImv2[i]+scale_im*tmp2;
    }

    /*----------end scale v2-----------------*/

    #pragma omp parallel for private(i)
    for (i = 0; i < komegaInfo->nz; i++){
      pi[i] = scale*pi[i];
    }

    //printf("scale %lg %lg \n", creal(scale), cimag(scale));

    scale = 1.0 / pi_old[komegaInfo->iz_seed];
    //printf("scale %lg %lg \n", creal(scale), cimag(scale));
    /*----------scale v3-----------------*/
    //#pragma omp parallel for private(i)
    //for (i = 0; i < komegaInfo->ndim; i++){
    //  v3[i] = scale*v3[i];
    //}
    scale_re = creal(scale);
    scale_im = cimag(scale);
    for(i=0;i<numCoeff;i++){
      tmp1 = creRev3[i];
      tmp2 = cimRev3[i];
      creRev3[i] = scale_re*creRev3[i]-scale_im*creImv3[i];
      creImv3[i] = scale_re*creImv3[i]+scale_im*tmp1;
      cimRev3[i] = scale_re*cimRev3[i]-scale_im*cimImv3[i];
      cimImv3[i] = scale_re*cimImv3[i]+scale_im*tmp2;
    } 
    /*----------end scale v3-----------------*/

    #pragma omp parallel for private(i)
    for (i = 0; i < komegaInfo->nz; i++){
      pi_old[i] = scale*pi_old[i];
    }


  }/* end if*/

  //printf("status %i %i %i \n", status[0], status[1], status[2] );

} /* END komega_COCG_seed_switch  */
/*==========================================================================*/


