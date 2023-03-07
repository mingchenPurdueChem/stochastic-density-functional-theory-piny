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

double complex fermi_fun(double complex x, double dmu, double beta, double epsilon) {
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
void solve_shifted_eqn_cocg( CP *cp, CLASS *class, GENERAL_DATA *general_data, int id, int is_mu_calc) { 
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
  KOMEGAINFO *komegaInfo;
  komegaInfo = (KOMEGAINFO*)malloc(sizeof(KOMEGAINFO));


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

  int itermax = rationalInfo->itermax;
  int ntgrid = rationalInfo->ntgrid;
  double threshold = rationalInfo->threshold;
  double complex *v12 = rationalInfo->v12;
  double complex *v2 = rationalInfo->v2;
  double complex *r_l = rationalInfo->r_l;
  double complex *zseed = rationalInfo->zseed;
  double complex *x = rationalInfo->x;
  double *rhs = rationalInfo->rhs;
/* =------------------------------------------------------------------------------------= */
  int status[3];
  double timeStart1, timeEnd1;
  double timeStart2, timeEnd2;
  double t1, t2, tot, tot1;
  int i, j, iter, jiter, ndim, nz;

  tot = 0.0;
  tot1 = 0.0;
  nz = 2*ntgrid;
  ndim = nfft2; // ndim;
  spinFlag = 0;
  
  printf("--- Starting solving shifted COCG eqn  ---- \n");
  timeStart1 = omp_get_wtime();
  
  genNoiseOrbitalRealRational(cp,cpcoeffs_pos, v2, id); 
  
  if(is_mu_calc == 1) {
    for (i=0; i < nfft2; i++){
      rhs[i] = creal(v2[i]);
    }
  }
  
  //printf("zdim %i dim %i ndim %i \n", zdim, dim, ndim);
  
  timeEnd1 = omp_get_wtime(); 
  timeStart2 = omp_get_wtime();
  //komega_cocg_init(&ndim, &ndim, &nz, x, zseed, &itermax, &threshold, NULL);
  komega_COCG_init(komegaInfo, ndim, ndim, nz, x, zseed, itermax, threshold);
  for (iter = 0; iter < itermax; iter++){
  
    for (i = 0; i < ndim; i++) {
      r_l[i] = v2[i];
    }
    t1 = omp_get_wtime(); 
    calcCoefForceWrapSCFReal(class,general_data,cp,cpcoeffs_pos,clatoms_pos,v2,v12,spinFlag);
    t2 = omp_get_wtime(); 
    tot = tot + (t2-t1); 
  
    t1 = omp_get_wtime(); 
    //komega_cocg_update(v12, v2, x, r_l, status);
    komega_COCG_update(komegaInfo, v12, v2, x, r_l, status);
    t2 = omp_get_wtime(); 
    tot1 = tot1 + (t2-t1); 
  
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
  komega_COCG_finalize(komegaInfo);
  
  timeEnd2 = omp_get_wtime(); 
  printf("--- Finished solving shifted COCG eqn  ---- %lg %lg %lg %lg \n", timeEnd1-timeStart1, timeEnd2-timeStart2, tot, tot1);
}
/*==========================================================================*/

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
  RATIONALINFO *rationalInfo = stodftInfo->rationalInfo;

  int expanType      = stodftInfo->expanType;
  int polynormLength = stodftInfo->polynormLength;
  int numChemPot     = stodftInfo->numChemPot;
  int smearOpt       = stodftInfo->smearOpt;
  int filterDiagFlag = stodftInfo->filterDiagFlag;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
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


  double *expanCoeff = (double*)stodftCoefPos->expanCoeff;
  double timeStart,timeEnd;
  double timeStart2,timeEnd2;
  double timeStart3,timeEnd3;
  double deltaTime = 0.0;
  double deltaTime2 = 0.0;
  double dot;


  int nfft          = cp_para_fft_pkg3d_lg->nfft;
  int nfft2         = nfft/2;
  double complex sum; 
  int ntgrid = rationalInfo->ntgrid;

  double complex *fun_p = rationalInfo->fun_p;
  double complex *fun_m = rationalInfo->fun_m;
  double complex *x = rationalInfo->x;
  //double preRat = rationalInfo->preRat;
  double *frhs = rationalInfo->frhs;

  printf("start tessssssssssssst RA \n");

  printf("myidState = %i %i %lg %lg %i %lg %lg %lg %lg %i %lg \n", 
   myidState,rationalInfo->ntgrid, rationalInfo->dmu, rationalInfo->threshold,
   rationalInfo->itermax, rationalInfo->large_dmu,
   rationalInfo->small_dmu, rationalInfo->maxmu, rationalInfo->epsilon, expanType,
   rationalInfo->init_mu);

/******************************************************************************/
printf("Starting Filtering with Rational Approximation\n");
/******************************************************************************/

init_zseed(cp, 0);

for(iState=0;iState<numStateUpProc;iState++){

  //printf("myidState %i %i %i \n", myidState, numStateUpProc, iState);
  solve_shifted_eqn_cocg( cp, class, general_data, iState, 0);

  for (int i = 0; i<nfft2;i++){
    sum = 0.0 + 0.0 *I ;
    for (int j =0; j < ntgrid; j++){
      sum = sum + fun_p[j] * x[j*nfft2 + i] + fun_m[j] * x[(j+ntgrid)*nfft2 + i];
    }
    frhs[i] = rationalInfo->preRat*cimag(sum);
  }

  rhsReal(class, general_data, cp, cpcoeffs_pos, clatoms_pos, frhs, iState);
}

/******************************************************************************/
printf("Finished Filtering with Rational Approximation\n");
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
init_zseed(cp, 1);

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
  solve_shifted_eqn_cocg( cp, class, general_data, iState, 1);

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
init_zseed(cp, 1);

dsum = 0.0;
for(iState=0;iState<numStateUpProc;iState++){

  //printf("myidState %i %i %i %i %i \n", myidState, numStateUpProc, iState, numProcStates, numStateStoUp);
  solve_shifted_eqn_cocg( cp, class, general_data, iState, 1);

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
init_zseed(cp, 1);

dsum = 0.0;
for(iState=0;iState<numStateUpProc;iState++){

  //printf("myidState %i %i %i %i %i \n", myidState, numStateUpProc, iState, numProcStates, numStateStoUp);
  solve_shifted_eqn_cocg( cp, class, general_data, iState, 1);

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
   fmu = -(NElecTot - numElecTrue)*dNdm; 

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
void init_zseed( CP *cp, int is_fun_sq) {

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


for (int i = 0; i < ntgrid; i++){
  if(is_fun_sq == 1) {
    fun_p[i] = fermi_fun(ksi_p[i], dmu, stodftInfo->beta, epsilon)*fermi_fun(ksi_p[i], dmu, stodftInfo->beta, epsilon);
    fun_m[i] = fermi_fun(ksi_m[i], dmu, stodftInfo->beta, epsilon)*fermi_fun(ksi_m[i], dmu, stodftInfo->beta, epsilon);
  }else{
    fun_p[i] = fermi_fun(ksi_p[i], dmu, stodftInfo->beta, epsilon);
    fun_m[i] = fermi_fun(ksi_m[i], dmu, stodftInfo->beta, epsilon);
  }
  fun_p[i] = fun_p[i]*cn_c[i]*dn_c[i] / ((kinv - sn_c[i])*(kinv - sn_c[i])) / ksi_p[i];
  fun_m[i] = fun_m[i]*cn_c[i]*dn_c[i] / ((kinv - sn_c[i])*(kinv - sn_c[i])) / ksi_m[i];
}

for (int j =0; j < ntgrid; j++){
  zseed[j] =  ksi_p[j];
  zseed[j + ntgrid] =  ksi_m[j];
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
void komega_COCG_init(KOMEGAINFO *komegaInfo, int ndim0, int nl0, int nz0, double complex *x, double complex *z0, int itermax0, double threshold0){

  int i;
  //KOMEGAINFO *komegaInfo;
  //komegaInfo = (KOMEGAINFO*)malloc(sizeof(KOMEGAINFO));
  
  komegaInfo->ndim = ndim0;
  komegaInfo->nl = nl0;
  komegaInfo->nz = nz0;
  komegaInfo->itermax = itermax0;
  komegaInfo->threshold = threshold0;
  komegaInfo->almost0 = 1.0E-50;

  printf("%i %i %i %i %lg\n", komegaInfo->ndim, komegaInfo->nl, komegaInfo->nz, komegaInfo->itermax, komegaInfo->threshold);

  komegaInfo->z = (double complex*)malloc((komegaInfo->nz)*sizeof(double complex));  
  komegaInfo->v3 = (double complex*)malloc((komegaInfo->ndim)*sizeof(double complex));  
  komegaInfo->pi = (double complex*)malloc((komegaInfo->nz)*sizeof(double complex));  
  komegaInfo->pi_old = (double complex*)malloc((komegaInfo->nz)*sizeof(double complex));  
  komegaInfo->p = (double complex*)malloc((komegaInfo->nz*komegaInfo->nl)*sizeof(double complex));  
  komegaInfo->lz_conv = (int*)malloc((komegaInfo->nz)*sizeof(int));  

  //printf("here 000 \n");

  double complex *z = komegaInfo->z;
  double complex *v3 = komegaInfo->v3;
  double complex *pi = komegaInfo->pi;
  double complex *pi_old = komegaInfo->pi_old;
  double complex *p = komegaInfo->p;
  int *lz_conv = komegaInfo->lz_conv; 
  

  for ( i = 0; i < komegaInfo->nz; i++){
    z[i] = z0[i];
  }

  //printf("inside %lg %lg \n", creal(z0[0]) , creal(z[0]));

  for ( i = 0; i < komegaInfo->ndim; i++){
    v3[i] = 0.0 + 0.0*I;
  }

  for ( i = 0; i < komegaInfo->nl*komegaInfo->nz; i++){
    p[i] = 0.0 + 0.0*I;
    x[i] = 0.0 + 0.0*I;
  }

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

  for ( i = 0; i < komegaInfo->nz; i++){
    lz_conv[i] = 0;
  }
  


} /* End komega_COCG_init */
/*==========================================================================*/

/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
void komega_COCG_finalize(KOMEGAINFO *komegaInfo){


  free(komegaInfo->z);
  free(komegaInfo->v3);
  free(komegaInfo->pi);
  free(komegaInfo->p);
  free(komegaInfo->lz_conv);
  //free(komegaInfo->pi_old);

} /* End komega_COCG_finalize */
/*==========================================================================*/

/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
void komega_COCG_update(KOMEGAINFO *komegaInfo, double complex *v12, double complex *v2, double complex *x, double complex *r_l, int *status){

  double complex rho_old, alpha_denom, cdotp;
  double complex conts1, conts2;
  int i,j;
  double complex *v3 = komegaInfo->v3; 
  double complex *pi = komegaInfo->pi; 
  int *lz_conv = komegaInfo->lz_conv; 

  komegaInfo->iter = komegaInfo->iter + 1;

  //printf("here 0 %i \n", komegaInfo->iter);
  for (i = 0; i < 3; i++){
    status[i] = 0;
  }

  //printf("here 1 \n");

  rho_old = komegaInfo->rho;

  komegaInfo->rho = 0.0 +0.0*I;
  for (i = 0; i < komegaInfo->ndim; i++){
    komegaInfo->rho = komegaInfo->rho + v2[i]*v2[i];
  }
  //printf("here 1 %lg %lg \n", creal(komegaInfo->rho), cimag(komegaInfo->rho));

  if(komegaInfo->iter == 1){
    komegaInfo->beta =  0.0 +0.0*I;
  }
  else{
    komegaInfo->beta = komegaInfo->rho / rho_old;
  }

  //printf("here 2 \n");

  for (i = 0; i < komegaInfo->ndim; i++){
    v12[i] = komegaInfo->z_seed * v2[i] - v12[i];
  }

  komegaInfo->alpha_old = komegaInfo->alpha;

  cdotp = 0.0 +0.0*I;
  for (i = 0; i < komegaInfo->ndim; i++){
    cdotp = cdotp + v2[i]*v12[i];
  }

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

  /* call Shifted equation */

  komega_COCG_shiftedeqn(komegaInfo, r_l, x);

  /* Update residual */

  conts1 = (1.0 + komegaInfo->alpha * komegaInfo->beta / komegaInfo->alpha_old);
  conts2 = komegaInfo->alpha * komegaInfo->beta / komegaInfo->alpha_old;

  for (i = 0; i < komegaInfo->ndim; i++){
    v12[i] = conts1*v2[i] - komegaInfo->alpha * v12[i] - conts2*v3[i];
  }

  for (i = 0; i < komegaInfo->ndim; i++){
    v3[i] = v2[i];
    v2[i] = v12[i];
  }

 // printf("v2 %lg %lg \n", creal(v2[0]), cimag(v2[0]));
 // printf("v2 %lg %lg \n", creal(v2[9]), cimag(v2[9]));
 // printf("v2 %lg %lg \n", creal(v2[99]), cimag(v2[99]));

 // printf("v3 %lg %lg \n", creal(v3[0]), cimag(v3[0]));
 // printf("v3 %lg %lg \n", creal(v3[9]), cimag(v3[9]));
 // printf("v3 %lg %lg \n", creal(v3[99]), cimag(v3[99]));
  /* Seed Switching  */
  komega_COCG_seed_switch(komegaInfo, v2,status);

  /* Convergence check  */
  cdotp = 0.0 + 0.0*I;
  for (i = 0; i < komegaInfo->ndim; i++){
   cdotp = cdotp + v2[i]* (creal(v2[i]) - cimag(v2[i])*I);
  }

    //printf(" %lg %lg \n", creal(cdotp), cimag(cdotp));

    v12[0] = sqrt(creal(cdotp)) + 0.0*I; //TODO
    komegaInfo->resnorm = creal(v12[0]);

  for (i = 0; i < komegaInfo->nz; i++){
    if( cabs(v12[0]/pi[i]) < komegaInfo->threshold ) lz_conv[i] = 1;
  }

  if( creal(v12[0]) < komegaInfo->threshold){
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


}/* End komega_COCG_update */
/*==========================================================================*/

/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
void komega_COCG_shiftedeqn(KOMEGAINFO *komegaInfo, double complex *r_l, double complex *x) {

int *lz_conv = komegaInfo->lz_conv; 
double complex *z = komegaInfo->z; 
double complex *p = komegaInfo->p; 
double complex *pi = komegaInfo->pi; 
double complex *pi_old = komegaInfo->pi_old; 
double complex pi_new;
double complex cons;
int iz, i;

  for ( iz = 0; iz < komegaInfo->nz; iz++){
    if(lz_conv[iz] == 1) continue;
  
    pi_new = (1.0 + komegaInfo->alpha * (z[iz] - komegaInfo->z_seed)) * pi[iz]
            - komegaInfo->alpha * komegaInfo->beta / komegaInfo->alpha_old * (pi_old[iz] - pi[iz]);
  
    //printf(" %i %lg %lg \n", iz, creal(pi_new), cimag(pi_new));
    for ( i = 0; i < komegaInfo->nl; i++){
      p[i*komegaInfo->nz + iz] = r_l[i] / pi[iz] + (pi_old[iz] / pi[iz])*(pi_old[iz] / pi[iz]) * komegaInfo->beta * p[i*komegaInfo->nz + iz];
    }
  
    cons = pi[iz]/ pi_new * komegaInfo->alpha;
    for ( i = 0; i < komegaInfo->nl; i++){
      x[i*komegaInfo->nz + iz] = x[i*komegaInfo->nz + iz] + cons * p[i*komegaInfo->nz + iz];
    }
  
    pi_old[iz] = pi[iz];
    pi[iz] = pi_new;
  
  }

  //printf(" %i %lg %lg \n", iz, creal(x[0]), cimag(x[0]));
  //printf(" %i %lg %lg \n", iz, creal(x[4999]), cimag(x[4999]));
  //printf(" %i %lg %lg \n", iz, creal(x[komegaInfo->nz*komegaInfo->nl-1]), cimag(x[komegaInfo->nz*komegaInfo->nl-1]));

} /* END komega_COCG_shiftedeqn */
/*==========================================================================*/


/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
void komega_COCG_seed_switch(KOMEGAINFO *komegaInfo, double complex *v2, int *status) {

double complex *z = komegaInfo->z; 
double complex *v3 = komegaInfo->v3; 
double complex *pi = komegaInfo->pi; 
double complex *pi_old = komegaInfo->pi_old; 
int *lz_conv = komegaInfo->lz_conv;
//int iz_seed = komegaInfo->iz_seed;
int i;
double complex scale;
double minimum;
int location;

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
    for (i = 0; i < komegaInfo->ndim; i++){
      v2[i] = scale*v2[i];
    }
    for (i = 0; i < komegaInfo->nz; i++){
      pi[i] = scale*pi[i];
    }

    //printf("scale %lg %lg \n", creal(scale), cimag(scale));

    scale = 1.0 / pi_old[komegaInfo->iz_seed];
    //printf("scale %lg %lg \n", creal(scale), cimag(scale));
    for (i = 0; i < komegaInfo->ndim; i++){
      v3[i] = scale*v3[i];
    }
    for (i = 0; i < komegaInfo->nz; i++){
      pi_old[i] = scale*pi_old[i];
    }


  }/* end if*/

  //printf("status %i %i %i \n", status[0], status[1], status[2] );

} /* END komega_COCG_seed_switch  */
/*==========================================================================*/


