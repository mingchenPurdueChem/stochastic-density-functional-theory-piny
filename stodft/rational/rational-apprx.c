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

  printf("myidState = %i %i %lg %lg %i %lg %lg \n", 
   myidState,rationalInfo->ntgrid, rationalInfo->dmu, rationalInfo->threshold,
   rationalInfo->itermax, rationalInfo->large_dmu,  rationalInfo->small_dmu);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/
