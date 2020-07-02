/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: filter.c                                       */
/*                                                                          */
/* This routine constructs filters of Fermi function.                       */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_stodft_local.h"

#define TIME_CP_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void filterNewtonPolyHerm(CP *cp,CLASS *class,GENERAL_DATA *general_data,
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
  int imu,iCoeff,iPoly,indexStart;
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


  double *expanCoeff = (double*)stodftCoefPos->expanCoeff;
  double timeStart,timeEnd; 
  double timeStart2,timeEnd2;
  double timeStart3,timeEnd3;
  double deltaTime = 0.0;
  double deltaTime2 = 0.0;

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
//debug
/*  
  int numPointTest = 100;
  int iPoint;
  double pointTest;
  double energyMinTest = -0.1150735;
  double energyMaxTest = 0.10104944;
  double deltPoint = (energyMaxTest-energyMinTest)/numPointTest;
  double pointScale;
  double funValue,prod;
  for(iPoint=0;iPoint<numPointTest;iPoint++){
    pointTest = energyMinTest+(iPoint+0.5)*deltPoint;
    pointScale = (pointTest-energyMean)*scale;
    funValue = expanCoeff[0];
    prod = 1.0;
    for(iPoly=1;iPoly<polynormLength;iPoly++){
      prod *= pointScale-sampPoint[iPoly-1];
      funValue += expanCoeff[iPoly]*prod;
    }
    printf("TestFunExpan %lg %lg %lg\n",pointTest,pointScale,funValue);
  }
  fflush(stdout);
  exit(0);
*/
//debug

  if(storeChebyMomentsFlag){
    stodftCoefPos->chebyMomentsUp = (double*)cmalloc((polynormLength+1)*sizeof(double));
    if(cpLsda==1&&numStateDnProc!=0){
      stodftCoefPos->chebyMomentsDn = (double*)cmalloc((polynormLength+1)*sizeof(double));
    }
  }

  omp_set_num_threads(numThreads);
 
/*==========================================================================*/
/* 0) Copy the initial stochastic orbital */
 
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
      entropyUpRe[iCoeff] = entropyExpanCoeff[imu]*cre_up[iCoeff];
      entropyUpIm[iCoeff] = entropyExpanCoeff[imu]*cim_up[iCoeff];
    }//endfor iCoeff
    if(cpLsda==1&&numStateDnProc!=0){
      #pragma omp parallel for private(iCoeff)
      for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
        entropyDnRe[iCoeff] = entropyExpanCoeff[imu]*cre_dn[iCoeff];
        entropyDnIm[iCoeff] = entropyExpanCoeff[imu]*cim_dn[iCoeff];
      }//endfor iCoeff
    }//endif cpLsda
  }//endif smearOpt

/*==========================================================================*/
/* 1) Loop over all polynomial terms (iPoly=0<=>polynomial order 1) */

  timeStart = omp_get_wtime();

  for(iPoly=1;iPoly<polynormLength;iPoly++){
    if(iPoly%1000==0&&myidState==0){
      printf("%lg%% ",iPoly*100.0/polynormLength);
      fflush(stdout);
    }
    timeStart3 = omp_get_wtime();  
    normHNewtonHerm(cp,class,general_data,
                 cpcoeffs_pos,clatoms_pos,sampPoint[iPoly-1]);
    timeEnd3 = omp_get_wtime();
    deltaTime2 += timeEnd3-timeStart3;
    timeStart2 = omp_get_wtime();
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
  Barrier(comm_states);

  if(myidState==0){
    printf("Process ID %i filter time %.8lg total NormH time %.8lg Accumulate P(H)|phi> time %.8lg\n",myidState,timeProc,deltaTime2,deltaTime);
    printf("Process ID %i Nlpp %.8lg Apply-KS-pot %.8lg Pack-fft %.8lg Unpack-fft %.8lg kinetic %.8lg scale-H|phi> %.8lg\n",myidState,stodftInfo->cputime_new[0],stodftInfo->cputime3,stodftInfo->cputime2,stodftInfo->cputime4,stodftInfo->cputime5,stodftInfo->cputime7);
    printf("Process ID Nlpp-part1 %.8lg Nlpp-part2 %.8lg\n",stodftInfo->cputime0,stodftInfo->cputime1);
    printf("Process ID %i FFTW3D time %.8lg FFTW3D-to-r-pre %.8lg FFTW3D-to-r-post %.8lg FFTW3D-to-g-pre %.8lg FFTW3D-to-g-post %.8lg\n",myidState,cp_sclr_fft_pkg3d_sm->cputime,cp_sclr_fft_pkg3d_sm->cputime1,cp_sclr_fft_pkg3d_sm->cputime2,cp_sclr_fft_pkg3d_sm->cputime3,cp_sclr_fft_pkg3d_sm->cputime4);
    printf("Process ID %i cpy-wf %.8lg %.8lg add force %.8lg %.8lg prepare-loc %.8lg zero %.8lg\n",myidState,stodftInfo->cputime_new[1],stodftInfo->cputime_new[3],stodftInfo->cputime_new[2],stodftInfo->cputime_new[4],stodftInfo->cputime_new[5],stodftInfo->cputime_new[6]);
   fflush(stdout);
  }

  Barrier(comm_states);

  //debug
  /*
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    printf("stowftesttest %.16lg %.16lg\n",stoWfUpRe[0][iCoeff],stoWfUpIm[0][iCoeff]);
  }
  */

/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void filterChebyPolyHerm(CP *cp,CLASS *class,GENERAL_DATA *general_data,
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
//debug
/*  
  int numPointTest = 100;
  int iPoint;
  double pointTest;
  double energyMinTest = -0.1150735;
  double energyMaxTest = 0.10104944;
  double deltPoint = (energyMaxTest-energyMinTest)/numPointTest;
  double pointScale;
  double funValue,prod;
  for(iPoint=0;iPoint<numPointTest;iPoint++){
    pointTest = energyMinTest+(iPoint+0.5)*deltPoint;
    pointScale = (pointTest-energyMean)*scale;
    funValue = expanCoeff[0];
    prod = 1.0;
    for(iPoly=1;iPoly<polynormLength;iPoly++){
      prod *= pointScale-sampPoint[iPoly-1];
      funValue += expanCoeff[iPoly]*prod;
    }
    printf("TestFunExpan %lg %lg %lg\n",pointTest,pointScale,funValue);
  }
  fflush(stdout);
  exit(0);
*/
//debug

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
      entropyUpRe[iCoeff] = entropyExpanCoeff[imu]*cre_up[iCoeff];
      entropyUpIm[iCoeff] = entropyExpanCoeff[imu]*cim_up[iCoeff];
    }//endfor iCoeff
    if(cpLsda==1&&numStateDnProc!=0){
      #pragma omp parallel for private(iCoeff)
      for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
        entropyDnRe[iCoeff] = entropyExpanCoeff[imu]*cre_dn[iCoeff];
        entropyDnIm[iCoeff] = entropyExpanCoeff[imu]*cim_dn[iCoeff];
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
  Barrier(comm_states);

  if(myidState==0){
    printf("Process ID %i filter time %.8lg total NormH time %.8lg Accumulate P(H)|phi> time %.8lg\n",myidState,timeProc,deltaTime2,deltaTime);
    printf("Process ID %i Nlpp %.8lg Apply-KS-pot %.8lg Pack-fft %.8lg Unpack-fft %.8lg kinetic %.8lg scale-H|phi> %.8lg\n",myidState,stodftInfo->cputime_new[0],stodftInfo->cputime3,stodftInfo->cputime2,stodftInfo->cputime4,stodftInfo->cputime5,stodftInfo->cputime7);
    printf("Process ID Nlpp-part1 %.8lg Nlpp-part2 %.8lg\n",stodftInfo->cputime0,stodftInfo->cputime1);
    printf("Process ID %i FFTW3D time %.8lg FFTW3D-to-r-pre %.8lg FFTW3D-to-r-post %.8lg FFTW3D-to-g-pre %.8lg FFTW3D-to-g-post %.8lg\n",myidState,cp_sclr_fft_pkg3d_sm->cputime,cp_sclr_fft_pkg3d_sm->cputime1,cp_sclr_fft_pkg3d_sm->cputime2,cp_sclr_fft_pkg3d_sm->cputime3,cp_sclr_fft_pkg3d_sm->cputime4);
    printf("Process ID %i cpy-wf %.8lg %.8lg add force %.8lg %.8lg prepare-loc %.8lg zero %.8lg\n",myidState,stodftInfo->cputime_new[1],stodftInfo->cputime_new[3],stodftInfo->cputime_new[2],stodftInfo->cputime_new[4],stodftInfo->cputime_new[5],stodftInfo->cputime_new[6]);
   fflush(stdout);
  }

  Barrier(comm_states);

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


#ifdef FAST_FILTER   

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void filterNewtonPolyHermFake(CP *cp,CLASS *class,GENERAL_DATA *general_data,
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
  int imu,iCoeff,iPoly,indexStart,iState,jState;
  int index1,index2,index3;
  int startIndex;

  int numStatePrintUp = stodftInfo->numStatePrintUp;
  int numStatePrintDn = stodftInfo->numStatePrintDn;

  MPI_Comm comm_states   =    communicate->comm_states;


  double energyDiff  = stodftInfo->energyDiff;
  double energyMin   = stodftInfo->energyMin;
  double energyMax   = stodftInfo->energyMax;
  double energyMean  = stodftInfo->energyMean;
  double scale       = newtonInfo->scale;
  double prefact;
  double polyCoeff;
  double prod,sum;
  double x,y,z;
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
  double *moUpRePrint = stodftCoefPos->moUpRePrint;
  double *moUpImPrint = stodftCoefPos->moUpImPrint;
  double *entropyUpRe = stodftCoefPos->entropyUpRe;
  double *entropyUpIm = stodftCoefPos->entropyUpIm;
  double *entropyDnRe = stodftCoefPos->entropyDnRe;
  double *entropyDnIm = stodftCoefPos->entropyDnIm;
  double *entropyExpanCoeff = stodftCoefPos->entropyExpanCoeff;


  double *expanCoeff = (double*)stodftCoefPos->expanCoeff;
  double timeStart,timeEnd; 
  double timeStart2,timeEnd2;
  double timeStart3,timeEnd3;
  double deltaTime = 0.0;
  double deltaTime2 = 0.0;

  double *energyLevel = stodftCoefPos->energyLevel;
  double *occupNumber;
  double *energyScale;
  double *entropyState;
  double *wfDot;


/*==========================================================================*/
/* i) Generate occupatation number */

  occupNumber = (double *)cmalloc(numChemPot*numStatePrintUp*sizeof(double));
  entropyState = (double *)cmalloc(numStatePrintUp*sizeof(double));
  energyScale = (double *)cmalloc(numStatePrintUp*sizeof(double));
  wfDot = (double *)cmalloc(numStateUpProc*numStatePrintUp*sizeof(double));

  for(iState=0;iState<numStatePrintUp;iState++){
    energyScale[iState] = (energyLevel[iState]-energyMean)*scale;
    //printf("iState %i %lg\n",iState,energyScale[iState]);
  }

  for(iState=0;iState<numStatePrintUp;iState++){
    prod = 1.0;
    for(imu=0;imu<numChemPot;imu++){
      occupNumber[imu*numStatePrintUp+iState] = expanCoeff[imu];
    }
    for(iPoly=1;iPoly<polynormLength;iPoly++){
      prod *= energyScale[iState]-sampPoint[iPoly-1];
      for(imu=0;imu<numChemPot;imu++){
        polyCoeff = expanCoeff[iPoly*numChemPot+imu];
        occupNumber[imu*numStatePrintUp+iState] += polyCoeff*prod;
      }//endfor imu
    }//endfor iPoly
    if(smearOpt>0&&filterDiagFlag==0){
      prod = 1.0;
      entropyState[iState] = entropyExpanCoeff[0];
      for(iPoly=1;iPoly<polynormLength;iPoly++){
        prod *= energyScale[iState]-sampPoint[iPoly-1];
        polyCoeff = entropyExpanCoeff[iPoly];
        entropyState[iState] += polyCoeff*prod;
      }//endfor iPoly
    }//endif smearOpt
  }//endfor iState

  /*
  for(imu=0;imu<numChemPot;imu++){
    for(iState=0;iState<numStatePrintUp;iState++){
      printf("imu %i iState %i %lg\n",imu,iState,occupNumber[imu*numStatePrintUp+iState]);
    }
  }
  */
  /*
  for(iPoly=1;iPoly<polynormLength;iPoly++){
    printf("coeff iPoly %lg\n",expanCoeff[iPoly*numChemPot]);
  }
  */

/*==========================================================================*/
/* ii) Filtering stochastic orbital */

  #pragma omp parallel for private(iState,jState,iCoeff,index1,index2,sum)
  for(iState=0;iState<numStateUpProc;iState++){
    index1 = iState*numCoeff;
    for(jState=0;jState<numStatePrintUp;jState++){
      index2 = jState*numCoeff;
      sum = 0.0;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
        sum += cre_up[index1+iCoeff]*moUpRePrint[index2+iCoeff]+
               cim_up[index1+iCoeff]*moUpImPrint[index2+iCoeff];
      }
      sum *= 2.0;
      sum += cre_up[index1+numCoeff]*moUpRePrint[index2+numCoeff];
      wfDot[iState*numStatePrintUp+jState] = sum;
    }
  }

  #pragma omp parallel for private(imu,iState,jState,iCoeff,x,y)
  for(imu=0;imu<numChemPot;imu++){
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      stoWfUpRe[imu][iCoeff] = 0.0;
      stoWfUpIm[imu][iCoeff] = 0.0;
    }
    for(iState=0;iState<numStateUpProc;iState++){
      for(jState=0;jState<numStatePrintUp;jState++){
        x = occupNumber[imu*numStatePrintUp+jState];
        y = wfDot[iState*numStatePrintUp+jState]*x;
        for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
          stoWfUpRe[imu][iState*numCoeff+iCoeff] += y*moUpRePrint[jState*numCoeff+iCoeff];
          stoWfUpIm[imu][iState*numCoeff+iCoeff] += y*moUpImPrint[jState*numCoeff+iCoeff];
        }
      }
    }
  }

  if(smearOpt>0&&filterDiagFlag==0){
    #pragma omp parallel for private(iCoeff)
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      entropyUpRe[iCoeff] = 0.0;
      entropyUpIm[iCoeff] = 0.0;
    }
    #pragma omp parallel for private(iState,jState,iCoeff,x,y)
    for(iState=0;iState<numStateUpProc;iState++){
      for(jState=0;jState<numStatePrintUp;jState++){
        x = entropyState[jState];
        y = wfDot[iState*numStatePrintUp+jState]*x;
        for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
          entropyUpRe[iState*numCoeff+iCoeff] += y*moUpRePrint[jState*numCoeff+iCoeff];
          entropyUpIm[iState*numCoeff+iCoeff] += y*moUpImPrint[jState*numCoeff+iCoeff];
        }//endfor iCoeff
      }//endfor jState
    }//endfor iState
    //printf("entropyUpRe %lg\n",entropyUpRe[1]);
  }//endif smearOpt


  free(occupNumber);
  free(energyScale);
  free(wfDot);
  free(entropyState);
 
  Barrier(comm_states);

/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void filterChebyPolyHermFake(CP *cp,CLASS *class,GENERAL_DATA *general_data,
			  int ip_now,int inverseFlag)
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

  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;

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
  int imu,iCoeff,iPoly,indexStart,iState,jState;
  int index1,index2,index3;
  int startIndex;
  int storeChebyMomentsFlag = stodftInfo->storeChebyMomentsFlag;
  int i,j,k;

  int numStatePrintUp = stodftInfo->numStatePrintUp;
  int numStatePrintDn = stodftInfo->numStatePrintDn;

  MPI_Comm comm_states   =    communicate->comm_states;

  double energyDiff  = stodftInfo->energyDiff;
  double energyMin   = stodftInfo->energyMin;
  double energyMax   = stodftInfo->energyMax;
  double energyMean  = stodftInfo->energyMean;
  double scale       = chebyshevInfo->scale;
  double prefact;
  double polyCoeff;
  double prod,sum;
  double x,y,z;
  double t0,t1,t2;
  double timeProc,timeTot;
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;
  double *moUpRePrint = stodftCoefPos->moUpRePrint;
  double *moUpImPrint = stodftCoefPos->moUpImPrint;
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

  double *energyLevel = stodftCoefPos->energyLevel;
  double *occupNumber;
  double *energyScale;
  double *entropyState;
  double *wfDot;

/*==========================================================================*/
/* i) Generate occupatation number */

  occupNumber = (double *)cmalloc(numChemPot*numStatePrintUp*sizeof(double));
  entropyState = (double *)cmalloc(numStatePrintUp*sizeof(double));
  energyScale = (double *)cmalloc(numStatePrintUp*sizeof(double));
  wfDot = (double *)cmalloc(numStateUpProc*numStatePrintUp*sizeof(double));

  for(iState=0;iState<numStatePrintUp;iState++){
    energyScale[iState] = (energyLevel[iState]-energyMean)*scale;
    //printf("iState %i %lg\n",iState,energyScale[iState]);
  }

  for(iState=0;iState<numStatePrintUp;iState++){
    //iPoly=0 and iPoly=1
    t0 = 1;
    t1 = energyScale[iState];
    for(imu=0;imu<numChemPot;imu++){
      occupNumber[imu*numStatePrintUp+iState] = expanCoeff[imu]
                                 +expanCoeff[numChemPot+imu]*t1;
    }
    if(smearOpt>0&&filterDiagFlag==0){
      entropyState[iState] = entropyExpanCoeff[0]+entropyExpanCoeff[1]*t1;
    }
    for(iPoly=2;iPoly<polynormLength;iPoly++){
      t2 = 2*energyScale[iState]*t1-t0;
      for(imu=0;imu<numChemPot;imu++){
        polyCoeff = expanCoeff[iPoly*numChemPot+imu];
        occupNumber[imu*numStatePrintUp+iState] += polyCoeff*t2;
      }//endfor imu
      if(smearOpt>0&&filterDiagFlag==0){
        polyCoeff = entropyExpanCoeff[iPoly];
        entropyState[iState] += polyCoeff*t2;
      }
      t0 = t1;
      t1 = t2;
    }//endfor iPoly
  }//endfor iState

  //DEBUG
  /*
  if(myidState==0){
    for(imu=0;imu<numChemPot;imu++){
      for(iState=0;iState<numStatePrintUp;iState++){
        printf("iiiiiiiimu %i iState %i %lg\n",imu,iState,occupNumber[imu*numStatePrintUp+iState]);
      }
    }
  }
  */
  
  
  /*
  for(iPoly=1;iPoly<polynormLength;iPoly++){
    printf("coeff iPoly %lg\n",expanCoeff[iPoly*numChemPot]);
  }
  */

/*==========================================================================*/
/* ii) Filtering stochastic orbital */

  #pragma omp parallel for private(iState,jState,iCoeff,index1,index2,sum)
  for(iState=0;iState<numStateUpProc;iState++){
    index1 = iState*numCoeff;
    for(jState=0;jState<numStatePrintUp;jState++){
      index2 = jState*numCoeff;
      sum = 0.0;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
        sum += cre_up[index1+iCoeff]*moUpRePrint[index2+iCoeff]+
               cim_up[index1+iCoeff]*moUpImPrint[index2+iCoeff];
      }
      sum *= 2.0;
      sum += cre_up[index1+numCoeff]*moUpRePrint[index2+numCoeff];
      wfDot[iState*numStatePrintUp+jState] = sum;
    }
  }

  #pragma omp parallel for private(imu,iState,jState,iCoeff,x,y)
  for(imu=0;imu<numChemPot;imu++){
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      stoWfUpRe[imu][iCoeff] = 0.0;
      stoWfUpIm[imu][iCoeff] = 0.0;
    }
    for(iState=0;iState<numStateUpProc;iState++){
      for(jState=0;jState<numStatePrintUp;jState++){
        x = occupNumber[imu*numStatePrintUp+jState];
        y = wfDot[iState*numStatePrintUp+jState]*x;
        for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
          stoWfUpRe[imu][iState*numCoeff+iCoeff] += y*moUpRePrint[jState*numCoeff+iCoeff];
          stoWfUpIm[imu][iState*numCoeff+iCoeff] += y*moUpImPrint[jState*numCoeff+iCoeff];
        }
      }
    }
  }

  if(inverseFlag==1){
    #pragma omp parallel for private(iState,iCoeff)
    for(iState=0;iState<numStateUpProc;iState++){
      for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
        stoWfUpRe[numChemPot-1][iState*numCoeff+iCoeff] = cre_up[iState*numCoeff+iCoeff]
                                        -stoWfUpRe[numChemPot-1][iState*numCoeff+iCoeff];
        stoWfUpIm[numChemPot-1][iState*numCoeff+iCoeff] = cim_up[iState*numCoeff+iCoeff]
                                        -stoWfUpIm[numChemPot-1][iState*numCoeff+iCoeff];
      }
    }
  }

  //test
  /*
  for(i=0;i<numProcStates;i++){
    if(myidState==i){
      for(iState=0;iState<numStateUpProc;iState++){
        for(imu=0;imu<numChemPot;imu++){
          printf("iState %i imu %i stowf %lg %lg\n",
                   iState,imu,stoWfUpRe[imu][iState*numCoeff+1],stoWfUpIm[imu][iState*numCoeff+1]);
        } 
      }
    }
  }
  */

  if(smearOpt>0&&filterDiagFlag==0){
    #pragma omp parallel for private(iCoeff)
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      entropyUpRe[iCoeff] = 0.0;
      entropyUpIm[iCoeff] = 0.0;
    }
    #pragma omp parallel for private(iState,jState,iCoeff,x,y)
    for(iState=0;iState<numStateUpProc;iState++){
      for(jState=0;jState<numStatePrintUp;jState++){
        x = entropyState[jState];
        y = wfDot[iState*numStatePrintUp+jState]*x;
        for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
          entropyUpRe[iState*numCoeff+iCoeff] += y*moUpRePrint[jState*numCoeff+iCoeff];
          entropyUpIm[iState*numCoeff+iCoeff] += y*moUpImPrint[jState*numCoeff+iCoeff];
        }//endfor iCoeff
      }//endfor jState
    }//endfor iState
    //printf("entropyUpRe %lg\n",entropyUpRe[1]);
  }//endif smearOpt

/*==========================================================================*/
/* iii) Calculate Chebyshev Momentum if necessary */

  if(storeChebyMomentsFlag==1){
    stodftCoefPos->chebyMomentsUp = (double*)cmalloc((polynormLength+1)*sizeof(double));
    //stodftCoefPos->chebyMomentsUp = (double*)crealloc(stodftCoefPos->chebyMomentsUp,
    //                                               (polynormLength+1)*sizeof(double));
    calcChebyMomentsFake(cp,class,general_data,ip_now);
    printf("Finish fake moment calculation\n");
  }

  free(occupNumber);
  free(energyScale);
  free(wfDot);
  free(entropyState);
 
  Barrier(comm_states);

/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/


#endif

