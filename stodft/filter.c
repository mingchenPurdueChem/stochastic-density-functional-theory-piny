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
    printf("Process ID %i cpy-wf %.8lg %.8lg add force %.8lg %.8lg prepare-loc %.8lg zero %.8lg",myidState,stodftInfo->cputime_new[1],stodftInfo->cputime_new[3],stodftInfo->cputime_new[2],stodftInfo->cputime_new[4],stodftInfo->cputime_new[5],stodftInfo->cputime_new[6]);
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

