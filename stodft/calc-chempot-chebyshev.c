/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: calc-chempot-chebyshev.c                       */
/*                                                                          */
/* This routine calculate the correct chemical potential via chebyshev      */
/* momentum.                                                                */
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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcChemPotCheby(CP *cp,CLASS *class,GENERAL_DATA *general_data,
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
  int numChebyGrid = polynormLength*numFFTGridMutpl;
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
  Barrier(comm_states);


/*==========================================================================*/
/* II) Solve N(mu)=Ne */

  if(myidState==0){
    printf("Start Calculating Chemical Potential\n");
    for(iPoly=0;iPoly<=polynormLength;iPoly++){
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
    
    //printf("2222 %.16lg\n",calcNumElecCheby(cp,0.4990160113690848,chebyCoeffs));
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

  free(chebyCoeffs);
  free(stodftCoefPos->chebyMomentsUp);
  if(cpLsda==1&&numStateDnProc!=0)free(stodftCoefPos->chebyMomentsDn);
  if(numProcStates>1&&myidState==0)free(chebyMomentsTemp);
  //exit(0);
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcChebyMoments(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                          int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Generate Chebyshev moments						 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CLATOMS_POS*  clatoms_pos     = &(class->clatoms_pos[ip_now]);
  COMMUNICATE *communicate      = &(cp->communicate);

  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;

  int expanType      = stodftInfo->expanType;
  int polynormLength = stodftInfo->polynormLength;
  int numChebyMoments = stodftInfo->numChebyMoments;
  int numChemPot     = stodftInfo->numChemPot;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int cpLsda         = cpopts->cp_lsda;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int myidState       = communicate->myid_state;
  int numThreads      = communicate->numThreads;
  int iState,iCoeff,iPoly,indexStart;
  int iOff1,iOff2;
  MPI_Comm comm_states   =    communicate->comm_states;


  double energyDiff  = stodftInfo->energyDiff;
  double energyMin   = stodftInfo->energyMin;
  double energyMax   = stodftInfo->energyMax;
  double energyMean  = stodftInfo->energyMean;
  double scale       = chebyshevInfo->scale;
  double polyCoeff;
  double timeProc,timeTot;
  double dot;
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;

  double *wfUpRe0,*wfUpIm0,*wfUpRe1,*wfUpIm1,*wfUpRe2,*wfUpIm2,*wfUpRe3,wfUpIm3;
  double *wfDnRe0,*wfDnIm0,*wfDnRe1,*wfDnIm1,*wfDnRe2,*wfDnIm2,*wfDnRe3,wfDnIm3;
  double *chebyMomentsUp = stodftCoefPos->chebyMomentsUp;
  double *chebyMomentsDn = stodftCoefPos->chebyMomentsDn;

  double *expanCoeff = (double*)stodftCoefPos->expanCoeff;
  double timeStart,timeEnd;

/*==========================================================================*/
/* 0) Allocate temp memories */

  stodftCoefPos->wfUpRe0 = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  stodftCoefPos->wfUpIm0 = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  stodftCoefPos->wfUpRe1 = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  stodftCoefPos->wfUpIm1 = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  wfUpRe0 = stodftCoefPos->wfUpRe0;
  wfUpIm0 = stodftCoefPos->wfUpIm0;
  wfUpRe1 = stodftCoefPos->wfUpRe1;
  wfUpIm1 = stodftCoefPos->wfUpIm1;
  if(cpLsda==1&&numStateDnProc!=0){
    stodftCoefPos->wfDnRe0 = (double*)cmalloc(numCoeffDnTotal*sizeof(double))-1;
    stodftCoefPos->wfDnIm0 = (double*)cmalloc(numCoeffDnTotal*sizeof(double))-1;
    stodftCoefPos->wfDnRe1 = (double*)cmalloc(numCoeffDnTotal*sizeof(double))-1;
    stodftCoefPos->wfDnIm1 = (double*)cmalloc(numCoeffDnTotal*sizeof(double))-1;
    wfDnRe0 = stodftCoefPos->wfDnRe0;
    wfDnIm0 = stodftCoefPos->wfDnIm0;
    wfDnRe1 = stodftCoefPos->wfDnRe1;
    wfDnIm1 = stodftCoefPos->wfDnIm1;
  }

  //double *numElecStoWf = stodftInfo->numElecStoWf;

/*==========================================================================*/
/* 1) Generate noise wave function */

  genNoiseOrbitalReal(cp,cpcoeffs_pos);

  //debug
  /*
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cre_up[iCoeff] = stodftCoefPos->creTest[iCoeff];
    cim_up[iCoeff] = stodftCoefPos->cimTest[iCoeff];
  }
  */

/*==========================================================================*/
/* 2) Copy the 0th order to wfUp(Dn)Re(Im)0 */

  /*
  for(iPoly=0;iPoly<=polynormLength;iPoly++){
    chebyMomentsUp[iPoly] = 0.0;
    if(cpLsda==1&&numStateDnProc!=0)chebyMomentsDn[iPoly] = 0.0;
  }
  */
  //double testsum = 0.0;
  //double testsum1;
  chebyMomentsUp[0] = 0.0;
  if(numThreads==1){
    memcpy(&wfUpRe0[1],&cre_up[1],numCoeffUpTotal*sizeof(double));
    memcpy(&wfUpIm0[1],&cim_up[1],numCoeffUpTotal*sizeof(double));
    memcpy(&wfUpRe1[1],&cre_up[1],numCoeffUpTotal*sizeof(double));
    memcpy(&wfUpIm1[1],&cim_up[1],numCoeffUpTotal*sizeof(double));
  }
  else{
    #pragma omp parallel for private(iCoeff)
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      wfUpRe0[iCoeff] = cre_up[iCoeff];
      wfUpIm0[iCoeff] = cim_up[iCoeff];
      wfUpRe1[iCoeff] = cre_up[iCoeff];
      wfUpIm1[iCoeff] = cim_up[iCoeff];
    }
  }
  for(iState=0;iState<numStateUpProc;iState++){
    //iOff1 = iState*polynormLength;
    iOff2 = iState*numCoeff;
    dot = 2.0*ddotBlasWrapperThreads(numCoeff,&wfUpRe0[iOff2+1],1,
	    &wfUpRe0[iOff2+1],1,numThreads)+
	    2.0*ddotBlasWrapperThreads(numCoeff,&wfUpIm0[iOff2+1],1,
	    &wfUpIm0[iOff2+1],1,numThreads)-
	    wfUpRe0[iOff2+numCoeff]*wfUpRe0[iOff2+numCoeff];   
    //numElecStoWf[iState] = dot;
    chebyMomentsUp[0] += dot;
    /*
    testsum1 = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      testsum1 += wfUpRe0[iOff2+iCoeff]*wfUpRe0[iOff2+iCoeff]+
		 wfUpIm0[iOff2+iCoeff]*wfUpIm0[iOff2+iCoeff];
    }
    testsum += 2.0*testsum1+wfUpRe0[iOff2+numCoeff]*wfUpRe0[iOff2+numCoeff];
    */
  }//endfor iState
  if(cpLsda==1&&numStateDnProc!=0){
    chebyMomentsDn[0] = 0.0;
    memcpy(&wfDnRe0[1],&cre_dn[1],numCoeffDnTotal*sizeof(double));
    memcpy(&wfDnIm0[1],&cim_dn[1],numCoeffDnTotal*sizeof(double));
    memcpy(&wfDnRe1[1],&cre_dn[1],numCoeffDnTotal*sizeof(double));
    memcpy(&wfDnIm1[1],&cim_dn[1],numCoeffDnTotal*sizeof(double));
    for(iState=0;iState<numStateDnProc;iState++){
      //iOff1 = iState*polynormLength;
      iOff2 = iState*numCoeff;
      dot = 2.0*ddotBlasWrapper(numCoeff,&wfDnRe0[iOff2+1],1,&wfDnRe0[iOff2+1],1)+
	      2.0*ddotBlasWrapper(numCoeff,&wfDnIm0[iOff2+1],1,&wfDnIm0[iOff2+1],1)-
	      wfDnRe0[iOff2+numCoeff]*wfDnRe0[iOff2+numCoeff];
      chebyMomentsDn[0] += dot;
    }//endfor iState
  }//endif

/*==========================================================================*/
/* 2) Calculate the 1st Order */
  
  normHCheby(cp,class,general_data,cpcoeffs_pos,clatoms_pos,1);
  chebyMomentsUp[1] = 0.0;
  for(iState=0;iState<numStateUpProc;iState++){
    //iOff1 = iState*polynormLength;
    iOff2 = iState*numCoeff;
    dot = 2.0*ddotBlasWrapperThreads(numCoeff,&wfUpRe0[iOff2+1],1,
	    &fcre_up[iOff2+1],1,numThreads)+
            2.0*ddotBlasWrapperThreads(numCoeff,&wfUpIm0[iOff2+1],1,
	    &fcim_up[iOff2+1],1,numThreads)-
            wfUpRe0[iOff2+numCoeff]*fcre_up[iOff2+numCoeff];
    //numElecStoWf[numStateUpProc+iState] = dot;
    chebyMomentsUp[1] += dot;
  }
  memcpy(&cre_up[1],&fcre_up[1],numCoeffUpTotal*sizeof(double));
  memcpy(&cim_up[1],&fcim_up[1],numCoeffUpTotal*sizeof(double));
  if(cpLsda==1&&numStateDnProc!=0){
    chebyMomentsDn[1] = 0.0;
    for(iState=0;iState<numStateDnProc;iState++){
      //iOff1 = iState*polynormLength;
      iOff2 = iState*numCoeff;
      dot = 2.0*ddotBlasWrapperThreads(numCoeff,&wfDnRe0[iOff2+1],1,
	      &fcre_dn[iOff2+1],1,numThreads)+
	      2.0*ddotBlasWrapperThreads(numCoeff,&wfDnIm0[iOff2+1],1,
	      &fcim_dn[iOff2+1],1,numThreads)-
	      wfDnRe0[iOff2+numCoeff]*fcre_dn[iOff2+numCoeff];
      chebyMomentsUp[1] += dot;
    }
    memcpy(&cre_dn[1],&fcre_dn[1],numCoeffDnTotal*sizeof(double));
    memcpy(&cim_dn[1],&fcim_dn[1],numCoeffDnTotal*sizeof(double));
  }

/*==========================================================================*/
/* 3) Calculate All the other Order */

  //printf("full %i half %i\n",polynormLength,2*numChebyMoments-2);
  //exit(0);
  
  for(iPoly=2;iPoly<numChebyMoments;iPoly++){
  //for(iPoly=2;iPoly<polynormLength;iPoly++){
    if(iPoly%1000==0&&myidState==0){
      printf("iPoly %i\n",iPoly);
    }
    chebyMomentsUp[iPoly] = 0.0;
    chebyMomentsUp[iPoly*2] = 0.0;
    chebyMomentsUp[iPoly*2-1] = 0.0;
    normHCheby(cp,class,general_data,cpcoeffs_pos,clatoms_pos,iPoly);
    for(iState=0;iState<numStateUpProc;iState++){
      //iOff1 = iState*polynormLength;
      iOff2 = iState*numCoeff;
      dot = 2.0*ddotBlasWrapperThreads(numCoeff,&wfUpRe0[iOff2+1],1,
	    &fcre_up[iOff2+1],1,numThreads)+
	    2.0*ddotBlasWrapperThreads(numCoeff,&wfUpIm0[iOff2+1],1,
	    &fcim_up[iOff2+1],1,numThreads)-
	    wfUpRe0[iOff2+numCoeff]*fcre_up[iOff2+numCoeff]; 
      //numElecStoWf[iPoly*numStateUpProc+iState] = dot;
      chebyMomentsUp[iPoly] += dot;
      /*
      chebyMomentsUp[iPoly*2-2] +=
              2.0*ddotBlasWrapper(numCoeff,&cre_up[iOff2+1],1,&cre_up[iOff2+1],1)+
              2.0*ddotBlasWrapper(numCoeff,&cim_up[iOff2+1],1,&cim_up[iOff2+1],1)-
              cre_up[iOff2+numCoeff]*cre_up[iOff2+numCoeff];         
      */
      dot = 2.0*ddotBlasWrapperThreads(numCoeff,&fcre_up[iOff2+1],1,
	    &cre_up[iOff2+1],1,numThreads)+
            2.0*ddotBlasWrapperThreads(numCoeff,&fcim_up[iOff2+1],1,
	    &cim_up[iOff2+1],1,numThreads)-
            fcre_up[iOff2+numCoeff]*cre_up[iOff2+numCoeff];
      //numElecStoWf[(2*iPoly-1)*numStateUpProc+iState] = dot;
      chebyMomentsUp[iPoly*2-1] += dot;
      dot = 2.0*ddotBlasWrapperThreads(numCoeff,&fcre_up[iOff2+1],1,
	    &fcre_up[iOff2+1],1,numThreads)+
            2.0*ddotBlasWrapperThreads(numCoeff,&fcim_up[iOff2+1],1,
	    &fcim_up[iOff2+1],1,numThreads)-
            fcre_up[iOff2+numCoeff]*fcre_up[iOff2+numCoeff];
      //numElecStoWf[(2*iPoly)*numStateUpProc+iState] = dot;
      chebyMomentsUp[iPoly*2] += dot;
      //numElecStoWf[(2*iPoly-1)*numStateUpProc+iState] = 
      //    2.0*numElecStoWf[(2*iPoly-1)*numStateUpProc+iState]-numElecStoWf[numStateUpProc+iState];
      //numElecStoWf[(2*iPoly)*numStateUpProc+iState] = 
      //    2.0*numElecStoWf[(2*iPoly)*numStateUpProc+iState]-numElecStoWf[iState];
    }//endfor iState
    chebyMomentsUp[iPoly*2] = 2.0*chebyMomentsUp[iPoly*2]-chebyMomentsUp[0];
    chebyMomentsUp[iPoly*2-1] = 2.0*chebyMomentsUp[iPoly*2-1]-chebyMomentsUp[1];
    /*
    if(myidState==0){
      printf("%i %lg %i %lg %i %lg\n",iPoly,chebyMomentsUp[iPoly],iPoly*2-1,chebyMomentsUp[iPoly*2-1],
  	    iPoly*2,chebyMomentsUp[iPoly*2]);
    }
    */
    if(numThreads==1){
      memcpy(&wfUpRe1[1],&cre_up[1],numCoeffUpTotal*sizeof(double));
      memcpy(&wfUpIm1[1],&cim_up[1],numCoeffUpTotal*sizeof(double));
      memcpy(&cre_up[1],&fcre_up[1],numCoeffUpTotal*sizeof(double));
      memcpy(&cim_up[1],&fcim_up[1],numCoeffUpTotal*sizeof(double));    
    }
    else{
      #pragma omp parallel for private(iCoeff)
      for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
	wfUpRe1[iCoeff] = cre_up[iCoeff];
	wfUpIm1[iCoeff] = cim_up[iCoeff];
	cre_up[iCoeff] = fcre_up[iCoeff];
	cim_up[iCoeff] = fcim_up[iCoeff];
      }
    }
    if(cpLsda==1&&numStateDnProc!=0){
      chebyMomentsDn[iPoly] = 0.0;
      chebyMomentsDn[iPoly*2] = 0.0;
      chebyMomentsDn[iPoly*2-1] = 0.0;
      for(iState=0;iState<numStateDnProc;iState++){
	//iOff1 = iState*polynormLength;
	iOff2 = iState*numCoeff;
	dot = 2.0*ddotBlasWrapper(numCoeff,&wfDnRe0[iOff2+1],1,&fcre_dn[iOff2+1],1)+
		2.0*ddotBlasWrapper(numCoeff,&wfDnIm0[iOff2+1],1,&fcim_dn[iOff2+1],1)-
		wfDnRe0[iOff2+numCoeff]*fcre_dn[iOff2+numCoeff];
	chebyMomentsDn[iPoly] += dot;
	/*
	chebyMomentsDn[iPoly*2-2] +=
		2.0*ddotBlasWrapper(numCoeff,&cre_dn[iOff2+1],1,&cre_dn[iOff2+1],1)+
		2.0*ddotBlasWrapper(numCoeff,&cim_dn[iOff2+1],1,&cim_dn[iOff2+1],1)-
		cre_dn[iOff2+numCoeff]*cre_dn[iOff2+numCoeff];
	*/
	dot = 2.0*ddotBlasWrapper(numCoeff,&fcre_dn[iOff2+1],1,&cre_dn[iOff2+1],1)+
		2.0*ddotBlasWrapper(numCoeff,&fcim_dn[iOff2+1],1,&cim_dn[iOff2+1],1)-
		fcre_dn[iOff2+numCoeff]*cre_dn[iOff2+numCoeff];
	chebyMomentsDn[iPoly*2-1] += dot;
        dot = 2.0*ddotBlasWrapper(numCoeff,&fcre_dn[iOff2+1],1,&fcre_dn[iOff2+1],1)+
                2.0*ddotBlasWrapper(numCoeff,&fcim_dn[iOff2+1],1,&fcim_dn[iOff2+1],1)-
                fcre_dn[iOff2+numCoeff]*fcre_dn[iOff2+numCoeff];
	chebyMomentsDn[iPoly*2] += dot;
      }//endfor iState
      chebyMomentsDn[iPoly*2] = 2.0*chebyMomentsDn[iPoly*2]-chebyMomentsDn[0];
      chebyMomentsDn[iPoly*2-1] = 2.0*chebyMomentsDn[iPoly*2-1]-chebyMomentsDn[1];
      memcpy(&wfDnRe1[1],&cre_dn[1],numCoeffDnTotal*sizeof(double));
      memcpy(&wfDnIm1[1],&cim_dn[1],numCoeffDnTotal*sizeof(double));
      memcpy(&cre_dn[1],&fcre_dn[1],numCoeffDnTotal*sizeof(double));
      memcpy(&cim_dn[1],&fcim_dn[1],numCoeffDnTotal*sizeof(double));
    }//endif
  }//endfor iPoly

  free(&wfUpRe0[1]);
  free(&wfUpIm0[1]);
  free(&wfUpRe1[1]);
  free(&wfUpIm1[1]);
  if(cpLsda==1&&numStateDnProc!=0){
    free(&wfDnRe0[1]);
    free(&wfDnIm0[1]);
    free(&wfDnRe1[1]);
    free(&wfDnIm1[1]);
  }


/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double calcNumElecCheby(CP *cp,double chemPot,double *chebyCoeffs)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
 
  int occNumber = stodftInfo->occNumber;
  int polynormLength = stodftInfo->polynormLength;
  int cpLsda         = cpopts->cp_lsda;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int iPoly,iState;
  int numThreads = cp->communicate.numThreads;
  
  double *chebyMomentsUp = stodftCoefPos->chebyMomentsUp;
  double *chebyMomentsDn = stodftCoefPos->chebyMomentsDn;
  double numElec = 0.0;
  
  calcChebyCoeff(stodftInfo,stodftCoefPos,chemPot,chebyCoeffs);
  // debug
  /*
  double numElecTest;
  double *numElecStoWf = stodftInfo->numElecStoWf;
  for(iState=0;iState<numStateUpProc;iState++){
    numElecTest = 0.0;
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      numElecTest += numElecStoWf[iPoly*numStateUpProc+iState]*chebyCoeffs[iPoly];
    }
    printf("iState %i numElec %.16lg\n",iState,numElecTest*occNumber);
  }
  */
  //printf("chebyCoeffs %lg %lg %lg\n",chebyCoeffs[0],chebyCoeffs[1],chebyCoeffs[2]);
  numElec += ddotBlasWrapperThreads(polynormLength,chebyCoeffs,1,chebyMomentsUp,1,numThreads);
  //debug
  /*
  for(iPoly=0;iPoly<polynormLength;iPoly++){
    numElec += chebyCoeffs[iPoly]*chebyMomentsUp[iPoly];
    printf("iPoly %i coeff %lg moments %lg numElec %lg\n",iPoly,chebyCoeffs[iPoly],chebyMomentsUp[iPoly],numElec);
  }
  */
  if(cpLsda==1&&numStateDnProc!=0){
    numElec += ddotBlasWrapperThreads(polynormLength,chebyCoeffs,1,chebyMomentsDn,1,numThreads);
  }
  
  //exit(0);
  return numElec*occNumber;

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

#ifdef FAST_FILTER   

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcChebyMomentsFake(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                          int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Generate Chebyshev moments                                            */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  CPCOEFFS_POS *cpcoeffs_pos = &(cp->cpcoeffs_pos[1]);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);


  int iSto,iState,iCoeff,iPoly;
  int index1,index2,index3;
  int polynormLength = stodftInfo->polynormLength;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numStatePrintUp = stodftInfo->numStatePrintUp;

  double Smin           = chebyshevInfo->Smin;
  double Smax           = chebyshevInfo->Smax;
  double scale          = chebyshevInfo->scale;
  double energyMean     = stodftInfo->energyMean;
  double energyDiff     = stodftInfo->energyDiff;
  double prefact        = -scale*energyMean;
  double scale1         = -scale*0.25;
  double scale2         = -scale*0.5;
  double x,y;
  double sum;

  double *stoDetDot;
  double *chebyEnergy;
  double *energyLevel = stodftCoefPos->energyLevel;
  double *chebyMomentsUp = stodftCoefPos->chebyMomentsUp;
  double *moUpRePrint = stodftCoefPos->moUpRePrint;
  double *moUpImPrint = stodftCoefPos->moUpImPrint;
  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;

/*--------------------------------------------------------------------------*/
/* i) dot product between deterministic MO and noise orbital */  

  stoDetDot = (double*)cmalloc(numStatePrintUp*sizeof(double));

  genNoiseOrbitalReal(cp,cpcoeffs_pos);

  for(iState=0;iState<numStatePrintUp;iState++){
    index2 = iState*numCoeff;
    stoDetDot[iState] = 0.0;
    for(iSto=0;iSto<numStateUpProc;iSto++){
      index1 = iSto*numCoeff;
      sum = 0.0;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
        sum += coeffReUp[index1+iCoeff]*moUpRePrint[index2+iCoeff]+
               coeffImUp[index1+iCoeff]*moUpImPrint[index2+iCoeff];
      }
      sum *= 2.0;
      sum += coeffReUp[index1+numCoeff]*moUpRePrint[index2+numCoeff];
      stoDetDot[iState] += sum*sum;
    }
  }

/*--------------------------------------------------------------------------*/
/* ii) Chebyshev poly for all energy levels */

  chebyEnergy = (double *)cmalloc(polynormLength*numStatePrintUp*sizeof(double));

  // Scale energy
  for(iState=0;iState<numStatePrintUp;iState++){
    chebyEnergy[iState] = 1.0;
    chebyEnergy[numStatePrintUp+iState] = (energyLevel[iState]-energyMean)*scale;   
  }

  for(iPoly=2;iPoly<polynormLength;iPoly++){
    index1 = iPoly*numStatePrintUp;
    index2 = (iPoly-1)*numStatePrintUp;
    index3 = (iPoly-2)*numStatePrintUp;
    for(iState=0;iState<numStatePrintUp;iState++){
      x = chebyEnergy[numStatePrintUp+iState];
      chebyEnergy[index1+iState] = 2*x*chebyEnergy[index2+iState]
                                   -chebyEnergy[index3+iState];
    }
  }

/*--------------------------------------------------------------------------*/
/* iii) Calculate Chebyshev Momentum */

  for(iPoly=0;iPoly<polynormLength;iPoly++){
    chebyMomentsUp[iPoly] = 0.0;
    for(iState=0;iState<numStatePrintUp;iState++){
      chebyMomentsUp[iPoly] += chebyEnergy[iPoly*numStatePrintUp+iState]*
                                 stoDetDot[iState];
    }
  }

  printf("11111111 %lg\n",chebyMomentsUp[0]);

  free(stoDetDot);
  free(chebyEnergy);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

#endif

