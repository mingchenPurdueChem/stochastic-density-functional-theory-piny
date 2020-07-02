/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: coeff.c                                        */
/*                                                                          */
/* This routine costruct a_n for f(z)=\sum a_nP_n(z)                        */
/* polynomial.                                                              */
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
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_stodft_local.h"

#define TIME_CP_OFF
#define SQRTFERMI

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genNewtonHermit(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This routine first guess a polynormial chain by an emperical formular */
/* The chain length is then increased to satisfied the fitting error for */
/* ALL different chemical potentials. I'd like to put the chemPot and	 */
/* numChemPot in the argument since we also need this when we use cheby  */
/* to calculate chemical potential. In this case we have numChemPot=1	 */
/* but I want to estimate polynormLength for two chem pots.		 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  NEWTONINFO *newtonInfo  = stodftInfo->newtonInfo;

  double fitErrTol = stodftInfo->fitErrTol;
  double fitErr;
  double beta = stodftInfo->beta;
  double energyDiff = stodftInfo->energyDiff;
  double *sampPoint;
  double *expanCoeff;
  double *sampPointUnscale;
  double *chemPot = stodftCoefPos->chemPot;

  int polynormLength = (int)(2.0*beta*energyDiff); //initial chain length, try drop the mutiplier 4.0 2.0->1.0
  int numChemPot     = stodftInfo->numChemPot;
  int smearOpt       = stodftInfo->smearOpt;
  int totalPoly      = polynormLength*numChemPot;  //iniital total polynormial

  FERMIFUNR fermiFunction = stodftInfo->fermiFunctionReal;
  
/*==========================================================================*/
/* I) Generate coeffcients for initial chain length  */
 
  printf("==============================================\n");
  printf("Start Calculating Polynormial Chain Length:\n");
  printf("totalPoly %i polynormLength %i numChemPot %i\n",totalPoly,polynormLength,numChemPot);
  printf("beta %lg energyDiff %lg\n",beta,energyDiff);
  stodftInfo->polynormLength = polynormLength;
  stodftCoefPos->expanCoeff = (double *)crealloc(stodftCoefPos->expanCoeff,
						totalPoly*sizeof(double));
  newtonInfo->sampPoint = (double *)crealloc(newtonInfo->sampPoint,
					     polynormLength*sizeof(double));
  newtonInfo->sampPointUnscale = (double *)crealloc(newtonInfo->sampPointUnscale,
						    polynormLength*sizeof(double));
  genSampNewtonHermit(stodftInfo,stodftCoefPos);  
  
  genCoeffNewtonHermit(stodftInfo,stodftCoefPos);
  
  fitErr = calcFitErrorNewton(stodftInfo,stodftCoefPos);
  
  printf("-----------------------------------------------\n");
  printf("Polynormial Length %i\n",polynormLength);
  printf("Fitting Error %lg Tolerance %lg\n",fitErr,fitErrTol);
  fflush(stdout);

/*==========================================================================*/
/* II) Increase chain length if polynomial expansion doesn't reach error    */
/*     torlerance */

  while(fitErr>fitErrTol){
    //free(stodftCoefPos->expanCoeff);
    //free(newtonInfo->sampPoint);
    //free(newtonInfo->sampPointUnscale);
    polynormLength += 1000;
    stodftInfo->polynormLength = polynormLength;
    totalPoly = polynormLength*numChemPot;
    stodftCoefPos->expanCoeff = (double *)crealloc(stodftCoefPos->expanCoeff,
						    totalPoly*sizeof(double));
    newtonInfo->sampPoint = (double *)crealloc(newtonInfo->sampPoint,
						polynormLength*sizeof(double));
    newtonInfo->sampPointUnscale = (double *)crealloc(newtonInfo->sampPointUnscale,
						    polynormLength*sizeof(double));
    
    genSampNewtonHermit(stodftInfo,stodftCoefPos);
    genCoeffNewtonHermit(stodftInfo,stodftCoefPos);
    fitErr = calcFitErrorNewton(stodftInfo,stodftCoefPos);  
    printf("-----------------------------------------------\n");
    printf("Polynormial Length %i\n",polynormLength);
    printf("Fitting Error %lg Tolerance %lg\n",fitErr,fitErrTol);
    fflush(stdout);
  }

  printf("Finish Calculating Polynormial Chain Length.\n");
  printf("The final Polynormial Chain Length is %i.\n",polynormLength);
  printf("==============================================\n");
  fflush(stdout);
  //exit(0);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genNewtonHermitTrueChemPot(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate expanCoeff for different chemical potentials                */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  NEWTONINFO *newtonInfo  = stodftInfo->newtonInfo;

  double fitErrTol = stodftInfo->fitErrTol;
  double fitErr;
  double beta = stodftInfo->beta;
  double energyDiff = stodftInfo->energyDiff;
  double *sampPoint;
  double *expanCoeff;
  double *sampPointUnscale;
  double *chemPot = stodftCoefPos->chemPot;

  int numChemPot = stodftInfo->numChemPot;
  int polynormLength = stodftInfo->polynormLength; //initial chain length, try drop the mutiplier 4.0
  int totalPoly      = polynormLength*numChemPot;  //iniital total polynormial
  int smearOpt       = stodftInfo->smearOpt;

  FERMIFUNR fermiFunction = stodftInfo->fermiFunctionReal;
/*==========================================================================*/
/* I) Generate coeffcients for initial chain length  */
  
  printf("chemPot %.16lg\n",chemPot[0]);

  printf("==============================================\n");
  printf("Start Calculating Polynormial Coeffcients w.r.t true Chem Pot:\n");
  //free(stodftCoefPos->expanCoeff);
  //free(newtonInfo->sampPoint);
  //free(newtonInfo->sampPointUnscale);

  stodftCoefPos->expanCoeff = (double *)crealloc(stodftCoefPos->expanCoeff,
                                                totalPoly*sizeof(double));
  newtonInfo->sampPoint = (double *)crealloc(newtonInfo->sampPoint,
                                             polynormLength*sizeof(double));
  newtonInfo->sampPointUnscale = (double *)crealloc(newtonInfo->sampPointUnscale,
                                                    polynormLength*sizeof(double));
  genSampNewtonHermit(stodftInfo,stodftCoefPos);

  genCoeffNewtonHermit(stodftInfo,stodftCoefPos);
  
  fitErr = calcFitErrorNewton(stodftInfo,stodftCoefPos);
  printf("fit error %.16lg\n",fitErr);

  if(smearOpt>0){
    stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));
    genCoeffNewtonHermitEntropy(stodftInfo,stodftCoefPos);
  }

  //exit(0);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genCoeffNewtonHermit(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate expanCoeff for different chemical potentials		 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  NEWTONINFO *newtonInfo  = stodftInfo->newtonInfo;
   

  int polynormLength     = stodftInfo->polynormLength;
  int numChemPot	 = stodftInfo->numChemPot;
  int energyWindowOn     = stodftInfo->energyWindowOn;
  int iPoly,jPoly,imu;
  int fragWindowFlag     = stodftInfo->fragWindowFlag;

  double Smin		   = newtonInfo->Smin;
  double Smax		   = newtonInfo->Smax;
  double scale		   = newtonInfo->scale;
  double chemPotTemp;
  double *sampPoint        = (double*)newtonInfo->sampPoint;
  double *expanCoeff       = (double*)stodftCoefPos->expanCoeff;
  double *sampPointUnscale = (double*)newtonInfo->sampPointUnscale;
  double *chemPot	   = stodftCoefPos->chemPot;

  double beta = stodftInfo->beta;
  double funValue,sum,prod;

  double timeStart,timeEnd;
  //FILE *fileCoeff = fopen("coeff","w");
  FERMIFUNR fermiFunction;
  FERMIFUNLR fermiFunctionLongDouble;

  if(energyWindowOn==0){
    fermiFunction = stodftInfo->fermiFunctionReal;
  }
  else{
    fermiFunctionLongDouble = stodftInfo->fermiFunctionLongDouble;
  }
 
  /*
  for(iPoly=0;iPoly<polynormLength;iPoly++){
    funValue = fermiFunction(sampPointUnscale[iPoly],chemPot[0],beta);
    printf("funValue %lg %lg\n",sampPointUnscale[iPoly],funValue);
  }

  fflush(stdout);
  exit(0);
  */
  cputime(&timeStart);  
 
  if(energyWindowOn==0){ //no energy window
    for(imu=0;imu<numChemPot;imu++){
      printf("chemPot %.16lg\n",chemPot[imu]);
      #ifdef SQRTFERMI
      funValue = sqrt(fermiFunction(sampPointUnscale[0],chemPot[imu],beta));
      #else
      funValue = fermiFunction(sampPointUnscale[0],chemPot[imu],beta);
      #endif

      expanCoeff[imu] = funValue;

      #ifdef SQRTFERMI
      funValue = sqrt(fermiFunction(sampPointUnscale[1],chemPot[imu],beta));
      #else
      funValue = fermiFunction(sampPointUnscale[1],chemPot[imu],beta);
      #endif

      expanCoeff[numChemPot+imu] = (funValue-expanCoeff[imu])/(sampPoint[1]-sampPoint[0]);
      for(iPoly=2;iPoly<polynormLength;iPoly++){

	#ifdef SQRTFERMI
	funValue = sqrt(fermiFunction(sampPointUnscale[iPoly],chemPot[imu],beta));
	#else
	funValue = fermiFunction(sampPointUnscale[iPoly],chemPot[imu],beta);
	#endif

	sum = funValue-expanCoeff[imu];
	prod = 1.0;
	for(jPoly=1;jPoly<iPoly;jPoly++){
	  prod *= (sampPoint[iPoly]-sampPoint[jPoly-1]);
	  sum -= expanCoeff[jPoly*numChemPot+imu]*prod;
	}//endfor jPoly
	prod *= (sampPoint[iPoly]-sampPoint[iPoly-1]);
	expanCoeff[iPoly*numChemPot+imu] = sum/prod;
      }//endfor iPoly
    }//endfor imu
  }
  else{
    // imu=0
    printf("chemPot %.16lg\n",chemPot[0]);
    funValue = sqrt(fermiFunctionLongDouble(sampPointUnscale[0],chemPot[0],beta));
    expanCoeff[0] = funValue;
    funValue = sqrt(fermiFunctionLongDouble(sampPointUnscale[1],chemPot[0],beta));
    expanCoeff[numChemPot] = (funValue-expanCoeff[0])/(sampPoint[1]-sampPoint[0]);
    for(iPoly=2;iPoly<polynormLength;iPoly++){

      funValue = sqrt(fermiFunctionLongDouble(sampPointUnscale[iPoly],chemPot[0],beta));

      sum = funValue-expanCoeff[0];
      prod = 1.0;
      for(jPoly=1;jPoly<iPoly;jPoly++){
	prod *= (sampPoint[iPoly]-sampPoint[jPoly-1]);
	sum -= expanCoeff[jPoly*numChemPot]*prod;
      }//endfor jPoly
      prod *= (sampPoint[iPoly]-sampPoint[iPoly-1]);
      expanCoeff[iPoly*numChemPot] = sum/prod;
    }//endfor iPoly


    // imu>=1
    if(fragWindowFlag==0){
      for(imu=1;imu<numChemPot;imu++){
        printf("chemPot %.16lg\n",chemPot[imu]);
        funValue = sqrt(fermiFunctionLongDouble(sampPointUnscale[0],chemPot[imu],beta)-
                        fermiFunctionLongDouble(sampPointUnscale[0],chemPot[imu-1],beta));

        expanCoeff[imu] = funValue;

        funValue = sqrt(fermiFunctionLongDouble(sampPointUnscale[1],chemPot[imu],beta)-
                        fermiFunctionLongDouble(sampPointUnscale[1],chemPot[imu-1],beta));

        expanCoeff[numChemPot+imu] = (funValue-expanCoeff[imu])/(sampPoint[1]-sampPoint[0]);
        for(iPoly=2;iPoly<polynormLength;iPoly++){

          funValue = sqrt(fermiFunctionLongDouble(sampPointUnscale[iPoly],chemPot[imu],beta)-
                          fermiFunctionLongDouble(sampPointUnscale[iPoly],chemPot[imu-1],beta));

          sum = funValue-expanCoeff[imu];
          prod = 1.0;
          for(jPoly=1;jPoly<iPoly;jPoly++){
            prod *= (sampPoint[iPoly]-sampPoint[jPoly-1]);
            sum -= expanCoeff[jPoly*numChemPot+imu]*prod;
          }//endfor jPoly
          prod *= (sampPoint[iPoly]-sampPoint[iPoly-1]);
          expanCoeff[iPoly*numChemPot+imu] = sum/prod;
        }//endfor iPoly
      }//endfor imu    
    }
    else{
      for(imu=1;imu<numChemPot-1;imu++){
        printf("imu %i chemPot %.16lg\n",imu,chemPot[imu]);
        funValue = sqrt(fermiFunctionLongDouble(sampPointUnscale[0],chemPot[imu],beta)-
                        fermiFunctionLongDouble(sampPointUnscale[0],chemPot[imu-1],beta));

        expanCoeff[imu] = funValue;

        funValue = sqrt(fermiFunctionLongDouble(sampPointUnscale[1],chemPot[imu],beta)-
                        fermiFunctionLongDouble(sampPointUnscale[1],chemPot[imu-1],beta));

        expanCoeff[numChemPot+imu] = (funValue-expanCoeff[imu])/(sampPoint[1]-sampPoint[0]);
        for(iPoly=2;iPoly<polynormLength;iPoly++){

          funValue = sqrt(fermiFunctionLongDouble(sampPointUnscale[iPoly],chemPot[imu],beta)-
                          fermiFunctionLongDouble(sampPointUnscale[iPoly],chemPot[imu-1],beta));

          sum = funValue-expanCoeff[imu];
          prod = 1.0;
          for(jPoly=1;jPoly<iPoly;jPoly++){
            prod *= (sampPoint[iPoly]-sampPoint[jPoly-1]);
            sum -= expanCoeff[jPoly*numChemPot+imu]*prod;
          }//endfor jPoly
          prod *= (sampPoint[iPoly]-sampPoint[iPoly-1]);
          expanCoeff[iPoly*numChemPot+imu] = sum/prod;
        }//endfor iPoly
      }//endfor imu
      imu = numChemPot-1;
      funValue = sqrt(1.0-fermiFunctionLongDouble(sampPointUnscale[0],chemPot[imu-1],beta));
      expanCoeff[imu] = funValue;
      funValue = sqrt(1.0-fermiFunctionLongDouble(sampPointUnscale[1],chemPot[imu-1],beta));
      expanCoeff[numChemPot+imu] = (funValue-expanCoeff[imu])/(sampPoint[1]-sampPoint[0]);
      for(iPoly=2;iPoly<polynormLength;iPoly++){
        funValue = sqrt(1.0-fermiFunctionLongDouble(sampPointUnscale[iPoly],chemPot[imu-1],beta));
        sum = funValue-expanCoeff[imu];
        prod = 1.0;
        for(jPoly=1;jPoly<iPoly;jPoly++){
          prod *= (sampPoint[iPoly]-sampPoint[jPoly-1]);
          sum -= expanCoeff[jPoly*numChemPot+imu]*prod;
        }//endfor jPoly
        prod *= (sampPoint[iPoly]-sampPoint[iPoly-1]);
        expanCoeff[iPoly*numChemPot+imu] = sum/prod;
      }//endfor iPoly       
    }//endif fragOpt
  }//endif 

  cputime(&timeEnd);

  printf("Coeff time %lg\n",timeEnd-timeStart);

  /*
  for(iPoly=1;iPoly<polynormLength;iPoly++){
    printf("coeff iPoly %lg\n",expanCoeff[iPoly*numChemPot]);
  }
  */

  //debug
  /*
  FILE *filecoeff = fopen("coeff-out","w");
  for(iPoly=0;iPoly<polynormLength;iPoly++){
    fprintf(filecoeff,"%i %.16lg\n",iPoly,expanCoeff[iPoly*numChemPot]);
  }
  fclose(filecoeff);
  */
  
  /*
  for(imu=0;imu<numChemPot;imu++){
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      fprintf(fileCoeff,"%.13lg\n",expanCoeff[iPoly*numChemPot+imu]);
    }
  }
  fclose(fileCoeff); 
  */
  //exit(0);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genSampNewtonHermit(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate expanCoeff for different chemical potentials                */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  NEWTONINFO *newtonInfo  = stodftInfo->newtonInfo;
  int polynormLength	     = stodftInfo->polynormLength;
  int numSampCand	     = 32*polynormLength;
  int iPoly,jPoly,iCand;
  int objMaxIndex;
  int numThreads = stodftInfo->numThreads;
  int iThread;
  int *objMaxIndexThreads = (int*)cmalloc(numThreads*sizeof(int));
  double Smin = newtonInfo->Smin;
  double Smax = newtonInfo->Smax;
  double scale = 1.0/newtonInfo->scale;
  double energyMin = stodftInfo->energyMin;
  double delta = (Smax-Smin)/(double)(numSampCand-1);
  double obj,objMax,diff;
  double *sampPoint = (double*)newtonInfo->sampPoint;
  double *sampPointUnscale = (double*)newtonInfo->sampPointUnscale;
  double *sampCand = (double*)cmalloc(numSampCand*sizeof(double));
  double *objValueArray = (double*)cmalloc(numSampCand*sizeof(double));
  double *objMaxThreads = (double*)cmalloc(numThreads*sizeof(double));

  double timeStart,timeEnd;
  double prod;

/*==========================================================================*/
/* 0) Generate sample candidates in range [Smin,Smax] */
  timeStart = omp_get_wtime();
 
  stodftInfo->polynormLength = polynormLength;
  for(iCand=0;iCand<numSampCand;iCand++)sampCand[iCand] = Smin+iCand*delta;

/*==========================================================================*/
/* 1) Select samples form sample candidates  */

  sampPoint[0] = sampCand[0];

  for(iCand=0;iCand<numSampCand;iCand++){
    diff = sampCand[iCand]-sampPoint[0];
    diff = diff*diff;
    if(diff<1.0e-20) objValueArray[iCand] = -1.0e30;
    else objValueArray[iCand] = log(diff);
  }

  omp_set_num_threads(numThreads);
  for(iPoly=1;iPoly<polynormLength;iPoly++){
    //if(iPoly%1000==0)printf("iPoly %i\n",iPoly);
    #pragma omp parallel private(iThread,iCand)
    {
      iThread = omp_get_thread_num();
      objMaxThreads[iThread] = -100000.0;
      #pragma omp for
      for(iCand=0;iCand<numSampCand;iCand++){
	if(objValueArray[iCand]>objMaxThreads[iThread]){
	  objMaxThreads[iThread] = objValueArray[iCand];
	  objMaxIndexThreads[iThread] = iCand;
	}
      }//endfor iCand
    }//omp
    objMax = -100000.0;
    for(iThread=0;iThread<numThreads;iThread++){
      if(objMaxThreads[iThread]>objMax){
	objMax = objMaxThreads[iThread];
	objMaxIndex = objMaxIndexThreads[iThread];
      }
    }
    sampPoint[iPoly] = sampCand[objMaxIndex];
    #pragma omp parallel for private(iCand,diff)
    for(iCand=0;iCand<numSampCand;iCand++){
      diff = sampCand[iCand]-sampPoint[iPoly];
      diff = diff*diff;
      if(diff<1.0e-20) objValueArray[iCand] += -1.0e30;
      else objValueArray[iCand] += log(diff);
    }//endfor iCand
  }//endfor iPoly

  for(iPoly=0;iPoly<polynormLength;iPoly++){
    sampPointUnscale[iPoly] = (sampPoint[iPoly]-Smin)*scale+energyMin;
  }

  //debug

  timeEnd = omp_get_wtime();
  printf("Sampling interpolation point time is %lgs.\n",timeEnd-timeStart);
  
  /*
  FILE *fileSampPoint = fopen("samp-point","w");
  for(iPoly=0;iPoly<polynormLength;iPoly++){
    fprintf(fileSampPoint,"%.16lg %.16lg\n",sampPoint[iPoly],sampPointUnscale[iPoly]);
  }
  fclose(fileSampPoint);
  */
  
  free(&sampCand[0]);
  free(&objValueArray[0]);
  free(&objMaxIndexThreads[0]);
  free(&objMaxThreads[0]);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genEigenOrb(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                 CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CELL *cell                    = &(general_data->cell);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  EWALD *ewald                  = &(general_data->ewald);
  EWD_SCR *ewd_scr              = &(class->ewd_scr);
  ATOMMAPS *atommaps            = &(class->atommaps);
  FOR_SCR *for_scr              = &(class->for_scr);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  PTENS *ptens                  = &(general_data->ptens);
  SIMOPTS *simopts              = &(general_data->simopts);

  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPEWALD *cpewald              = &(cp->cpewald);
  CPSCR *cpscr                  = &(cp->cpscr);
  PSEUDO *pseudo                = &(cp->pseudo);
  COMMUNICATE *communicate      = &(cp->communicate);

  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm;
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_sm            = &(cp->cp_para_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_dens_cp_box   = &(cp->cp_sclr_fft_pkg3d_dens_cp_box);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box   = &(cp->cp_para_fft_pkg3d_dens_cp_box);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_lg;
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg;
  CP_COMM_STATE_PKG *cp_comm_state_pkg_up          = &(cp->cp_comm_state_pkg_up);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_dn          = &(cp->cp_comm_state_pkg_dn);

  int numStateUpProc  = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc  = cpcoeffs_info->nstate_dn_proc;
  int numCoeff        = cpcoeffs_info->ncoef;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int cpLsda          = cpopts->cp_lsda;
  int realSparseOpt   = cpopts->realSparseOpt;
  int numCoeffUpTot   = numStateUpProc*numCoeff;
  int numCoeffDnTot   = numStateUpProc*numCoeff;
  int iCoeff;
  int cpWaveMin     = simopts->cp_wave_min;
  int cpMin         = simopts->cp_min;
  int cpWaveMinPimd = simopts->cp_wave_min_pimd;
  int cpMinOn = cpWaveMin + cpMin + cpWaveMinPimd;


  double *kseig_vals = (double*)cmalloc(numStateUpProc*sizeof(double))-1;
  double *kseig_vecs = (double*)cmalloc(numStateUpProc*numStateUpProc*sizeof(double))-1;
  double *ksmat_test = (double*)cmalloc(numStateUpProc*numStateUpProc*sizeof(double))-1;
  double *ks_scr = (double*)cmalloc(numStateUpProc*numStateUpProc*sizeof(double))-1;
  double *rs_scr1 = (double*)cmalloc(numStateUpProc*numStateUpProc*sizeof(double))-1;
  double *rs_scr2 = (double*)cmalloc(numStateUpProc*numStateUpProc*sizeof(double))-1;
  int *icoef_orth_up    = &(cpcoeffs_pos->icoef_orth_up);
  int *icoef_form_up    = &(cpcoeffs_pos->icoef_form_up);
  int *ifcoef_orth_up   = &(cpcoeffs_pos->ifcoef_orth_up);
  int *ifcoef_form_up   = &(cpcoeffs_pos->ifcoef_form_up);
  int *ioff_upt      = cpcoeffs_info->ioff_upt;
  double kseig_sum;

  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *cre_temp = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  double *cim_temp = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;

  //if(realSparseOpt==0){
  cp_sclr_fft_pkg3d_sm = &(cp->cp_sclr_fft_pkg3d_sm);
  cp_sclr_fft_pkg3d_lg = &(cp->cp_sclr_fft_pkg3d_lg);
  cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  //}
  /*
  else{
    cp_sclr_fft_pkg3d_sm = &(cp->cp_sclr_fft_pkg3d_sparse);
    cp_sclr_fft_pkg3d_lg = &(cp->cp_sclr_fft_pkg3d_sparse);
    cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_sparse);
  }
  */

  calcCoefForceWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

  cp_condiag_ksmat(cre_up,cim_up,*icoef_form_up,*icoef_orth_up,fcre_up,fcim_up,
                 *ifcoef_form_up,*ifcoef_orth_up,kseig_vals,kseig_vecs,
                  ksmat_test,ks_scr,rs_scr1,rs_scr2,ioff_upt,
                  &(cp->cp_comm_state_pkg_up),&kseig_sum);

  cp_rotate_vector(cre_up,cim_up,*icoef_form_up,
                      kseig_vecs,ioff_upt,cre_temp,cim_temp,
                      &(cp->cp_comm_state_pkg_up));

  free(&kseig_vals[1]);
  free(&kseig_vecs[1]);
  free(&ksmat_test[1]);
  free(&ks_scr[1]);
  free(&rs_scr1[1]);
  free(&rs_scr2[1]);
  free(&cre_temp[1]);
  free(&cim_temp[1]);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double calcFitErrorNewton(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
   
  NEWTONINFO *newtonInfo  = stodftInfo->newtonInfo;

  int numPointTest = 1000;
  int numChemPot = stodftInfo->numChemPot;
  int polynormLength = stodftInfo->polynormLength;
  int energyWindowOn     = stodftInfo->energyWindowOn;
  int fragWindowFlag     = stodftInfo->fragWindowFlag;
  int iPoint,iChem,iPoly;
  double pointTest;
  double beta       = stodftInfo->beta;
  double energyMin  = stodftInfo->energyMin;
  double energyMax  = stodftInfo->energyMax;
  double energyMean = stodftInfo->energyMean;
  double deltPoint  = (energyMax-energyMin)/numPointTest;
  double scale	    = newtonInfo->scale;
  double pointScale;
  double funValue,prod,funBM;
  double fitErr = -1.0e30;
  double diff;

  double *chemPot          = stodftCoefPos->chemPot;
  double *sampPoint = newtonInfo->sampPoint;
  double *expanCoeff = stodftCoefPos->expanCoeff;
  double *sampPointUnscale = newtonInfo->sampPointUnscale;

  //FILE *fileTestFun = fopen("fermi-fun","w");

  FERMIFUNR fermiFunction;
  FERMIFUNLR fermiFunctionLongDouble;

  if(energyWindowOn==0){
    fermiFunction = stodftInfo->fermiFunctionReal;
  }
  else{
    fermiFunctionLongDouble = stodftInfo->fermiFunctionLongDouble;
  }
 
  /*
  for(iPoly=1;iPoly<polynormLength;iPoly++){
    printf("coeff iPoly %lg\n",expanCoeff[iPoly*numChemPot]);
  }
  */

  if(energyWindowOn==0){
    fitErr = -1.0e30;
    for(iChem=0;iChem<numChemPot;iChem++){
      for(iPoint=0;iPoint<numPointTest;iPoint++){
	pointTest = energyMin+(iPoint+0.5)*deltPoint;
	pointScale = (pointTest-energyMean)*scale;
	funValue = expanCoeff[iChem];
	prod = 1.0;
	for(iPoly=1;iPoly<polynormLength;iPoly++){
	  prod *= pointScale-sampPoint[iPoly-1];
	  funValue += expanCoeff[iPoly*numChemPot+iChem]*prod;
	}
	#ifdef SQRTFERMI
	funBM = sqrt(fermiFunction(pointTest,chemPot[iChem],beta));
	#else
	funBM = fermiFunction(pointTest,chemPot[iChem],beta);
	#endif
	//funBM = fermiFunction(pointTest,chemPot[iChem],beta);
	//funBM = sqrt(fermiFunction(pointTest,chemPot[iChem],beta));
	diff = fabs(funValue-funBM);
	if(diff>fitErr)fitErr = diff;
	//printf("TestFunExpan %lg %lg %lg %lg\n",pointTest,pointScale,funValue,funBM);
	//fprintf(fileTestFun,"%lg %lg %lg\n",pointTest,funValue,funBM);
      }
    }
  }
  else{ //energy window
    fitErr = -1.0e30;
    for(iPoint=0;iPoint<numPointTest;iPoint++){
      pointTest = energyMin+(iPoint+0.5)*deltPoint;
      pointScale = (pointTest-energyMean)*scale;
      funValue = expanCoeff[0];
      prod = 1.0;
      for(iPoly=1;iPoly<polynormLength;iPoly++){
	prod *= pointScale-sampPoint[iPoly-1];
	funValue += expanCoeff[iPoly*numChemPot]*prod;
      }
      funBM = sqrt(fermiFunctionLongDouble(pointTest,chemPot[0],beta));
      //funBM = fermiFunction(pointTest,chemPot[iChem],beta);
      //funBM = sqrt(fermiFunction(pointTest,chemPot[iChem],beta));
      diff = fabs(funValue-funBM);
      if(diff>fitErr)fitErr = diff;
      //printf("TestFunExpan %lg %lg %lg %lg\n",pointTest,pointScale,funValue,funBM);
      //fprintf(fileTestFun,"%lg %lg %lg\n",pointTest,funValue,funBM);
    }
    if(fragWindowFlag==0){
      fitErr = -1.0e30;
      for(iChem=1;iChem<numChemPot;iChem++){
        for(iPoint=0;iPoint<numPointTest;iPoint++){
          pointTest = energyMin+(iPoint+0.5)*deltPoint;
          pointScale = (pointTest-energyMean)*scale;
          funValue = expanCoeff[iChem];
          prod = 1.0;
          for(iPoly=1;iPoly<polynormLength;iPoly++){
            prod *= pointScale-sampPoint[iPoly-1];
            funValue += expanCoeff[iPoly*numChemPot+iChem]*prod;
          }
          funBM = sqrt(fermiFunctionLongDouble(pointTest,chemPot[iChem],beta)-
                       fermiFunctionLongDouble(pointTest,chemPot[iChem-1],beta));
          //funBM = fermiFunction(pointTest,chemPot[iChem],beta);
          //funBM = sqrt(fermiFunction(pointTest,chemPot[iChem],beta));
          diff = fabs(funValue-funBM);
          if(diff>fitErr)fitErr = diff;
          //printf("TestFunExpan %lg %lg %lg %lg\n",pointTest,pointScale,funValue,funBM);
          //fprintf(fileTestFun,"%lg %lg %lg\n",pointTest,funValue,funBM);
        }
      }
    }
    else{
      for(iChem=1;iChem<numChemPot-1;iChem++){
        for(iPoint=0;iPoint<numPointTest;iPoint++){
          pointTest = energyMin+(iPoint+0.5)*deltPoint;
          pointScale = (pointTest-energyMean)*scale;
          funValue = expanCoeff[iChem];
          prod = 1.0;
          for(iPoly=1;iPoly<polynormLength;iPoly++){
            prod *= pointScale-sampPoint[iPoly-1];
            funValue += expanCoeff[iPoly*numChemPot+iChem]*prod;
          }
          funBM = sqrt(fermiFunctionLongDouble(pointTest,chemPot[iChem],beta)-
                       fermiFunctionLongDouble(pointTest,chemPot[iChem-1],beta));
          //funBM = fermiFunction(pointTest,chemPot[iChem],beta);
          //funBM = sqrt(fermiFunction(pointTest,chemPot[iChem],beta));
          diff = fabs(funValue-funBM);
          if(diff>fitErr)fitErr = diff;
          //printf("TestFunExpan %lg %lg %lg %lg\n",pointTest,pointScale,funValue,funBM);
          //fprintf(fileTestFun,"%lg %lg %lg\n",pointTest,funValue,funBM);
        }//endfor iPoint
      }//endfor iChem
      for(iPoint=0;iPoint<numPointTest;iPoint++){
        pointTest = energyMin+(iPoint+0.5)*deltPoint;
        pointScale = (pointTest-energyMean)*scale;
        funValue = expanCoeff[iChem];
        prod = 1.0;
        for(iPoly=1;iPoly<polynormLength;iPoly++){
          prod *= pointScale-sampPoint[iPoly-1];
          funValue += expanCoeff[iPoly*numChemPot+iChem]*prod;
        }
        funBM = sqrt(1.0-fermiFunctionLongDouble(pointTest,chemPot[numChemPot-2],beta));
        //funBM = fermiFunction(pointTest,chemPot[iChem],beta);
        //funBM = sqrt(fermiFunction(pointTest,chemPot[iChem],beta));
        diff = fabs(funValue-funBM);
        if(diff>fitErr)fitErr = diff;
        //printf("TestFunExpan %lg %lg %lg %lg\n",pointTest,pointScale,funValue,funBM);
        //fprintf(fileTestFun,"%lg %lg %lg\n",pointTest,funValue,funBM);
      }//endfor iPoint
    }//endif fragOpt
  }//endif energyWindowOn

  //fclose(fileTestFun);
  return fitErr;

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genChebyHermit(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos,int filterFlag)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This routine is the Chebyshev polynormial version of genNewtonHermit. */
/* The flag filterFlag determine the filter type. It could be Fermi      */
/* function or */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;

  double fitErrTol = stodftInfo->fitErrTol;
  double fitErr;
  double beta = stodftInfo->beta;
  double energyDiff = stodftInfo->energyDiff;
  double *sampPoint;
  double *expanCoeff,*expanCoeffLocal;
  double *sampPointUnscale;
  double *chemPot = stodftCoefPos->chemPot;

  int iChem,iPoly;
  int polynormLength = (int)(2.0*beta*energyDiff); //initial chain length, try drop the mutiplier 4.0 2.0->1.0
  int numChemPot     = stodftInfo->numChemPot;
  int smearOpt       = stodftInfo->smearOpt;
  int totalPoly      = polynormLength*numChemPot;  //iniital total polynormial

/*==========================================================================*/
/* I) Generate coeffcients for initial chain length  */
  
  printf("==============================================\n");
  printf("Start Calculating Polynormial Chain Length:\n");
  printf("totalPoly %i polynormLength %i numChemPot %i\n",totalPoly,polynormLength,numChemPot);
  printf("beta %lg energyDiff %lg\n",beta,energyDiff);

  stodftInfo->polynormLength = polynormLength;
  stodftCoefPos->expanCoeff = (double *)crealloc(stodftCoefPos->expanCoeff,
                                                totalPoly*sizeof(double));

  printf("ssssssssssssssssssssss\n");
  calcChebyCoeffWrapper(stodftInfo,stodftCoefPos,filterFlag);
  
  fitErr = calcFitErrorCheby(stodftInfo,stodftCoefPos,filterFlag);

  printf("-----------------------------------------------\n");
  printf("Polynormial Length %i\n",polynormLength);
  printf("Fitting Error %lg Tolerance %lg\n",fitErr,fitErrTol);
  fflush(stdout);

  while(fitErr>fitErrTol){
    polynormLength += 1000;
    stodftInfo->polynormLength = polynormLength;
    totalPoly = polynormLength*numChemPot;
    stodftCoefPos->expanCoeff = (double *)crealloc(stodftCoefPos->expanCoeff,
                                                totalPoly*sizeof(double));
    calcChebyCoeffWrapper(stodftInfo,stodftCoefPos,filterFlag);
    fitErr = calcFitErrorCheby(stodftInfo,stodftCoefPos,filterFlag);
    printf("-----------------------------------------------\n");
    printf("Polynormial Length %i\n",polynormLength);
    printf("Fitting Error %lg Tolerance %lg\n",fitErr,fitErrTol);
    fflush(stdout);
  }

  expanCoeffLocal = (double*)cmalloc(totalPoly*sizeof(double));
  expanCoeff = (double*)stodftCoefPos->expanCoeff;
  memcpy(expanCoeffLocal,expanCoeff,totalPoly*sizeof(double));
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      expanCoeff[iPoly*numChemPot+iChem] = 
                      expanCoeffLocal[iChem*polynormLength+iPoly];
    }
  }
  free(expanCoeffLocal);

  printf("Finish Calculating Polynormial Chain Length.\n");
  printf("The final Polynormial Chain Length is %i.\n",polynormLength);
  printf("==============================================\n");
  fflush(stdout);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genChebyHermitTrueChemPot(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos,
                               int filterFlag)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This routine is the Chebyshev polynormial version of genNewtonHermit  */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;

  double fitErrTol = stodftInfo->fitErrTol;
  double fitErr;
  double *expanCoeff,*expanCoeffLocal;
  double *chemPot = stodftCoefPos->chemPot;

  int iChem,iPoly;
  int polynormLength = stodftInfo->polynormLength; 
  int numChemPot     = stodftInfo->numChemPot;
  int smearOpt       = stodftInfo->smearOpt;
  int totalPoly      = polynormLength*numChemPot;  //iniital total polynormial

/*==========================================================================*/
/* I) Generate coeffcients for initial chain length  */

  printf("==============================================\n");
  printf("Start Calculating Polynormial Coeffcients w.r.t true Chem Pot:\n");

  stodftCoefPos->expanCoeff = (double *)crealloc(stodftCoefPos->expanCoeff,
                                                totalPoly*sizeof(double));
  expanCoeffLocal = (double *)cmalloc(totalPoly*sizeof(double));

  calcChebyCoeffWrapper(stodftInfo,stodftCoefPos,filterFlag);
  
  fitErr = calcFitErrorCheby(stodftInfo,stodftCoefPos,filterFlag);

  // Transpose
  expanCoeff = (double*)stodftCoefPos->expanCoeff;
  memcpy(expanCoeffLocal,expanCoeff,totalPoly*sizeof(double));
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      expanCoeff[iPoly*numChemPot+iChem] = 
                                 expanCoeffLocal[iChem*polynormLength+iPoly];
      //printf("iPoly %i iChem %i coeff %lg\n",iPoly,iChem,expanCoeff[iPoly*numChemPot+iChem]);
    }
  }

  if(smearOpt>0){
    stodftCoefPos->entropyExpanCoeff = (double*)cmalloc(polynormLength*sizeof(double));
    calcChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[0],
                   &(stodftCoefPos->entropyExpanCoeff[0]),13,stodftInfo->chemPotTrue);
    //memcpy(expanCoeffLocal,stodftCoefPos->entropyExpanCoeff,totalPoly*sizeof(double));    
  }


  free(expanCoeffLocal);
  printf("-----------------------------------------------\n");
  printf("Fitting Error %lg Tolerance %lg\n",fitErr,fitErrTol);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcChebyCoeffWrapper(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos,
                           int filterFlag)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Wrapper routine for calculating Cheby coeff with different filterFlag.*/
/* filterFlag = 0: Fermi functions for all chemical potentials           */
/* filterFlag = 1: Energy window without fragmentation, no window for    */
/*                 unoccupied orbitals.                                  */
/* filterFlag = 2: Energy window with fragmentation, the first time      */
/*                 filtering without multiplying F                       */
/* filterFlag = 3: Energy window with fragmentation, the second time     */
/*                 filtering with multiplying F                          */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int polynormLength     = stodftInfo->polynormLength;
  int numChemPot         = stodftInfo->numChemPot;
  int energyWindowOn     = stodftInfo->energyWindowOn;
  int iPoly,jPoly,iChem;
  int fragWindowFlag     = stodftInfo->fragWindowFlag;
  int numChemPotTemp;
  int numFFTGridMutpl = 32;
  int numChebyGridInit = polynormLength*numFFTGridMutpl;
  int numChebyGrid;

  double *expanCoeff       = (double*)stodftCoefPos->expanCoeff;
  double *chemPot          = stodftCoefPos->chemPot;
  
  double chemPotTrue       = stodftInfo->chemPotTrue;
  double funValue,sum,prod;

  fftw_complex *chebyCoeffsFFT,*funValGridFFT;

  numChebyGrid = roundFFT(numChebyGridInit);
  printf("numChebyGridInit %i numChebyGrid %i\n",numChebyGridInit,numChebyGrid);
  

  stodftInfo->numChebyGrid = numChebyGrid;
  printf("ttt111\n");
  stodftCoefPos->chebyCoeffsFFT = fftw_malloc(numChebyGrid*sizeof(fftw_complex));
  printf("ttt222\n");
  stodftCoefPos->funValGridFFT = fftw_malloc(numChebyGrid*sizeof(fftw_complex));
  printf("ttt333\n");
  chebyCoeffsFFT = stodftCoefPos->chebyCoeffsFFT;
  printf("ttt444\n");
  funValGridFFT = stodftCoefPos->funValGridFFT;
  printf("ttt555\n");
  stodftInfo->fftwPlanForward = fftw_plan_dft_1d(numChebyGrid,funValGridFFT,chebyCoeffsFFT,
                                  FFTW_FORWARD,FFTW_MEASURE);
  printf("ttt666\n");

  for(iPoly=0;iPoly<polynormLength*numChemPot;iPoly++)expanCoeff[iPoly] = 0.0;

  //for(iChem=0;iChem<numChemPot;iChem++)printf("11111111 chemPot %lg\n",chemPot[iChem]);

  switch(filterFlag){
    case 0:
      for(iChem=0;iChem<numChemPot;iChem++){
        calcChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[iChem],
                       &expanCoeff[iChem*polynormLength],2,0.0);
      }
      break;
    case 1:
      calcChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[0],&expanCoeff[0],2,0.0);
      for(iChem=1;iChem<numChemPot;iChem++){
        calcChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[iChem-1],
                                   &expanCoeff[iChem*polynormLength],10,0.0);
      }
      break;
    case 2: // 1st time filter
      calcChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[0],&expanCoeff[0],2,0.0);
      for(iChem=1;iChem<numChemPot-1;iChem++){
        calcChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[iChem-1],
                                   &expanCoeff[iChem*polynormLength],10,0.0);
      }
      calcChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[numChemPot-2],
                     &expanCoeff[(numChemPot-1)*polynormLength],12,0.0);
      break;
    case 3: // 2nd time filter with true chemical potential
      calcChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[0],&expanCoeff[0],4,chemPotTrue);
      for(iChem=1;iChem<numChemPot-1;iChem++){
        calcChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[iChem-1],
                                   &expanCoeff[iChem*polynormLength],6,chemPotTrue);
      }
      calcChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[numChemPot-2],
                     &expanCoeff[(numChemPot-1)*polynormLength],8,chemPotTrue);
      break;
  }

  /*
  // Transpose 
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      expanCoeff[iPoly*numChemPot+iChem] = coeffTemp[iChem*polynormLength+iPoly];
    }
  }
  */
  
  fftw_destroy_plan(stodftInfo->fftwPlanForward);
  fftw_free(stodftCoefPos->chebyCoeffsFFT);
  fftw_free(stodftCoefPos->funValGridFFT);
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcChebyCoeff(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos,
		    double *chemPot,double *chebyCoeffs,int funFlag,
                    double chemPotTest)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* T*/
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;

  int iGrid,iCoeff,iChem;
  int polynormLength  = stodftInfo->polynormLength;
  int numChebyMoments = stodftInfo->numChebyMoments;
  int numChebyGrid    = stodftInfo->numChebyGrid;
  int energyWindowOn  = stodftInfo->energyWindowOn;
  
  double energyDiff = stodftInfo->energyDiff*0.5;
  double energyMean = stodftInfo->energyMean;
  double beta = stodftInfo->beta;
  double x;
  double pre = 2.0*M_PI;
  double temp;
  long double y1,y2;
  long double templd;

  fftw_plan fftwPlanForward = stodftInfo->fftwPlanForward;
  fftw_complex *chebyCoeffsFFT = stodftCoefPos->chebyCoeffsFFT;
  fftw_complex *funValGridFFT = stodftCoefPos->funValGridFFT;

  FERMIFUNR fermiFunction;
  FERMIFUNLR fermiFunctionLongDouble;

  fermiFunction = stodftInfo->fermiFunctionReal;
  if(energyWindowOn==1)fermiFunctionLongDouble = stodftInfo->fermiFunctionLongDouble;

  for(iGrid=0;iGrid<numChebyGrid;iGrid++){
    x = energyDiff*cos(pre*(double)iGrid/(double)numChebyGrid)+energyMean;
    // P_N = I-\sum_i P_i
    switch(funFlag){
      case 1: // F
        funValGridFFT[iGrid] = fermiFunction(x,chemPot[0],beta);
        break;
      case 2: // \sqrt{F}
        funValGridFFT[iGrid] = sqrt(fermiFunction(x,chemPot[0],beta));
        break;
      case 3: // F*P_0
        funValGridFFT[iGrid] = fermiFunction(x,chemPotTest,beta)*
                               fermiFunction(x,chemPot[0],beta);
        break;
      case 4: // \sqrt{F*P_0}
        temp = fermiFunction(x,chemPotTest,beta)*
               fermiFunction(x,chemPot[0],beta);
        funValGridFFT[iGrid] = sqrt(temp);
        break;          
      case 5: // F*P_i
        temp = (double)(fermiFunctionLongDouble(x,chemPot[1],beta)-
                        fermiFunctionLongDouble(x,chemPot[0],beta));
        funValGridFFT[iGrid] = temp*fermiFunction(x,chemPotTest,beta);
        break;
      case 6: // \sqrt{F*P_i}
        temp = (double)(fermiFunctionLongDouble(x,chemPot[1],beta)-
                        fermiFunctionLongDouble(x,chemPot[0],beta));
        //y1 = fermiFunctionLongDouble(x,chemPot[1],beta);
        //y2 = fermiFunctionLongDouble(x,chemPot[0],beta);
        //if(temp<0.0)printf("Negative filter function %lg x %lg beta %lg chemPot %lg %lg %lg %lg\n",temp,x,beta,chemPot[1],chemPot[0],(double)(y1),(double)(y2));
        funValGridFFT[iGrid] = sqrt(temp*fermiFunction(x,chemPotTest,beta));  
        break;
      case 7: // F*P_N
        temp = (double)(1.0-fermiFunctionLongDouble(x,chemPot[0],beta));
        funValGridFFT[iGrid] = temp*fermiFunction(x,chemPotTest,beta);
        break;
      case 8: // \sqrt{F*P_N}
        temp = (double)(1.0-fermiFunctionLongDouble(x,chemPot[0],beta));
        funValGridFFT[iGrid] = sqrt(temp*fermiFunction(x,chemPotTest,beta));
        break;
      case 9: // P_i
        funValGridFFT[iGrid] = (double)(fermiFunctionLongDouble(x,chemPot[1],beta)-
                                        fermiFunctionLongDouble(x,chemPot[0],beta));
        break;
      case 10: // \sqrt{P_i}
        temp = (double)(fermiFunctionLongDouble(x,chemPot[1],beta)-
                                        fermiFunctionLongDouble(x,chemPot[0],beta)); 
        funValGridFFT[iGrid] = sqrt(temp);
        break;
      case 11: // P_N
        funValGridFFT[iGrid] = (double)(1.0-fermiFunctionLongDouble(x,chemPot[0],beta));
        break;
      case 12: // \sqrt{P_N}
#ifdef FAST_FILTER
        // We can't use deterministic orbitals in the window for unoccupied orbitals.
        // Therefre we use |chi>-(I-\sqrt{P_N})|chi> 
        temp = (double)(1.0-fermiFunctionLongDouble(x,chemPot[0],beta));
        funValGridFFT[iGrid] = 1.0-sqrt(temp);
#else
        temp = (double)(1.0-fermiFunctionLongDouble(x,chemPot[0],beta));
        funValGridFFT[iGrid] = sqrt(temp);
#endif
        break;
      case 13: // Entropy
        funValGridFFT[iGrid] = sqrt(-entropyReal(x,chemPotTest,beta));
        break;
      default:
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Unsupported filter function Flag %i!\n",funFlag);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(0);
    }
    //if(iGrid%1000==0)printf("x %lg funVal %lg\n",x,funValGridFFT[iGrid]);
  }

  fftw_execute(fftwPlanForward);
  
  for(iCoeff=1;iCoeff<polynormLength;iCoeff++){
    //printf("real %lg\n",creal(chebyCoeffsFFT[iCoeff]));
    chebyCoeffs[iCoeff] = 2.0*creal(chebyCoeffsFFT[iCoeff])/(double)numChebyGrid;
  }
  chebyCoeffs[0] = creal(chebyCoeffsFFT[0])/(double)numChebyGrid;

  //testChebyCoeff(stodftInfo,stodftCoefPos,chemPot,chebyCoeffs);
  //exit(0);
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double calcFitErrorCheby(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos,
                       int filterFlag)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This routine estimates the fitting error of the Chebyshev polynormial */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;

  int numPointTest = 1000;
  int numChemPot = stodftInfo->numChemPot;
  int polynormLength = stodftInfo->polynormLength;
  int energyWindowOn     = stodftInfo->energyWindowOn;
  int fragWindowFlag     = stodftInfo->fragWindowFlag;
  int iChem,iPoly;
  
  double testTemp;
  double testMax = -100000.0;
  double chemPotTrue       = stodftInfo->chemPotTrue;
  double *chemPot          = stodftCoefPos->chemPot;

  double *expanCoeff       = (double*)stodftCoefPos->expanCoeff;


  switch(filterFlag){
    case 0:
      for(iChem=0;iChem<numChemPot;iChem++){
        testTemp = testChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[iChem],
                                &expanCoeff[iChem*polynormLength],2,0.0);
        if(testTemp>testMax)testMax = testTemp;
      }
      break;
    case 1:
      testTemp = testChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[0],
                                &expanCoeff[0],2,0.0);
      if(testTemp>testMax)testMax = testTemp;
      for(iChem=1;iChem<numChemPot;iChem++){
        testTemp = testChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[iChem-1],
                                  &expanCoeff[iChem*polynormLength],10,0.0);
        if(testTemp>testMax)testMax = testTemp;
      }
      break;
    case 2:
      testTemp = testChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[0],
                                &expanCoeff[0],2,0.0);
      if(testTemp>testMax)testMax = testTemp;
      for(iChem=1;iChem<numChemPot-1;iChem++){
        testTemp = testChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[iChem-1],
                                  &expanCoeff[iChem*polynormLength],10,0.0); 
        if(testTemp>testMax)testMax = testTemp;
      }
      testTemp = testChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[numChemPot-2],
                                &expanCoeff[(numChemPot-1)*polynormLength],12,0.0);
      if(testTemp>testMax)testMax = testTemp;
      break;
    case 3:
      testTemp = testChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[0],
                                &expanCoeff[0],4,chemPotTrue);
      if(testTemp>testMax)testMax = testTemp;
      for(iChem=1;iChem<numChemPot-1;iChem++){
        testTemp = testChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[iChem-1],
                                  &expanCoeff[iChem*polynormLength],6,chemPotTrue); 
        if(testTemp>testMax)testMax = testTemp;
      }
      testTemp = testChebyCoeff(stodftInfo,stodftCoefPos,&chemPot[numChemPot-2],
                                &expanCoeff[(numChemPot-1)*polynormLength],8,chemPotTrue);
      if(testTemp>testMax)testMax = testTemp;
      break;
    default:
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Unsupported Filter Flag!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);    
  }
   

  return testMax;

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double testChebyCoeff(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos,
                     double *chemPot,double *chebyCoeffs,int funFlag,
                     double chemPotTest)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  CHEBYSHEVINFO *chebyshevInfo = stodftInfo->chebyshevInfo;
  int numPointTest = 1000;
  int polynormLength = stodftInfo->polynormLength;
  int iGrid,iCoeff;
  int energyWindowOn     = stodftInfo->energyWindowOn;
  
  double beta       = stodftInfo->beta;
  double energyMin  = stodftInfo->energyMin;
  double energyMax  = stodftInfo->energyMax;
  double energyMean = stodftInfo->energyMean;
  double deltPoint  = (energyMax-energyMin)/numPointTest;
  double deltaE = 2.0/(energyMax-energyMin);
  double point,pointScale;
  double diff,funTrue;
  double diffMax = -100000.0;
  double temp;

  FERMIFUNR fermiFunction = stodftInfo->fermiFunctionReal;
  FERMIFUNLR fermiFunctionLongDouble;
  double *x = (double*)cmalloc(numPointTest*sizeof(double));
  double *xScale = (double*)cmalloc(numPointTest*sizeof(double));
  double *T1 = (double*)cmalloc(numPointTest*sizeof(double));
  double *T2 = (double*)cmalloc(numPointTest*sizeof(double));
  double *T3 = (double*)cmalloc(numPointTest*sizeof(double));
  double *funVal = (double*)cmalloc(numPointTest*sizeof(double));

  if(energyWindowOn==1){
    fermiFunctionLongDouble = stodftInfo->fermiFunctionLongDouble;
  }

  for(iGrid=0;iGrid<numPointTest;iGrid++){
    x[iGrid] = energyMin+(iGrid+0.5)*deltPoint;
    xScale[iGrid] = (x[iGrid]-energyMean)*deltaE;
  }
  for(iGrid=0;iGrid<numPointTest;iGrid++){
    funVal[iGrid] = chebyCoeffs[0]+chebyCoeffs[1]*xScale[iGrid];
    T1[iGrid] = 1.0;
    T2[iGrid] = xScale[iGrid];
  }
  for(iCoeff=2;iCoeff<polynormLength;iCoeff++){
    for(iGrid=0;iGrid<numPointTest;iGrid++){
      T3[iGrid] = 2.0*xScale[iGrid]*T2[iGrid]-T1[iGrid];
      funVal[iGrid] += chebyCoeffs[iCoeff]*T3[iGrid];
      T1[iGrid] = T2[iGrid];
      T2[iGrid] = T3[iGrid];
    }
  }
  /*
  for(iCoeff=0;iCoeff<polynormLength;iCoeff++){
    printf("iCoeff %i coeff %lg\n",iCoeff,chebyCoeffs[iCoeff]);
  }
  */
  for(iGrid=0;iGrid<numPointTest;iGrid++){
    switch(funFlag){
      case 1: // F
        funTrue = fermiFunction(x[iGrid],chemPot[0],beta);
        break;
      case 2: // \sqrt{F}
        funTrue = sqrt(fermiFunction(x[iGrid],chemPot[0],beta));
        break;
      case 3: // F*P_0
        funTrue = fermiFunction(x[iGrid],chemPotTest,beta)*
                               fermiFunction(x[iGrid],chemPot[0],beta);
        break;
      case 4: // \sqrt{F*P_0}
        temp = fermiFunction(x[iGrid],chemPotTest,beta)*
               fermiFunction(x[iGrid],chemPot[0],beta);
        funTrue = sqrt(temp);
        break;          
      case 5: // F*P_i
        temp = (double)(fermiFunctionLongDouble(x[iGrid],chemPot[1],beta)-
                        fermiFunctionLongDouble(x[iGrid],chemPot[0],beta));
        funTrue = temp*fermiFunction(x[iGrid],chemPotTest,beta);
        break;
      case 6: // \sqrt{F*P_i}
        temp = (double)(fermiFunctionLongDouble(x[iGrid],chemPot[1],beta)-
                        fermiFunctionLongDouble(x[iGrid],chemPot[0],beta));
        funTrue = sqrt(temp*fermiFunction(x[iGrid],chemPotTest,beta));  
        break;
      case 7: // F*P_N
        temp = (double)(1.0-fermiFunctionLongDouble(x[iGrid],chemPot[0],beta));
        funTrue = temp*fermiFunction(x[iGrid],chemPotTest,beta);
        break;
      case 8: // \sqrt{F*P_N}
        temp = (double)(1.0-fermiFunctionLongDouble(x[iGrid],chemPot[0],beta));
        funTrue = sqrt(temp*fermiFunction(x[iGrid],chemPotTest,beta));
        break;
      case 9: // P_i
        funTrue = (double)(fermiFunctionLongDouble(x[iGrid],chemPot[1],beta)-
                        fermiFunctionLongDouble(x[iGrid],chemPot[0],beta));
        break;
      case 10: // \sqrt{P_i}
        temp = (double)(fermiFunctionLongDouble(x[iGrid],chemPot[1],beta)-
                        fermiFunctionLongDouble(x[iGrid],chemPot[0],beta));
        funTrue = sqrt(temp);
        break;
      case 11: // P_N
        funTrue = (double)(1.0-fermiFunctionLongDouble(x[iGrid],chemPot[0],beta));
        break;
      case 12: // \sqrt{P_N}
#ifdef FAST_FILTER
        temp = (double)(1.0-fermiFunctionLongDouble(x[iGrid],chemPot[0],beta));
        funTrue = 1.0-sqrt(temp);
#else
        temp = (double)(1.0-fermiFunctionLongDouble(x[iGrid],chemPot[0],beta));
        funTrue = sqrt(temp);
#endif
        break;      
      default:
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Unsupported filter function Flag!\n");
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(0);
    }

    diff = fabs(funTrue-funVal[iGrid]);
    //printf("x %lg xScale %lg funTrue %lg funVal %lg\n",x[iGrid],xScale[iGrid],funTrue,funVal[iGrid]);
    if(diff>diffMax)diffMax = diff;
  }
  //printf("diff %.10lg\n",diffMax);  

  free(x);
  free(xScale);
  free(T1);
  free(T2);
  free(T3);
  free(funVal);

  return diffMax;

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genCoeffNewtonHermitEntropy(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate expensions for entropy calculation. ONLY WORKS WITH CHEBY	 */
/* chemical potential option                                             */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  NEWTONINFO *newtonInfo  = stodftInfo->newtonInfo;

  int polynormLength     = stodftInfo->polynormLength;
  int iPoly,jPoly,imu;

  double Smin		   = newtonInfo->Smin;
  double Smax		   = newtonInfo->Smax;
  double scale		   = newtonInfo->scale;
  double *sampPoint        = (double*)newtonInfo->sampPoint;
  double *sampPointUnscale = (double*)newtonInfo->sampPointUnscale;
  double *entropyExpanCoeff = stodftCoefPos->entropyExpanCoeff;

  double beta = stodftInfo->beta;
  double chemPotTrue = stodftInfo->chemPotTrue;
  double funValue,sum,prod;

  double timeStart,timeEnd;
  //FILE *fileCoeff = fopen("coeff","w");
 
  /*
  for(iPoly=0;iPoly<polynormLength;iPoly++){
    funValue = fermiFunction(sampPointUnscale[iPoly],chemPot[0],beta);
    printf("funValue %lg %lg\n",sampPointUnscale[iPoly],funValue);
  }

  fflush(stdout);
  exit(0);
  */
  funValue = sqrt(-entropyReal(sampPointUnscale[0],chemPotTrue,beta));

  entropyExpanCoeff[0] = funValue;

  funValue = sqrt(-entropyReal(sampPointUnscale[1],chemPotTrue,beta));

  entropyExpanCoeff[1] = (funValue-entropyExpanCoeff[0])/(sampPoint[1]-sampPoint[0]);
  for(iPoly=2;iPoly<polynormLength;iPoly++){

    funValue = sqrt(-entropyReal(sampPointUnscale[iPoly],chemPotTrue,beta));

    sum = funValue-entropyExpanCoeff[0];
    prod = 1.0;
    for(jPoly=1;jPoly<iPoly;jPoly++){
      prod *= (sampPoint[iPoly]-sampPoint[jPoly-1]);
      sum -= entropyExpanCoeff[jPoly]*prod;
    }//endfor jPoly
    prod *= (sampPoint[iPoly]-sampPoint[iPoly-1]);
    entropyExpanCoeff[iPoly] = sum/prod;
  }//endfor iPoly

  /*
  for(iPoly=1;iPoly<polynormLength;iPoly++){
    printf("coeff iPoly %lg\n",expanCoeff[iPoly*numChemPot]);
  }
  */

  //debug
  /*
  FILE *filecoeff = fopen("coeff-out","w");
  for(iPoly=0;iPoly<polynormLength;iPoly++){
    fprintf(filecoeff,"%i %.16lg\n",iPoly,expanCoeff[iPoly*numChemPot]);
  }
  fclose(filecoeff);
  */
  
  /*
  for(imu=0;imu<numChemPot;imu++){
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      fprintf(fileCoeff,"%.13lg\n",expanCoeff[iPoly*numChemPot+imu]);
    }
  }
  fclose(fileCoeff); 
  */
  //exit(0);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
int roundFFT(int numGridIn)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* With Chebyshev Polynormial we want the numChebyGrid to be 2^n*3^m     */
/* m, n are choosen to be closest to 32*polynormLength (numGridIn).      */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
/*==========================================================================*/
  int m,n;
  int nmax,mnow;
  int i,j,k;
  int numGridOut = 1;
  double ln2 = log(2.0);
  double ln3 = log(3.0);
  double lnn = log(numGridIn);
  double min;
  double x,y,z;
  double diff;

  nmax = (int)(lnn/ln2)+1;
  min = fabs(pow(2.0,nmax)-numGridIn);
  n = nmax;
  m = 0;

  for(i=0;i<nmax;i++){
    x = fmax((lnn-i*ln2)/ln3,0.0);
    mnow = (int)(x)+1;
    diff = fabs(pow(2.0,i)*pow(3.0,mnow)-numGridIn);
    //printf("difffffff %lg min %lg\n",diff,min);
    if(diff<min){
      m = mnow;
      n = i;
      min = diff;
    }
  }

  //printf("mmmmmmm %i nnnnnnn %i\n",m,n);
  for(i=0;i<n;i++)numGridOut *= 2;
  for(i=0;i<m;i++)numGridOut *= 3;

  return numGridOut;
  
}/*end Routine*/
/*==========================================================================*/

