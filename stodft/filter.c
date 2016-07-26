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
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CLATOMS_POS*  clatoms_pos     = &(class->clatoms_pos[ip_now]);

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
  int imu,iCoeff,iPoly,indexStart;
  int startIndex;

  double energyDiff  = stodftInfo->energyDiff;
  double energyMin   = stodftInfo->energyMin;
  double energyMax   = stodftInfo->energyMax;
  double energyMean  = stodftInfo->energyMean;
  double scale       = newtonInfo->scale;
  double polyCoeff;
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
 
/*==========================================================================*/
/* 0) Copy the initial stochastic orbital */
 
  for(imu=0;imu<numChemPot;imu++){
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      stoWfUpRe[imu][iCoeff] = expanCoeff[imu]*cre_up[iCoeff];
      stoWfUpIm[imu][iCoeff] = expanCoeff[imu]*cim_up[iCoeff];
    }//endfor iCoeff
    if(cpLsda==1&&numStateDnProc!=0){
      for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
	stoWfDnRe[imu][iCoeff] = expanCoeff[imu]*cre_dn[iCoeff];
	stoWfDnIm[imu][iCoeff] = expanCoeff[imu]*cim_dn[iCoeff];
      }//endfor iCoeff      
    }//endif 
  }//endfor imu

/*==========================================================================*/
/* 1) Loop over all polynomial terms (iPoly=0<=>polynomial order 1) */

  cputime(&timeStart);  

  for(iPoly=1;iPoly<polynormLength;iPoly++){
    if(iPoly%1000==0){
      printf("%lg%% ",iPoly*100.0/polynormLength);
      fflush(stdout);
    }  
    normHNewtonHerm(cp,class,general_data,
                 cpcoeffs_pos,clatoms_pos,sampPoint[iPoly-1]);
    for(imu=0;imu<numChemPot;imu++){
      polyCoeff = expanCoeff[iPoly*numChemPot+imu];
      for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
	stoWfUpRe[imu][iCoeff] += polyCoeff*cre_up[iCoeff];	
        stoWfUpIm[imu][iCoeff] += polyCoeff*cim_up[iCoeff];                       

      }//endfor iCoeff
      if(cpLsda==1&&numStateDnProc!=0){
        for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
	  stoWfDnRe[imu][iCoeff] += polyCoeff*cre_dn[iCoeff];                     
	  stoWfDnIm[imu][iCoeff] += polyCoeff*cim_dn[iCoeff];
        }//endfor iCoeff        
      }//endif 
    }//endfor imu
  }//endfor iPoly
  cputime(&timeEnd);
  printf("\n");
  printf("Filter time is %lg\n",timeEnd-timeStart);

/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genEnergyMax(CP *cp,CLASS *class,GENERAL_DATA *general_data,
	  CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate expanCoeff for different chemical potentials                */
/* Call this function only for the first process.	     */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  COMMUNICATE *communicate      = &(cp->communicate);

  double randMin    = -1.0;
  double randMax    = 1.0;
  double length	    = 0.0;
  double energyConv     = 1.0e-6;
  double energy	    = 0.0;
  double energyOld;

  int numStateUpProc  = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc  = cpcoeffs_info->nstate_dn_proc;
  int numCoeff        = cpcoeffs_info->ncoef;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int cpLsda          = cpopts->cp_lsda;
  int numCoeffUpTot   = numStateUpProc*numCoeff;
  int numCoeffDnTot   = numStateDnProc*numCoeff;

  int numIteration    = 100;
  int iIter;
  int iState,iCoeff,iCoeffStart,index1,index2;

  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *coeffReUpBackup = (double*)cmalloc((numCoeff+1)*sizeof(double));
  double *coeffImUpBackup = (double*)cmalloc((numCoeff+1)*sizeof(double));
  double *randTrail = (double *)cmalloc((2*numCoeff+1)*sizeof(double));

/*==========================================================================*/
/* I) Set parameters and backup			        */

  printf("==============================================\n");
  printf("Estimate Energy Upperbound: \n");
  fflush(stdout);

  cpcoeffs_info->nstate_up_proc = 1;
  cpcoeffs_info->nstate_dn_proc = 1;

  /*
  for(iCoeff=1;iCoeff<=numCoeffUpTot;iCoeff++){
    fcre_up[iCoeff] = 0.0;
    fcim_up[iCoeff] = 0.0;
  }
  */

  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    coeffReUpBackup[iCoeff] = cre_up[iCoeff];
    coeffImUpBackup[iCoeff] = cim_up[iCoeff];
  }//endfor iCoeff
  /*
  length = 0.0;
  for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
    length += cre_up[iCoeff]*cre_up[iCoeff]+cim_up[iCoeff]*cim_up[iCoeff];
  }
  length *= 2.0;
  length += cre_up[numCoeff]*cre_up[numCoeff];
  length = sqrt(length);
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cre_up[iCoeff] /= length;
    cim_up[iCoeff] /= length;
  }
  */ 

/*==========================================================================*/
/* II) Prepare a random initial wave function                               */
/*     We can't use the readin coeff since it is an eigenfunction of H      */
/*     Pick a random orbital and normalize it.		            */

#ifdef MKL_RANDOM
  VSLStreamStatePtr stream;
  int errcode;
  int seed = 1;
  errcode = vslNewStream(&stream,VSL_BRNG_MCG31,seed);
  errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,2*numCoeff,randTrail,randMin,randMax);
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cre_up[iCoeff] = randTrail[iCoeff-1];
  }
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cim_up[iCoeff] = randTrail[iCoeff-1+numCoeff];
  }
  cim_up[numCoeff] = 0.0;//Keep everything real
#endif
#ifndef MKL_RANDOM
  //whatever random number is good, I'm using Gaussian in this case
  double seed = 15.0;
  int iseed;
  gaussran(2*numCoeff,&iseed,&iseed,&seed,randTrail);
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cre_up[iCoeff] = randTrail[iCoeff];
  }
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cim_up[iCoeff] = randTrail[iCoeff+numCoeff];
  }
  cim_up[numCoeff] = 0.0;
#endif
   
  //Normalize the trail wave function
  for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
    length += cre_up[iCoeff]*cre_up[iCoeff]+cim_up[iCoeff]*cim_up[iCoeff];
  }
  length *= 2.0;
  length += cre_up[numCoeff]*cre_up[numCoeff];
  length = sqrt(length);
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cre_up[iCoeff] /= length;
    cim_up[iCoeff] /= length;
  }

/*==========================================================================*/
/* III) Loop over iteration numbers                                           */

  for(iIter=0;iIter<numIteration;iIter++){

/*--------------------------------------------------------------------------*/
/* ii) Calcualte H|phi>				    */

    //calcCoefForceWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    calcCoefForceWrapReduce(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

/*--------------------------------------------------------------------------*/
/* iii) Calcluate <phi|H|phi>	                                            */
    
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      fcre_up[iCoeff] *= -0.5;
      fcim_up[iCoeff] *= -0.5; 
    }
    fcre_up[numCoeff] *= -1;
    
    energyOld = energy;
    energy = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      energy += fcre_up[iCoeff]*cre_up[iCoeff]+fcim_up[iCoeff]*cim_up[iCoeff];
    }
    energy *= 2.0;
    energy += fcre_up[numCoeff]*cre_up[numCoeff];
    // We already normalize wf to 1.0, so we don't have scaling 0.5 here
    if(iIter%100==0){
      printf("iStep %i Energy %lg\n",iIter,energy);
      fflush(stdout);
    }

/*--------------------------------------------------------------------------*/
/* iv) Copy H|phi> to |phi> and normalize it                                */

    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      cre_up[iCoeff] = fcre_up[iCoeff];
      cim_up[iCoeff] = fcim_up[iCoeff];
    }//endfor iCoeff    
    
    length = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      length += cre_up[iCoeff]*cre_up[iCoeff]+cim_up[iCoeff]*cim_up[iCoeff];
    }
    length *= 2.0;
    length += cre_up[numCoeff]*cre_up[numCoeff];
    length = sqrt(length);
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      cre_up[iCoeff] /= length;
      cim_up[iCoeff] /= length;
    }

/*--------------------------------------------------------------------------*/
/* v) Check Convergence		                                    */
    if(fabs(energyOld-energy)<energyConv)break;
	
  }//endfor iIter

/*==========================================================================*/
/* IV) Restore flags and clean up                                           */

  stodftInfo->energyMax = energy*1.1;
  cpcoeffs_info->nstate_up_proc = numStateUpProc;
  cpcoeffs_info->nstate_dn_proc = numStateDnProc;
 
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cre_up[iCoeff] = coeffReUpBackup[iCoeff];
    cim_up[iCoeff] = coeffImUpBackup[iCoeff];
  } 

  free(coeffReUpBackup);
  free(coeffImUpBackup);
  free(randTrail);
  //fflush(stdout);
  //exit(0);

  printf("Finish estimating energy upperbound.\n");
  printf("The energy upperbound is %lg.\n",stodftInfo->energyMax);
  printf("==============================================\n");
  fflush(stdout);


/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genEnergyMin(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                  CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate expanCoeff for different chemical potentials                */
/* Call this function only for the first process.                        */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  COMMUNICATE *communicate      = &(cp->communicate);

  double randMin        = -1.0;
  double randMax        = 1.0;
  double length         = 0.0;
  double energyConv     = 1.0e-6;
  double energyMax  = stodftInfo->energyMax;
  double energy,energyOld;

  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numIteration   = 1000;
  int iIter;
  int iState,iCoeff,iCoeffStart,index1,index2;

  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *coeffReUpBackup = (double*)cmalloc((numCoeff+1)*sizeof(double));
  double *coeffImUpBackup = (double*)cmalloc((numCoeff+1)*sizeof(double));
  double *randTrail = (double *)cmalloc((2*numCoeff+1)*sizeof(double));

/*==========================================================================*/
/* I) Set parameters and backup                                             */

  printf("==============================================\n");
  printf("Estimate Energy Lowerbound:\n");

  cpcoeffs_info->nstate_up_proc = 1;
  cpcoeffs_info->nstate_dn_proc = 1;

  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    coeffReUpBackup[iCoeff] = cre_up[iCoeff];
    coeffImUpBackup[iCoeff] = cim_up[iCoeff];
  }//endfor iCoeff

/*==========================================================================*/
/* II) Prepare a random initial wave function                               */
/*     We can't use the readin coeff since it is an eigenfunction of H      */
/*     Pick a random orbital and normalize it.                              */
#ifdef MKL_RANDOM
  VSLStreamStatePtr stream;
  int errcode;
  int seed = 1;
  errcode = vslNewStream(&stream,VSL_BRNG_MCG31,seed);
  errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,2*numCoeff,randTrail,randMin,randMax);
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cre_up[iCoeff] = randTrail[iCoeff-1];
  }
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cim_up[iCoeff] = randTrail[iCoeff-1+numCoeff];
  }
  cim_up[numCoeff] = 0.0;//Keep everything real
#endif
#ifndef MKL_RANDOM
  //whatever random number is good, I'm using Gaussian in this case
  double seed = 1.0;
  int iseed;
  gaussran(2*numCoeff,&iseed,&iseed,&seed,randTrail);
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cre_up[iCoeff] = randTrail[iCoeff];
  }
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cim_up[iCoeff] = randTrail[iCoeff+numCoeff];
  }
  cim_up[numCoeff] = 0.0;
#endif

  //Normalize the trail wave function
  for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
    length += cre_up[iCoeff]*cre_up[iCoeff]+cim_up[iCoeff]*cim_up[iCoeff];
  }
  length *= 2.0;
  length += cre_up[numCoeff]*cre_up[numCoeff];
  length = sqrt(length);
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cre_up[iCoeff] /= length;
    cim_up[iCoeff] /= length;
  }


/*==========================================================================*/
/* III) Loop over iteration numbers                                           */

  for(iIter=0;iIter<numIteration;iIter++){

/*--------------------------------------------------------------------------*/
/* ii) Calcualte H|phi>                                                     */

    //calcCoefForceWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    calcCoefForceWrapReduce(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

/*--------------------------------------------------------------------------*/
/* iii) Calcluate <phi|H|phi>                                               */

    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      fcre_up[iCoeff] *= -0.5;
      fcim_up[iCoeff] *= -0.5;
    }
    fcre_up[numCoeff] *= -1;

    energyOld = energy;
    energy = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      energy += fcre_up[iCoeff]*cre_up[iCoeff]+fcim_up[iCoeff]*cim_up[iCoeff];
    }
    energy *= 2.0;
    energy += fcre_up[numCoeff]*cre_up[numCoeff];
    // We already normalize wf to 1.0, so we don't have scaling 0.5 here
    if(iIter%100==0){
      printf("iStep %i Energy %lg\n",iIter,energy);
      fflush(stdout);
    }

/*--------------------------------------------------------------------------*/
/* iv) Calculate (Emax-H)|phi> and normalize it                             */

    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      cre_up[iCoeff] = cre_up[iCoeff]*energyMax-fcre_up[iCoeff];
      cim_up[iCoeff] = cim_up[iCoeff]*energyMax-fcim_up[iCoeff];
    }//endfor iCoeff    

    length = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      length += cre_up[iCoeff]*cre_up[iCoeff]+cim_up[iCoeff]*cim_up[iCoeff];
    }
    length *= 2.0;
    length += cre_up[numCoeff]*cre_up[numCoeff];
    length = sqrt(length);
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      cre_up[iCoeff] /= length;
      cim_up[iCoeff] /= length;
    }

/*--------------------------------------------------------------------------*/
/* v) Check Convergence                                                     */
    if(fabs(energyOld-energy)<energyConv)break;

  }//endfor iIter

/*==========================================================================*/
/* IV) Restore flags and clean up                                           */

  if(energy>0.0)stodftInfo->energyMin = energy*0.9;
  else stodftInfo->energyMin = energy*1.1;

  cpcoeffs_info->nstate_up_proc = numStateUpProc;
  cpcoeffs_info->nstate_dn_proc = numStateDnProc;

  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cre_up[iCoeff] = coeffReUpBackup[iCoeff];
    cim_up[iCoeff] = coeffImUpBackup[iCoeff];
  }

  free(coeffReUpBackup);
  free(coeffImUpBackup);
  free(randTrail);


  printf("Finish estimating energy lowerbound.\n");
  printf("The energy lowerbound is %lg.\n",stodftInfo->energyMin);
  printf("==============================================\n");
  fflush(stdout);
  //exit(0);
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcEnergyChemPot(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                       CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate average energy for each chemical potential, then            */
/* Interpolate the energy for the correct chemical potential             */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  COMMUNICATE *communicate      = &(cp->communicate);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  CPEWALD *cpewald              = &(cp->cpewald);

  int cpLsda         = cpopts->cp_lsda;
  int numStateStoUp  = stodftInfo->numStateStoUp;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numChemPot = stodftInfo->numChemPot;
  int occNumber = stodftInfo->occNumber;
  int myidState         = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int iState,iCoeff,iChem;
  int ioff,iis;

  double tpi = 2.0*M_PI;
  double eke,ekeDn;
  double chemPotTrue = stodftInfo->chemPotTrue;
  double energyKineticTemp,energyNLTemp;

  double *energyKNL = stodftInfo->energyKNL;
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  double *ak2_sm  =  cpewald->ak2_sm;
  double *chemPot = stodftCoefPos->chemPot;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  MPI_Comm commStates = communicate->comm_states;

/*--------------------------------------------------------------------------*/
/* I) Generate kinetic energy and nonlocal pseudopotential energy for       */
/*    each chemical potential.                                              */

  for(iChem=0;iChem<numChemPot;iChem++){
    stat_avg->vrecip = 0.0;
    stat_avg->cp_enl = 0.0;
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      cre_up[iCoeff] = stoWfUpRe[iChem][iCoeff];
      cim_up[iCoeff] = stoWfUpIm[iChem][iCoeff];
      fcre_up[iCoeff] = 0.0;
      fcim_up[iCoeff] = 0.0;
    }//endfor iCoeff
    if(cpLsda==1&&numStateDnProc!=0){
      for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
        cre_dn[iCoeff] = stoWfDnRe[iChem][iCoeff];
        cim_dn[iCoeff] = stoWfDnIm[iChem][iCoeff];
        fcre_dn[iCoeff] = 0.0;
        fcim_dn[iCoeff] = 0.0;
      }//endfor iCoeff
    }//endif cpLsda

    eke = 0.0;
    for(iState=0;iState<numStateUpProc;iState++){
      ioff = iState*numCoeff;
      for(iCoeff=1;iCoeff<=numCoeff-1;iCoeff++){
        iis = ioff+iCoeff;
        eke += 2.0*ak2_sm[iCoeff]*(cre_up[iis]*cre_up[iis]+cim_up[iis]*cim_up[iis]);
      }//endfor iCoeff
    }//endfor iState
    eke *= occNumber*0.5;
    if(cpLsda==1&&numStateDnProc!=0){
      ekeDn = 0.0;
      for(iState=0;iState<numStateDnProc;iState++){
        ioff = iState*numCoeff;
        for(iCoeff=1;iCoeff<=numCoeff-1;iCoeff++){
          iis = ioff+iCoeff;
          ekeDn += 2.0*ak2_sm[iCoeff]*(cre_dn[iis]*cre_dn[iis]+cim_dn[iis]*cim_dn[iis]);
        }//endfor iCoeff
      }//endfor iState
      eke += ekeDn*occNumber*0.5;
    }//endif cpLsda
    stat_avg->cp_eke = eke;

    calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    calcCoefForceExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    stat_avg->cp_enl *= occNumber;

    if(numProcStates>1){
      Reduce(&(stat_avg->cp_enl),&energyNLTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);
      Reduce(&(stat_avg->cp_eke),&energyKineticTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);
    }
    else{
      energyNLTemp = stat_avg->cp_enl;
      energyKineticTemp = stat_avg->cp_eke;
    }

    if(myidState==0){
      energyNLTemp /= numStateStoUp;
      energyKineticTemp /= numStateStoUp;
      energyKNL[iChem] = energyNLTemp+energyKineticTemp;
      //printf("iChem %i chemPot %lg K %lg NL %lg\n",iChem,chemPot[iChem],energyKineticTemp,energyNLTemp);
    }
  }//endfor iChem
  /*
  // debug
  printf("Hartree Energy: %.6lg\n",stat_avg->cp_ehart);
  printf("External Potential Energy(Local Pseudopotential): %.6lg\n",stat_avg->cp_eext);
  printf("Exchange-Correlation Energy: %.6lg\n",stat_avg->cp_exc);
  */


/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcTotEnergy(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                  CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate average energy for each chemical potential, then            */
/* Interpolate the energy for the correct chemical potential		 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  COMMUNICATE *communicate      = &(cp->communicate);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);

  int numChemPot = stodftInfo->numChemPot;
  int myidState         = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int iState,iCoeff,iChem;

  double chemPotTrue = stodftInfo->chemPotTrue;
  double energyKineticTemp,energyNLTemp;
  double energyHartTemp,energyExtTemp,energyExcTemp;
  double energyTrue,energyTotElec;

  double *chemPot = stodftCoefPos->chemPot;
  double *energyKNL = stodftInfo->energyKNL;
  double *lagFunValue = (double*)cmalloc(numChemPot*sizeof(double));

  MPI_Comm commStates = communicate->comm_states; 

/*--------------------------------------------------------------------------*/
/* II) Interpolate the correct energy				            */
    
  if(myidState==0){
    //debug
    /*
    printf("chemPotTrue %lg\n",chemPotTrue);
    for(iChem=0;iChem<numChemPot;iChem++){
      printf("%lg %lg\n",chemPot[iChem],energyKNL[iChem]);
    }
    */
    energyTrue = calcLagrangeInterpFun(numChemPot,chemPotTrue,chemPot,energyKNL,lagFunValue);
  }
  
/*--------------------------------------------------------------------------*/
/* III) Reduce all the other energy terms calculated from density           */
 
  if(numProcStates>1){
    Reduce(&(stat_avg->cp_ehart),&energyHartTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);   
    Reduce(&(stat_avg->cp_eext),&energyExtTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&(stat_avg->cp_exc),&energyExcTemp,1,MPI_DOUBLE,MPI_SUM,0,commStates);
  }
  else{
    energyHartTemp = stat_avg->cp_ehart;
    energyExtTemp = stat_avg->cp_eext;
    energyExcTemp = stat_avg->cp_exc;
  }

/*--------------------------------------------------------------------------*/
/* IV) Output the energy Term					            */

  if(myidState==0){
    energyTotElec = energyTrue+energyHartTemp+energyExtTemp+energyExcTemp;
    printf("==============================================\n");
    printf("Output Energy\n");
    printf("==============================================\n");
    printf("Kinetic Energy+NLPP: %.6lg\n",energyTrue);
    printf("Hartree Energy:      %.6lg\n",energyHartTemp);
    printf("Ext Energy:          %.6lg\n",energyExtTemp); 
    printf("Ex-Cor Energy:       %.6lg\n",energyExcTemp); 
    printf("Total Elec Energy    %.6lg\n",energyTotElec);
    printf("==============================================\n");
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

