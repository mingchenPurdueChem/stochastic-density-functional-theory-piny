/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: filter.c                                       */
/*                                                                          */
/* This routine calculate the spectral range of KS Hamiltonian.             */
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
void genEnergyMax(CP *cp,CLASS *class,GENERAL_DATA *general_data,
      CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate expanCoeff for different chemical potentials                */
/* Call this function only for the first process.        */
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
  double energyConv     = 1.0e-5;
  double energy	    = 0.0;
  double energyOld;

  int numStateUpProc  = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc  = cpcoeffs_info->nstate_dn_proc;
  int numCoeff        = cpcoeffs_info->ncoef;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int cpLsda          = cpopts->cp_lsda;
  int numCoeffUpTot   = numStateUpProc*numCoeff;
  int numCoeffDnTot   = numStateDnProc*numCoeff;
  int myidState = communicate->myid_state;
  int iScf = stodftInfo->iScf;

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
  double *testWfMaxRe = stodftCoefPos->testWfMaxRe;
  double *testWfMaxIm = stodftCoefPos->testWfMaxIm;

/*==========================================================================*/
/* I) Set parameters and backup		        */

  if(myidState==0){
    printf("==============================================\n");
    printf("Estimate Energy Upperbound: \n");
    fflush(stdout);
  }

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
  
  //debug
  /*
  length = 0.0;
  for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
    length += cre_up[iCoeff]*cre_up[iCoeff]+cim_up[iCoeff]*cim_up[iCoeff];
  }
  length *= 2.0;
  length += cre_up[numCoeff]*cre_up[numCoeff];
  length = sqrt(length);
  printf("length %lg\n",length);
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    cre_up[iCoeff] /= length;
    cim_up[iCoeff] /= length;
  }
  */ 

/*==========================================================================*/
/* II) Prepare a random initial wave function                               */
/*     We can't use the readin coeff since it is an eigenfunction of H      */
/*     Pick a random orbital and normalize it.	                */

  if(iScf==1){
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
  }
  else{
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      cre_up[iCoeff] = testWfMaxRe[iCoeff-1];
      cim_up[iCoeff] = testWfMaxIm[iCoeff-1];
    }
  }
   
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

  energyOld = 10000000.0;
  energy = 100000.0;
  iIter = 0;

  while(fabs(energyOld-energy)>energyConv){
  //for(iIter=0;iIter<numIteration;iIter++){

/*--------------------------------------------------------------------------*/
/* ii) Calcualte H|phi>		        */

    //calcCoefForceWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    /*
    calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    calcCoefForceWrapReduce(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    */
    calcCoefForceWrapSCF(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    //printf("1111111 cp_enl %lg\n",general_data->stat_avg.cp_enl);

/*--------------------------------------------------------------------------*/
/* iii) Calcluate <phi|H|phi>	                                            */
    
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      fcre_up[iCoeff] *= -0.25;
      fcim_up[iCoeff] *= -0.25; 
    }
    fcre_up[numCoeff] *= -0.5;
    
    energyOld = energy;
    energy = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      energy += fcre_up[iCoeff]*cre_up[iCoeff]+fcim_up[iCoeff]*cim_up[iCoeff];
    }
    energy *= 2.0;
    energy += fcre_up[numCoeff]*cre_up[numCoeff];
    // We already normalize wf to 1.0, so we don't have scaling 0.5 here
    if(myidState==0){
      //if(iIter%100==0){
      if(iIter%100==0){
	printf("iStep %i Energy %lg\n",iIter,energy);
	fflush(stdout);
      }
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
/* v) Check Convergence	                                        */
    iIter += 1;
    //if(fabs(energyOld-energy)<energyConv)break;
    
  }//endfor iIter

/*==========================================================================*/
/* IV) Restore flags and clean up                                           */

  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    testWfMaxRe[iCoeff-1] = cre_up[iCoeff];
    testWfMaxIm[iCoeff-1] = cim_up[iCoeff];
  }

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

  if(myidState==0){
    printf("Finish estimating energy upperbound.\n");
    printf("The energy upperbound is %lg.\n",stodftInfo->energyMax);
    printf("==============================================\n");
    fflush(stdout);
  }


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
  double energyConv     = 1.0e-5;
  double energyMax  = stodftInfo->energyMax;
  double energy,energyOld;

  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int myidState = communicate->myid_state;
  int numIteration   = 1000;
  int iIter;
  int iState,iCoeff,iCoeffStart,index1,index2;
  int iScf = stodftInfo->iScf;

  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *coeffReUpBackup = (double*)cmalloc((numCoeff+1)*sizeof(double));
  double *coeffImUpBackup = (double*)cmalloc((numCoeff+1)*sizeof(double));
  double *randTrail = (double *)cmalloc((2*numCoeff+1)*sizeof(double));
  double *testWfMinRe = stodftCoefPos->testWfMinRe;
  double *testWfMinIm = stodftCoefPos->testWfMinIm;

/*==========================================================================*/
/* I) Set parameters and backup                                             */

  if(myidState==0){
    printf("==============================================\n");
    printf("Estimate Energy Lowerbound:\n");
  }

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
  if(iScf==1){
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
  }
  else{
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      cre_up[iCoeff] = testWfMinRe[iCoeff-1];
      cim_up[iCoeff] = testWfMinIm[iCoeff-1];
    }
  }

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
  energyOld = 10000000.0;
  energy = 100000.0;
  iIter = 0;

  while(fabs(energyOld-energy)>energyConv){
  //for(iIter=0;iIter<numIteration;iIter++){

/*--------------------------------------------------------------------------*/
/* ii) Calcualte H|phi>                                                     */

    //calcCoefForceWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    /*
    calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    calcCoefForceWrapReduce(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    */
    calcCoefForceWrapSCF(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

/*--------------------------------------------------------------------------*/
/* iii) Calcluate <phi|H|phi>                                               */

    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      fcre_up[iCoeff] *= -0.25;
      fcim_up[iCoeff] *= -0.25;
    }
    fcre_up[numCoeff] *= -0.5;

    energyOld = energy;
    energy = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      energy += fcre_up[iCoeff]*cre_up[iCoeff]+fcim_up[iCoeff]*cim_up[iCoeff];
    }
    energy *= 2.0;
    energy += fcre_up[numCoeff]*cre_up[numCoeff];
    // We already normalize wf to 1.0, so we don't have scaling 0.5 here
    if(myidState==0){
      if(iIter%100==0){
    printf("iStep %i Energy %lg\n",iIter,energy);
    fflush(stdout);
      }
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
    iIter += 1;
    //if(fabs(energyOld-energy)<energyConv)break;

  }//endfor iIter

/*==========================================================================*/
/* IV) Restore flags and clean up                                           */
  
  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    testWfMinRe[iCoeff-1] = cre_up[iCoeff];
    testWfMinIm[iCoeff-1] = cim_up[iCoeff];
  }
  stodftInfo->eigValMin = energy;

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

  if(myidState==0){
    printf("Finish estimating energy lowerbound.\n");
    printf("The energy lowerbound is %lg.\n",stodftInfo->energyMin);
    printf("==============================================\n");
    fflush(stdout);
  }
  //exit(0);
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


