/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: gen-noise.c                                  */
/*                                                                          */
/* Subroutine to generate noise orbitals                                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_stodft_local.h"
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genNoiseOrbital(CP *cp,CPCOEFFS_POS *cpcoeffs_pos)
/*========================================================================*/
{/*begin routine*/
/**************************************************************************/
/* We shall direct work on reciprocal space, instead of generating noise  */
/* on real space since the high frequency(could be important) part will   */
/* discard when reduce the reciprocal lattice from a cubic to a ball.     */
/* To satisfies the completeness the coeffecients are required to be	  */
/* a equally random choice between 1/sqrt(2) and -1/sqrt(2) except k=0    */
/* k=0 term should be a random number of +/-1.		      */
/**************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
/*------------------------------------------------------------------------*/
#include "../typ_defs/typ_mask.h"
  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR        *cpscr        = &(cp->cpscr);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  COMMUNICATE   *communicate      = &(cp->communicate);
 
  int iStat,iCoeff,iOff,iOff2,iProc;
  int cpLsda = cpopts->cp_lsda;
  int numRandNum;
  int numStatUpProc = cpcoeffs_info->nstate_up_proc;
  int numStatDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff = cpcoeffs_info->ncoef;
  int numStatUpTot = numStatUpProc*numCoeff;
  int numStatDnTot = numStatDnProc*numCoeff;
  int numProcStates             = communicate->np_states;
  int myidState                 = communicate->myid_state;
  MPI_Comm comm_states   =    communicate->comm_states;
 
  double ranValue = 1.0/sqrt(2.0);
  double *randNumSeedTot = stodftInfo->randSeedTot;
  double *randNum;
  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;

  numRandNum = numStatUpTot*2;
  if(cpLsda==1)numRandNum += numStatDnTot*2;
  randNum = (double*)cmalloc(numRandNum*sizeof(double));

  // Generate the random number seeds from the given random number seed
  if(myidState==0){  
#ifdef MKL_RANDOM
    VSLStreamStatePtr stream;
    int errcode;
    int seed = (int)(stodftInfo->randSeed);
    errcode = vslNewStream(&stream,VSL_BRNG_MCG31,seed);
    errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,
			    numProcStates,randNumSeedTot,0.0,100000.0);
#endif
#ifndef MKL_RANDOM
    double seed = stodftInfo->randSeed;
    //printf("seed!!!!!!!! %lg\n",seed);
    //whatever random number is good, I'm using Gaussian in this case
    //double seed = 8.3;
    //double seed = 2.5;
    int iseed;
    //printf("numRandTot %i numStateUpTot %i numCoeff %i\n",numRandTot,numStatUpTot,numCoeff);
    //fflush(stdout);
    gaussran2(numProcStates,&iseed,&iseed,&seed,randNumSeedTot);
    for(iProc=0;iProc<numProcStates;iProc++){
      randNumSeedTot[iProc] = randNumSeedTot[iProc]*randNumSeedTot[iProc]*10000.0;
    }
#endif
  }
  // Bcast the random number seeds
  if(numProcStates>1){
    Bcast(randNumSeedTot,numProcStates,MPI_DOUBLE,0,comm_states);
    Barrier(comm_states);
  }
  
  // Generate random orbital
#ifdef MKL_RANDOM
  VSLStreamStatePtr streamNew;
  int errcodeNew;
  int seedNew = (int)randNumSeedTot[myidState];
  errcode = vslNewStream(&streamNew,VSL_BRNG_MCG31,seedNew);
  errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,streamNew,
			  numRandNum,randNum,-1.0,1.0);
#endif
#ifndef MKL_RANDOM
  double seedNew = randNumSeedTot[myidState];
  int iseedNew;
  //printf("seedNew %p\n",randNumSeedTot);
  gaussran2(numRandNum,&iseedNew,&iseedNew,&seedNew,randNum);
#endif

  for(iStat=0;iStat<numStatUpProc;iStat++){
    iOff = iStat*numCoeff;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      if(randNum[2*(iOff+iCoeff-1)]>0)coeffReUp[iOff+iCoeff] = ranValue;
      else coeffReUp[iOff+iCoeff] = -ranValue;
      if(randNum[2*(iOff+iCoeff-1)+1]>0)coeffImUp[iOff+iCoeff] = ranValue;
      else coeffImUp[iOff+iCoeff] = -ranValue;
    }
    if(randNum[2*(iStat*numCoeff+numCoeff-1)]>0)coeffReUp[iStat*numCoeff+numCoeff] = 1.0;
    else coeffReUp[iStat*numCoeff+numCoeff] = -1.0;
    coeffImUp[iStat*numCoeff+numCoeff] = 0.0;
  }
  if(cpLsda==1){
    for(iStat=0;iStat<numStatDnProc;iStat++){
      iOff = iStat*numCoeff;
      iOff2 = numStatUpProc*numCoeff+iStat*numCoeff;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	if(randNum[2*(iOff2+iCoeff-1)]>0)coeffReDn[iOff+iCoeff] = ranValue;
	else coeffReDn[iOff+iCoeff] = -ranValue;
	if(randNum[2*(iOff2+iCoeff-1)+1]>0)coeffImDn[iOff+iCoeff] = ranValue;
	else coeffImDn[iOff+iCoeff] = -ranValue;
      }
      if(randNum[2*(iOff2+numCoeff-1)]>0)coeffReDn[iStat*numCoeff+numCoeff] = 1.0;
      else coeffReDn[iStat*numCoeff+numCoeff] = -1.0;
      coeffImDn[iStat*numCoeff+numCoeff] = 0.0;
    }
  }
  
  free(randNum);

/*--------------------------------------------------------------------------*/  
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genNoiseOrbitalReal(CP *cp,CPCOEFFS_POS *cpcoeffs_pos)
/*========================================================================*/
{/*begin routine*/
/**************************************************************************/
/* Let's try to work on real space (put +-1) and FFT transform to k space */
/**************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
/*------------------------------------------------------------------------*/
#include "../typ_defs/typ_mask.h"
  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR        *cpscr        = &(cp->cpscr);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  COMMUNICATE   *communicate      = &(cp->communicate);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm = &(cp->cp_sclr_fft_pkg3d_sm);
 
  int iStat,iGrid,iOff,iOff2,iProc,iCoeff;
  int cpLsda = cpopts->cp_lsda;
  int numRandNum;
  int numStatUpProc = cpcoeffs_info->nstate_up_proc;
  int numStatDnProc = cpcoeffs_info->nstate_dn_proc;
  int nfft	    = cp_para_fft_pkg3d_lg->nfft;
  int nfft2	    = nfft/2;
  int numCoeff = cpcoeffs_info->ncoef;
  int numStatUpTot = numStatUpProc*numCoeff;
  int numStatDnTot = numStatDnProc*numCoeff;
  int numProcStates             = communicate->np_states;
  int myidState                 = communicate->myid_state;
  MPI_Comm comm_states   =    communicate->comm_states;
 
  double ranValue = sqrt(nfft2);
  double *randNumSeedTot = stodftInfo->randSeedTot;
  double *randNum;
  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;
  double *zfft      = cpscr->cpscr_wave.zfft;
  double *zfft_temp = cpscr->cpscr_wave.zfft_tmp;

  numRandNum = numStatUpProc*nfft2;
  if(cpLsda==1)numRandNum += numStatDnProc*nfft2;
  randNum = (double*)cmalloc(numRandNum*sizeof(double));

  // Generate the random number seeds from the given random number seed
  if(myidState==0){  
#ifdef MKL_RANDOM
    VSLStreamStatePtr stream;
    int errcode;
    int seed = (int)(seed = stodftInfo->randSeed);
    errcode = vslNewStream(&stream,VSL_BRNG_MCG31,seed);
    errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,
			    numProcStates,randNumSeedTot,0.0,100000.0);
#endif
#ifndef MKL_RANDOM
    double seed = stodftInfo->randSeed;
    //printf("seed!!!!!!!! %lg\n",seed);
    //whatever random number is good, I'm using Gaussian in this case
    //double seed = 8.3;
    //double seed = 2.5;
    int iseed;
    //printf("numRandTot %i numStateUpTot %i numCoeff %i\n",numRandTot,numStatUpTot,numCoeff);
    //fflush(stdout);
    double x;
    gaussran2(numProcStates,&iseed,&iseed,&seed,randNumSeedTot);
    for(iProc=0;iProc<numProcStates;iProc++){
      x = randNumSeedTot[iProc]*randNumSeedTot[iProc];
      if(x>=1.0)randNumSeedTot[iProc] = x*100.0;
      else randNumSeedTot[iProc] = 100.0/x;
    }
#endif
  }
  // Bcast the random number seeds
  if(numProcStates>1){
    Bcast(randNumSeedTot,numProcStates,MPI_DOUBLE,0,comm_states);
    Barrier(comm_states);
  }
  
  // Generate random orbital
#ifdef MKL_RANDOM
  VSLStreamStatePtr streamNew;
  int errcodeNew;
  int seedNew = (int)randNumSeedTot[myidState];
  errcode = vslNewStream(&streamNew,VSL_BRNG_MCG31,seedNew);
  errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,streamNew,
			  numRandNum,randNum,-1.0,1.0);
#endif
#ifndef MKL_RANDOM
  double seedNew = randNumSeedTot[myidState];
  //printf("proc %i seed %lg\n",myidState,seedNew);
  int iseedNew;
  gaussran2(numRandNum,&iseedNew,&iseedNew,&seedNew,randNum);
  //printf("randNum[1] %lg\n",randNum[1]);
#endif
  for(iStat=1;iStat<=numStatUpTot;iStat++){
    coeffReUp[iStat] = 0.0;
    coeffImUp[iStat] = 0.0;
  }
  //debug
  /*
  char fileNameRand[100];
  FILE *fileRand;
  sprintf(fileNameRand,"rand-%i",myidState);
  fileRand = fopen(fileNameRand,"w");
  for(iStat=0;iStat<numStatUpProc;iStat++){
    for(iGrid=0;iGrid<nfft2;iGrid++){
      if(randNum[iStat*nfft2+iGrid]<0.0)fprintf(fileRand,"-1.0\n");
      else fprintf(fileRand,"1.0\n");
    }
  }
  fclose(fileRand);
  */
  //debug

  
  /*
  char fileNameRand[100];
  FILE *fileRand;
  printf("I'm reading noise orbital\n");
  sprintf(fileNameRand,"rand-%i",myidState);
  fileRand = fopen(fileNameRand,"r");
  for(iStat=0;iStat<numStatUpProc;iStat++){
    for(iGrid=0;iGrid<nfft2;iGrid++){
      fscanf(fileRand,"%lg",&randNum[iStat*nfft2+iGrid]);
      //if(randNum[iStat*nfft2+iGrid]<0.0)fprintf(fileRand,"-1.0\n");
      //else fprintf(fileRand,"1.0\n");
    }
  }
  fclose(fileRand);
  */
  
  
  for(iStat=0;iStat<numStatUpProc;iStat++){
    for(iGrid=0;iGrid<nfft2;iGrid++){
      if(randNum[iStat*nfft2+iGrid]<0.0) zfft[iGrid*2+1] = -ranValue;
      else zfft[iGrid*2+1] = ranValue;
      zfft[iGrid*2+2] = 0.0;
    }
    iOff = iStat*numCoeff;
    para_fft_gen3d_bck_to_g(zfft,zfft_temp,cp_sclr_fft_pkg3d_sm);
    sngl_upack_coef_sum(&coeffReUp[iOff],&coeffImUp[iOff],zfft,
                            cp_sclr_fft_pkg3d_sm);
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      coeffReUp[iOff+iCoeff] *= 0.25;
      coeffImUp[iOff+iCoeff] *= 0.25;
    }
    coeffReUp[iOff+numCoeff] *= 0.5;
    coeffImUp[iOff+numCoeff] = 0.0;
    //printf("myid %i iStat %i %.16lg\n",myidState,iStat,coeffReUp[iStat*numCoeff+1]);
  }
  if(cpLsda==1){
    for(iStat=1;iStat<=numStatDnTot;iStat++){
      coeffReDn[iStat] = 0.0;
      coeffImDn[iStat] = 0.0;
    }
    for(iStat=0;iStat<numStatDnProc;iStat++){
      iOff = numStatUpProc*nfft2+iStat*nfft2;
      for(iGrid=0;iGrid<nfft2;iGrid++){
	if(randNum[iOff+iGrid]<0.0)zfft[iGrid*2+1] = -ranValue;
	else zfft[iGrid*2+1] = ranValue;
	zfft[iGrid*2+2] = 0.0;
      }
      iOff = iStat*numCoeff;
      para_fft_gen3d_bck_to_g(zfft,zfft_temp,cp_sclr_fft_pkg3d_sm);
      sngl_upack_coef_sum(&coeffReDn[iOff],&coeffImDn[iOff],zfft,
			      cp_sclr_fft_pkg3d_sm);
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	coeffReUp[iOff+iCoeff] *= 0.25;
	coeffImUp[iOff+iCoeff] *= 0.25;
      }
      coeffReUp[iOff+numCoeff] *= 0.5;
      coeffImUp[iOff+numCoeff] = 0.0;
    }
  }
  /* 
  char fileNameRand[100];
  FILE *fileRand;
  double test;
  //sprintf(fileNameRand,"rand-%i",myidState);
  for(iStat=1;iStat<=numStatUpTot;iStat++){
    coeffReUp[iStat] = 0.0;
    coeffImUp[iStat] = 0.0;
  }
  fileRand = fopen("rand-all","r");
  for(iStat=0;iStat<numStatUpProc;iStat++){
    for(iGrid=0;iGrid<nfft2;iGrid++){
      fscanf(fileRand,"%lg",&test);
      if(test<0.0) zfft[iGrid*2+1] = -ranValue;
      else zfft[iGrid*2+1] = ranValue;
      zfft[iGrid*2+2] = 0.0;
    }
    iOff = iStat*numCoeff;
    para_fft_gen3d_bck_to_g(zfft,zfft_temp,cp_sclr_fft_pkg3d_sm);
    sngl_upack_coef_sum(&coeffReUp[iOff],&coeffImUp[iOff],zfft,
                        cp_sclr_fft_pkg3d_sm);
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      coeffReUp[iOff+iCoeff] *= 0.25;
      coeffImUp[iOff+iCoeff] *= 0.25;
    }
    coeffReUp[iOff+numCoeff] *= 0.5;
    coeffImUp[iOff+numCoeff] = 0.0;
  }
  */
  free(randNum);
  if(numProcStates>1)Barrier(comm_states);
  //fflush(stdout);
  //exit(0);

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

