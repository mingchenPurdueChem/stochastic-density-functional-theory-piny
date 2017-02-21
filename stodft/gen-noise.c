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
 
  int iStat,iCoeff,iOff;
  int cpLsda = cpopts->cp_lsda;
  int numStatUpProc = cpcoeffs_info->nstate_up_proc;
  int numStatDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff = cpcoeffs_info->ncoef;
  int numStatUpTot = numStatUpProc*numCoeff;
  int numStatDnTot = numStatDnProc*numCoeff;
  int numProcStates             = communicate->np_states;
  int myidState                 = communicate->myid_state;
  int numRandTot = stodftInfo->numRandTot;
  int *noiseSendCounts = stodftInfo->noiseSendCounts;
  int *noiseDispls = stodftInfo->noiseDispls;
  MPI_Comm comm_states   =    communicate->comm_states;
  
  double ranValue = 1.0/sqrt(2.0);
  double *randNumTot,*randNum;
  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;

  if(myidState==0)randNumTot = (double*)cmalloc(numRandTot*sizeof(double));
  if(numProcStates>1){
    randNum = (double*)cmalloc(noiseSendCounts[myidState]*sizeof(double));
  }
  else{
    randNum = randNumTot;
  }

  if(myidState==0){  
#ifdef MKL_RANDOM
    VSLStreamStatePtr stream;
    int errcode;
    int seed = 1;
    errcode = vslNewStream(&stream,VSL_BRNG_MCG31,seed);
    errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,numRandTot,randNumTot,-1.0,1.0);
#endif
#ifndef MKL_RANDOM
    //whatever random number is good, I'm using Gaussian in this case
    double seed = 8.3;
    //double seed = 2.5;
    int iseed;
    printf("numRandTot %i numStateUpTot %i numCoeff %i\n",numRandTot,numStatUpTot,numCoeff);
    fflush(stdout);
    gaussran2(numRandTot,&iseed,&iseed,&seed,randNumTot);
#endif
  }
  if(numProcStates>1){
    Scatterv(randNumTot,noiseSendCounts,noiseDispls,MPI_DOUBLE,
             randNum,noiseSendCounts[myidState],MPI_DOUBLE,0,comm_states);
    Barrier(comm_states);
    if(myidState==0)free(randNumTot);
  }

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
  
  free(randNum);

/*--------------------------------------------------------------------------*/  
}/*end routine*/
/*==========================================================================*/
