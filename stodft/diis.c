/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: diis.c	                                    */
/*                                                                          */
/* This routine calculate diis density update.                              */
/*	                                                                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_stodft_local.h"

#include "complex.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genDensityMix(CP *cp,int iScf)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the routine to calculate the diis density                     */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPSCR *cpscr		        = &(cp->cpscr);
  COMMUNICATE *communicate      = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  
  int iGrid;
  int rhoRealGridNum = stodftInfo->rhoRealGridNum;
  int numDiis = stodftInfo->numDiis;
  int numStepMix = stodftInfo->numStepMix;
  int cpLsda = cpopts->cp_lsda;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;

  double mixRatio1 = stodftInfo->mixRatio;
  double mixRatio2 = 1.0-mixRatio1;

  double *rhoUpCorrect = stodftCoefPos->rhoUpCorrect;
  double *rhoDnCorrect = stodftCoefPos->rhoDnCorrect;
  double *rhoUp = cpscr->cpscr_rho.rho_up;
  double *rhoDn = cpscr->cpscr_rho.rho_dn;

  double **rhoUpBank = stodftCoefPos->rhoUpBank;
  double **rhoDnBank = stodftCoefPos->rhoDnBank;
  double **rhoUpErr  = stodftCoefPos->rhoUpErr;
  double **rhoDnErr  = stodftCoefPos->rhoDnErr;

  //updateBank(rhoUpCorrect,rhoUpBank,iScf);
  //updateErr(rhoUpBank,rhoUpErr,iScf);

  if(iScf==0){//Initial Step
    updateBank(todftInfo,stodftCoefPos,rhoUpCorrect,rhoUpBank);
  }
  else if(iScf<=numStepMix){//Mixing Step
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
      rhoUp[iGrid+1] = rhoUpCorrect[iGrid]*mixRatio1+rhoUpBank[0][iGrid]*mixRatio2;
    }
    updateErr(stodftInfo,stodftCoefPos,rhoUpCorrect,rhoUpErr,rhoUpBank);
    updateBank(stodftInfo,stodftCoefPos,rhoUpCorrect,rhoUpBank);
  }
  else{//diis
    updateErr(stodftInfo,stodftCoefPos,rhoUpCorrect,rhoUpErr,rhoUpBank);
    calcDensityDiis(cp,rhoUpBank,rhoUpErr);   
    updateBank(stodftInfo,stodftCoefPos,rhoUpCorrect,rhoUpBank);
  }
  if(cpLsda==1&&numStateDnProc>0){
    //updateBank(rhoDnCorrect,rhoDnBank,iScf);
    //updateErr(rhoDnBank,rhoDnErr,iScf);

    if(iScf==0){//Initial Step
      updateBank(stodftInfo,stodftCoefPos,rhoDnCorrect,rhoDnBank);
    }
    else if(iScf<=numStepMix){//Mixing Step
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
	rhoDn[iGrid+1] = rhoDnBank[0][iGrid]*mixRatio1+rhoDnBank[1][iGrid]*mixRatio2;
      }
      updateErr(stodftInfo,stodftCoefPos,rhoDnCorrect,rhoDnErr,rhoDnBank);
      updateBank(stodftInfo,stodftCoefPos,rhoDnCorrect,rhoDnBank);
    }
    else{//diis
      updateErr(stodftInfo,stodftCoefPos,rhoDnCorrect,rhoDnErr,rhoDnBank);
      calcDensityDiis(cp,rhoDnBank,rhoDnErr);      
      updateBank(stodftInfo,stodftCoefPos,rhoDnCorrect,rhoDnBank);
    }
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void updateBank(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos,
		double *rho,double **rhoBank)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the routine to calculate the diis density                     */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
  int numDiis = stodftInfo->numDiis;
  int rhoRealGridNum = stodftInfo->rhoRealGridNum;
  int iDiis;
  
  if(rhoBank[numDiis-1]!=NULL)free(rhoBank[numDiis-1]);
  for(iDiis=numDiis-2;iDiis>-1;iDiis--)rhoBank[iDiis+1] = rhoBank[iDiis];
  rhoBank[0] = (double*)cmalloc(rhoRealGridNum*sizeof(double));
  mencpy(rhoBank[0],rho,rhoRealGridNum*sizeof(double));

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void updateErr(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos,
                double *rho,double **rhoErr,double **rhoBank)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the routine to calculate the diis density                     */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int numDiis = stodftInfo->numDiis;
  int rhoRealGridNum = stodftInfo->rhoRealGridNum;
  int iDiis,iGrid;
  
  if(rhoErr[numDiis-1]!=NULL)free(rhoErr[numDiis-1]);
  for(iDiis=numDiis-2;iDiis>-1;iDiis--)rhoErr[iDiis+1] = rhoErr[iDiis];
  rhoErr[0] = (double*)cmalloc(rhoRealGridNum*sizeof(double));
  memcpy(rhoErr[0],rho,rhoRealGridNum*sizeof(double));
  
  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoErr[0][iGrid] -= rhoBank[0][iGrid];
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcDensityDiis(CP *cp,double **rhoBank,double **rhoErr)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the routine to calculate the diis density                     */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPSCR *cpscr                  = &(cp->cpscr);
  COMMUNICATE *communicate      = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);

  int numDiis = stodftInfo->numDiis;
  int myidState = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int rhoRealGridNum = stodftInfo->rhoRealGridNum;
  int diisMatrixCalcFullFlag = stodftInfo->diisMatrixCalcFullFlag;
  int iDiis,jDiis,index1,index2;
  int incx = 1;
  int incy = 1;

  double dotProc;
  double dotTot;   
  double *diisMatrixTemp = (double*)cmalloc((numDiis-1)*(numDiis-1)*sizeof(double));
  double *diisMatrix = stodftCoefPos->diisMatrix;
  double *svdLinSol = (double*)cmalloc((numDiis+1)*sizeof(double));
  double *b = (double*)cmalloc((numDiis+1)*sizeof(double));
  double *diisCoeff = stodftCoefPos->diisCoeff;
  double *rhoUp = cpscr->cpscr_rho.rho_up;
  double *rhoDn = cpscr->cpscr_rho.rho_dn;
  double *rhoUpCorrect = stodftCoefPos->rhoUpCorrect;
  double *rhoDnCorrect = stodftCoefPos->rhoDnCorrect;

  MPI_Comm commStates = communicate->comm_states;

/*==========================================================================*/
/* I) For the first time, calculate the whole diisMatrix		    */

  if(diisMatrixCalcFullFlag==1){
    // Calculate the dot prod part
    for(iDiis=0;iDiis<numDiis;iDiis++){
      for(jDiis=iDiis;jDiis<numDiis;jDiis++){
	dotProc = DDOT(&rhoRealGridNum,rhoErr[iDiis],&incx,rhoErr[jDiis],&incy);
	if(numProcStates>1)Reduce(&dotProc,&dotTot,1,MPI_DOUBLE,MPI_SUM,0,commStates);
	else dotTot = dotProc;
	if(myidState==0){
	  index1 = iDiis*(numDiis+1)+jDiis;
	  index2 = jDiis*(numDiis+1)+iDiis;
	  diisMatrix[index1] = dotTot;
	  diisMatrix[index2] = dotTot;
	}//endif myidState
      }//endfor jDiis
    }//endfor iDiis
    // Calculate the Lagrange Multiplier part
    for(iDiis=0;iDiis<numDiis;iDiis++){
      diisMatrix[numDiis*(numDiis+1)+iDiis] = -1.0;
      diisMatrix[iDiis*(numDiis+1)+numDiis] = -1.0;
    }
    diisMatrix[numDiis*(numDiis+1)+numDiis] = 0.0;
    // reset the flag so that we don't do it again
    stodftInfo->diisMatrixCalcFullFlag = 0;
  }
/*==========================================================================*/
/* II) Update partially the diisMatrix		                            */
  else{
/*--------------------------------------------------------------------------*/
/* i) Transfer the [0,numDiis-2]*[0,numDiis-2] block to			    */
/*    to [1,numDiis-1]*[i,numDiis-1] block				    */

    if(myidState==0){
      for(iDiis=0;iDiis<numDiis-1;iDiis++){
	for(jDiis=0;jDiis<numDiis-1;jDiis++){
	  index1 = iDiis*(numDiis+1)+jDiis; //index in diisMatrix
	  index2 = iDiis*(numDiis-1)+jDiis; //index in diisMatrixTemp
	  diisMatrixTemp[index2] = diisMatrix[index1];
	}//endfor jDiis
      }//endfor iDiis

      for(iDiis=0;iDiis<numDiis-1;iDiis++){
	for(jDiis=0;jDiis<numDiis-1;jDiis++){
	  index1 = iDiis*(numDiis-1)+jDiis; //index in diisMatrixTemp
	  index2 = (iDiis+1)*(numDiis+1)+jDiis+1 //index in new diisMatrix
	  diisMatrix[index2] = diisMatrixTemp[index1];
	}//endfor jDiis
      }//endfor iDiis
    }//endif myidState

/*--------------------------------------------------------------------------*/
/* ii) Calculate the dot product of rhoErr[0] with rhoErr[i], 0<=i<numDiis  */
  
    for(iDiis=0;iDiis<numDiis;iDiis++){
      dotProc = DDOT(&rhoRealGridNum,rhoErr[0],&incx,rhoErr[iDiis],&incy);
      if(numProcStates>1)Reduce(&dotProc,&dotTot,1,MPI_DOUBLE,MPI_SUM,0,commStates);
      else dotTot = dotProc;
      if(myidState==0){
	diisMatrix[iDiis] = dotTot;
	diisMatrix[iDiis*(numDiis+1)] = dotTot;
      }//endif myidState
    }//endfor iDiis
  }//endif diisMatrixCalcFullFlag
/*==========================================================================*/
/* III) Safely reverse the diisMatrix with SVD				    */

  if(myidState==0){
    for(iDiis=0;iDiis<numDiis;iDiis++)b[iDiis] = 0.0;
    b[numDiis] = -1.0;

    matrixInvSVD(diisMatrix,b,svdLinSol,numDiis+1);       
  }

  if(numProcStates>1){
    Bcast(svdLinSol,numDiis+1,MPI_DOUBLE,0,commStates);
  }

  for(iDiis=0;iDiis<numDiis;iDiis++){
    diisCoeff[iDiis] = svdLinSol[iDiis];    
  }
  stodftInfo->lambdaDiis = svdLinSol[numDiis];

/*==========================================================================*/
/* IV) Linear Combinination of all densities                                */
  
  // I may change this 
  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
    rhoUp[iGrid+1] = diisCoeff[0]*rhoUpCorrect[iGrid];
    for(iDiis=1;iDiis<numDiis;iDiis++){
      rhoUp[iGrid+1] = diisCoeff[iDiis]*rhoBank[iDiis-1][iGrid];
    }
  } 

  free(diisMatrixTemp);
  free(svdLinSol);
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void matrixInvSVD(double *mat,double *b,double *x,int ndim)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This routine solve the linear equation via SVD                        */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  char jobu = 'A';
  char jobvt = 'A';
  char trans;
  int m = ndim;
  int n = ndim;
  int lda = ndim;
  int ldu = ndim;
  int ldvt = ndim;
  int lwork = 10*ndim;
  int info;
  int incx,incy;
  int i,j,k;
  double alpha;
  double beta;

  double *A = (double*)cmalloc(ndim*ndim*sizeof(double));
  double *s = (double*)cmalloc(ndim*sizeof(double));
  double *u = (double*)cmalloc(ndim*ndim*sizeof(double));
  double *vt = (double*)cmalloc(ndim*ndim*sizeof(double));
  double *work = (double*)cmalloc(lwork*sizeof(double));
  double *u_b = (double*)cmalloc(ndim*sizeof(double));

  
  memcpy(A,mat,ndim*ndim*sizeof(double));

  DGESVD(&jobu,&jobvt,&m,&n,A,&lda,s,u,&ldu,vt,&ldvt,work,&lwork,&info);
  
  if(info==0){
    trans = 't';
    alpha = 1.0;
    beta = 0.0;
    incx = 1;
    incy = 1;
    DGEMV(&trans,&ndim,&ndim,&alpha,u,&lda,b,&incx,&beta,u_b,&incy);
    for(i=0;i<ndim;i++)u_b[i] /= s[i];
    DGEMV(&trans,&ndim,&ndim,&alpha,vt,&lda,u_b,&incx,&beta,x,&incy);
  }
  else if(info<0){
    printf("The %ith parameter is illegal.\n",info);
  }
  else{
    printf("Bidiagonal form B fails to converge.\n");
  }
  free(A);
  free(s);
  free(u);
  free(vt);
  free(work); 
  free(u_b);
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

