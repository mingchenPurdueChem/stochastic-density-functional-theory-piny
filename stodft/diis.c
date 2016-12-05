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
  
  int iGrid,iDiis;
  int rhoRealGridNum = stodftInfo->rhoRealGridNum;
  int numDiis = stodftInfo->numDiis;
  int numDiisNow = MIN(iScf,numDiis);
  int numStepMix = stodftInfo->numStepMix;
  int cpLsda = cpopts->cp_lsda;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;

  double mixRatio1 = stodftInfo->mixRatioBig;
  double mixRatio2 = 1.0-mixRatio1;

  double *rhoUpCorrect = stodftCoefPos->rhoUpCorrect;
  double *rhoDnCorrect = stodftCoefPos->rhoDnCorrect;
  double *rhoUp = cpscr->cpscr_rho.rho_up;
  double *rhoDn = cpscr->cpscr_rho.rho_dn;
  double *rhoUpOld = stodftCoefPos->rhoUpOld;
  double *rhoDnOld = stodftCoefPos->rhoDnOld;

  double **rhoUpBank = stodftCoefPos->rhoUpBank;
  double **rhoDnBank = stodftCoefPos->rhoDnBank;
  double **rhoUpErr  = stodftCoefPos->rhoUpErr;
  double **rhoDnErr  = stodftCoefPos->rhoDnErr;

  //updateBank(rhoUpCorrect,rhoUpBank,iScf);
  //updateErr(rhoUpBank,rhoUpErr,iScf);
  stodftInfo->numDiisNow = numDiisNow;

  if(iScf==0){//Initial Step
    updateBank(stodftInfo,stodftCoefPos,rhoUpCorrect,rhoUpBank);
    memcpy(rhoUpOld,rhoUpCorrect,rhoRealGridNum*sizeof(double));
  }
  else if(iScf<=numStepMix){//Mixing Step
    updateErr(stodftInfo,stodftCoefPos,rhoUpCorrect,rhoUpOld,rhoUpErr);
    updateBank(stodftInfo,stodftCoefPos,rhoUpCorrect,rhoUpBank);
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
      rhoUp[iGrid+1] = rhoUpCorrect[iGrid]*mixRatio2+rhoUpBank[0][iGrid]*mixRatio1;
    }
    //updateErr(stodftInfo,stodftCoefPos,&rhoUp[1],rhoUpErr,rhoUpBank);
    //updateBank(stodftInfo,stodftCoefPos,&rhoUp[1],rhoUpBank);
    memcpy(rhoUpOld,&rhoUp[1],rhoRealGridNum*sizeof(double));
  }
  else{//diis
    updateErr(stodftInfo,stodftCoefPos,rhoUpCorrect,rhoUpOld,rhoUpErr);
    updateBank(stodftInfo,stodftCoefPos,rhoUpCorrect,rhoUpBank);
    calcDensityDiis(cp,rhoUpBank,rhoUpErr);   
    //updateBank(stodftInfo,stodftCoefPos,&rhoUp[1],rhoUpBank);
    memcpy(rhoUpOld,&rhoUp[1],rhoRealGridNum*sizeof(double));
  }
  if(cpLsda==1&&numStateDnProc>0){
    //updateBank(rhoDnCorrect,rhoDnBank,iScf);
    //updateErr(rhoDnBank,rhoDnErr,iScf);

    if(iScf==1){//Initial Step
      updateBank(stodftInfo,stodftCoefPos,rhoDnCorrect,rhoDnBank);
      memcpy(rhoDnOld,rhoDnCorrect,rhoRealGridNum*sizeof(double));
    }
    else if(iScf<=numStepMix){//Mixing Step
      updateErr(stodftInfo,stodftCoefPos,rhoDnCorrect,rhoDnOld,rhoDnErr);
      updateBank(stodftInfo,stodftCoefPos,rhoDnCorrect,rhoDnBank);
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
	rhoDn[iGrid+1] = rhoDnCorrect[iGrid]*mixRatio2+rhoDnBank[0][iGrid]*mixRatio1;
      }
      //updateErr(stodftInfo,stodftCoefPos,&rhoDn[1],rhoDnErr,rhoDnBank);
      //updateBank(stodftInfo,stodftCoefPos,&rhoDn[1],rhoDnBank);
      memcpy(rhoDnOld,&rhoDn[1],rhoRealGridNum*sizeof(double));
    }
    else{//diis
      updateErr(stodftInfo,stodftCoefPos,rhoDnCorrect,rhoDnOld,rhoDnErr);
      updateBank(stodftInfo,stodftCoefPos,rhoDnCorrect,rhoDnBank);
      calcDensityDiis(cp,rhoDnBank,rhoDnErr);      
      //updateBank(stodftInfo,stodftCoefPos,&rhoDn[1],rhoDnBank);
      memcpy(rhoDnOld,&rhoDn[1],rhoRealGridNum*sizeof(double));
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
  
  if(rhoBank[numDiis]!=NULL)free(rhoBank[numDiis]);
  for(iDiis=numDiis-1;iDiis>-1;iDiis--)rhoBank[iDiis+1] = rhoBank[iDiis];
  rhoBank[0] = (double*)cmalloc(rhoRealGridNum*sizeof(double));
  memcpy(rhoBank[0],rho,rhoRealGridNum*sizeof(double));

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void updateErr(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos,
                double *rho,double *rhoOld,double **rhoErr)
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
  
  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoErr[0][iGrid] -= rhoOld[iGrid];
  
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
  int numDiisNow = stodftInfo->numDiisNow;
  int myidState = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int rhoRealGridNum = stodftInfo->rhoRealGridNum;
  int diisMatrixCalcFullFlag = stodftInfo->diisMatrixCalcFullFlag;
  int iDiis,jDiis,index1,index2;
  int iGrid;
  int incx = 1;
  int incy = 1;

  double dotProc;
  double dotTot;   
  double matrixElemMax = -1000.0;
  double matrixElemNorm;
  double *diisMatrix = (double*)cmalloc((numDiisNow+1)*(numDiisNow+1)*sizeof(double));
  double *svdLinSol = (double*)cmalloc((numDiisNow+1)*sizeof(double));
  double *b = (double*)cmalloc((numDiisNow+1)*sizeof(double));
  double *diisCoeff = stodftCoefPos->diisCoeff;
  double *rhoUp = cpscr->cpscr_rho.rho_up;
  double *rhoDn = cpscr->cpscr_rho.rho_dn;
  double *rhoUpCorrect = stodftCoefPos->rhoUpCorrect;
  double *rhoDnCorrect = stodftCoefPos->rhoDnCorrect;

  MPI_Comm commStates = communicate->comm_states;

/*==========================================================================*/
/* I) For the first time, calculate the whole diisMatrix		    */

  // Calculate the dot prod part
  if(myidState==0){
    printf("Mixing: DIIS Step\n");
  }

  for(iDiis=0;iDiis<numDiisNow;iDiis++){
    for(jDiis=iDiis;jDiis<numDiisNow;jDiis++){
      dotProc = DDOT(&rhoRealGridNum,rhoErr[iDiis],&incx,rhoErr[jDiis],&incy);
      if(numProcStates>1)Reduce(&dotProc,&dotTot,1,MPI_DOUBLE,MPI_SUM,0,commStates);
      else dotTot = dotProc;
      if(myidState==0){
	index1 = iDiis*(numDiisNow+1)+jDiis;
	index2 = jDiis*(numDiisNow+1)+iDiis;
	diisMatrix[index1] = dotTot;
	diisMatrix[index2] = dotTot;
      }//endif myidState
    }//endfor jDiis
  }//endfor iDiis

  if(myidState==0){
    //Get the largest fabs(matrix element)
    for(iDiis=0;iDiis<numDiisNow;iDiis++){
      for(jDiis=iDiis;jDiis<numDiisNow;jDiis++){
	matrixElemNorm = fabs(diisMatrix[iDiis*(numDiisNow+1)+jDiis]);
	if(matrixElemNorm>matrixElemMax){
	  matrixElemMax = matrixElemNorm;
	}
      }
    }
    

    // Calculate the Lagrange Multiplier part
    for(iDiis=0;iDiis<numDiisNow;iDiis++){
      diisMatrix[numDiisNow*(numDiisNow+1)+iDiis] = -matrixElemMax;
      diisMatrix[iDiis*(numDiisNow+1)+numDiisNow] = -matrixElemMax;
    }
    diisMatrix[numDiisNow*(numDiisNow+1)+numDiisNow] = 0.0;

    //debug
    for(iDiis=0;iDiis<numDiisNow+1;iDiis++){
      for(jDiis=0;jDiis<numDiisNow+1;jDiis++){
	printf("%.16lg,",diisMatrix[iDiis*(numDiisNow+1)+jDiis]);
      }
      printf("\n");
    }
  // reset the flag so that we don't do it again
/*==========================================================================*/
/* II) Safely reverse the diisMatrix with SVD				    */

    for(iDiis=0;iDiis<numDiisNow;iDiis++)b[iDiis] = 0.0;
    b[numDiisNow] = -matrixElemMax;

    matrixInvSVD(diisMatrix,b,svdLinSol,numDiisNow+1);       

    for(iDiis=0;iDiis<numDiisNow+1;iDiis++){
      printf("%.16lg ",svdLinSol[iDiis]);
    }
    printf("\n");


  }//endif myidState


  if(numProcStates>1){
    Bcast(svdLinSol,numDiisNow+1,MPI_DOUBLE,0,commStates);
  }
  
  for(iDiis=0;iDiis<numDiisNow;iDiis++){
    diisCoeff[iDiis] = svdLinSol[iDiis];    
  }
  stodftInfo->lambdaDiis = svdLinSol[numDiisNow];

/*==========================================================================*/
/* III) Linear Combinination of all densities                               */
  
  // I may change this 
  /*
  //mix density then push stack
  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
    rhoUp[iGrid+1] = diisCoeff[0]*rhoUpCorrect[iGrid];
    for(iDiis=1;iDiis<numDiisNow;iDiis++){
      rhoUp[iGrid+1] += diisCoeff[iDiis]*rhoBank[iDiis-1][iGrid];
    }
  }
  */
  
  // push stack then mix density
  int iScf = stodftInfo->iScf;
  double mixRatioSM = 0.1; //alpha
  double pre = mixRatioSM-1.0;
  if(iScf<=5){
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
      rhoUp[iGrid+1] = 0.0;
      for(iDiis=0;iDiis<numDiisNow;iDiis++){
	rhoUp[iGrid+1] += diisCoeff[iDiis]*(rhoBank[iDiis][iGrid]+
			    pre*rhoErr[iDiis][iGrid]);
      }//endfor iDiis
    }//endfor iGrid
  }//endif iScf
  else{
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
      rhoUp[iGrid+1] = 0.0;
      for(iDiis=0;iDiis<numDiisNow;iDiis++)rhoUp[iGrid+1] += diisCoeff[iDiis]*rhoBank[iDiis][iGrid];
    }
  }

  free(svdLinSol);
  free(diisMatrix);
  free(b);
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
  double smax = -1.0e30;
  double smin = 1.0e30;
  double sratio;

  double *A = (double*)cmalloc(ndim*ndim*sizeof(double));
  double *s = (double*)cmalloc(ndim*sizeof(double));
  double *sinv = (double*)cmalloc(ndim*sizeof(double));
  double *u = (double*)cmalloc(ndim*ndim*sizeof(double));
  double *vt = (double*)cmalloc(ndim*ndim*sizeof(double));
  double *work = (double*)cmalloc(lwork*sizeof(double));
  double *u_b = (double*)cmalloc(ndim*sizeof(double));

  
  memcpy(A,mat,ndim*ndim*sizeof(double));

  DGESVD(&jobu,&jobvt,&m,&n,A,&lda,s,u,&ldu,vt,&ldvt,work,&lwork,&info);

  for(i=0;i<ndim;i++){
    if(s[i]>smax)smax = s[i];
    if(s[i]<smin)smin = s[i];
  }
  if(smin>0.0)sratio = smax/smin;
  else sratio = 1.0e30;
  printf("Singular Value range is from %.6lg to %.6lg\n",smin,smax);
  for(i=0;i<ndim;i++){
    if(s[i]/smax>1.0e-10)sinv[i] = 1.0/s[i];
    else sinv[i] = 0.0;
  }
  
  if(info==0){
    trans = 't';
    alpha = 1.0;
    beta = 0.0;
    incx = 1;
    incy = 1;
    DGEMV(&trans,&ndim,&ndim,&alpha,u,&lda,b,&incx,&beta,u_b,&incy);
    
    //for(i=0;i<ndim;i++)u_b[i] /= s[i];
    for(i=0;i<ndim;i++)u_b[i] *= sinv[i];
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
  free(sinv);
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

