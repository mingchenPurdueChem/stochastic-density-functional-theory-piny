/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: calc-chempot.c                                 */
/*                                                                          */
/*  This routine calculate the correct chemical potential, either from	    */
/*  interpolation or from half polynomial expension(Chebshev only)	    */
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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcChemPotInterp(CP *cp)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* I decide to use Lagrange polynormial for interpolation b/c it is easy */
/* to get the linear coeffcients for density. The Lagrange polynormial   */
/* interpolation shares the formular : f(x)=\sum a_i l_i(x). and l_i(x)  */
/* =\Pi_{j\neq i}(x-x_j)/(x_i-x_j). We need to evaluate l_i' for         */
/* optimization. l_i'=l_i*g_i where g_i=(\sum_{k\neq i}1/(x-x_k)).	 */
/* First we will interpolate Ne(mu). Then we will solve Ne(mu*)=tot # of */
/* electrons. Then we shall get rhoUp(mu*) and rhoDn(mu*) if necessary.  */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  STODFTINFO *stodftInfo = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos = cp->stodftCoefPos;
  CPOPTS *cpopts = &(cp->cpopts);
  CPSCR *cpscr = &(cp->cpscr);  
  COMMUNICATE *commCP = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);

  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

  int i,j,k;
  int iChem,iGrid;
  int chemPotIndex;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numChemPot = stodftInfo->numChemPot;
  int numChemProc = stodftInfo->numChemProc;
  int cpLsda = cpopts->cp_lsda;
  int cpParaOpt      = cpopts->cp_para_opt;
  int rhoRealGridNum = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot = stodftInfo->rhoRealGridTot;
  int numFFTProc = cp_para_fft_pkg3d_lg->nfft_proc;
  int numFFT2Proc = numFFTProc/2; 
  int numFFT = cp_para_fft_pkg3d_lg->nfft;
  int numFFT2 = numFFT/2;
  int myidState = commCP->myid_state;
  int numProcStates = commCP->np_states;
  MPI_Comm comm_states = commCP->comm_states;
  
  int *densityMap = stodftInfo->densityMap;
  int *indexChemProc = stodftInfo->indexChemProc;
  int *chemProcIndexInv = stodftInfo->chemProcIndexInv;
  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;
  
  double numElecTrue = stodftInfo->numElecTrue;
  double chemPotTrue;
  double numElecDiff;
  double *numElectron = stodftCoefPos->numElectron;
  double *chemPot = stodftCoefPos->chemPot;
  double *interpCoef = (double*)cmalloc(numChemPot*sizeof(double)); //Interpolation Coeffcients
  double *rhoTemp = (double*)cmalloc(numChemPot*rhoRealGridNum*sizeof(double));
  double **rhoUpChemPot = stodftCoefPos->rhoUpChemPot;
  double **rhoDnChemPot = stodftCoefPos->rhoDnChemPot;
  double *rhoUpCorrect = stodftCoefPos->rhoUpCorrect;
  double *rhoDnCorrect = stodftCoefPos->rhoDnCorrect;
  double *rhoUp = cpscr->cpscr_rho.rho_up;
  double *rhoDn = cpscr->cpscr_rho.rho_dn;

  if(myidState==0){
    chemPotTrue = solveLagrangePolyInterp(numChemPot,chemPot,numElectron,numElecTrue,interpCoef,&numElecDiff);
    printf("Chem Pot %.6lg DNe %.6lg\n",chemPotTrue,numElecDiff);
  }
  if(numProcStates>1){
    Bcast(&chemPotTrue,1,MPI_DOUBLE,0,comm_states);
    Bcast(interpCoef,numChemPot,MPI_DOUBLE,0,comm_states);
  }
  stodftInfo->chemPotTrue = chemPotTrue;
  printf("%p %p\n",rhoRealSendCounts,rhoRealDispls);

  if(myidState==0){
    for(i=0;i<4;i++)printf("rhoRealSendCounts %i rhoRealDispls %i\n",rhoRealSendCounts[i],rhoRealDispls[i]);
  }
  
  if(cpParaOpt==0){
    for(iChem=0;iChem<numChemProc;iChem++){
      chemPotIndex = chemProcIndexInv[iChem];
      if(numProcStates>1){
	Scatterv(rhoUpChemPot[iChem],rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
		&rhoTemp[chemPotIndex*rhoRealGridNum],rhoRealGridNum,MPI_DOUBLE,myidState,comm_states);
      }
      else{
	memcpy(&rhoTemp[chemPotIndex*rhoRealGridNum],rhoUpChemPot[iChem],rhoRealGridNum*sizeof(double));
	//memcpy(rhoUpChemPot[iChem],&rhoTemp[chemPotIndex*rhoRealGridNum],rhoRealGridNum*sizeof(double));
	//for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoTemp[chemPotIndex*rhoRealGridNum+iGrid] = rhoUpChemPot[iChem][iGrid];
      }//endif       
    }//endfor iChem
  }

  printf("coef %lg\n",interpCoef[0]);

  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
    rhoUpCorrect[iGrid] = 0.0;
    for(iChem=0;iChem<numChemPot;iChem++){
      rhoUpCorrect[iGrid] += interpCoef[iChem]*rhoTemp[iChem*rhoRealGridNum+iGrid];
    }//endfor iChem
    //printf("iGrid %i rhoCorrect %lg\n",iGrid,rhoUpCorrect[iGrid]*0.0009250463018013585);
  }//endfor iGrid  
  if(cpLsda==1&&numStateDnProc>0){
    if(cpParaOpt==0){
      for(iChem=0;iChem<numChemProc;iChem++){
	chemPotIndex = chemProcIndexInv[iChem];
	if(numProcStates>1){
	  Scatterv(rhoDnChemPot[iChem],rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
		  &rhoTemp[chemPotIndex*rhoRealGridNum],rhoRealGridNum,MPI_DOUBLE,myidState,comm_states);
	}
	else{
	  memcpy(&rhoTemp[chemPotIndex*rhoRealGridNum],rhoDnChemPot[iChem],rhoRealGridNum*sizeof(double));
	}//endif numProcStates
      }//endfor iChem
    }//endif cpParaOpt
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
      rhoDnCorrect[iGrid] = 0.0;
      for(iChem=0;iChem<numChemPot;iChem++){
	rhoDnCorrect[iGrid] += interpCoef[iChem]*rhoTemp[iChem*rhoRealGridNum+iGrid];
      }//endfor iChem
    }//endfor iGrid  
  }

  //debug
  /*
  double testNumElec = 0.0;
  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
    testNumElec += rhoUpCorrect[iGrid];
  }
  testNumElec /= rhoRealGridTot;
  printf("tot number of electron after interp %.16lg\n",testNumElec);
  */

  
  free(interpCoef);
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genChemPotInterpPoints(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This function generate initial chemical potential The initial chem	 */
/* potential and the band gap is provided by the user. This is a easy    */
/* version with static chemical potential and static gap (I may need to  */
/* up date this in the future for dynamic chem pot and gap.)             */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int numChemPot = stodftInfo->numChemPot;
  int iNode;
  double chemPotInit = stodftInfo->chemPotInit;
  double gapInit = stodftInfo->gapInit;
  double factor = M_PI*0.5/numChemPot;
  double *chemPot = stodftCoefPos->chemPot;
  double *chebyNode = (double*)cmalloc(numChemPot*sizeof(double));
  
  // Generate Chebyshev nodes
  for(iNode=0;iNode<numChemPot;iNode++){
    chebyNode[numChemPot-iNode-1] = cos((2.0*iNode+1.0)*factor);
  }
  // Scale the nodes to correct chem pot and gap
  
  for(iNode=0;iNode<numChemPot;iNode++){
    chemPot[iNode] = chebyNode[iNode]*0.5*gapInit+chemPotInit;
    //printf("iNode %i chemPot %lg\n",iNode,chemPot[iNode]);
  }
  
  free(chebyNode);
  //exit(0);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double solveLagrangePolyInterp(int numSamp,double *x, double *y,double target,
				double *interpCoef,double *numElecDiff)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This function solve \sum a_i l_i(x)=target, given {x} and {y} as	 */
/* interpolating points. x is in ascent order and y[x] should be 	 */
/* monotonic. The interpolation coefficients are actually Lagrange	 */
/* polynormial values.							 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int iSamp;

  double ymin = y[0];
  double ymax = y[numSamp-1];
  double xopt,x1,x2,x1old,x2old,y1,y2,y1old,y2old;
  double tol = 1.0e-7;
  double tolnow = 1.0;
  double value;
  double deriv;

/*=======================================================================*/
/* I. Find ymin and ymax, make sure ymin<target<ymax			 */
  
  if(target<=ymin||target>=ymax){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("You are choosing chemical potential range that is\n");
    printf("far from true number of electron. The interpolation\n");
    printf("could be very unstable. I'm gonna to kill the process.\n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(0);
  }
/*=======================================================================*/
/* II. Get the initial guess of chemical potential.                      */

  iSamp = 0;
  while(y[iSamp]<target)iSamp += 1;
  xopt = x[iSamp-1]+(target-y[iSamp-1])*(x[iSamp]-x[iSamp-1])/(y[iSamp]-y[iSamp-1]);

  value = calcLagrangeInterpFun(numSamp,xopt,x,y,interpCoef);
  //deriv  = calcLagrangeInterpDrv(numSamp,xopt,x,y,interpCoef);
  tolnow = fabs(value-target);
  //printf("initial %.6lg\n",xopt);
  if(value>target){
    x1 = x[iSamp-1];
    x2 = xopt;
    y1 = y[iSamp-1];
    y2 = value;
  }
  else{
    x1 = xopt;
    x2 = x[iSamp];
    y1 = value;
    y2 = y[iSamp];
  }

  while(tolnow>tol){
    x1old = x1;
    x2old = x2;
    y1old = y1;
    y2old = y2;
    xopt = x1+(target-y1)*(x2-x1)/(y2-y1);
    value = calcLagrangeInterpFun(numSamp,xopt,x,y,interpCoef);
    if(value<target){
      x1 = x1old;
      x2 = xopt;
      y1 = y1old;
      y2 = value;
    }
    else{
      x1 = xopt;
      x2 = x2old;
      y1 = value;
      y2 = y2old;
    }
    tolnow = fabs(value-target);
    /*
    xopt -= value/deriv;
    value = calcLagrangeInterpFun(numSamp,xopt,x,y,interpCoef)-target;
    deriv = calcLagrangeInterpDrv(numSamp,xopt,x,y,interpCoef);
    tolnow = fabs(value);
    */
    //printf("opt-tol %.6lg %.6lg %.6lg %.6lg %.6lg\n",x1,x2,xopt,value,tolnow);
  }
  *numElecDiff = tolnow;
  
  return xopt;
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double calcLagrangeInterpFun(int numSamp,double xopt,double *x,double *y,
			     double *lagFunValue)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This function calculate Lagrange interpolation value at xopt          */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int iTerm,jTerm;

  int sign1,sign2;
  int sampflag = 0;
  int sampIndex;
  double sum = 0.0;
  double prod,prod2;
  double diff1,diff2;

  for(iTerm=0;iTerm<numSamp;iTerm++){
    diff1 = fabs(xopt-x[iTerm]);
    if(diff1<1.0e-16){
      sampflag += 1;
      sampIndex = iTerm;
    }
  }

  if(sampflag>0)sum = y[sampIndex];
  else{
    for(iTerm=0;iTerm<numSamp;iTerm++){
      prod = 0.0;
      prod2 = 0.0;
      sign1 = 1;
      sign2 = 1;
      for(jTerm=0;jTerm<numSamp;jTerm++){
        if(jTerm!=iTerm){
          diff1 = xopt-x[jTerm];
          diff2 = x[iTerm]-x[jTerm];
          if(diff1<0)sign1 *= -1;
          if(diff2<0)sign2 *= -1;
          prod += 0.5*log(diff1*diff1);
          prod2 += 0.5*log(diff2*diff2);
          //prod *= xopt-x[jTerm];
          //prod2 *= x[iTerm]-x[jTerm];
        }
      }
      if(fabs(prod2)<1e-50)printf("iTerm %i prod2 %lg\n",iTerm,prod2);
      //lagFunValue[iTerm] = prod/prod2;
      lagFunValue[iTerm] = sign1*sign2*exp(prod-prod2);
      sum += y[iTerm]*lagFunValue[iTerm];
    }
  }
  return sum;

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double calcLagrangeInterpDrv(int numSamp,double xopt,double *x,double *y,
                             double *lagFunValue)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This function solve \sum a_i l_i(x)=target, given {x} and {y} as      */
/* interpolating points. x is in ascent order and y[x] should be         */
/* monotonic.                                                            */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int iTerm,jTerm;
  double sum;
  double drv = 0.0;
  
  for(iTerm=0;iTerm<numSamp;iTerm++){
    sum = 0.0;
    for(jTerm=0;jTerm<numSamp;jTerm++){
      if(jTerm!=iTerm)sum += 1.0/(xopt-x[jTerm]);
    }
    drv += sum*lagFunValue[iTerm]*y[iTerm];
  }
  return drv;

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/




