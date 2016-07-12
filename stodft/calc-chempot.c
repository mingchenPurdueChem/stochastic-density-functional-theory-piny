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
  STODFTINFO *stodftInfo = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos = cp->stodftCoefPos;
  CPOPTS *cpOpts = &(cp->cpopts);
  COMMUNICATE *commCP = &(cp->communicate);

  int iChem,iGrid;
  int chemPotIndex;
  int numChemPot = stodftInfo->numChemPot;
  int cpLsda = cpopts->cp_lsda;
  int numFFTProc = cp_para_fft_pkg3d_lg->nfft_proc;
  int numFFT2Proc = numFFTProc/2; 
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
  double *numElectron = stodftCoefPos->numElectron;
  double *chemPot = stodftCoefPos->chemPot;
  double *interpCoef = (double*)cmalloc(numChemPot*sizeof(double)); //Interpolation Coeffcients
  double *rhoTemp = (double*)cmalloc(numChemPot*numFFT2Proc);
  double **rhoUp = stodftCoefPos->rhoUp;
  double **rhoDn = stodftCoefPos->rhoDn;
  double *rhoUpCorrect = stodftCoefPos->rhoUpCorrect;
  double *rhoDnCorrect = stodftCoefPos->rhoDnCorrect;

  if(myidState==0){
    printf("==============================================\n");
    printf("Start interpolating number of electrons.\n");
    chemPotTrue = solveLagrangePolyInterp(numChemPot,chemPot,numElectron,numElecTrue,interpCoef);
  }
  Bcast(&chemPotTrue,1,MPI_DOUBLE,0,commCP);
  Bcast(interpCoef,numChemPot,MPI_DOUBLE,0,commCP);

  for(iChem=0;iChem<numChemProc;iChem++){
    chemPotIndex = chemProcIndexInv[iChem];
    Scatterv(rhoUp[iChem],rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	    &rhoTemp[chemPotIndex*numFFT2Proc],numFFT2Proc,MPI_DOUBLE,myidState,comm_states);
    
  }
  for(iGrid=0;iGrid<numFFT2Proc;iGrid++){
    rhoUpCorrect[iGrid] = 0.0;
    for(iChem=0;iChem<numChemPot;iChem++){
      rhoUpCorrect[iGrid] += interpCoef[iChem]*rhoUp[iChem*numFFT2Proc+iGrid];
    }//endfor iChem
  }//endfor iGrid  
  if(cpLsda==1){
    for(iChem=0;iChem<numChemProc;iChem++){
      chemPotIndex = chemProcIndexInv[iChem];
      Scatterv(rhoDn[iChem],rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	      &rhoTemp[chemPotIndex*numFFT2Proc],numFFT2Proc,MPI_DOUBLE,myidState,comm_states);

    }
    for(iGrid=0;iGrid<numFFT2Proc;iGrid++){
      rhoDnCorrect[iGrid] = 0.0;
      for(iChem=0;iChem<numChemPot;iChem++){
	rhoDnCorrect[iGrid] += interpCoef[iChem]*rhoDn[iChem*numFFT2Proc+iGrid];
      }//endfor iChem
    }//endfor iGrid  
  }

  
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
  double gapInit = stodftInfo->gapInput;
  double factor = M_PI*0.5/numChemPot;
  double *chemPot = stodftCoefPos->chemPot;
  double *chebyNode = (double*)cmalloc(numChemPot*sizeof(double));
  
  // Generate Chebyshev nodes
  for(iNode=0;iNode<numChemPot;iNode++){
    chebyNode[iNode] = cos(2.0*iNode+1*factor);
  }
  // Scale the nodes to correct chem pot and gap
  
  for(iNode=0;iNode<numChemPot;iNode++){
    chemPot[iNode] = chebyNode[iNode]*0.5*gapInit+chemPotInit;
  }
  
  free(chebyNode);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double solveLagrangePolyInterp(int numSamp,double *x, double *y,double target,
				double *interpCoef)
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
  double ymax = y[iSamp-1];
  double xopt;
  double tol = 1.0e-7;
  double tolnow = 1.0;
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

  tolnow = calcLagrangeInterpFun(numSamp,xopt,x,y,target,lagFunValue);
  deriv  = calcLagrangeInterpDrv(numSamp,xopt,x,y,target,lagFunValue);

  while(torlnow>torl){
    xopt -= tolnow/deriv;
    tolnow = calcLagrangeInterpFun(numSamp,xopt,x,y,target,interpCoef);
    deriv = calcLagrangeInterpDrv(numSamp,xopt,x,y,target,interpCoef);
  }

  return xopt;
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double calcLagrangeInterpFun(int numSamp,double xopt,double *x,double *y,
			     double target,double *lagFunValue)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This function calculate Lagrange interpolation value at xopt          */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int iTerm,jTerm;
  
  double sum = 0.0;
  double prod,prod2;

  for(iTerm=0;iTerm<numSamp;iTerm++){
    prod = 1.0;
    prod2 = 1.0;
    for(jTerm=0;jTerm<numSamp;jTerm++){
      if(jTerm!=iTerm){
	prod *= xopt-x[jTerm];
	prod2 *= x[iTerm]-x[jTerm];
      }
    }
    lagFunValue[iTerm] = prod/prod2;
    sum += y[iTerm]*lagFunValue[iTerm];
  }
  return sum;

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double calcLagrangeInterpDrv(int numSamp,double xopt,double *x,double *y,
                             double target,double *lagFunValue)
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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genInterpCoef(double xtrue,int numSamp,double x, double y,double *interpCoef)
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
  double prod,prod2;
  
  for(iTerm=0;iTerm<numSamp;iTerm++){
    prod = 1.0;
    prod2 = 1.0;
    for(jTerm=0;jTerm<numSamp;jTerm++){
      if(jTerm!=iTerm){
	prod *= xtrue-x[jTerm];
	prod2 *= x[iTerm]-x[jTerm];
      }
    }
    interpCoef = prod/prod2;
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/




