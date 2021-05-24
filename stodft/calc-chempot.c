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

  int i,j,k;
  int iChem,iGrid,iProc;
  int chemPotIndex,chemPotIndexProc;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numChemPot = stodftInfo->numChemPot;
  int numChemProc = stodftInfo->numChemProc;
  int cpLsda = cpopts->cp_lsda;
  int cpParaOpt      = cpopts->cp_para_opt;
  int rhoRealGridNum = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot = stodftInfo->rhoRealGridTot;
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
    printf("Chem Pot %.16lg DNe %.16lg\n",chemPotTrue,numElecDiff);
  }
  if(numProcStates>1){
    Bcast(&chemPotTrue,1,MPI_DOUBLE,0,comm_states);
    Bcast(interpCoef,numChemPot,MPI_DOUBLE,0,comm_states);
  }
  stodftInfo->chemPotTrue = chemPotTrue;

  /*
  if(myidState==3){
    for(i=0;i<4;i++)printf("rhoRealSendCounts %i rhoRealDispls %i\n",rhoRealSendCounts[i],rhoRealDispls[i]);
  }
  */
  if(cpParaOpt==0){
    for(iChem=0;iChem<numChemPot;iChem++){
      if(numProcStates>1){
        Scatterv(rhoUpChemPot[indexChemProc[iChem]],rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
                 &rhoTemp[iChem*rhoRealGridNum],rhoRealGridNum,MPI_DOUBLE,densityMap[iChem],comm_states);  
      }
      else{
	memcpy(&rhoTemp[iChem*rhoRealGridNum],rhoUpChemPot[iChem],rhoRealGridNum*sizeof(double));
      }
    }
    /*
    for(iChem=0;iChem<numChemProc;iChem++){
      chemPotIndex = chemProcIndexInv[iChem];
      //printf("iChem %i chemPotIndex %i\n",iChem,chemPotIndex);
      if(numProcStates>1){
	for(iProc=0;iProc<numProcStates;iProc++){
	  if(myidState==iProc)chemPotIndexProc = chemPotIndex;
	  Bcast(&chemPotIndexProc,1,MPI_INT,iProc,comm_states);
	  printf("myid %i iChem %i iProc %i chemPotIndexProc %i\n",myidState,iChem,iProc,chemPotIndexProc);
	  Scatterv(rhoUpChemPot[iChem],rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
		   &rhoTemp[chemPotIndexProc*rhoRealGridNum],rhoRealGridNum,MPI_DOUBLE,iProc,comm_states);
	}
      }
      else{
	memcpy(&rhoTemp[chemPotIndex*rhoRealGridNum],rhoUpChemPot[iChem],rhoRealGridNum*sizeof(double));
	//memcpy(rhoUpChemPot[iChem],&rhoTemp[chemPotIndex*rhoRealGridNum],rhoRealGridNum*sizeof(double));
	//for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoTemp[chemPotIndex*rhoRealGridNum+iGrid] = rhoUpChemPot[iChem][iGrid];
      }//endif       
    }//endfor iChem
    */
  }

  //printf("coef %lg\n",interpCoef[0]);

  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
    rhoUpCorrect[iGrid] = 0.0;
    for(iChem=0;iChem<numChemPot;iChem++){
      rhoUpCorrect[iGrid] += interpCoef[iChem]*rhoTemp[iChem*rhoRealGridNum+iGrid];
    }//endfor iChem
    if(rhoUpCorrect[iGrid]<0.0){
      printf("Negative Density Value!\n");
    }
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
  int filterDiagFlag = stodftInfo->filterDiagFlag;
  int energyWindowOn = stodftInfo->energyWindowOn;
  int fragWindowFlag = stodftInfo->fragWindowFlag;
  int iScf = stodftInfo->iScf;
  double chemPotInit = stodftInfo->chemPotInit;
  double gapInit = stodftInfo->gapInit;
  double factor = M_PI*0.5/numChemPot;
  double *chemPot = stodftCoefPos->chemPot;
  double *chemPotBackUp = stodftCoefPos->chemPotBackUp;
  double *chebyNode = (double*)cmalloc(numChemPot*sizeof(double));
  double energyMin = chemPotInit-gapInit*0.5;
  double dE = gapInit/numChemPot;
  double chemPotTrue;
  FILE *fileTest;
  
  // Generate Chebyshev nodes
  if(filterDiagFlag==0){
    if(energyWindowOn==0){
      for(iNode=0;iNode<numChemPot;iNode++){
	chebyNode[numChemPot-iNode-1] = cos((2.0*iNode+1.0)*factor);
      }
      // Scale the nodes to correct chem pot and gap
      for(iNode=0;iNode<numChemPot;iNode++){
	chemPot[iNode] = chebyNode[iNode]*0.5*gapInit+chemPotInit;
	//printf("iNode %i chemPot %lg\n",iNode,chemPot[iNode]);
      }
      chemPot[0] = chemPotInit-gapInit*0.5;
      chemPot[1] = chemPotInit+gapInit*0.5;
    }
    else{
      if(fragWindowFlag==0){ //regular energy window
        if(iScf<1000000){
          energyMin = stodftInfo->energyMin;
          chemPotTrue = stodftInfo->chemPotTrue;
          dE = (chemPotTrue-energyMin)/numChemPot; 
          for(iNode=0;iNode<numChemPot;iNode++){
            chemPot[iNode] = energyMin+(iNode+1)*dE;
            //printf("0000000 iNode %i energyMin %lg chemPotTrue %lg dE %lg chemPot %lg\n",
            //        iNode,energyMin,chemPotTrue,dE,chemPot[iNode]);
          }
          // double check the last chemical potential is correct...
          // TEST READ chempot
          fileTest = fopen("chempot","r");
          for(iNode=0;iNode<numChemPot;iNode++)fscanf(fileTest,"%lg",&chemPot[iNode]);
          fclose(fileTest);
          chemPot[numChemPot-1] = chemPotTrue;
          for(iNode=0;iNode<numChemPot;iNode++)chemPotBackUp[iNode] = chemPot[iNode];
        }
        else{
          // I have problems in scf convergence, let's test fixing other energy window locations
          for(iNode=0;iNode<numChemPot;iNode++)chemPot[iNode] = chemPotBackUp[iNode];
          chemPot[numChemPot-1] = chemPotTrue;
        }
      }
      else{ 
        // With fragment density, the chemical potential does not dependent on 
        // the CURRENT chemical potential. In the first step, chemPot[numChemPot-2] 
        // is the chemPotInit+0.5*gapInit. In the first 2 to N steps, chemPot[numChemPot-2]
        // is the chemical potential in the previous step. After N steps, chemPot[numChemPot-2] 
        // is fixed. 
        if(iScf<10){
          energyMin = stodftInfo->energyMin;
          // The first SCF step: use chemPotInit
          // Otherwise: Use the chemical potential from the 
          // previous SCF step
          if(iScf==1)chemPotTrue = chemPotInit+0.5*gapInit;
          else chemPotTrue = stodftInfo->chemPotTrue;
          // DEBUG
          //chemPotTrue = chemPotInit+0.5*gapInit;
          dE = (chemPotTrue-energyMin)/(numChemPot-1.0); 
          for(iNode=0;iNode<numChemPot-1;iNode++){
            chemPot[iNode] = energyMin+(iNode+1)*dE;
          }
          // DEBUG
          fileTest = fopen("chempot","r");
          for(iNode=0;iNode<numChemPot-1;iNode++)fscanf(fileTest,"%lg",&chemPot[iNode]);
          fclose(fileTest);
          chemPot[numChemPot-1] = 100000.0;
          for(iNode=0;iNode<numChemPot;iNode++)chemPotBackUp[iNode] = chemPot[iNode];
          //for(iNode=0;iNode<numChemPot;iNode++)printf("ccccccccchempot %i %lg\n",iNode,chemPot[iNode]);
        }
        else{
          for(iNode=0;iNode<numChemPot;iNode++)chemPot[iNode] = chemPotBackUp[iNode];
        }
      }//endif fragWindowFlag
    }//endif energyWindowOn    
  }//endif filterDiagFlag
  else{
    for(iNode=0;iNode<numChemPot;iNode++){
      chemPot[iNode] = energyMin+(iNode+1)*dE;
    }
  }

  //debug
  /*
  double deltChemPot = gapInit/numChemPot;
  double chemPotMin = chemPotInit-0.5*gapInit;
  for(iNode=0;iNode<numChemPot;iNode++){
    chemPot[iNode] = chemPotMin+(iNode+0.5)*deltChemPot;
  }
  
  free(chebyNode);
  */
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

  //debug
  /*
  double dx = (x[numSamp-1]-x[0])/1000.0;
  int i;
  double xi,test;
  for(i=0;i<1000;i++){
    xi = x[0]+(i+0.5)*dx;
    test = calcLagrangeInterpFun(numSamp,xi,x,y,interpCoef);
    printf("testinterp %.16lg %.16lg\n",xi,test);
  }
  for(iSamp=0;iSamp<numSamp;iSamp++){
    printf("testinterp %.16lg %.16lg\n",x[iSamp],y[iSamp]);
  }
  */
  
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
    //printf("x1 %.10lg x2 %.10lg y1 %.10lg y2 %.10lg xopt %.10lg\n",x1,x2,y1,y2,xopt);
    value = calcLagrangeInterpFun(numSamp,xopt,x,y,interpCoef);
    //printf("value %.10lg\n",value);
    if(value<target){
      x1 = xopt;
      x2 = x2old;
      y1 = value;
      y2 = y2old;
    }
    else{
      x1 = x1old;
      x2 = xopt;
      y1 = y1old;
      y2 = value;
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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void updateChemPot(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos)
/*========================================================================*/
{/*begin routine*/
/**************************************************************************/
/* We may want to dynamically change the chemical ponteital interpolation */
/* range. As I have tested, more dense the interpolation points, better   */
/* convergence we will have. We can not affort very dense interpolation   */
/* points covering from <HOMO to >LUMO. Therefore, we shall dynamically   */
/* shrink the interpolation range. Considering the chemical potential can */
/* jump between HOMO and LUMO at the beginning of a SCF. The interp range */
/* shall cover the whole gap until the chemical potential pinning around  */
/* HOMO or LUMO for 5 steps. Then we use the last step chemical potential */
/* u_n as the center and 2*max_{i<=5}{|u_{n-i}-u_n|} as range.            */
/**************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
/*------------------------------------------------------------------------*/

  int iScf = stodftInfo->iScf;
  int iPot,jPot;
  int numChemPot = stodftInfo->numChemPot;
  int iNode;

  double factor = M_PI*0.5/numChemPot;
  double chemPotTrue = stodftInfo->chemPotTrue;
  double diff;
  double rangeMax = -100000.0;
  double *chemPotHistory = stodftInfo->chemPotHistory;
  double *chemPot = stodftCoefPos->chemPot;
  double *chebyNode = (double*)cmalloc(numChemPot*sizeof(double));

/*======================================================================*/
/* I) Record the true chemical potential in this step                   */

  chemPotHistory[iScf-1] = chemPotTrue;

/*======================================================================*/
/* II) Record the true chemical potential in this step                  */

  if(iScf>=5){
    for(iPot=iScf-5;iPot<iScf;iPot++){
      for(jPot=iPot+1;jPot<iScf;jPot++){
        diff = fabs(chemPotHistory[iPot]-chemPotHistory[jPot]);
        if(diff>rangeMax)rangeMax = diff;
      }//endfor jPot
    }//endfor iPot

    // Generate Chebyshev nodes
    for(iNode=0;iNode<numChemPot;iNode++){
      chebyNode[numChemPot-iNode-1] = cos((2.0*iNode+1.0)*factor);
    }
    // Scale the nodes to chemPotTrue+-rangeMax

    for(iNode=0;iNode<numChemPot;iNode++){
      chemPot[iNode] = chebyNode[iNode]*rangeMax+chemPotTrue;
      //printf("iNode %i chemPot %lg\n",iNode,chemPot[iNode]);
    }
  }//endif iScf

  free(chebyNode);
/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void adjChemPot(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos)
/*========================================================================*/
{/*begin routine*/
/**************************************************************************/
/* Double the chem pot range if it is too small				  */
/**************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
/*------------------------------------------------------------------------*/

  int numChemPot = stodftInfo->numChemPot;
  int iNode;

  double factor = M_PI*0.5/numChemPot;
  double numElecTrue = stodftInfo->numElecTrue;

  double range,center;
  double *chemPot = stodftCoefPos->chemPot;
  double *chebyNode = (double*)cmalloc(numChemPot*sizeof(double));
  double *numElectron = stodftCoefPos->numElectron;


/*======================================================================*/
/* I) Get the center and range			                        */

  center = 0.5*(chemPot[0]+chemPot[numChemPot-1]);
  range = chemPot[numChemPot-1]-chemPot[0];


/*======================================================================*/
/* II) Shift the center						        */

  if(numElectron[0]>numElecTrue)center -= range;
  else if(numElectron[numChemPot-1]<numElecTrue)center += range;

/*======================================================================*/
/* II) Record the true chemical potential in this step                  */

  // Generate Chebyshev nodes
  for(iNode=0;iNode<numChemPot;iNode++){
    chebyNode[numChemPot-iNode-1] = cos((2.0*iNode+1.0)*factor);
  }
  // Scale the nodes to chemPotTrue+-rangeMax

  for(iNode=0;iNode<numChemPot;iNode++){
    //[-range+center,range+center]
    chemPot[iNode] = chebyNode[iNode]*0.5*range+center;
    //printf("iNode %i chemPot %lg\n",iNode,chemPot[iNode]);
  }

  free(chebyNode);
/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcChemPotMetal(CP *cp,double *numOccDetProc)
/*========================================================================*/
{/*begin routine*/
/**************************************************************************/
/* For filter diag, we need to determine the chemical potential in the    */
/* smeared occupation number                                              */
/**************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
/*------------------------------------------------------------------------*/
  STODFTINFO *stodftInfo = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos = cp->stodftCoefPos;
  CPOPTS *cpopts = &(cp->cpopts);
  CPSCR *cpscr = &(cp->cpscr);
  COMMUNICATE *communicate = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);

  int myidState = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int smearOpt = stodftInfo->smearOpt;
  int iState;
  int maxIters = 1000;
  int cpLsda = cpopts->cp_lsda;
  int numChemPot = stodftInfo->numChemPot;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numStateUpTotal = stodftInfo->numStateStoUp;
  int numStateUpProcTotal = numStateUpProc*numChemPot;
  int numStateUpAllProc = numStateUpTotal*numChemPot;
  int numStateUpIdp = stodftInfo->numStateUpIdp;
  int numStatePrintUp = stodftInfo->numStatePrintUp;
  int numStatePrintDn = stodftInfo->numStatePrintDn;

  int *numStates = stodftInfo->numStates;
  int *dsplStates = stodftInfo->dsplStates;
  int *numStates2 = stodftInfo->numStates2;
  int *dsplStates2 = stodftInfo->dsplStates2;
  MPI_Comm comm_states = communicate->comm_states;


  double smearTemperature = stodftInfo->smearTemperature;
  double invTemp = 1.0/smearTemperature;
  double chemPotUp,chemPotDn;
  double numElecTrueUp = stodftInfo->numElecTrueUp;
  double numElecTrueDn = stodftInfo->numElecTrueDn;
  double numElecTrue = stodftInfo->numElecTrue;
  double numElecUpNow,numElecDnNow,numElecDiff;
  double numElecTol = 1.0e-11*numElecTrue;
  double chemPotMin,chemPotMax;
  double numElecMin,numElecMax;
  double chemPotNew,numElecNew;
  double energyMin = stodftInfo->energyMin;
  double pre = sqrt(2.0);

  double *energyLevel = stodftCoefPos->energyLevel;
  double *numOccDetAll;

  
  /*
  // For closed shell system, we will use half the number of 
  // electrons to determine the chemical potential so that 
  // we don't need *2 factor in smearing function. 
  if(cpLsda=0){
    numElecTrueUp = numElecTrue;
  }
  else{
    numElecTrueUp = numElecTrueUp;
  }
  */
  if(myidState==0){
    numOccDetAll = (double*)cmalloc(numStateUpIdp*sizeof(double)); 
    for(iState=0;iState<numStateUpIdp;iState++){
      numOccDetAll[iState] = 0.0;
    }
    // Now I only have spin-unpolarize filter diag
    numElecDiff = 1000;
    chemPotMin = energyLevel[0];
    chemPotMax = energyLevel[numStatePrintUp-1];
    numElecMin = calcNumElecSmear(smearOpt,smearTemperature,chemPotMin,
                             energyLevel,numStatePrintUp);
    numElecMax = calcNumElecSmear(smearOpt,smearTemperature,chemPotMax,
                             energyLevel,numStatePrintUp);
    if(numElecMin>numElecTrueUp||numElecMax<numElecTrueUp){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("The number of electron %lg is out of range %lg %lg\n",
              numElecTrueUp,numElecMin,numElecMax);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }
    chemPotNew = (numElecTrueUp-numElecMin)*(chemPotMax-chemPotMin)/(numElecMax-numElecMin)+
                 chemPotMin;
    numElecNew = calcNumElecSmear(smearOpt,smearTemperature,chemPotNew,
                             energyLevel,numStatePrintUp);
    //int icount = 0;
    while(fabs(numElecNew-numElecTrueUp)>numElecTol){
      //printf("chemPotMin %.8lg chemPotMax %.8lg chemPotNew %.8lg numElecMin %.8lg numElecMax %.8lg numElecNew %.8lg diff %.8lg\n",chemPotMin,chemPotMax,chemPotNew,numElecMin,numElecMax,numElecNew,fabs(numElecTrueUp-numElecNew));
      if(numElecNew>numElecTrueUp){
        chemPotMax = chemPotNew;
        numElecMax = numElecNew;
      }
      if(numElecNew<numElecTrueUp){
        chemPotMin = chemPotNew;
        numElecMin = numElecNew;
      }
      chemPotNew = (numElecTrueUp-numElecMin)*(chemPotMax-chemPotMin)/(numElecMax-numElecMin)+
                   chemPotMin;
      numElecNew = calcNumElecSmear(smearOpt,smearTemperature,chemPotNew,
                             energyLevel,numStatePrintUp);
      //icount += 1;
      //if(icount==5)exit(0);
    }


    // I'll do it later
    //if(cpLsda=1){
    //}
    
    printf("Finish Calculating Chemical Potential for Smearing\n");  
    printf("The correct chemical potential is %.16lg Ne %.16lg DNe %.16lg\n",chemPotNew,numElecNew,
              fabs(numElecNew-numElecTrueUp));
    switch(smearOpt){
      case 1:
        for(iState=0;iState<numStatePrintUp;iState++){
          numOccDetAll[iState] = 1.0/(1.0+exp(invTemp*(energyLevel[iState]-chemPotNew)));
        }
        break;
      case 2:
        for(iState=0;iState<numStatePrintUp;iState++){
          numOccDetAll[iState] = 0.5*(1.0-erf(invTemp*(energyLevel[iState]-chemPotNew)));
        }
        break;
    }//endswitch
    printf("Orbital energies and occupation numbers\n");
    for(iState=0;iState<numStatePrintUp;iState++){
      printf("%.8lg %.8lg\n",energyLevel[iState],numOccDetAll[iState]);
    }
  }//endif myidState

  if(numProcStates>1){
    Barrier(comm_states);
    Scatterv(numOccDetAll,numStates2,dsplStates2,MPI_DOUBLE,numOccDetProc,
             numStateUpProc,MPI_DOUBLE,0,comm_states);
  }
  else{
    memcpy(numOccDetProc,numOccDetAll,numStateUpIdp*sizeof(double));
  }
  for(iState=0;iState<numStateUpProc;iState++){
    numOccDetProc[iState] = sqrt(numOccDetProc[iState])*pre;
  }

  if(myidState==0){
    free(numOccDetAll);
  }

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genChemPotEnergyWindows(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This function generate chemical potentials for each energy window.    */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  


/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double calcNumElecSmear(int smearOpt,double smearTemperature,double chemPot,
                 double *orbital,int numStates)
/*========================================================================*/
{/*begin routine*/
/**************************************************************************/
/* For filter diag, we need to determine the chemical potential in the    */
/* smeared occupation number                                              */
/**************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
/*------------------------------------------------------------------------*/
  int iState;
  double numElec = 0.0;
  double x;
  double invTemp = 1.0/smearTemperature;

  switch(smearOpt){
    case 0:
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("Have to choose a smearing scheme when you work with metallic\n");
      printf("system. \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
      break;
    case 1:
      for(iState=0;iState<numStates;iState++){
        x = 1.0/(1.0+exp(invTemp*(orbital[iState]-chemPot)));
        //printf("%i %lg %lg %lg\n",iState,chemPot,orbital[iState],x);
        numElec += x;
      }//endfor
      break;
    case 2:
      for(iState=0;iState<numStates;iState++){
        x = 0.5*(1.0-erf(invTemp*(orbital[iState]-chemPot)));
        numElec += x;
      }//endfor
      break;
  }//endswitch

  return numElec;
  
/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



