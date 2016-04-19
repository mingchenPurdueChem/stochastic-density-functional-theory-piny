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
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define TIME_CP_OFF

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
  int iPoly,jPoly,imu;

  double Smin		   = newtonInfo->Smin;
  double Smax		   = newtonInfo->Smax;
  double scale		   = newtonInfo->scale;
  double *sampPoint        = (double*)newtonInfo->sampPoint;
  double *expanCoeff       = (double*)stodftCoefPos->expanCoeff;
  double *sampPointUnscale = (double*)newtonInfo->sampPointUnscale;
  double *chemPot	   = stodftCoefPos->chemPot;

  double beta = stodftInfo->beta;
  double funValue,sum,prod;

  FERMIFUNC fermiFunction = stodftInfo->fermiFunction;

  for(imu=0;imu<numChemPot;imu++){
    funValue = fermiFunction(sampPointUnscale[0],chemPot[imu],beta);
    expanCoeff[imu] = funValue;
    funValue = fermiFunction(sampPointUnscale[1],chemPot[imu],beta);
    expanCoeff[numChemPot+imu] = (funValue-expanCoeff[imu])/(sampPoint[1]-sampPoint[0]);
    for(iPoly=2;iPoly<polynormLength;iPoly++){
      funValue = fermiFunction(sampPointUnscale[1],chemPot[imu],beta);
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

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genCoeffNewtonNoHermit(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate expanCoeff for different chemical potentials                */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  NEWTONINFO *newtonInfo  = stodftInfo->newtonInfo;


  int polynormLength     = stodftInfo->polynormLength;
  int numChemPot         = stodftInfo->numChemPot;
  int iPoly,jPoly,imu;

  double Smin              = newtonInfo->Smin;
  double Smax              = newtonInfo->Smax;
  double scale             = newtonInfo->scale;
  double complex *sampPoint        = (double complex*)newtonInfo->sampPoint;
  double complex *expanCoeff       = (double complex*)stodftCoefPos->expanCoeff;
  double complex *sampPointUnscale = (double complex*)newtonInfo->sampPointUnscale;
  double *chemPot          = stodftCoefPos->chemPot;

  double beta = stodftInfo->beta;
  double funValue,sum,prod;

  FERMIFUNC fermiFunction = stodftInfo->fermiFunction;

  for(imu=0;imu<numChemPot;imu++){
    funValue = fermiFunction(sampPointUnscale[0],chemPot[imu],beta);
    expanCoeff[imu] = funValue;
    funValue = fermiFunction(sampPointUnscale[1],chemPot[imu],beta);
    expanCoeff[numChemPot+imu] = (funValue-expanCoeff[imu])/(sampPoint[1]-sampPoint[0]);
    for(iPoly=2;iPoly<polynormLength;iPoly++){
      funValue = fermiFunction(sampPointUnscale[iPoly],chemPot[imu],beta);
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
  double Smin = stodftInfo->Smin;
  double Smax = stodftInfo->Smax;
  double scale = 1.0/newtonInfo->scale;
  double energyMin = stodftInfo->energyMin;
  double delta = (Smax-Smin)/(double)(numSampCand-1);
  double obj,objMax,diff;
  double *sampPoint = (double*)newtonInfo->sampPoint;
  double *sampPointUnscale = (double*)newtonInfo->sampPointUnscale;
  double *sampCand = (double*)cmalloc(numSampCand*sizeof(double));

/*==========================================================================*/
/* 0) Generate sample candidates in range [Smin,Smax] */
  
  for(iPoly=0;iPoly<polynormLength;iPoly++)sampCand[iPoly] = Smin+iPoly*delta;

/*==========================================================================*/
/* 1) Select samples form sample candidates  */

  sampPoint[0] = sampCand[0];
  for(iPoly=1;iPoly<polynormLength;iPoly++){
    objMax = -100000.0
    for(iCand=0;iCand<numSampCand;iCand++){
      obj = 0.0;
      for(jPoly=0;jPoly<iPoly;jPoly++){
	diff = sampCand[iCand]-sampPoint[jPoly];
	if(diff<1.0e-10)obj += -1.0e30;
	else obj += log(diff*diff);
      }//endfor jPoly
      if(obj>objMax){
	objMax = obj;
	objMaxIndex = iCand;
      }//endif
    }//endfor iCand
    sampPoint[iPoly] = sampCand[objMaxIndex];
  }//endfor iPoly

/*==========================================================================*/
/* 2) Rescale the sample points to energy space  */
  
  for(iPoly=0;iPoly<polynormLength;iPoly++){
    sampPointUnscale[iPoly] = (sampPoint[iPoly]-Smin)*scale+energyMin;
  }
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

