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
void genNewtonHermit(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This routine first guess a polynormial chain by an emperical formular */
/* The chain length is then increased to satisfied the fitting error for */
/* ALL different chemical potentials					 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  NEWTONINFO *newtonInfo  = stodftInfo->newtonInfo;

  double fitErrTol = stodftInfo->fitErrTol;
  double fitErr;
  double beta = stodftInfo->beta;
  double energyDiff = stodftInfo->energyDiff;
  double *chemPot          = stodftCoefPos->chemPot;
  double *sampPoint;
  double *expanCoeff;
  double *sampPointUnscale;

  int numChemPot     = stodftInfo->numChemPot;
  int polynormLength = (int)(4.0*beta*energyDiff); //initial chain length
  int totalPoly      = polynormLength*numChemPot;  //iniital total polynormial

  FERMIFUNR fermiFunction = stodftInfo->fermiFunctionReal;
  
/*==========================================================================*/
/* I) Generate coeffcients for initial chain length  */
 
  printf("polynormLength %i\n",polynormLength); 
  stodftInfo->polynormLength = polynormLength;
  stodftCoefPos->expanCoeff = (double *)cmalloc(totalPoly*sizeof(double));
  newtonInfo->sampPoint = (double *)cmalloc(polynormLength*sizeof(double));
  newtonInfo->sampPointUnscale = (double *)cmalloc(polynormLength*sizeof(double));
 
  genSampNewtonHermit(stodftInfo,stodftCoefPos);  
  
  genCoeffNewtonHermit(stodftInfo,stodftCoefPos);
  
  fitErr = calcFitError(stodftInfo,stodftCoefPos);

/*==========================================================================*/
/* II) Increase chain length if polynomial expansion doesn't reach error    */
/*     torlerance */

  while(fitErr>fitErrTol){
    free(stodftCoefPos->expanCoeff);
    free(newtonInfo->sampPoint);
    free(newtonInfo->sampPointUnscale);
    polynormLength += 1000;
    stodftInfo->polynormLength = polynormLength;
    totalPoly = polynormLength*numChemPot;
    stodftCoefPos->expanCoeff = (double *)cmalloc(totalPoly*sizeof(double));
    newtonInfo->sampPoint = (double *)cmalloc(polynormLength*sizeof(double));
    newtonInfo->sampPointUnscale = (double *)cmalloc(polynormLength*sizeof(double));
    
    genSampNewtonHermit(stodftInfo,stodftCoefPos);
    genCoeffNewtonHermit(stodftInfo,stodftCoefPos);
    fitErr = calcFitError(stodftInfo,stodftCoefPos);  
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/



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

  double timeStart,timeEnd;
  //FILE *fileCoeff = fopen("coeff","w");
  FERMIFUNR fermiFunction = stodftInfo->fermiFunctionReal;
 
  /*
  for(iPoly=0;iPoly<polynormLength;iPoly++){
    funValue = fermiFunction(sampPointUnscale[iPoly],chemPot[0],beta);
    printf("funValue %lg %lg\n",sampPointUnscale[iPoly],funValue);
  }

  fflush(stdout);
  exit(0);
  */
  cputime(&timeStart);  
 
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

  cputime(&timeEnd);

  printf("Coeff time %lg\n",timeEnd-timeStart);
  //debug
  /*
  for(imu=0;imu<numChemPot;imu++){
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      printf("mu %lg expanCoeff %lg\n",chemPot[imu],expanCoeff[iPoly]);
    }
  }
  */
  /*
  for(imu=0;imu<numChemPot;imu++){
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      fprintf(fileCoeff,"%.13lg\n",expanCoeff[iPoly*numChemPot+imu]);
    }
  }
  fclose(fileCoeff); 
  */

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
  double Smin = newtonInfo->Smin;
  double Smax = newtonInfo->Smax;
  double scale = 1.0/newtonInfo->scale;
  double energyMin = stodftInfo->energyMin;
  double delta = (Smax-Smin)/(double)(numSampCand-1);
  double obj,objMax,diff;
  double *sampPoint = (double*)newtonInfo->sampPoint;
  double *sampPointUnscale = (double*)newtonInfo->sampPointUnscale;
  double *sampCand = (double*)cmalloc(numSampCand*sizeof(double));
  double *objValueArray = (double*)cmalloc(numSampCand*sizeof(double));

  double timeStart,timeEnd;
  double prod;
  FILE *fileSampPoint = fopen("samp-point","w");

/*==========================================================================*/
/* 0) Generate sample candidates in range [Smin,Smax] */
  cputime(&timeStart); 
 
  for(iCand=0;iCand<numSampCand;iCand++)sampCand[iCand] = Smin+iCand*delta;

/*==========================================================================*/
/* 1) Select samples form sample candidates  */

  sampPoint[0] = sampCand[0];

  for(iCand=0;iCand<numSampCand;iCand++){
    diff = sampCand[iCand]-sampPoint[0];
    diff = diff*diff;
    if(diff<1.0e-20) objValueArray[iCand] = -1.0e30;
    else objValueArray[iCand] = log(diff);
  }

  for(iPoly=1;iPoly<polynormLength;iPoly++){
    if(iPoly%1000==0)printf("iPoly %i\n",iPoly);
    objMax = -100000.0;
    for(iCand=0;iCand<numSampCand;iCand++){
      if(objValueArray[iCand]>objMax){
	objMax = objValueArray[iCand];
	objMaxIndex = iCand;
      }
    }//endfor iCand
    sampPoint[iPoly] = sampCand[objMaxIndex];
    for(iCand=0;iCand<numSampCand;iCand++){
      diff = sampCand[iCand]-sampPoint[iPoly];
      diff = diff*diff;
      if(diff<1.0e-20) objValueArray[iCand] = -1.0e30;
      else objValueArray[iCand] = log(diff);
    }//endfor iCand
  }//endfor iPoly

  //debug
  /*
  for(iPoly=0;iPoly<polynormLength;iPoly++){
    printf("iPoly %i samp %lg\n",iPoly,sampPoint[iPoly]);
  }
  */
  cputime(&timeEnd);
  printf("Samp time %lg\n",timeEnd-timeStart);
  for(iPoly=0;iPoly<polynormLength;iPoly++){
    fprintf(fileSampPoint,"%.13lg\n",sampPoint[iPoly]);
  }
  fclose(fileSampPoint);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genEigenOrb(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                 CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CELL *cell                    = &(general_data->cell);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  EWALD *ewald                  = &(general_data->ewald);
  EWD_SCR *ewd_scr              = &(class->ewd_scr);
  ATOMMAPS *atommaps            = &(class->atommaps);
  FOR_SCR *for_scr              = &(class->for_scr);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  PTENS *ptens                  = &(general_data->ptens);
  SIMOPTS *simopts              = &(general_data->simopts);

  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPEWALD *cpewald              = &(cp->cpewald);
  CPSCR *cpscr                  = &(cp->cpscr);
  PSEUDO *pseudo                = &(cp->pseudo);
  COMMUNICATE *communicate      = &(cp->communicate);

  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm            = &(cp->cp_sclr_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_sm            = &(cp->cp_para_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_dens_cp_box   = &(cp->cp_sclr_fft_pkg3d_dens_cp_box);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box   = &(cp->cp_para_fft_pkg3d_dens_cp_box);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_lg             = &(cp->cp_sclr_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg             = &(cp->cp_para_fft_pkg3d_lg);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_up          = &(cp->cp_comm_state_pkg_up);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_dn          = &(cp->cp_comm_state_pkg_dn);

  int numStateUpProc  = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc  = cpcoeffs_info->nstate_dn_proc;
  int numCoeff        = cpcoeffs_info->ncoef;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int cpLsda          = cpopts->cp_lsda;
  int numCoeffUpTot   = numStateUpProc*numCoeff;
  int numCoeffDnTot   = numStateUpProc*numCoeff;
  int iCoeff;
  int cpWaveMin     = simopts->cp_wave_min;
  int cpMin         = simopts->cp_min;
  int cpWaveMinPimd = simopts->cp_wave_min_pimd;
  int cpMinOn = cpWaveMin + cpMin + cpWaveMinPimd;


  double *kseig_vals = (double*)cmalloc(numStateUpProc*sizeof(double))-1;
  double *kseig_vecs = (double*)cmalloc(numStateUpProc*numStateUpProc*sizeof(double))-1;
  double *ksmat_test = (double*)cmalloc(numStateUpProc*numStateUpProc*sizeof(double))-1;
  double *ks_scr = (double*)cmalloc(numStateUpProc*numStateUpProc*sizeof(double))-1;
  double *rs_scr1 = (double*)cmalloc(numStateUpProc*numStateUpProc*sizeof(double))-1;
  double *rs_scr2 = (double*)cmalloc(numStateUpProc*numStateUpProc*sizeof(double))-1;
  int *icoef_orth_up    = &(cpcoeffs_pos->icoef_orth_up);
  int *icoef_form_up    = &(cpcoeffs_pos->icoef_form_up);
  int *ifcoef_orth_up   = &(cpcoeffs_pos->ifcoef_orth_up);
  int *ifcoef_form_up   = &(cpcoeffs_pos->ifcoef_form_up);
  int *ioff_upt      = cpcoeffs_info->ioff_upt;
  double kseig_sum;

  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *cre_temp = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;
  double *cim_temp = (double*)cmalloc(numCoeffUpTot*sizeof(double))-1;

  calcCoefForceWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);

  cp_condiag_ksmat(cre_up,cim_up,*icoef_form_up,*icoef_orth_up,fcre_up,fcim_up,
                 *ifcoef_form_up,*ifcoef_orth_up,kseig_vals,kseig_vecs,
                  ksmat_test,ks_scr,rs_scr1,rs_scr2,ioff_upt,
                  &(cp->cp_comm_state_pkg_up),&kseig_sum);

  cp_rotate_vector(cre_up,cim_up,*icoef_form_up,
                      kseig_vecs,ioff_upt,cre_temp,cim_temp,
                      &(cp->cp_comm_state_pkg_up));

  free(&kseig_vals[1]);
  free(&kseig_vecs[1]);
  free(&ksmat_test[1]);
  free(&ks_scr[1]);
  free(&rs_scr1[1]);
  free(&rs_scr2[1]);
  free(&cre_temp[1]);
  free(&cim_temp[1]);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double calcFitError(STODFTINFO *stodftInfo,STODFTCOEFPOS *stodftCoefPos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
   
  NEWTONINFO *newtonInfo  = stodftInfo->newtonInfo;

  int numPointTest = 100;
  int numChemPot = stodftInfo->numChemPot;
  int polynormLength = stodftInfo->polynormLength;
  int iPoint,iChem,iPoly;
  double pointTest;
  double beta       = stodftInfo->beta;
  double energyMin  = stodftInfo->energyMin;
  double energyMax  = stodftInfo->energyMax;
  double energyMean = stodftInfo->energyMean;
  double deltPoint  = (energyMax-energyMin)/numPointTest;
  double scale	    = newtonInfo->scale;
  double pointScale;
  double funValue,prod,funBM;
  double fitErr = -1.0e30;
  double diff;

  double *chemPot          = stodftCoefPos->chemPot;
  double *sampPoint = newtonInfo->sampPoint;
  double *expanCoeff = stodftCoefPos->expanCoeff;
  double *sampPointUnscale = newtonInfo->sampPointUnscale;


  FERMIFUNR fermiFunction = stodftInfo->fermiFunctionReal;


  for(iChem=0;iChem<numChemPot;iChem++){
    for(iPoint=0;iPoint<numPointTest;iPoint++){
      pointTest = energyMin+(iPoint+0.5)*deltPoint;
      pointScale = (pointTest-energyMean)*scale;
      funValue = expanCoeff[0];
      prod = 1.0;
      for(iPoly=1;iPoly<polynormLength;iPoly++){
	prod *= pointScale-sampPoint[iPoly-1];
	funValue += expanCoeff[iPoly]*prod;
      }
      funBM = fermiFunction(pointScale,chemPot[iChem],beta);
      diff = fabs(funValue-funBM);
      if(diff>fitErr)fitErr = diff;
      printf("TestFunExpan %lg %lg %lg\n",pointTest,pointScale,funValue);
    }
  }

  return fitErr;

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

