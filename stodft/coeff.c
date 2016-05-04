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
#include "../proto_defs/proto_stodft_local.h"

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

  FERMIFUNR fermiFunction = stodftInfo->fermiFunctionReal;

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

  FERMIFUNC fermiFunction = stodftInfo->fermiFunctionComplex;

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
  double Smin = newtonInfo->Smin;
  double Smax = newtonInfo->Smax;
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
    objMax = -100000.0;
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
/* Call this function only for the first process.			 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CELL *cell			= &(general_data->cell);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  EWALD *ewald		        = &(general_data->ewald);
  EWD_SCR *ewd_scr		= &(class->ewd_scr);
  ATOMMAPS *atommaps		= &(class->atommaps);
  FOR_SCR *for_scr		= &(class->for_scr);
  STAT_AVG *stat_avg		= &(general_data->stat_avg);
  PTENS *ptens			= &(general_data->ptens);
  SIMOPTS *simopts		= &(general_data->simopts);

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


  double randMin	= -1.0;
  double randMax	= 1.0;
  double length		= 0.0;
  double energyConv     = 1.0e-6;
  double energy		= 0.0;
  double energyOld;

  int numStateUpProc  = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc  = cpcoeffs_info->nstate_dn_proc;
  int numCoeff        = cpcoeffs_info->ncoef;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numIteration    = 100;
  int iIter;
  int iState,iCoeff,iCoeffStart,index1,index2;
  int cpWaveMin     = simopts->cp_wave_min;
  int cpMin         = simopts->cp_min;
  int cpWaveMinPimd = simopts->cp_wave_min_pimd;
  int cpMinOn = cpWaveMin + cpMin + cpWaveMinPimd;

  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *coeffReUpBackup = (double*)cmalloc((numCoeff+1)*sizeof(double));
  double *coeffImUpBackup = (double*)cmalloc((numCoeff+1)*sizeof(double));
  double *randTrail = (double *)cmalloc((2*numCoeff+1)*sizeof(double));

/*==========================================================================*/
/* I) Set parameters and backup						    */

  printf("==============================================");
  printf("Estimate Largest Energy\n");

  cpcoeffs_info->nstate_up_proc = 1;
  cpcoeffs_info->nstate_dn_proc = 1;

  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    fcre_up[iCoeff] = 0.0;
    fcim_up[iCoeff] = 0.0;
    coeffReUpBackup[iCoeff] = cre_up[iCoeff];
    coeffImUpBackup[iCoeff] = cim_up[iCoeff];
  }//endfor iCoeff

/*==========================================================================*/
/* II) Prepare a random initial wave function                               */
/*     We can't use the readin coeff since it is an eigenfunction of H      */
/*     Pick a random orbital and normalize it.			            */
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
  double seed = 2.0;
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

  for(iIter=0;iIter<numIteration;iIter++){

/*--------------------------------------------------------------------------*/
/* i) Reset force                                                           */
    
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      fcre_up[iCoeff] = 0.0;
      fcim_up[iCoeff] = 0.0;
    }

/*--------------------------------------------------------------------------*/
/* ii) Calcualte H|phi>							    */

    control_cp_eext_recip(clatoms_info,clatoms_pos,cpcoeffs_info,
			 cpcoeffs_pos,cpewald,cpscr,cpopts,pseudo,
			 ewd_scr,atommaps,cell,ewald,ptens,&(stat_avg->vrecip),
			 &(stat_avg->cp_enl),communicate,for_scr,cpDualGridOptOn,
			 cp_para_fft_pkg3d_lg);

    coef_force_control(cpopts,cpcoeffs_info,cpcoeffs_pos,cpscr,ewald,cpewald,
		      cell,stat_avg,pseudo->vxc_typ,ptens->pvten_tmp,pseudo->gga_cut,
		      pseudo->alpha_conv_dual,pseudo->n_interp_pme_dual,cpMinOn,
		      communicate,cp_comm_state_pkg_up,
		       cp_comm_state_pkg_dn,cp_para_fft_pkg3d_lg,cp_sclr_fft_pkg3d_lg,
		       cp_para_fft_pkg3d_dens_cp_box,cp_sclr_fft_pkg3d_dens_cp_box,
		       cp_para_fft_pkg3d_sm,cp_sclr_fft_pkg3d_sm,cpDualGridOptOn);

/*--------------------------------------------------------------------------*/
/* iii) Calcluate <phi|H|phi>	                                            */
    
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      fcre_up[iCoeff] *= -0.5;
      fcim_up[iCoeff] *= -0.5; 
    }
    fcre_up[numCoeff] *= -1;
    
    energyOld = energy;
    energy = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      energy += fcre_up[iCoeff]*cre_up[iCoeff]+fcim_up[iCoeff]+cim_up[iCoeff];
    }
    energy *= 2.0;
    energy += fcre_up[numCoeff]*cre_up[numCoeff];
    // We already normalize wf to 1.0, so we don't have scaling 0.5 here
    printf("iStep %i Energy %lg\n",iIter,energy);

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
/* v) Check Convergence			                                    */
    if(fabs(energyOld-energy)<energyConv)break;
    	
  }//endfor iIter

/*==========================================================================*/
/* IV) Restore flags and clean up                                           */

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
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPEWALD *cpewald              = &(cp->cpewald);
  CPSCR *cpscr                  = &(cp->cpscr);
  PSEUDO *pseudo                = &(cp->pseudo);
  COMMUNICATE *communicate      = &(cp->communicate);

  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm            = &(cp->cp_sclr_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_sm            = &(cp->cp_para_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_dens_cp_box   = &(cp->cp_sclr_fft_pkg3d_dens_cp_box);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box   = &(cp->cp_para_fft_pkg3d_dens_cp_box);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_lg            = &(cp->cp_sclr_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg            = &(cp->cp_para_fft_pkg3d_lg);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_up         = &(cp->cp_comm_state_pkg_up);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_dn         = &(cp->cp_comm_state_pkg_dn);


  double randMin        = -1.0;
  double randMax        = 1.0;
  double length         = 0.0;
  double energyConv     = 1.0e-6;
  double energyMax	= stodftInfo->energyMax;
  double energy,energyOld;

  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numIteration   = 100;
  int iIter;
  int iState,iCoeff,iCoeffStart,index1,index2;
  int cpWaveMin     = simopts->cp_wave_min;
  int cpMin         = simopts->cp_min;
  int cpWaveMinPimd = simopts->cp_wave_min_pimd;
  int cpMinOn = cpWaveMin + cpMin + cpWaveMinPimd;

  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *coeffReUpBackup = (double*)cmalloc((numCoeff+1)*sizeof(double));
  double *coeffImUpBackup = (double*)cmalloc((numCoeff+1)*sizeof(double));
  double *randTrail = (double *)cmalloc((2*numCoeff+1)*sizeof(double));

/*==========================================================================*/
/* I) Set parameters and backup                                             */

  printf("==============================================");
  printf("Estimate Largest Energy\n");

  cpcoeffs_info->nstate_up_proc = 1;
  cpcoeffs_info->nstate_dn_proc = 1;

  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    fcre_up[iCoeff] = 0.0;
    fcim_up[iCoeff] = 0.0;
    coeffReUpBackup[iCoeff] = cre_up[iCoeff];
    coeffImUpBackup[iCoeff] = cim_up[iCoeff];
  }//endfor iCoeff

/*==========================================================================*/
/* II) Prepare a random initial wave function                               */
/*     We can't use the readin coeff since it is an eigenfunction of H      */
/*     Pick a random orbital and normalize it.                              */
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

  for(iIter=0;iIter<numIteration;iIter++){

/*--------------------------------------------------------------------------*/
/* i) Reset force                                                           */

    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      fcre_up[iCoeff] = 0.0;
      fcim_up[iCoeff] = 0.0;
    }

/*--------------------------------------------------------------------------*/
/* ii) Calcualte H|phi>                                                     */

    control_cp_eext_recip(clatoms_info,clatoms_pos,cpcoeffs_info,
                         cpcoeffs_pos,cpewald,cpscr,cpopts,pseudo,
                         ewd_scr,atommaps,cell,ewald,ptens,&(stat_avg->vrecip),
                         &(stat_avg->cp_enl),communicate,for_scr,cpDualGridOptOn,
                         cp_para_fft_pkg3d_lg);

    coef_force_control(cpopts,cpcoeffs_info,cpcoeffs_pos,cpscr,ewald,cpewald,
                      cell,stat_avg,pseudo->vxc_typ,ptens->pvten_tmp,pseudo->gga_cut,
                      pseudo->alpha_conv_dual,pseudo->n_interp_pme_dual,cpMinOn,
                      communicate,cp_comm_state_pkg_up,
                       cp_comm_state_pkg_dn,cp_para_fft_pkg3d_lg,cp_sclr_fft_pkg3d_lg,
                       cp_para_fft_pkg3d_dens_cp_box,cp_sclr_fft_pkg3d_dens_cp_box,
                       cp_para_fft_pkg3d_sm,cp_sclr_fft_pkg3d_sm,cpDualGridOptOn);
/*--------------------------------------------------------------------------*/
/* iii) Calcluate <phi|H|phi>                                               */

    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      fcre_up[iCoeff] *= -0.5;
      fcim_up[iCoeff] *= -0.5;
    }
    fcre_up[numCoeff] *= -1;

    energyOld = energy;
    energy = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      energy += fcre_up[iCoeff]*cre_up[iCoeff]+fcim_up[iCoeff]+cim_up[iCoeff];
    }
    energy *= 2.0;
    energy += fcre_up[numCoeff]*cre_up[numCoeff];
    // We already normalize wf to 1.0, so we don't have scaling 0.5 here
    printf("iStep %i Energy %lg\n",iIter,energy);

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
    if(fabs(energyOld-energy)<energyConv)break;

  }//endfor iIter

/*==========================================================================*/
/* IV) Restore flags and clean up                                           */

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

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


