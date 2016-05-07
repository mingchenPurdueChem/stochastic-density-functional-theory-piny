/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: filter.c                                       */
/*                                                                          */
/* This routine constructs filters of Fermi function.                       */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_stodft_local.h"

#define TIME_CP_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void filterNewtonPolyHerm(CP *cp,CLASS *class,GENERAL_DATA *general_data,
			  int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CLATOMS_POS*  clatoms_pos     = &(class->clatoms_pos[ip_now]);

  NEWTONINFO *newtonInfo = stodftInfo->newtonInfo;

  int expanType	     = stodftInfo->expanType;
  int polynormLength = stodftInfo->polynormLength;
  int numChemPot     = stodftInfo->numChemPot;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int cpLsda	     = cpopts->cp_lsda;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int imu,iCoeff,iPoly,indexStart;
  int startIndex;

  double energyDiff  = stodftInfo->energyDiff;
  double energyMin   = stodftInfo->energyMin;
  double energyMax   = stodftInfo->energyMax;
  double energyMean  = stodftInfo->energyMean;
  double scale       = newtonInfo->scale;
  double polyCoeff;
  double *sampPoint = (double*)newtonInfo->sampPoint;
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;


  double *expanCoeff = (double*)stodftCoefPos->expanCoeff;
  

//debug
/*  
  int numPointTest = 100;
  int iPoint;
  double pointTest;
  double energyMinTest = -0.1150735;
  double energyMaxTest = 0.10104944;
  double deltPoint = (energyMaxTest-energyMinTest)/numPointTest;
  double pointScale;
  double funValue,prod;
  for(iPoint=0;iPoint<numPointTest;iPoint++){
    pointTest = energyMinTest+(iPoint+0.5)*deltPoint;
    pointScale = (pointTest-energyMean)*scale;
    funValue = expanCoeff[0];
    prod = 1.0;
    for(iPoly=1;iPoly<polynormLength;iPoly++){
      prod *= pointScale-sampPoint[iPoly-1];
      funValue += expanCoeff[iPoly]*prod;
    }
    printf("TestFunExpan %lg %lg %lg\n",pointTest,pointScale,funValue);
  }
  fflush(stdout);
  exit(0);
*/
//debug
 
/*==========================================================================*/
/* 0) Copy the initial stochastic orbital */
 
  for(imu=0;imu<numChemPot;imu++){
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      stoWfUpRe[imu][iCoeff] = expanCoeff[imu]*cre_up[iCoeff];
      stoWfUpIm[imu][iCoeff] = expanCoeff[imu]*cim_up[iCoeff];
    }//endfor iCoeff
    if(cpLsda==1&&numStateDnProc!=0){
      for(iCoeff=1;iCoeff<=numCoeffDnTotal;iCoeff++){
	stoWfDnRe[imu][iCoeff] = expanCoeff[imu]*cre_dn[iCoeff];
	stoWfDnIm[imu][iCoeff] = expanCoeff[imu]*cim_dn[iCoeff];
      }//endfor iCoeff      
    }//endif 
  }//endfor imu

/*==========================================================================*/
/* 1) Loop over all polynomial terms (iPoly=0<=>polynomial order 1) */
  
  for(iPoly=1;iPoly<polynormLength;iPoly++){
    printf("iPoly %i\n",iPoly);
    normHNewtonHerm(cp,class,general_data,
                 cpcoeffs_pos,clatoms_pos,sampPoint[iPoly-1]);
    for(imu=0;imu<numChemPot;imu++){
      polyCoeff = expanCoeff[iPoly*numChemPot+imu];
      for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
	stoWfUpRe[imu][iCoeff] += polyCoeff*cre_up[iCoeff];	
        stoWfUpIm[imu][iCoeff] += polyCoeff*cim_up[iCoeff];                       

      }//endfor iCoeff
      if(cpLsda==1&&numStateDnProc!=0){
        for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
	  stoWfDnRe[imu][iCoeff] += polyCoeff*cre_dn[iCoeff];                     
	  stoWfDnIm[imu][iCoeff] += polyCoeff*cim_dn[iCoeff];
        }//endfor iCoeff        
      }//endif 
    }//endfor imu
  }//endfor iPoly

/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/

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
/* Call this function only for the first process.	     */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CELL *cell		= &(general_data->cell);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  EWALD *ewald	            = &(general_data->ewald);
  EWD_SCR *ewd_scr	= &(class->ewd_scr);
  ATOMMAPS *atommaps	    = &(class->atommaps);
  FOR_SCR *for_scr	= &(class->for_scr);
  STAT_AVG *stat_avg	    = &(general_data->stat_avg);
  PTENS *ptens		= &(general_data->ptens);
  SIMOPTS *simopts	= &(general_data->simopts);

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


  double randMin    = -1.0;
  double randMax    = 1.0;
  double length	    = 0.0;
  double energyConv     = 1.0e-6;
  double energy	    = 0.0;
  double energyOld;

  int numStateUpProc  = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc  = cpcoeffs_info->nstate_dn_proc;
  int numCoeff        = cpcoeffs_info->ncoef;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int cpLsda          = cpopts->cp_lsda;
  int numCoeffUpTot   = numStateUpProc*numCoeff;
  int numCoeffDnTot   = numStateDnProc*numCoeff;

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
/* I) Set parameters and backup			        */

  printf("==============================================\n");
  printf("Estimate Largest Energy\n");

  cpcoeffs_info->nstate_up_proc = 1;
  cpcoeffs_info->nstate_dn_proc = 1;

  /*
  for(iCoeff=1;iCoeff<=numCoeffUpTot;iCoeff++){
    fcre_up[iCoeff] = 0.0;
    fcim_up[iCoeff] = 0.0;
  }
  */

  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
    coeffReUpBackup[iCoeff] = cre_up[iCoeff];
    coeffImUpBackup[iCoeff] = cim_up[iCoeff];
  }//endfor iCoeff
  /*
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
  */ 

/*==========================================================================*/
/* II) Prepare a random initial wave function                               */
/*     We can't use the readin coeff since it is an eigenfunction of H      */
/*     Pick a random orbital and normalize it.		            */

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
  double seed = 15.0;
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
    
    /*
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      fcre_up[iCoeff] = 0.0;
      fcim_up[iCoeff] = 0.0;
    }
    */
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      fcre_up[iCoeff] = 0.0;
      fcim_up[iCoeff] = 0.0;
    }

/*--------------------------------------------------------------------------*/
/* ii) Calcualte H|phi>				    */

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
      energy += fcre_up[iCoeff]*cre_up[iCoeff]+fcim_up[iCoeff]*cim_up[iCoeff];
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
/* v) Check Convergence		                                    */
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
  //fflush(stdout);
  //exit(0);

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
  double energyMax  = stodftInfo->energyMax;
  double energy,energyOld;

  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numIteration   = 1000;
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

  printf("==============================================\n");
  printf("Estimate Smallest Energy\n");

  cpcoeffs_info->nstate_up_proc = 1;
  cpcoeffs_info->nstate_dn_proc = 1;

  for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
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
      energy += fcre_up[iCoeff]*cre_up[iCoeff]+fcim_up[iCoeff]*cim_up[iCoeff];
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

  //fflush(stdout);
  //exit(0);
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

