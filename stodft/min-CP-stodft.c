/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: min-CP-stodft.c                              */
/*                                                                          */
/* Subroutine for SCF calculation			                    */
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
void scfStodft(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                    CP *cp,int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  CELL *cell			 = &(general_data->cell);  
  PTENS *ptens			 = &(general_data->ptens);
  
  STODFTINFO *stodftInfo         = cp->stodftInfo;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CPEWALD *cpewald              = &(cp->cpewald);
  CPSCR *cpscr                  = &(cp->cpscr);
  PSEUDO *pseudo                = &(cp->pseudo);
  COMMUNICATE *communicate      = &(cp->communicate);

  FOR_SCR      *for_scr         = &(class->for_scr);
  EWD_SCR      *ewd_scr         = &(class->ewd_scr);
  CLATOMS_INFO *clatoms_info	= &(class->clatoms_info);
  CLATOMS_POS *clatoms_pos	= &(class->clatoms_pos[ip_now]);
  ATOMMAPS *atommaps		= &(class->atommaps);

  int iperd            		= cell->iperd;
  int iScf,iCell;
  int numScf			= stodftInfo->numScf; //Need claim this in cp
  int cpLsda 			= cpopts->cp_lsda;
  int checkPerdSize 		= cpopts->icheck_perd_size;
  int checkDualSize 		= cpopts->icheck_dual_size;
  int cpDualGridOptOn 	 	= cpopts->cp_dual_grid_opt;
  int numProcStates 		= communicate->np_states; 
  int myidState 		= communicate->myid_state;
  int coefFormUp 		= cpcoeffs_pos->icoef_form_up;
  int coefOrthUp                = cpcoeffs_pos->icoef_orth_up;
  int forceCoefFormUp 		= cpcoeffs_pos->ifcoef_form_up;
  int forceCoefOrthUp           = cpcoeffs_pos->ifcoef_orth_up;
  int coefFormDn                = cpcoeffs_pos->icoef_form_dn;
  int coefOrthDn                = cpcoeffs_pos->icoef_orth_dn;
  int forceCoefFormDn           = cpcoeffs_pos->ifcoef_form_dn;
  int forceCoefOrthDn           = cpcoeffs_pos->ifcoef_orth_dn;
  int numStateUp		= cpcoeffs_info->nstate_up_proc;
  int numStateDn                = cpcoeffs_info->nstate_dn_proc;

  int *pcoefFormUp		     = &(cpcoeffs_pos->icoef_form_up);
  int *pcoefOrthUp		     = &(cpcoeffs_pos->icoef_orth_up);
  int *pforceCoefFormUp              = &(cpcoeffs_pos->ifcoef_form_up);
  int *pforceCoefOrthUp              = &(cpcoeffs_pos->ifcoef_orth_up);
  int *pcoefFormDn		     = &(cpcoeffs_pos->icoef_form_dn);
  int *pcoefOrthDn		     = &(cpcoeffs_pos->icoef_orth_dn);
  int *pforceCoefFormDn              = &(cpcoeffs_pos->icoef_form_dn);
  int *pforceCoefOrthDn              = &(cpcoeffs_pos->icoef_orth_dn);


  double tolEdgeDist 		= cpopts->tol_edge_dist;
  double energyDiff		= -1.0;
  double energyTol		= 1.0;

  double *coeffReUp        = cpcoeffs_pos->cre_up;
  double *coeffImUp        = cpcoeffs_pos->cim_up;
  double *forceCoeffReUp   = cpcoeffs_pos->fcre_up;
  double *forceCoeffImUp   = cpcoeffs_pos->fcre_up;
  double *coeffReDn        = cpcoeffs_pos->cre_dn;
  double *coeffImDn        = cpcoeffs_pos->cim_dn;
  double *forceCoeffReDn   = cpcoeffs_pos->fcre_dn;
  double *forceCoeffImDn   = cpcoeffs_pos->fcre_dn;
  double *cpScrCoeffReUp   = cpscr->cpscr_wave.cre_up;
  double *cpScrCoeffImUp   = cpscr->cpscr_wave.cim_up;
  double *cpScrCoeffReDn   = cpscr->cpscr_wave.cre_dn;
  double *cpScrCoeffImDn   = cpscr->cpscr_wave.cim_dn;
  double *ptensPvtenTmp    = ptens->pvten_tmp;
  

/*======================================================================*/
/* 0) Check the forms                                                   */

  if(numProcStates>1){
    if((coefFormUp+forceCoefFormUp)!=2){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP vectors are not in transposed form \n");
     printf("on state processor %d in min_STD_cp \n",myidState);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cpLsda==1){
     if((coefFormDn+forceCoefFormDn)!=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("on state processor %d in min_STD_cp \n",myidState);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/
  

/*======================================================================*/
/* 0.05) Check the approximations in the methods                        */


  if((iperd<3||iperd==4)&&checkPerdSize==1){
    cp_boundary_check(cell,clatoms_info,clatoms_pos,tolEdgeDist);
  }/*endif*/
  if(cpDualGridOptOn>=1&&checkDualSize==1){
    cp_dual_check(cell,clatoms_info,clatoms_pos,
                  atommaps->cp_atm_lst,tolEdgeDist);
  }/*endif*/

/*======================================================================*/
/* I) Orthogonalize the coefs if norbing                                */

/*======================================================================*/
/* II) In parallel, transpose coefs back to normal form                 */

  if(numProcStates>1){
    cp_transpose_bck(coeffReUp,coeffImUp,pcoefFormUp,
                    cpScrCoeffReUp,cpScrCoeffImUp,&(cp->cp_comm_state_pkg_up));
    if(cpLsda==1&&numStateDn>0){
      cp_transpose_bck(coeffReDn,coeffImDn,pcoefFormDn,
                     cpScrCoeffReDn,cpScrCoeffImDn,&(cp->cp_comm_state_pkg_dn));
    }/*endif*/
  }/*endif*/

/*======================================================================*/
/* III) Initialize forces, pressure tensor, inverse hmat                */

  cpcoeffs_pos->ifcoef_form_up = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1;
  if(cpLsda==1&&numStateDn>0){
    cpcoeffs_pos->ifcoef_form_dn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }/*endif*/

  for(iCell=1;iCell<=9;iCell++){ptensPvtenTmp[iCell] = 0.0;}
  gethinv(cell->hmat_cp,cell->hmati_cp,&(cell->vol_cp),iperd);
  gethinv(cell->hmat,cell->hmati,&(cell->vol),iperd);

/*======================================================================*/
/* IV) Get the total density, (spin densities too if lsda)              */
/*       and necessary gradients of density for GGA calculations        */

  //debug only
  //calcRhoDeterm(class,bonded,general_data,cp,cpcoeffs_pos);

  printf("Finish generating density\n");

/*==========================================================================*/
/* V) Calculate the non-local pseudopotential list                          */

  if(stodftInfo->vpsAtomListFlag==0||cpDualGridOptOn>= 1){
    control_vps_atm_list(pseudo,cell,clatoms_pos,clatoms_info,
                         atommaps,ewd_scr,for_scr,cpDualGridOptOn,
                         stodftInfo->vpsAtomListFlag);
    stodftInfo->vpsAtomListFlag = 1;
  }
  printf("Finish generating Pseudopotential list.\n");

/*======================================================================*/
/* V) Calculate Emin and Emax				                */
/*       and necessary gradients of density for GGA calculations        */


  //debug only
  //genStoOrbital(class,bonded,general_data,cp,ip_now);

  
  //exit(0);
/*======================================================================*/
/* V) SCF loop						                */

  for(iScf=0;iScf<numScf;iScf++){

/*----------------------------------------------------------------------*/
/* i) Generate stochastic WF for different chemical potentials          */

    genStoOrbital(class,bonded,general_data,cp,ip_now);

/*----------------------------------------------------------------------*/
/* 2)  Get the total density, for each chemical potential and get       */
/*     total electron number for each chemical potential	        */
    calcRhoSto(class,bonded,general_data,cp,cpcoeffs_pos);

/*----------------------------------------------------------------------*/
/* 3)  Interpolate the chemical potential w.r.t			        */
/*     Correct electron number.						*/

/*----------------------------------------------------------------------*/
/* ii) Use the interpolation results to generate correct density        */
/*     (or pick a closest one) and update the density			*/

  }//endfor iScf

/*======================================================================*/
/* VI) In parallel, transpose coefs and coef forces fwd                 */

  if(numProcStates>1){
    cp_transpose_fwd(coeffReUp,coeffImUp,pcoefFormUp,
                    cpScrCoeffReUp,cpScrCoeffImUp,&(cp->cp_comm_state_pkg_up));
    cp_transpose_fwd(forceCoeffReUp,forceCoeffImUp,pforceCoefFormUp,
                    cpScrCoeffReUp,cpScrCoeffImUp,&(cp->cp_comm_state_pkg_up));
    if(cpLsda==1&&numStateDn>0){
    cp_transpose_fwd(coeffReDn,coeffImDn,pcoefFormDn,
                    cpScrCoeffReDn,cpScrCoeffImDn,&(cp->cp_comm_state_pkg_dn));
    cp_transpose_fwd(forceCoeffReDn,forceCoeffImDn,pforceCoefFormDn,
                    cpScrCoeffReDn,cpScrCoeffImDn,&(cp->cp_comm_state_pkg_dn));
    }/*endif*/
  }/*endif*/



/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genStoOrbital(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
		    CP *cp,int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS  *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR	       *cpscr	     = &(cp->cpscr);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  NEWTONINFO   *newtonInfo      = stodftInfo->newtonInfo;
  COMMUNICATE   *commCP	        = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS  *cpcoeffs_pos   = &(cp->cpcoeffs_pos[ip_now]);
  MPI_Comm commStates = commCP->comm_states;
  

  int iPoly,iState,iState2,iCoeff;
  int polynormLength = stodftInfo->polynormLength;
  int expanType = stodftInfo->expanType;
  int numProcStates = commCP->np_states;
  int myidState = commCP->myid_state;
  int cpLsda = cpopts->cp_lsda;
  int cpGGA  = cpopts->cp_gga;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *forceFormUp  = &(cpcoeffs_pos->ifcoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *forceFormDn  = &(cpcoeffs_pos->ifcoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *forceOrthUp  = &(cpcoeffs_pos->ifcoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);
  int *forceOrthDn  = &(cpcoeffs_pos->ifcoef_orth_dn);
  int numCoeffUpTot   = numStateUpProc*numCoeff;
  int numCoeffDnTot   = numStateDnProc*numCoeff;
  
  double energyMin,energyMax;
  double Smin,Smax; 
  double energyDiff;
  double scale;
  double length;

  double *sampPoint = newtonInfo->sampPoint;
  double *sampPointUnscale = newtonInfo->sampPointUnscale;
  double *coeffReUp = cpcoeffs_pos->cre_up;
  double *coeffImUp = cpcoeffs_pos->cim_up;
  double *coeffReDn = cpcoeffs_pos->cre_dn;
  double *coeffImDn = cpcoeffs_pos->cim_dn;
  double *scrCoeffReUp   = cpscr->cpscr_wave.cre_up;
  double *scrCoeffImUp   = cpscr->cpscr_wave.cim_up;
  double *rhoCoeffReUp   = cpscr->cpscr_rho.rhocr_up;
  double *rhoCoeffImUp   = cpscr->cpscr_rho.rhoci_up;
  double *rhoUp          = cpscr->cpscr_rho.rho_up;
  double *rhoCoeffReUpDensCpBox = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
  double *rhoCoeffImUpDensCpBox = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
  double *divRhoxUp       = cpscr->cpscr_grho.d_rhox_up;
  double *divRhoyUp       = cpscr->cpscr_grho.d_rhoy_up;
  double *divRhozUp       = cpscr->cpscr_grho.d_rhoz_up;
  double *d2RhoUp        = cpscr->cpscr_grho.d2_rho_up;
  double *scrCoeffReDn   = cpscr->cpscr_wave.cre_dn;
  double *scrCoeffImDn   = cpscr->cpscr_wave.cim_dn;
  double *rhoCoeffReDn   = cpscr->cpscr_rho.rhocr_dn;
  double *rhoCoeffImDn   = cpscr->cpscr_rho.rhoci_dn;
  double *rhoDn          = cpscr->cpscr_rho.rho_dn;
  double *rhoCoeffReDnDensCpBox = cpscr->cpscr_rho.rhocr_dn_dens_cp_box;
  double *rhoCoeffImDnDensCpBox = cpscr->cpscr_rho.rhoci_dn_dens_cp_box;
  double *divRhoxDn       = cpscr->cpscr_grho.d_rhox_dn;
  double *divRhoyDn       = cpscr->cpscr_grho.d_rhoy_dn;
  double *divRhozDn       = cpscr->cpscr_grho.d_rhoz_dn;
  double *d2RhoDn         = cpscr->cpscr_grho.d2_rho_dn;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  
/*==========================================================================*/
/* 0) Check the forms							    */

//debug
  genEigenOrb(cp,class,general_data,cpcoeffs_pos,clatoms_pos);

  double *coeffReUpBackup = (double*)cmalloc((numCoeffUpTot+1)*sizeof(double));
  double *coeffImUpBackup = (double*)cmalloc((numCoeffUpTot+1)*sizeof(double));
  double *randTrail = (double *)cmalloc((2*numCoeffUpTot+1)*sizeof(double));

  for(iState=0;iState<numStateUpProc;iState++){
    length = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      length += coeffReUp[iState*numCoeff+iCoeff]*coeffReUp[iState*numCoeff+iCoeff]+
		coeffImUp[iState*numCoeff+iCoeff]*coeffImUp[iState*numCoeff+iCoeff];
    }
    length *= 2.0;
    length += coeffReUp[iState*numCoeff+numCoeff]*coeffReUp[iState*numCoeff+numCoeff];
    length = sqrt(length);
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      coeffReUp[iState*numCoeff+iCoeff] /= length;
      coeffImUp[iState*numCoeff+iCoeff] /= length;
    }
  }

  for(iCoeff=1;iCoeff<=numCoeffUpTot;iCoeff++){
    coeffReUpBackup[iCoeff] = coeffReUp[iCoeff];
    coeffImUpBackup[iCoeff] = coeffImUp[iCoeff];
  }
#ifdef MKL_RANDOM
  VSLStreamStatePtr stream;
  int errcode;
  int seed = 1;
  errcode = vslNewStream(&stream,VSL_BRNG_MCG31,seed);
  errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,2*numCoeffUpTot,randTrail,randMin,randMax);
  for(iCoeff=1;iCoeff<=numCoeffTot;iCoeff++){
    coeffReUp[iCoeff] = randTrail[iCoeff-1];
  }
  for(iCoeff=1;iCoeff<=numCoeffTot;iCoeff++){
    coeffImUp[iCoeff] = randTrail[iCoeff-1+numCoeffUpTot];
  }
  for(iState=0;iState<numStateUpProc;iState++)coeffImUp[iState*numCoeff+numCoeff] = 0.0;//Keep everything real
#endif
#ifndef MKL_RANDOM
  //whatever random number is good, I'm using Gaussian in this case
  double seed = 45.154;
  int iseed;
  gaussran(2*numCoeffUpTot,&iseed,&iseed,&seed,randTrail);
  for(iCoeff=1;iCoeff<=numCoeffUpTot;iCoeff++){
    coeffReUp[iCoeff] = randTrail[iCoeff-1];
  }
  for(iCoeff=1;iCoeff<=numCoeffUpTot;iCoeff++){
    coeffImUp[iCoeff] = randTrail[iCoeff-1+numCoeffUpTot];
  }
  for(iState=0;iState<numStateUpProc;iState++)coeffImUp[iState*numCoeff+numCoeff] = 0.0;//Keep everything real
#endif

  //Normalize the trail wave function
  for(iState=0;iState<numStateUpProc;iState++){
    length = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      length += coeffReUp[iState*numCoeff+iCoeff]*coeffReUp[iState*numCoeff+iCoeff]+
		coeffImUp[iState*numCoeff+iCoeff]*coeffImUp[iState*numCoeff+iCoeff];
    }
    length *= 2.0;
    length += coeffReUp[iState*numCoeff+numCoeff]*coeffReUp[iState*numCoeff+numCoeff];
    length = sqrt(length);
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      coeffReUp[iState*numCoeff+iCoeff] /= length;
      coeffImUp[iState*numCoeff+iCoeff] /= length;
    }
  }
//enddebug

  if(numProcStates>1){
    if(*coefFormUp+*forceFormUp!=2){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP vectors are not in transposed form \n");
     printf("on state processor %d in min_CG_cp \n",myidState);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cpLsda==1){
     if(*coefFormDn+*forceFormDn!=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("on state processor %d in min_CG_cp \n",myidState);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/

/*======================================================================*/
/* 0.05) Check the approximations in the methods                        */

/*======================================================================*/
/* II) In parallel, transpose coefs back to normal form                 */

  if(numProcStates>1){
   cp_transpose_bck(coeffReUp,coeffImUp,coefFormUp,
                    scrCoeffReUp,scrCoeffImUp,&(cp->cp_comm_state_pkg_up));
   if(cpLsda==1&&numStateDnProc>0){
    cp_transpose_bck(coeffReDn,coeffImDn,coefFormDn,
                     scrCoeffReDn,scrCoeffImDn,&(cp->cp_comm_state_pkg_dn));
   }/*endif*/

  }/*endif*/

/*======================================================================*/
/* III) Set flags							*/

  *forceFormUp = 0;
  cpcoeffs_pos->ifcoef_orth_up = 1; // temp keep this flag for debug filter
  if(cpLsda==1&&numStateDnProc!=0){
    *forceFormDn = 0;
    cpcoeffs_pos->ifcoef_orth_dn = 1;
  }

/*======================================================================*/
/* IV) Calculate Emax and Emin                                          */

  if(myidState==0){
    genEnergyMax(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
    genEnergyMin(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
  }
  if(numProcStates>1){
    Bcast(&(stodftInfo->energyMin),1,MPI_DOUBLE,1,commStates);
    Bcast(&(stodftInfo->energyMax),1,MPI_DOUBLE,1,commStates);
  }
  energyMin = stodftInfo->energyMin;
  energyMax = stodftInfo->energyMax;
  printf("energyMax %lg energyMin %lg\n",energyMax,energyMin);
  stodftInfo->energyDiff = energyMax-energyMin;
  energyDiff = stodftInfo->energyDiff;
  stodftInfo->energyMean = 0.5*(energyMin+energyMax);
  
  if(expanType==2){
    Smin = newtonInfo->Smin;
    Smax = newtonInfo->Smax;
    scale = (Smax-Smin)/energyDiff;
    newtonInfo->scale = scale;
    for(iPoly=0;iPoly<polynormLength;iPoly++){
      sampPointUnscale[iPoly] = (sampPoint[iPoly]-Smin)/scale+energyMin;
      //printf("sampunscale %lg\n",sampPointUnscale[iPoly]);
    }
    genCoeffNewtonHermit(stodftInfo,stodftCoefPos);
  }

/*======================================================================*/
/* V) Filter the stochastic orbitals					*/

  //debug 

  switch(expanType){
    case 2:
      filterNewtonPolyHerm(cp,class,general_data,ip_now);
      break;
  }

//debug
  for(iState=0;iState<numStateUpProc;iState++){
    length = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      length += stoWfUpRe[0][iState*numCoeff+iCoeff]*stoWfUpRe[0][iState*numCoeff+iCoeff]+
                stoWfUpIm[0][iState*numCoeff+iCoeff]*stoWfUpIm[0][iState*numCoeff+iCoeff];
    }
    length *= 2.0;
    length += stoWfUpRe[0][iState*numCoeff+numCoeff]*stoWfUpRe[0][iState*numCoeff+numCoeff];
    length = sqrt(length);
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      stoWfUpRe[0][iState*numCoeff+iCoeff] /= length;
      stoWfUpIm[0][iState*numCoeff+iCoeff] /= length;
    }
  }
  double dot;
  for(iState=0;iState<numStateUpProc;iState++){
    printf("iState %i ",iState);
    for(iState2=5;iState2<6;iState2++){
      dot = 0.0;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	dot += stoWfUpRe[0][iState*numCoeff+iCoeff]*coeffReUpBackup[iState2*numCoeff+iCoeff]+
	       stoWfUpIm[0][iState*numCoeff+iCoeff]*coeffImUpBackup[iState2*numCoeff+iCoeff];
      }
      dot *= 2.0;
      dot += stoWfUpRe[0][iState*numCoeff+numCoeff]*coeffReUpBackup[iState2*numCoeff+numCoeff];
      printf("%.10lg ",1.0-fabs(dot));
    }
    printf("\n");
  }
//enddebug


/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genNoiseOrbital(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                    CP *cp,int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
/*-----------------------------------------------------------------------*/
  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR        *cpscr        = &(cp->cpscr);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
 
  double *zfft           =    cpscr->cpscr_wave.zfft;
  double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
  
  
 


}/*end routine*/
/*==========================================================================*/

