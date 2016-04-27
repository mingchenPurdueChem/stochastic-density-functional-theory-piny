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
  STODFTCOEFFPOS *stodftCoeffPos = cp->stodftCoeffPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  CPEWALD *cpewald              = &(cp->cpewald);
  CPSCR *cpscr                  = &(cp->cpscr);
  PSEUDO *pseudo                = &(cp->pseudo);
  COMMUNICATE *communicate      = &(cp->communicate);
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

  double tolEdgeDist 		= cpopts->tol_edge_dist;

  double *coeffReUp        = cpcoeffs_pos->cre_up;
  double *coeffImUp        = cpcoeffs_pos->cim_up;
  double *coeffReDn        = cpcoeffs_pos->cre_dn;
  double *coeffImDn        = cpcoeffs_pos->cim_dn;
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
    cp_transpose_bck(coeffReUp,coeffImUp,coefFormUp,
                    cpScrCoeffReUp,cpScrCoeffImUp,&(cp->cp_comm_state_pkg_up));
    if(cpLsda==1&&numStateDn>0){
      cp_transpose_bck(coeffReDn,coeffImDn,coefFormDn,
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
  calcRhoDeterm(class,bonded,general_data,cp,cpcoeffs_pos);

  //debug only
  genStoOrbital(class,bonded,general_data,cp,ipNow);

  
  //exit(0);
/*======================================================================*/
/* V) SCF loop						                */

  for(iScf=0;iScf<numScf;iScf++){

/*----------------------------------------------------------------------*/
/* i) converge chemical potential					*/
    while (/* Chemical potential does not converge*/){

/*----------------------------------------------------------------------*/
/* 1) converge chemical potential                                       */

      genStoOrbital(class,bonded,general_data,cp,ipNow);

/*----------------------------------------------------------------------*/
/* 2)  Get the total density, for each chemical potential and get       */
/*     total electron number for each chemical potential	        */

/*----------------------------------------------------------------------*/
/* 3)  Using interpolation, calculate the chemical potential w.r.t      */
/*     Correct electron number.						*/

    }//endwhile
/*----------------------------------------------------------------------*/
/* ii) Generate stochastic orbital for the correct chemical potential   */
/*     (or pick a closest one) and update the density			*/

  }//endfor iScf

/*======================================================================*/
/* VI) In parallel, transpose coefs and coef forces fwd                 */

  if(np_states>1){
    cp_transpose_fwd(cre_up,cim_up,icoef_form_up,
                    cpscr_cre_up,cpscr_cim_up,&(cp->cp_comm_state_pkg_up));
    cp_transpose_fwd(fcre_up,fcim_up,ifcoef_form_up,
                    cpscr_cre_up,cpscr_cim_up,&(cp->cp_comm_state_pkg_up));
    if( (cp_lsda==1) && (nstate_dn > 0) ){
     cp_transpose_fwd(cre_dn,cim_dn,icoef_form_dn,
                     cpscr_cre_dn,cpscr_cim_dn,&(cp->cp_comm_state_pkg_dn));
     cp_transpose_fwd(fcre_dn,fcim_dn,ifcoef_form_dn,
                     cpscr_cre_dn,cpscr_cim_dn,&(cp->cp_comm_state_pkg_dn));
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
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS  *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  ATOMMAPS     *atommaps     = &(class->atommaps); 
  FOR_SCR      *for_scr      = &(class->for_scr);
  EWD_SCR      *ewd_scr      = &(class->ewd_scr);
  EWALD        *ewald        = &(general_data->ewald);
  CELL         *cell         = &(general_data->cell);
  STAT_AVG     *stat_avg     = &(general_data->stat_avg);
  PTENS        *ptens        = &(general_data->ptens);
  SIMOPTS      *simopts      = &(general_data->simopts);
  CPOPTS       *cpopts       = &(cp->cpopts);
  PSEUDO       *pseudo       = &(cp->pseudo);
  CPEWALD      *cpewald      = &(cp->cpewald)
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoeffPos = cp->stodftCoeffPos;
  COMMUNICATE   *commCP	        = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS  *cpcoeffs_pos   = &(cp->cpcoeffs_pos[ip_now]);
  

  int expanType = stodftInfo->expanType;
  int numProcStates = commCP->np_states;
  int myidStates = commCP->myidStates;
  int cpLsda = cpopts->cp_lsda;
  int cpGGA  = cpopts->cp_gga;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int numInterpPmeDual = pseudo->n_interp_pme_dual;
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

  char *ggaxType  = pseudo->ggax_typ;
  char *ggacType  = pseudo->ggac_typ; 

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
  double *d2RhoDn        = cpscr->cpscr_grho.d2_rho_dn;
  
/*==========================================================================*/
/* 0) Check the forms							    */


  if(numProcStates>1){
    if(*coefFormUp+*forceFormUp!=2){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP vectors are not in transposed form \n");
     printf("on state processor %d in min_CG_cp \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cpLsda==1){
     if(*coefFormDn+*forceFormDn!=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("on state processor %d in min_CG_cp \n",myid_state);
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
/* V) Filter the stochastic orbitals					*/

  switch(expanType){
    case 2:
      filterNewtonPolyHerm(cp,ip_now,ewald,ewd_scr,cell,clatoms_info,
		 clatoms_pos,atommaps,stat_avg,ptens,simopts,for_scr);
      break;
    case 3:
      filterNewtonPolyNoHerm(cp,ip_now,ewald,ewd_scr,cell,clatoms_info,
                 clatoms_pos,atommaps,stat_avg,ptens,simopts,for_scr);
      break;
      
  }


/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

