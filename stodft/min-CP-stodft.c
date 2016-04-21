/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NVE                                      */
/*                                                                          */
/* This subprogram integrates the system using Vel Verlet                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_integrate_cpmin_entry.h"
#include "../proto_defs/proto_integrate_cpmin_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_energy_ctrl_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void minStodft(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp,
               int ip_now)
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
/* 0) Check the forms                                                       */

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
/* 0.1) Set the GGA flags                                               */
  cpopts->cp_becke=0;
  cpopts->cp_pw91x=0;
  cpopts->cp_fila_1x=0;
  cpopts->cp_fila_2x=0;
  cpopts->cp_pbe_x=0;
  cpopts->cp_revpbe_x=0;
  cpopts->cp_rpbe_x=0;
  cpopts->cp_xpbe_x=0;
  cpopts->cp_brx89=0;
  cpopts->cp_brx2k=0;
  cpopts->cp_lyp=0;
  cpopts->cp_lypm1=0;
  cpopts->cp_pw91c=0;
  cpopts->cp_pbe_c=0;
  cpopts->cp_xpbe_c=0;
  cpopts->cp_tau1_c=0;
  cpopts->cp_debug_xc=0;
  if(cpGGA == 1){
    if(strcasecmp(ggaxType,"becke"   )==0){cpopts->cp_becke=1;}
    if(strcasecmp(ggaxType,"pw91x"   )==0){cpopts->cp_pw91x=1;}
    if(strcasecmp(ggaxType,"fila_1x" )==0){cpopts->cp_fila_1x=1;}
    if(strcasecmp(ggaxType,"fila_2x" )==0){cpopts->cp_fila_2x=1;}
    if(strcasecmp(ggaxType,"pbe_x"   )==0){cpopts->cp_pbe_x=1;}
    if(strcasecmp(ggaxType,"revpbe_x")==0){cpopts->cp_revpbe_x=1;}
    if(strcasecmp(ggaxType,"rpbe_x"  )==0){cpopts->cp_rpbe_x=1;}
    if(strcasecmp(ggaxType,"xpbe_x"  )==0){cpopts->cp_xpbe_x=1;}
    if(strcasecmp(ggaxType,"brx89"   )==0){cpopts->cp_brx89=1;}
    if(strcasecmp(ggaxType,"brx2k"   )==0){cpopts->cp_brx2k=1;}
    if(strcasecmp(ggacType,"lyp"     )==0){cpopts->cp_lyp=1;  }
    if(strcasecmp(ggacType,"lypm1"   )==0){cpopts->cp_lypm1=1;  }
    if(strcasecmp(ggacType,"pw91c"   )==0){cpopts->cp_pw91c=1;}
    if(strcasecmp(ggacType,"pbe_c"   )==0){cpopts->cp_pbe_c=1;}
    if(strcasecmp(ggacType,"xpbe_c"  )==0){cpopts->cp_xpbe_c=1;}
    if(strcasecmp(ggacType,"tau1_c"  )==0){cpopts->cp_tau1_c=1;}
    if(strcasecmp(ggacType,"debug97x")==0){cpopts->cp_debug_xc=1;}
  }/*endif*/

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
/* IV) Get the total density, (spin densities too if lsda)              */
/*       and necessary gradients of density for GGA calculations        */


  cp_rho_calc_hybrid(cpewald,cpscr,cpcoeffs_info,
                     ewald,cell,coeffReUp,coeffImUp,*coefFormUp,*coefOrthUp,
                     rhoCoeffReUp,rhoCoeffImUp,rhoUp,rhoCoeffReUpDensCpBox,
		     rhoCoeffImUpDensCpBox,divRhoxUp,divRhoyUp,
                     divRhozUp,d2RhoUp,numStateUpProc,numCoeff,
                     cpGGA,cpDualGridOptOn,numInterpPmeDual,communicate,
                     &(cp->cp_para_fft_pkg3d_lg),&(cp->cp_sclr_fft_pkg3d_lg),
                     &(cp->cp_para_fft_pkg3d_dens_cp_box),
                     &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
                     &(cp->cp_sclr_fft_pkg3d_sm));

  if(cpLsda==1&&numStateDnProc>0){
  cp_rho_calc_hybrid(cpewald,cpscr,cpcoeffs_info,
                     ewald,cell,coeffReDn,coeffImUp,*coefFormDn,*coefOrthDn,
                     rhoCoeffReDn,rhoCoeffImDn,rhoDn,rhoCoeffReUpDensCpBox,
		     rhoCoeffImDnDensCpBox,divRhoxDn,divRhoyDn,
		     divRhozDn,d2RhoDn,numStateDnProc,numCoeff,
		     cpGGA,cpDualGridOptOn,numInterpPmeDual,communicate,
                     &(cp->communicate),&(cp->cp_para_fft_pkg3d_lg),
                     &(cp->cp_sclr_fft_pkg3d_lg),
                     &(cp->cp_para_fft_pkg3d_dens_cp_box),
                     &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
                     &(cp->cp_sclr_fft_pkg3d_sm));
    for(i=1;i <= ncoef_l_proc;i++) {
      rhocr_up[i] += rhocr_dn[i];
      rhoci_up[i] += rhoci_dn[i];
    }/* endfor */
    if(cp_dual_grid_opt_on >= 1){
      for(i=1;i<= ncoef_l_dens_cp_box; i++){
        rhocr_up_dens_cp_box[i] += rhocr_dn_dens_cp_box[i];
        rhoci_up_dens_cp_box[i] += rhoci_dn_dens_cp_box[i];
      }/* endfor */
    } /* endif */
  }/* endif */

/*======================================================================*/
/* V) Calculate the non-local pseudopotential list			*/
  
  if(stodftInfo->vpsAtomListFlag==0||cp_dual_grid_opt_on >= 1){
    control_vps_atm_list(pseudo,cell,clatoms_pos,clatoms_info,
                         atommaps,ewd_scr,for_scr,cp_dual_grid_opt_on,
			 stodftInfo->vpsAtomListFlag);    
    stodftInfo->vpsAtomListFlag = 1;
  }

/*======================================================================*/
/* VI) Filter the stochastic orbitals					*/

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

