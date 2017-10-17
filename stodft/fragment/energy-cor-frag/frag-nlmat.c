/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: frag-nlmat.c                                 */
/*                                                                          */
/* This routine calculate the kinectic energy, non-local pseudo-potential   */
/* energy and nuclei force correction.                                      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"

#include "../proto_defs/proto_frag_local.h"
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcNonLocalMatrix(CP *cp, CP *cpMini, CLASS *classMini, 
			GENERAL_DATA *generalDataMini)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
/**************************************************************************/
/* This function calculate non-local pseudopotential matrix w.r.t. frag-  */
/* -ment MO, as well as the force component.				  */
/* The mojority part of the code is copied from control_cp_eext_recip	  */
/* Right now only KB form is included.					  */
/**************************************************************************/
/*======================================================================*/
/* Local Variable declaration                                           */
#include "../typ_defs/typ_mask.h"
  
  CLATOMS_INFO *clatoms_info = &(classMini->clatoms_info);
  CLATOMS_POS *clatoms_pos = &(classMini->clatoms_pos[1]);
  CPCOEFFS_INFO *cpcoeffs_info = &(cpMini->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos = &(cpMini->cpcoeffs_pos[1]);
  CPEWALD *cpewald = &(cpMini->cpewald);
  CPSCR *cpscr = &(cpMini->cpscr);
  CPOPTS *cpopts = &(cpMini->cpopts);
  PSEUDO *pseudo = &(cpMini->pseudo);
  EWD_SCR *ewd_scr = &(classMini->ewd_scr);
  ATOMMAPS *atommaps = &(classMini->atommaps);
  CELL *cell = &(generalDataMini->cell);
  EWALD *ewald = &(generalDataMini->ewald);
  PTENS *ptens = &(generalDataMini->ptens);
  COMMUNICATE *communicate = &(cp->communicate);
  FOR_SCR *for_scr = &(classMini->for_scr);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cpMini->cp_para_fft_pkg3d_lg);
  STAT_AVG *stat_avg = &(generalDataMini->stat_avg);
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;

  int cp_dual_grid_opt = cpopts->cp_dual_grid_opt;
  int idual_switch;
  int i,j,iii,igh;
  int iState,jState;
  int iFrag = fragInfo->iFrag;
  int nlmtot,ntot_up,ntot_dn;
  int nl_max_kb,np_nlmax_kb,nl_max_gh,np_nlmax_gh,nl_max_all,np_nlmax_all;
  int numAtomCalc = fragInfo->numAtomFragVnlCalc[iFrag];
  double vrecip,cp_enl,cp_enl_gh;
  double cpu1,cpu2;
/*----------------------------------------------------------------------*/
/*         Local Pointer declarations                                   */

  int iperd     = cell->iperd;

 /*-------------------------*/
 /* Pressure local pointers */
  int cp_ptens  = cpopts->cp_ptens_calc;
  double *pvten = ptens->pvten_tmp;

 /*------------------------------------------*/
 /* Non-local pseudopotential local pointers */

  int *np_nl            = pseudo->np_nl;
  int *np_nl_gh         = pseudo->np_nl_gh;

  int n_ang_max         = pseudo->n_ang_max;
  int n_ang_max_kb      = pseudo->n_ang_max_kb;
  int n_ang_max_gh      = pseudo->n_ang_max_gh;

  int n_rad_max         = pseudo->n_rad_max;
  double *vnlreal_up    = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up    = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn    = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn    = cpscr->cpscr_nonloc.vnlim_dn;

  double *dvnlreal_x_up = cpscr->cpscr_nonloc.dvnlre_x_up;
  double *dvnlreal_y_up = cpscr->cpscr_nonloc.dvnlre_y_up;
  double *dvnlreal_z_up = cpscr->cpscr_nonloc.dvnlre_z_up;
  double *dvnlimag_x_up = cpscr->cpscr_nonloc.dvnlim_x_up;
  double *dvnlimag_y_up = cpscr->cpscr_nonloc.dvnlim_y_up;
  double *dvnlimag_z_up = cpscr->cpscr_nonloc.dvnlim_z_up;
  double *dvnlreal_x_dn = cpscr->cpscr_nonloc.dvnlre_x_dn;
  double *dvnlreal_y_dn = cpscr->cpscr_nonloc.dvnlre_y_dn;
  double *dvnlreal_z_dn = cpscr->cpscr_nonloc.dvnlre_z_dn;
  double *dvnlimag_x_dn = cpscr->cpscr_nonloc.dvnlim_x_dn;
  double *dvnlimag_y_dn = cpscr->cpscr_nonloc.dvnlim_y_dn;
  double *dvnlimag_z_dn = cpscr->cpscr_nonloc.dvnlim_z_dn;

  double *dvnlreal_gxgx_up = cpscr->cpscr_nonloc.dvnlre_gxgx_up;
  double *dvnlimag_gxgx_up = cpscr->cpscr_nonloc.dvnlim_gxgx_up;
  double *dvnlreal_gygy_up = cpscr->cpscr_nonloc.dvnlre_gygy_up;
  double *dvnlimag_gygy_up = cpscr->cpscr_nonloc.dvnlim_gygy_up;
  double *dvnlreal_gzgz_up = cpscr->cpscr_nonloc.dvnlre_gzgz_up;
  double *dvnlimag_gzgz_up = cpscr->cpscr_nonloc.dvnlim_gzgz_up;
  double *dvnlreal_gxgy_up = cpscr->cpscr_nonloc.dvnlre_gxgy_up;
  double *dvnlimag_gxgy_up = cpscr->cpscr_nonloc.dvnlim_gxgy_up;
  double *dvnlreal_gygz_up = cpscr->cpscr_nonloc.dvnlre_gygz_up;
  double *dvnlimag_gygz_up = cpscr->cpscr_nonloc.dvnlim_gygz_up;
  double *dvnlreal_gxgz_up = cpscr->cpscr_nonloc.dvnlre_gxgz_up;
  double *dvnlimag_gxgz_up = cpscr->cpscr_nonloc.dvnlim_gxgz_up;
  double *dvnlreal_gxgx_dn = cpscr->cpscr_nonloc.dvnlre_gxgx_dn;
  double *dvnlimag_gxgx_dn = cpscr->cpscr_nonloc.dvnlim_gxgx_dn;
  double *dvnlreal_gygy_dn = cpscr->cpscr_nonloc.dvnlre_gygy_dn;
  double *dvnlimag_gygy_dn = cpscr->cpscr_nonloc.dvnlim_gygy_dn;
  double *dvnlreal_gzgz_dn = cpscr->cpscr_nonloc.dvnlre_gzgz_dn;
  double *dvnlimag_gzgz_dn = cpscr->cpscr_nonloc.dvnlim_gzgz_dn;
  double *dvnlreal_gxgy_dn = cpscr->cpscr_nonloc.dvnlre_gxgy_dn;
  double *dvnlimag_gxgy_dn = cpscr->cpscr_nonloc.dvnlim_gxgy_dn;
  double *dvnlreal_gygz_dn = cpscr->cpscr_nonloc.dvnlre_gygz_dn;
  double *dvnlimag_gygz_dn = cpscr->cpscr_nonloc.dvnlim_gygz_dn;
  double *dvnlreal_gxgz_dn = cpscr->cpscr_nonloc.dvnlre_gxgz_dn;
  double *dvnlimag_gxgz_dn = cpscr->cpscr_nonloc.dvnlim_gxgz_dn;

 /*----------------------*/
 /* Atom local pointers */
  int natm_tot             = clatoms_info->natm_tot;
  double *fx_tmp           = ewd_scr->fx;
  double *fy_tmp           = ewd_scr->fy;
  double *fz_tmp           = ewd_scr->fz;
  double *fx               = clatoms_pos->fx;
  double *fy               = clatoms_pos->fy;
  double *fz               = clatoms_pos->fz;
  double *x                = clatoms_pos->x;
  double *y                = clatoms_pos->y;
  double *z                = clatoms_pos->z;
  double *q                = clatoms_info->q;

  int hess_size;
  int hess_calc = clatoms_info->hess_calc;

  double *hess_xx         = clatoms_pos->hess_xx;
  double *hess_xy         = clatoms_pos->hess_xy;
  double *hess_xz         = clatoms_pos->hess_xz;
  double *hess_yy         = clatoms_pos->hess_yy;
  double *hess_yz         = clatoms_pos->hess_yz;
  double *hess_zz         = clatoms_pos->hess_zz;

 /*-----------------------------------*/
 /* Cp Option and form local pointers */
  int cp_lsda   = cpopts->cp_lsda;
  int icoef_orth_up  = cpcoeffs_pos->icoef_orth_up;
  int icoef_form_up  = cpcoeffs_pos->icoef_form_up;
  int ifcoef_form_up = cpcoeffs_pos->ifcoef_form_up;
  int icoef_orth_dn  = cpcoeffs_pos->icoef_orth_dn;
  int icoef_form_dn  = cpcoeffs_pos->icoef_form_dn;
  int ifcoef_form_dn = cpcoeffs_pos->ifcoef_form_dn;

 /*------------------------------*/
 /* Wave function local pointers */
  int nstate_up = cpcoeffs_info->nstate_up_proc;
  int nstate_dn = cpcoeffs_info->nstate_dn_proc;

 /*------------------------------*/
 /* Communciation local pointers */
  MPI_Comm comm_states = communicate->comm_states;
  int myid_state  = communicate->myid_state;
  int np_states   = communicate->np_states;
  int np_forc     = communicate->np_forc;

  double *vrecip_ret = &(stat_avg->vrecip);
  double *cp_enl_ret = &(stat_avg->cp_enl);

/*======================================================================*/
/* 0) Check the forms                                                   */

  if(icoef_orth_up!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The Up coefficients must be in orthogonal form    \n");
    printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }//endif icoef_orth_up
  if(cp_lsda==1 && nstate_dn != 0){
    if(icoef_orth_dn!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The dn coefficients must be in orthogonal form    \n");
      printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }//endif icoef_orth_dn
  }//endif cp_lsda

  if(np_states>1){
    if((icoef_form_up+ifcoef_form_up)!=0){
      //printf("icoef_form_up %i ifcoef_form_up %i\n",icoef_form_up,ifcoef_form_up);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The up coefs and coef forces must not be in transposed form\n");
      printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }//endif icoef_form_up+ifcoef_form_up
    if(cp_lsda==1 && nstate_dn != 0){
      if((icoef_form_dn+ifcoef_form_dn)!=0){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("The dn coefs and coef forces must not be in transposed form\n");
        printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);exit(1);
      }//endif icoef_form_dn+ifcoef_form_dn
    }//endif cp_lsda
  }//endif np_stastes

/*======================================================================*/
/* II) Malloc a hessian scratch vector if necessary                     */

  if(hess_calc==3&&np_states>1)hess_size = natm_tot*natm_tot;

/*======================================================================*/
/* III) Determine the maximum open non-local angular momentum channel   */
/*      for Kleinman-Bylander and Goedecker pseudo potentials           */

  nl_max_kb = -1;
  for(i=1;i<=(n_ang_max_kb+1);i++){
    if(np_nl[i]>0){nl_max_kb=i-1;}
  }//endfor i
  nl_max_gh = -1;
  for(i=1;i<=(n_ang_max_gh+1);i++){
    if(np_nl_gh[i]>0){nl_max_gh=i-1;}
  }//endfor i

/*======================================================================*/
/* IV) Determine the maximum number of atoms in any                     */
/*       open angular momentum channel                                  */

  np_nlmax_kb = 1;
  for(i=1;i<=(nl_max_kb+1);i++){
   np_nlmax_kb = MAX(np_nlmax_kb,np_nl[i]);
  }//endfor i
  np_nlmax_gh = 1;
  for(i=1;i<=(nl_max_gh+1);i++){
   np_nlmax_gh = MAX(np_nlmax_gh,np_nl_gh[i]);
  }//endfor
  np_nlmax_all = (np_nlmax_gh>np_nlmax_kb ? np_nlmax_gh:np_nlmax_kb);

/*======================================================================*/
/* V) Zero the non-local tensors                                        */

  nl_max_all = (nl_max_gh>nl_max_kb ? nl_max_gh:nl_max_kb);
  nlmtot = (nl_max_all+1)*(nl_max_all+1);
  ntot_up = nstate_up*np_nlmax_all*nlmtot*n_rad_max;
  ntot_dn = 0;

  if(nl_max_all>=0){
    for(i=1;i<=ntot_up;i++){
      vnlreal_up[i] = 0.0;
      vnlimag_up[i] = 0.0;
      dvnlreal_x_up[i] = 0.0;
      dvnlreal_y_up[i] = 0.0;
      dvnlreal_z_up[i] = 0.0;
      dvnlimag_x_up[i] = 0.0;
      dvnlimag_y_up[i] = 0.0;
      dvnlimag_z_up[i] = 0.0;
    }//endfor
    if(cp_ptens==1||hess_calc==3){
      for(i=1;i<=ntot_up;i++){
        dvnlreal_gxgx_up[i] = 0.0;
        dvnlreal_gzgz_up[i] = 0.0;
        dvnlreal_gygy_up[i] = 0.0;
        dvnlreal_gxgy_up[i] = 0.0;
        dvnlreal_gxgz_up[i] = 0.0;
        dvnlreal_gygz_up[i] = 0.0;
        dvnlimag_gxgx_up[i] = 0.0;
        dvnlimag_gxgy_up[i] = 0.0;
        dvnlimag_gygy_up[i] = 0.0;
        dvnlimag_gxgz_up[i] = 0.0;
        dvnlimag_gygz_up[i] = 0.0;
        dvnlimag_gzgz_up[i] = 0.0;
      }//endfor i
    }//endif:ptens
    if(cp_lsda==1){
      ntot_dn = nstate_dn*np_nlmax_all*nlmtot*n_rad_max;
      for(i=1;i<=ntot_dn;i++){
        vnlreal_dn[i]    = 0.0;
        vnlimag_dn[i]    = 0.0;
        dvnlreal_x_dn[i] = 0.0;
        dvnlreal_y_dn[i] = 0.0;
        dvnlreal_z_dn[i] = 0.0;
        dvnlimag_x_dn[i] = 0.0;
        dvnlimag_y_dn[i] = 0.0;
        dvnlimag_z_dn[i] = 0.0;
      }//endfor
      if(cp_ptens==1 || hess_calc == 3){
        for(i=1;i<=ntot_dn;i++){
 	  dvnlreal_gxgx_dn[i] = 0.0;
 	  dvnlreal_gxgy_dn[i] = 0.0;
	  dvnlreal_gxgz_dn[i] = 0.0;
	  dvnlreal_gygy_dn[i] = 0.0;
	  dvnlreal_gygz_dn[i] = 0.0;
	  dvnlreal_gzgz_dn[i] = 0.0;
	  dvnlimag_gxgx_dn[i] = 0.0;
	  dvnlimag_gxgy_dn[i] = 0.0;
	  dvnlimag_gxgz_dn[i] = 0.0;
	  dvnlimag_gygy_dn[i] = 0.0;
	  dvnlimag_gygz_dn[i] = 0.0;
	  dvnlimag_gzgz_dn[i] = 0.0;
        }//endfor i
      }//endif:ptens
    }//endif:lsda
  }//endif : non-local potential on

  // initialize matrix and force
  for(iState=0;iState<nstate_up;iState++){
    for(jState=0;jState<nstate_up;jState++){
      fragInfo->vnlMatrixUp[iFrag][iState*nstate_up+jState] = 0.0;
    }
  }
  if(cp_lsda==1 && nstate_dn != 0){
    for(iState=0;iState<nstate_dn;iState++){
      for(jState=0;jState<nstate_up;jState++){
	fragInfo->vnlMatrixDn[iFrag][iState*nstate_dn+jState] = 0.0;
      }
    }
  }

  for(i=0;i<numAtomCalc;i++){
    for(iState=0;iState<nstate_up;iState++){
      for(jState=0;jState<nstate_up;jState++){
	fragInfo->vnlFxMatrixUp[iFrag][i*nstate_up*nstate_up+iState*nstate_up+jState] = 0.0;
	fragInfo->vnlFyMatrixUp[iFrag][i*nstate_up*nstate_up+iState*nstate_up+jState] = 0.0;
	fragInfo->vnlFzMatrixUp[iFrag][i*nstate_up*nstate_up+iState*nstate_up+jState] = 0.0;
      }
    }
    if(cp_lsda==1 && nstate_dn != 0){
      for(i=0;i<numAtomCalc;i++){
	for(iState=0;iState<nstate_dn;iState++){
	  for(jState=0;jState<nstate_dn;jState++){
	    fragInfo->vnlFxMatrixDn[iFrag][i*nstate_dn*nstate_dn+iState*nstate_dn+jState] = 0.0;
	    fragInfo->vnlFyMatrixDn[iFrag][i*nstate_dn*nstate_dn+iState*nstate_dn+jState] = 0.0;
	    fragInfo->vnlFzMatrixDn[iFrag][i*nstate_dn*nstate_dn+iState*nstate_dn+jState] = 0.0;
	  }
	}
      }
    }
  }

  for(i=1;i<=natm_tot;i++){
    fx[i] = 0.0;
    fy[i] = 0.0;
    fz[i] = 0.0;
  }

/*======================================================================*/
/* VI) Perform the ewald sum/ cp local pseudopotential calculation      */
/*     I don't know whether this parted is needed or not. Keep it just  */
/*     to be safe.						        */

  // We currently don't need this for nl calculation

  /*
  idual_switch = 0; // cp_dual_grid_opt<=1 : get vrecip vext on dense grid 
                    // cp_dual_grid_opt==2 : get vrecip vext on sparse grid
  vrecip = 0.0;
  if(cpscr->cpscr_atom_pme.pme_on==1&&cp_dual_grid_opt==2){
    control_ewd_loc_pme(clatoms_info,clatoms_pos,cell,ptens,ewald,cpewald,
                        cpscr,pseudo,ewd_scr,cpopts,atommaps,
                        &vrecip, &(cpcoeffs_info->pseud_hess_loc),communicate,
                        for_scr,cp_dual_grid_opt,idual_switch,
                        cp_para_fft_pkg3d_lg);
  }
  else{
    control_ewd_loc(clatoms_info,clatoms_pos,cell,ptens,ewald,cpewald,
                  cpscr,pseudo,ewd_scr,cpopts,atommaps,
                  &vrecip,&(cpcoeffs_info->pseud_hess_loc),communicate,
                  for_scr,cp_dual_grid_opt,idual_switch);
  }//endif
  if(cp_dual_grid_opt==2){
    idual_switch = 1; //get vext on small dense grid
    control_ewd_loc(clatoms_info,clatoms_pos,cell,ptens,ewald,cpewald,
                    cpscr,pseudo,ewd_scr,cpopts,atommaps,
                    &vrecip,&(cpcoeffs_info->pseud_hess_loc),communicate,
                    for_scr,cp_dual_grid_opt,idual_switch);
  }//endif cp_dual_grid_opt
  */

/*======================================================================*/
/* VII) Get the nl pe, pvten and particle forces then the coef forces   */
  cp_enl    = 0.0;
  cp_enl_gh = 0.0;

  if(nl_max_all>=0&&n_rad_max>1){
     non_loc_chng_ord(clatoms_pos,clatoms_info,atommaps, pseudo,ewd_scr,for_scr,1);
  }//endif nl_max_all n_rad_max
/*-------------------------------------------------------------------------*/
/* A) KB/Goedecker NLs */

  if(nl_max_kb>=0&&(ntot_up+ntot_dn)>0){
    control_ewd_non_loc(clatoms_info,clatoms_pos,cpcoeffs_info,cpcoeffs_pos,
                        cell,ptens,cpewald,cpscr,pseudo,ewd_scr,
                        cpopts,atommaps,communicate,for_scr);
  }
  else get_ak2_sm(cpewald,cell);//endif nl_max_kb ntot_up+ntot_dn

  //printf("nl_max_kb %i ntot_up %i pseudo->np_nonloc_cp_box_kb %i\n",nl_max_kb,ntot_up,pseudo->np_nonloc_cp_box_kb);

  if((nl_max_kb >= 0)&&((ntot_up+ntot_dn)>0)&&(pseudo->np_nonloc_cp_box_kb>0) ){
    getnlPotPvFatmFrag(clatoms_info,clatoms_pos,cell,cpcoeffs_info,
    		       cpscr,ewd_scr,cpopts,pseudo,atommaps,fragInfo,&cp_enl,
    		       np_nlmax_kb,pvten);
    /*
    getnl_fcoef(clatoms_info,clatoms_pos,cpcoeffs_info,cpcoeffs_pos,
                  cpscr,ewd_scr,cpopts,pseudo,cpewald,atommaps,
                  cell,np_nlmax_kb,pvten,for_scr);
    */
  }//endif
  //printf("111111 fx[1] %lg fx[2] %lg fx[3] %lg\n",fx[1],fx[2],fx[3]);

/*======================================================================*/
/* VIII) Assign the potential energy                                    */

  *vrecip_ret += vrecip;
  *cp_enl_ret += cp_enl;

/*======================================================================*/
/* IX) Reduce particle forces if necessary                              */

  /*
  if(np_states>1 && np_forc == 1){
    reduce_cp_atm_forc(natm_tot,fx,fy,fz,fx_tmp,fy_tmp,fz_tmp,
                       comm_states,myid_state);
  }// endif:npstates 
  */

/*======================================================================*/
/* X) Reduce particle hessian if necessary                              */

  /*
  if(hess_calc == 3 && np_states>1){
    reduce_cp_hess_stuff(hess_xx,hess_yy,hess_zz,
                         hess_xy,hess_xz,hess_yz,hess_size,
                         myid_state,comm_states);
  }//endif
  */

/*======================================================================*/
    }/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void getnlPotPvFatmFrag(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
			CELL *cell,CPCOEFFS_INFO *cpcoeffs_info,CPSCR *cpscr,
			EWD_SCR *ewd_scr,CPOPTS *cpopts,
			PSEUDO *pseudo,ATOMMAPS *atommaps,FRAGINFO *fragInfo,
			double *cp_enl_ret, int np_nlmax,double *pvten)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
  int i,l,m,i_shift,ipart,is,iii,ioff,lp1;
  int iatm;
  int ind_loc,nl_max,irad,jrad;
  int nl_chan_max;
  double rvol_cp,vol_cp,cp_enl;
  double p11,p22,p33,p12,p13,p23;

/* Local pointers */
  int npart                = clatoms_info->natm_tot;
  double *fx               = clatoms_pos->fx;
  double *fy               = clatoms_pos->fy;
  double *fz               = clatoms_pos->fz;
  double *hess_xx          = clatoms_pos->hess_xx;
  double *hess_xy          = clatoms_pos->hess_xy;
  double *hess_xz          = clatoms_pos->hess_xz;
  double *hess_yy          = clatoms_pos->hess_yy;
  double *hess_yz          = clatoms_pos->hess_yz;
  double *hess_zz          = clatoms_pos->hess_zz;
  int natm_typ             = atommaps->natm_typ;
  int *iatm_typ            = atommaps->iatm_atm_typ;
  int *iatm_typ_nl         = atommaps->iatm_atm_typ_nl;
  double *hmat_cp          = cell->hmat_cp;

  int cp_ptens             = cpopts->cp_ptens_calc;
  int cp_lsda              = cpopts->cp_lsda;
  int atm_hess_calc        = clatoms_info->hess_calc;
  int nstate_up            = cpcoeffs_info->nstate_up_proc;
  int nstate_dn            = cpcoeffs_info->nstate_dn_proc;

  int np_nonloc_cp_box_kb  = pseudo->np_nonloc_cp_box_kb;
  double *vpsnorm          = pseudo->vpsnorm;
  int n_ang_max            = pseudo->n_ang_max;
  int n_ang_max_kb         = pseudo->n_ang_max_kb;
  int n_rad_max            = pseudo->n_rad_max;
  int *loc_opt             = pseudo->loc_opt;
  int *ip_nl               = pseudo->ip_nl;
  int *ip_nl_rev           = pseudo->ip_nl_rev;
  int *np_nl               = pseudo->np_nl;
  int *nrad_max_l          = pseudo->nrad_max_l;
  int **np_nl_rad_str      = pseudo->np_nl_rad_str;

  double *vnlreal_up       = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up       = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn       = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn       = cpscr->cpscr_nonloc.vnlim_dn;
  double *dvnlreal_x_up    = cpscr->cpscr_nonloc.dvnlre_x_up;
  double *dvnlreal_y_up    = cpscr->cpscr_nonloc.dvnlre_y_up;
  double *dvnlreal_z_up    = cpscr->cpscr_nonloc.dvnlre_z_up;
  double *dvnlimag_x_up    = cpscr->cpscr_nonloc.dvnlim_x_up;
  double *dvnlimag_y_up    = cpscr->cpscr_nonloc.dvnlim_y_up;
  double *dvnlimag_z_up    = cpscr->cpscr_nonloc.dvnlim_z_up;
  double *dvnlreal_x_dn    = cpscr->cpscr_nonloc.dvnlre_x_dn;
  double *dvnlreal_y_dn    = cpscr->cpscr_nonloc.dvnlre_y_dn;
  double *dvnlreal_z_dn    = cpscr->cpscr_nonloc.dvnlre_z_dn;
  double *dvnlimag_x_dn    = cpscr->cpscr_nonloc.dvnlim_x_dn;
  double *dvnlimag_y_dn    = cpscr->cpscr_nonloc.dvnlim_y_dn;
  double *dvnlimag_z_dn    = cpscr->cpscr_nonloc.dvnlim_z_dn;
  double *dvnlreal_gxgx_up = cpscr->cpscr_nonloc.dvnlre_gxgx_up;
  double *dvnlimag_gxgx_up = cpscr->cpscr_nonloc.dvnlim_gxgx_up;
  double *dvnlreal_gygy_up = cpscr->cpscr_nonloc.dvnlre_gygy_up;
  double *dvnlimag_gygy_up = cpscr->cpscr_nonloc.dvnlim_gygy_up;
  double *dvnlreal_gzgz_up = cpscr->cpscr_nonloc.dvnlre_gzgz_up;
  double *dvnlimag_gzgz_up = cpscr->cpscr_nonloc.dvnlim_gzgz_up;
  double *dvnlreal_gxgy_up = cpscr->cpscr_nonloc.dvnlre_gxgy_up;
  double *dvnlimag_gxgy_up = cpscr->cpscr_nonloc.dvnlim_gxgy_up;
  double *dvnlreal_gygz_up = cpscr->cpscr_nonloc.dvnlre_gygz_up;
  double *dvnlimag_gygz_up = cpscr->cpscr_nonloc.dvnlim_gygz_up;
  double *dvnlreal_gxgz_up = cpscr->cpscr_nonloc.dvnlre_gxgz_up;
  double *dvnlimag_gxgz_up = cpscr->cpscr_nonloc.dvnlim_gxgz_up;
  double *dvnlreal_gxgx_dn = cpscr->cpscr_nonloc.dvnlre_gxgx_dn;
  double *dvnlimag_gxgx_dn = cpscr->cpscr_nonloc.dvnlim_gxgx_dn;
  double *dvnlreal_gygy_dn = cpscr->cpscr_nonloc.dvnlre_gygy_dn;
  double *dvnlimag_gygy_dn = cpscr->cpscr_nonloc.dvnlim_gygy_dn;
  double *dvnlreal_gzgz_dn = cpscr->cpscr_nonloc.dvnlre_gzgz_dn;
  double *dvnlimag_gzgz_dn = cpscr->cpscr_nonloc.dvnlim_gzgz_dn;
  double *dvnlreal_gxgy_dn = cpscr->cpscr_nonloc.dvnlre_gxgy_dn;
  double *dvnlimag_gxgy_dn = cpscr->cpscr_nonloc.dvnlim_gxgy_dn;
  double *dvnlreal_gygz_dn = cpscr->cpscr_nonloc.dvnlre_gygz_dn;
  double *dvnlimag_gygz_dn = cpscr->cpscr_nonloc.dvnlim_gygz_dn;
  double *dvnlreal_gxgz_dn = cpscr->cpscr_nonloc.dvnlre_gxgz_dn;
  double *dvnlimag_gxgz_dn = cpscr->cpscr_nonloc.dvnlim_gxgz_dn;

  double *fxtemp           = ewd_scr->fx;
  double *fytemp           = ewd_scr->fy;
  double *fztemp           = ewd_scr->fz;
  double *vscr             = ewd_scr->fx2;
  double *vnorm            = ewd_scr->fy2;
  double *vnorm_now        = ewd_scr->fz2;

  // Fragment related
  int iFrag = fragInfo->iFrag;
  int *atomVnlCalcMapInv = fragInfo->atomFragVnlCalcMapInv[iFrag];
  int *atomVnlCalcFlag = fragInfo->atomFragVnlCalcFlag[iFrag];

  double *vnlFxMatrixUp,*vnlFxMatrixDn;
  double *vnlFyMatrixUp,*vnlFyMatrixDn;
  double *vnlFzMatrixUp,*vnlFzMatrixDn;
  double *vnlMatrixUp,*vnlMatrixDn;

  vnlFxMatrixUp = fragInfo->vnlFxMatrixUp[iFrag];
  vnlFyMatrixUp = fragInfo->vnlFyMatrixUp[iFrag];
  vnlFzMatrixUp = fragInfo->vnlFzMatrixUp[iFrag];
  vnlMatrixUp = fragInfo->vnlMatrixUp[iFrag];
  if(cp_lsda==1){
    vnlFxMatrixDn = fragInfo->vnlFxMatrixDn[iFrag];
    vnlFyMatrixDn = fragInfo->vnlFyMatrixDn[iFrag];
    vnlFzMatrixDn = fragInfo->vnlFzMatrixDn[iFrag];
    vnlMatrixDn = fragInfo->vnlMatrixDn[iFrag];
  }//endif
  

/*======================================================================*/
/* I) Useful constants                                                  */

  vol_cp   = getdeth(hmat_cp);
  rvol_cp  = 1.0/vol_cp;
  cp_enl   = 0.0;
  nl_max   = -1;
  for(i=1;i<=(n_ang_max_kb+1);i++){
    if(np_nl[i]>0){nl_max=i-1;}
  }/*endfor*/
  nl_chan_max = (nl_max+1)*(nl_max+1);

/*======================================================================*/
/* II) Loop over the open channels, the states and get the nl potent,   */
/*     pvten and  particle forces                                       */

  for(l=0;l<=nl_max;l++){
    lp1 = l+1;
    if(np_nl[lp1]>0){
      for(irad=1;irad<=nrad_max_l[lp1];irad++){
        for(jrad=irad;jrad<=nrad_max_l[lp1];jrad++){
/*-----------------------------------------------------------------------*/
/* i) Get the normalization scaled by the volume                        */
	  get_vpsnorm(vscr,vpsnorm,vnorm,iatm_typ_nl,natm_typ,np_nonloc_cp_box_kb,l,
		      n_ang_max,irad,jrad,n_rad_max);
	  i_shift = l*npart;
	  for(ipart=np_nl_rad_str[lp1][jrad];ipart<=np_nl[lp1];ipart++){
	    vnorm_now[ipart] = vnorm[ip_nl_rev[(ipart+i_shift)]]*rvol_cp;
	  }/*endfor*/
/*-----------------------------------------------------------------------*/
/* ii) Sum the contributions over the 2l+1 directions and the states    */
	   sumnlPotPvFatmHessFrag(npart,nstate_up,np_nlmax,nl_chan_max,
				  np_nl[lp1],l,np_nl_rad_str[lp1][jrad],
				  irad,jrad,
				  ip_nl,vnorm_now,vnlreal_up,vnlimag_up,
				  dvnlreal_gxgx_up,dvnlimag_gxgx_up,
				  dvnlreal_gygy_up,dvnlimag_gygy_up,
				  dvnlreal_gzgz_up,dvnlimag_gzgz_up,
				  dvnlreal_gxgy_up,dvnlimag_gxgy_up,
				  dvnlreal_gxgz_up,dvnlimag_gxgz_up,
				  dvnlreal_gygz_up,dvnlimag_gygz_up,
				  dvnlreal_x_up,dvnlimag_x_up,
				  dvnlreal_y_up,dvnlimag_y_up,
				  dvnlreal_z_up,dvnlimag_z_up,
				  fx,fy,fz,fxtemp,fytemp,fztemp,
				  hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
				  atm_hess_calc,cp_ptens,pvten,&cp_enl,
				  vnlFxMatrixUp,vnlFyMatrixUp,vnlFzMatrixUp,vnlMatrixUp,
				  atomVnlCalcMapInv,atomVnlCalcFlag);
	   if(cp_lsda==1){
	     sumnlPotPvFatmHessFrag(npart,nstate_dn,np_nlmax,nl_chan_max,
				    np_nl[lp1],l,np_nl_rad_str[lp1][jrad],
				    irad,jrad,
				    ip_nl,vnorm_now,vnlreal_dn,vnlimag_dn,
				    dvnlreal_gxgx_dn,dvnlimag_gxgx_dn,
				    dvnlreal_gygy_dn,dvnlimag_gygy_dn,
				    dvnlreal_gzgz_dn,dvnlimag_gzgz_dn,
				    dvnlreal_gxgy_dn,dvnlimag_gxgy_dn,
				    dvnlreal_gxgz_dn,dvnlimag_gxgz_dn,
				    dvnlreal_gygz_dn,dvnlimag_gygz_dn,
				    dvnlreal_x_dn,dvnlimag_x_dn,
				    dvnlreal_y_dn,dvnlimag_y_dn,
				    dvnlreal_z_dn,dvnlimag_z_dn,
				    fx,fy,fz,fxtemp,fytemp,fztemp,
				    hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
				    atm_hess_calc,cp_ptens,pvten,&cp_enl,
				    vnlFxMatrixDn,vnlFyMatrixDn,vnlFzMatrixDn,vnlMatrixDn,
				    atomVnlCalcMapInv,atomVnlCalcFlag);
	   }//endif
	}//endfor jrad
      }//endfor: radial channels irad
    }//endif: l channel open
  }//endfor: l channels

/*======================================================================*/
/* III) Assign the non-local energy  and add it to the pvten            */

  *cp_enl_ret = cp_enl;
  if(cp_ptens==1){
    pvten[1] += cp_enl;
    pvten[5] += cp_enl;
    pvten[9] += cp_enl;
  }//endif

/*======================================================================*/
  }/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void sumnlPotPvFatmHessFrag(int npart,int nstate,int np_nlmax,
				int nl_chan_max,int np_nl,int l,int np_nl_rad_str,
				int irad, int jrad,
				int *ip_nl,double *vnorm_now,
				double *vnlreal,double *vnlimag,
				double *dvnlreal_gxgx,double *dvnlimag_gxgx,
				double *dvnlreal_gygy,double *dvnlimag_gygy,
				double *dvnlreal_gzgz,double *dvnlimag_gzgz,
				double *dvnlreal_gxgy,double *dvnlimag_gxgy,
				double *dvnlreal_gxgz,double *dvnlimag_gxgz,
				double *dvnlreal_gygz,double *dvnlimag_gygz,
				double *dvnlreal_x,double *dvnlimag_x,
				double *dvnlreal_y,double *dvnlimag_y,
				double *dvnlreal_z,double *dvnlimag_z,
				double *fx,double *fy,double *fz,
				double *fxtemp,double *fytemp,double *fztemp,
				double *hess_xx,double *hess_xy,double *hess_xz,
				double *hess_yy,double *hess_yz,double *hess_zz,
				int atm_hess_calc,
				int cp_ptens,double *pvten, double *cp_enl_ret,
				double *vnlFxMatrix,double *vnlFyMatrix,
				double *vnlFzMatrix,double *vnlMatrix,
				int *atomVnlCalcMapInv,int *atomVnlCalcFlag)
/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
  int i,m,i_shift,ipart,is,js,iii,ioff,ltemp;
  int ind_loc,ind_lm;
  int ioff_i,ioff_j,ind_loc_i,ind_loc_j;
  int ioff_i_is,ioff_j_is,ioff_i_js,ioff_j_js;
  int ind_loc_i_is,ind_loc_i_js,ind_loc_j_is,ind_loc_j_js;
  int hess_ind;
  int ind_force_mat_1,ind_force_mat_2;
  int atomInd,countAtom;
  double p11,p22,p33,p12,p13,p23,cp_enl_now;
  double cp_enl = *cp_enl_ret;
/*==========================================================================*/
/* I) Loop over the 2*l+1 directions and sum the nl contributions           */
    for(m = 1;m<=(2*l+1);m++){
      ind_lm = m + l*l;
      for(is=1;is<=nstate;is++){
	for(js=is;js<=nstate;js++){
/*----------------------------------------------------------------------*/
/*  i) Get the contrib to the non-local energy                          */
	  cp_enl_now = 0.0;
	  ioff_i_is = (is-1)*np_nlmax+(ind_lm-1)*nstate*np_nlmax
                   +(irad-1)*nl_chan_max*nstate*np_nlmax;
          ioff_j_is = (is-1)*np_nlmax+(ind_lm-1)*nstate*np_nlmax
                   +(jrad-1)*nl_chan_max*nstate*np_nlmax;
          ioff_i_js = (js-1)*np_nlmax+(ind_lm-1)*nstate*np_nlmax
                   +(irad-1)*nl_chan_max*nstate*np_nlmax;
          ioff_j_js = (js-1)*np_nlmax+(ind_lm-1)*nstate*np_nlmax
                   +(jrad-1)*nl_chan_max*nstate*np_nlmax;
	  if(irad==jrad){
	    i_shift = l*npart;
            for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
	      atomInd = ip_nl[(ipart+i_shift)]-1;
	      //printf("atomInd %i\n",atomInd);
	      if(atomVnlCalcFlag[atomInd]==1){
		//printf("atomInd %i\n",atomInd);
		ind_loc_i_is = ipart + ioff_i_is;
		ind_loc_i_js = ipart + ioff_i_js;
		cp_enl_now = (vnlreal[ind_loc_i_is]*vnlreal[ind_loc_i_js]
			     +vnlimag[ind_loc_i_is]*vnlimag[ind_loc_i_js])*vnorm_now[ipart];
		vnlMatrix[(is-1)*nstate+js-1] += cp_enl_now;
                //printf("vnlreal %lg vnlimag %lg\n",vnlreal[ind_loc_i_is],vnlimag[ind_loc_i_is]);

		if(is==js){
		  cp_enl += cp_enl_now;
		  //if(ipart>=16&&atomInd<24)cp_enl += cp_enl_now;
		  //printf("cp_enl_now %lg\n",cp_enl_now);
		}
		else vnlMatrix[(js-1)*nstate+is-1] += cp_enl_now;
		//printf("atomInd %i cp_enl %lg is %i js %i vnlMatrix %lg %lg\n",atomInd,cp_enl,is,js,vnlMatrix[(is-1)*nstate+js-1],vnlMatrix[(js-1)*nstate+is-1]);
	      }//endif
            }//endfor
	  }//endif
	  else{
	    i_shift = l*npart;
            for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
	      atomInd = ip_nl[(ipart+i_shift)]-1;
	      //printf("atomInd %i\n",atomInd);
	      if(atomVnlCalcFlag[atomInd]==1){
		ind_loc_i_is = ipart + ioff_i_is;
		ind_loc_i_js = ipart + ioff_i_js;
		ind_loc_j_is = ipart + ioff_j_is;
		ind_loc_j_js = ipart + ioff_j_js;
		cp_enl_now = 2.0*(vnlreal[ind_loc_i_is]*vnlreal[ind_loc_j_js]
			      +vnlimag[ind_loc_i_is]*vnlimag[ind_loc_j_js])
			      *vnorm_now[ipart];
		vnlMatrix[(is-1)*nstate+js-1] += cp_enl_now;
		if(is==js){
		  cp_enl += cp_enl_now;
		  //printf("vnlreal %lg vnlimag %lg\n",vnlreal[ind_loc_i_is],vnlimag[ind_loc_i_is]);
		}
		else vnlMatrix[(js-1)*nstate+is-1] += cp_enl_now;
		//printf("atomInd %i cp_enl %lg is %i js %i vnlMatrix %i\n",atomInd,cp_enl,is,js,vnlMatrix[(js-1)*nstate+is-1]);
	      }//endif
            }//endfor
	  }
	  /*
	  cp_enl_now = 0.0;
	  ioff_i = (is-1)*np_nlmax+(ind_lm-1)*nstate*np_nlmax
		 +(irad-1)*nl_chan_max*nstate*np_nlmax;
	  ioff_j = (is-1)*np_nlmax+(ind_lm-1)*nstate*np_nlmax
		 +(jrad-1)*nl_chan_max*nstate*np_nlmax;
	  if(irad==jrad){
	    for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
	      ind_loc = ipart + ioff_i;
	      cp_enl_now = (vnlreal[ind_loc]*vnlreal[ind_loc]
			   +vnlimag[ind_loc]*vnlimag[ind_loc])*vnorm_now[ipart];
	      cp_enl += cp_enl_now;
	    }//endfor
	  }
	  else{
	    for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
	      ind_loc_i = ipart + ioff_i;
	      ind_loc_j = ipart + ioff_j;
	      cp_enl += 2.0*(vnlreal[ind_loc_i]*vnlreal[ind_loc_j]
			   +vnlimag[ind_loc_i]*vnlimag[ind_loc_j])
			   *vnorm_now[ipart];

	    }//endfor
	  }//endif
	  */

/*----------------------------------------------------------------------*/
/* ii) Get the contrib to non-local piece of the pressure tensor        */
	/* //Turn on ptens in the future
        if(cp_ptens==1){
         if(irad==jrad){
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc = ipart + ioff_i;
            p11 = 2.0*(dvnlreal_gxgx[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gxgx[ind_loc]*vnlimag[ind_loc])
                     *vnorm_now[ipart];
            p22 = 2.0*(dvnlreal_gygy[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gygy[ind_loc]*vnlimag[ind_loc])
                     *vnorm_now[ipart];
            p33 = 2.0*(dvnlreal_gzgz[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gzgz[ind_loc]*vnlimag[ind_loc])
                      *vnorm_now[ipart];
            p12 = 2.0*(dvnlreal_gxgy[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gxgy[ind_loc]*vnlimag[ind_loc])
                      *vnorm_now[ipart];
            p13 = 2.0*(dvnlreal_gxgz[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gxgz[ind_loc]*vnlimag[ind_loc])
                      *vnorm_now[ipart];
            p23 = 2.0*(dvnlreal_gygz[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gygz[ind_loc]*vnlimag[ind_loc])
                     *vnorm_now[ipart];
            pvten[1] += p11;
            pvten[5] += p22;
            pvten[9] += p33;
            pvten[4] += p12;
            pvten[2] += p12;
            pvten[7] += p13;
            pvten[3] += p13;
            pvten[8] += p23;
            pvten[6] += p23;
          }//endfor:ipart
         }else{
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc_i = ipart + ioff_i;
            ind_loc_j = ipart + ioff_j;
            p11 = 2.0*(dvnlreal_gxgx[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gxgx[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gxgx[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gxgx[ind_loc_j]*vnlimag[ind_loc_i])
                     *vnorm_now[ipart];
            p22 = 2.0*(dvnlreal_gygy[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gygy[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gygy[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gygy[ind_loc_j]*vnlimag[ind_loc_i])
                     *vnorm_now[ipart];
            p33 = 2.0*(dvnlreal_gzgz[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gzgz[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gzgz[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gzgz[ind_loc_j]*vnlimag[ind_loc_i])
                      *vnorm_now[ipart];
            p12 = 2.0*(dvnlreal_gxgy[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gxgy[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gxgy[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gxgy[ind_loc_j]*vnlimag[ind_loc_i])
                      *vnorm_now[ipart];
            p13 = 2.0*(dvnlreal_gxgz[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gxgz[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gxgz[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gxgz[ind_loc_j]*vnlimag[ind_loc_i])
                      *vnorm_now[ipart];
            p23 = 2.0*(dvnlreal_gygz[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gygz[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gygz[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gygz[ind_loc_j]*vnlimag[ind_loc_i])
                     *vnorm_now[ipart];
            pvten[1] += p11;
            pvten[5] += p22;
            pvten[9] += p33;
            pvten[4] += p12;
            pvten[2] += p12;
            pvten[7] += p13;
            pvten[3] += p13;
            pvten[8] += p23;
            pvten[6] += p23;
          }//endfor:ipart
         }//endif
        }//endif:cp_ptens on
	*/
/*----------------------------------------------------------------------*/
/* iii) Sum the non-local particle force                                  */

// All imaginary part of val and dvnl are 0
          if(irad==jrad){
            for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
              ind_loc_i_is = ipart + ioff_i_is;
              ind_loc_i_js = ipart + ioff_i_js;
              fxtemp[ipart] = -(dvnlreal_x[ind_loc_i_is]*vnlreal[ind_loc_i_js]
                               +dvnlreal_x[ind_loc_i_js]*vnlreal[ind_loc_i_is])
                               *vnorm_now[ipart];
              fytemp[ipart] = -(dvnlreal_y[ind_loc_i_is]*vnlreal[ind_loc_i_js]
                               +dvnlreal_y[ind_loc_i_js]*vnlreal[ind_loc_i_is])
                               *vnorm_now[ipart];
              fztemp[ipart] = -(dvnlreal_z[ind_loc_i_is]*vnlreal[ind_loc_i_js]
                               +dvnlreal_z[ind_loc_i_js]*vnlreal[ind_loc_i_is])
                               *vnorm_now[ipart];
	      //debug
	      /*
	      if(is==js){
		printf("dvnlreal_x %lg dvnlimag_x %lg dvnlreal_y %lg dvnlimag_y %lg dvnlreal_z %lg dvnlimag_z %lg\n",
			dvnlreal_x[ind_loc_i_is],dvnlimag_x[ind_loc_i_is],
			dvnlreal_y[ind_loc_i_is],dvnlimag_y[ind_loc_i_is],dvnlreal_z[ind_loc_i_is],
			dvnlimag_z[ind_loc_i_is]);
	      }
	      */
            }//endfor
          }
          else{
            for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
              ind_loc_i_is = ipart + ioff_i_is;
              ind_loc_i_js = ipart + ioff_i_js;
              ind_loc_j_is = ipart + ioff_j_is;
              ind_loc_j_js = ipart + ioff_j_js;
              fxtemp[ipart] = -(dvnlreal_x[ind_loc_i_is]*vnlreal[ind_loc_j_js]
                                +dvnlreal_x[ind_loc_j_is]*vnlreal[ind_loc_i_js]
                                +dvnlreal_x[ind_loc_i_js]*vnlreal[ind_loc_j_is]
                                +dvnlreal_x[ind_loc_j_js]*vnlreal[ind_loc_i_is])
                                *vnorm_now[ipart];
              fytemp[ipart] = -(dvnlreal_y[ind_loc_i_is]*vnlreal[ind_loc_j_js]
                                +dvnlreal_y[ind_loc_j_is]*vnlreal[ind_loc_i_js]
                                +dvnlreal_y[ind_loc_i_js]*vnlreal[ind_loc_j_is]
                                +dvnlreal_y[ind_loc_j_js]*vnlreal[ind_loc_i_is])
                                *vnorm_now[ipart];
              fztemp[ipart] = -(dvnlreal_z[ind_loc_i_is]*vnlreal[ind_loc_j_js]
                                +dvnlreal_z[ind_loc_j_is]*vnlreal[ind_loc_i_js]
                                +dvnlreal_z[ind_loc_i_js]*vnlreal[ind_loc_j_is]
                                +dvnlreal_z[ind_loc_j_js]*vnlreal[ind_loc_i_is])
                               *vnorm_now[ipart];
            }//endfor
          }//endif
          i_shift = l*npart;
	  countAtom = 0;
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            atomInd = ip_nl[(ipart+i_shift)]-1;
	    if(atomVnlCalcFlag[atomInd]==1){
	      //printf("countAtom %i\n",countAtom);
	      ind_force_mat_1 = countAtom*nstate*nstate+(is-1)*nstate+js-1;
	      ind_force_mat_2 = countAtom*nstate*nstate+(js-1)*nstate+is-1;
	      vnlFxMatrix[ind_force_mat_1] += fxtemp[ipart];
	      vnlFyMatrix[ind_force_mat_1] += fytemp[ipart];
	      vnlFzMatrix[ind_force_mat_1] += fztemp[ipart];
	      if(js==is){
		fx[countAtom+1] += fxtemp[ipart];
		fy[countAtom+1] += fytemp[ipart];
		fz[countAtom+1] += fztemp[ipart];
	      }
	      else{
		vnlFxMatrix[ind_force_mat_2] += fxtemp[ipart];
		vnlFyMatrix[ind_force_mat_2] += fytemp[ipart];
		vnlFzMatrix[ind_force_mat_2] += fztemp[ipart];
	      }
	      countAtom += 1;
	      //printf("atomInd %i vnlFxMatrix %lg vnlFyMatrix %lg vnlFzMatrix %lg\n",atomInd,vnlFxMatrix[ind_force_mat_1],vnlFyMatrix[ind_force_mat_1],vnlFzMatrix[ind_force_mat_1]);
	    }//endif atomVnlCalcFlag
          }/*endfor:atomic forces*/

	  /*     
	  if(irad==jrad){
            for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
	      ind_loc = ipart + ioff_i;
	      fxtemp[ipart] = -2.0*(dvnlreal_x[ind_loc]*vnlreal[ind_loc]
                               +dvnlimag_x[ind_loc]*vnlimag[ind_loc])
                               *vnorm_now[ipart];
	      fytemp[ipart] = -2.0*(dvnlreal_y[ind_loc]*vnlreal[ind_loc]
                               +dvnlimag_y[ind_loc]*vnlimag[ind_loc])
                               *vnorm_now[ipart];
	      fztemp[ipart] = -2.0*(dvnlreal_z[ind_loc]*vnlreal[ind_loc]
                               +dvnlimag_z[ind_loc]*vnlimag[ind_loc])
                               *vnorm_now[ipart];
            }//endfor
	  }
	  else{
            for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
              ind_loc_i = ipart + ioff_i;
              ind_loc_j = ipart + ioff_j;
	      fxtemp[ipart] = -2.0*(dvnlreal_x[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_x[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_x[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_x[ind_loc_j]*vnlimag[ind_loc_i])
                                *vnorm_now[ipart];
	      fytemp[ipart] = -2.0*(dvnlreal_y[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_y[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_y[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_y[ind_loc_j]*vnlimag[ind_loc_i])
                                *vnorm_now[ipart];
	      fztemp[ipart] = -2.0*(dvnlreal_z[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_z[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_z[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_z[ind_loc_j]*vnlimag[ind_loc_i])
                               *vnorm_now[ipart];
            }//endfor
          }//endif
          i_shift = l*npart;
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ltemp = ip_nl[(ipart+i_shift)];
            fx[ltemp] += fxtemp[ipart];
            fy[ltemp] += fytemp[ipart];
            fz[ltemp] += fztemp[ipart];
          }//endfor:atomic forces
	  */

/*----------------------------------------------------------------------*/
/* iv) Sum the non-local atomic hessian                                  */

        /* // Hessian comes later
        if(atm_hess_calc == 3){
         if(irad==jrad){
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc = ipart + ioff_i;
            fxtemp[ipart] = 2.0*(dvnlreal_gxgx[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gxgx[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_x[ind_loc]*dvnlreal_x[ind_loc]
                                +dvnlimag_x[ind_loc]*dvnlimag_x[ind_loc])
                                *vnorm_now[ipart];
            fytemp[ipart] = 2.0*(dvnlreal_gxgy[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gxgy[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_x[ind_loc]*dvnlreal_y[ind_loc]
                                +dvnlimag_x[ind_loc]*dvnlimag_y[ind_loc])
                                *vnorm_now[ipart];
            fztemp[ipart] = 2.0*(dvnlreal_gxgz[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gxgz[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_x[ind_loc]*dvnlreal_z[ind_loc]
                                +dvnlimag_x[ind_loc]*dvnlimag_z[ind_loc])
                                *vnorm_now[ipart];
          }// endfor ipart
         } else {
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc_i = ipart + ioff_i;
            ind_loc_j = ipart + ioff_j;
            fxtemp[ipart] = 2.0*(dvnlreal_gxgx[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gxgx[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gxgx[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gxgx[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_x[ind_loc_i]*dvnlreal_x[ind_loc_j]
                                     +dvnlimag_x[ind_loc_i]*dvnlimag_x[ind_loc_j]))
                                *vnorm_now[ipart];
            fytemp[ipart] = 2.0*(dvnlreal_gxgy[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gxgy[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gxgy[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gxgy[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_x[ind_loc_i]*dvnlreal_y[ind_loc_j]
                                     +dvnlimag_x[ind_loc_i]*dvnlimag_y[ind_loc_j]))
                                *vnorm_now[ipart];
            fztemp[ipart] = 2.0*(dvnlreal_gxgz[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gxgz[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gxgz[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gxgz[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_x[ind_loc_i]*dvnlreal_z[ind_loc_j]
                                     +dvnlimag_x[ind_loc_i]*dvnlimag_z[ind_loc_j]))
                                *vnorm_now[ipart];
          }// endfor
         }// endif
         i_shift = l*npart;
         for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
           ltemp = ip_nl[(ipart+i_shift)];
           hess_ind = (ltemp-1)*npart + ltemp;
           hess_xx[hess_ind] += fxtemp[ipart];
           hess_xy[hess_ind] += fytemp[ipart];
           hess_xz[hess_ind] += fztemp[ipart];
         }//endfor:atomic forces

         if(irad==jrad){
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc = ipart + ioff_i;
            fxtemp[ipart] = 2.0*(dvnlreal_gygy[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gygy[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_y[ind_loc]*dvnlreal_y[ind_loc]
                                +dvnlimag_y[ind_loc]*dvnlimag_y[ind_loc])
                                *vnorm_now[ipart];
            fytemp[ipart] = 2.0*(dvnlreal_gygz[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gygz[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_y[ind_loc]*dvnlreal_z[ind_loc]
                                +dvnlimag_y[ind_loc]*dvnlimag_z[ind_loc])
                                *vnorm_now[ipart];
            fztemp[ipart] = 2.0*(dvnlreal_gzgz[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gzgz[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_z[ind_loc]*dvnlreal_z[ind_loc]
                                +dvnlimag_z[ind_loc]*dvnlimag_z[ind_loc])
                                *vnorm_now[ipart];
          }// endfor ipart
         } else {
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc_i = ipart + ioff_i;
            ind_loc_j = ipart + ioff_j;
            fxtemp[ipart] = 2.0*(dvnlreal_gygy[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gygy[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gygy[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gygy[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_y[ind_loc_i]*dvnlreal_y[ind_loc_j]
                                     +dvnlimag_y[ind_loc_i]*dvnlimag_y[ind_loc_j]))
                                *vnorm_now[ipart];
            fytemp[ipart] = 2.0*(dvnlreal_gygz[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gygz[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gygz[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gygz[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_y[ind_loc_i]*dvnlreal_z[ind_loc_j]
                                     +dvnlimag_y[ind_loc_i]*dvnlimag_z[ind_loc_j]))
                                *vnorm_now[ipart];
            fztemp[ipart] = 2.0*(dvnlreal_gzgz[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gzgz[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gzgz[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gzgz[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_z[ind_loc_i]*dvnlreal_z[ind_loc_j]
                                     +dvnlimag_z[ind_loc_i]*dvnlimag_z[ind_loc_j]))
                                *vnorm_now[ipart];
          }// endfor ipart
         }// endif
         i_shift = l*npart;
         for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
           ltemp = ip_nl[(ipart+i_shift)];
           hess_ind = (ltemp-1)*npart + ltemp;
           hess_yy[hess_ind] += fxtemp[ipart];
           hess_yz[hess_ind] += fytemp[ipart];
           hess_zz[hess_ind] += fztemp[ipart];
         }//endfor:atomic forces 

        }// endif: atm hess calc
	*/
      }//endfor loop over j states
    }//endfor:loop over i states
  }//endfor: loop over the m channels

/*==========================================================================*/
/* II) Set the return values                                               */

 printf("cp_enl %lg\n",cp_enl);

 *cp_enl_ret = cp_enl;

/*======================================================================*/
  }/*end routine*/
/*======================================================================*/

