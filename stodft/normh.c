/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: normh.c                                        */
/*                                                                          */
/* This routine costruct P_N(H)|phi> where P_N is some                      */
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
void normHNewtonHerm(CP *cp,CPCOEFFS_POS *cpcoeffs_pos,CPCOEFFS_INFO *cpcoeffs_info,
		 CELL *cell,CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
		 EWALD *ewald,EWD_SCR *ewd_scr,ATOMMAPS *atommaps,FOR_SCR *for_scr,
		 STAT_AVG *stat_avg,PTENS *ptens,double zn)


/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Let's first build a simple version. We will first apply this on       */
/* Hermitian operator. So all the sample points will locate on the real  */
/* axis. zn is the interpolation point. For details please read		 */
/* Ashkenazi et.al. JChemPhys 103, 10005(1995).				 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPEWALD *cpewald              = &(cp->cpewald);
  CPSCR *cpscr                  = &(cp->cpscr);
  CPOPTS *cpopts                = &(cp->cpopts);
  PSEUDO *pseudo                = &(cp->pseudo);
  COMMUNICATE *communicate      = &(cp->communicate);

  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm;            = &(cp->cp_sclr_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_sm;		   = &(cp->cp_para_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_dens_cp_box;   = &(cp->cp_sclr_fft_pkg3d_dens_cp_box);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box;   = &(cp->cp_para_fft_pkg3d_dens_cp_box);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_lg;		   = &(cp->cp_sclr_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg;		   = &(cp->cp_para_fft_pkg3d_lg);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_up          = &(cp_comm_state_pkg_up);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_dn          = &(cp_comm_state_pkg_dn);

  NWETONINFO *newtonInfo = stodftInfo->newtonInfo;
  double Smin		= newtonInfo->Smin;
  double Smax		= newtonInfo->Smax;
  double scale		= newtonInfo->scale;
  double energyMean	= stodftInfo->energyMean;
  double energyDiff	= stodftInfo->energyDiff;
  double prefact	= -scale*energyMean-zn;
  double scale1		= -scale*0.5;
  double scale2		= -scale;

  int cpLsda         = cpopts->cp_lsda;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int cp_dual_grid_opt_on  = cpopts->cp_dual_grid_opt;
  int numCoeffM1     = numCoeff-1;
  int incx = 1;
  int incy = 1;
  int iState,iCoeff,iCoeffStart,index1,index2;
  
  double *expanCoeff = (double*)stodftCoefPos->expanCoeff;
  double complex *wfInUp   = stodftCoefPos->wfInUp;
  double complex *wfInDn   = stodftCoefPos->wfInDn;
  double complex *wfOutUp  = stodftCoefPos->wfOutUp;
  double complex *wfOutDn  = stodftCoefPos->wfOutDn;
  
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;

/*==========================================================================*/
/* 0) Copy the input wave function to CP coeff and zero the force */
 
  memcpy(&(cre_up[1]),wfInReUp,numCoeffUpTotal);
  memcpy(&(cim_up[1]),wfInImUp,numCoeffUpTotal);
  for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
    fcre_up[iCoeff] = 0.0;
    fcim_up[iCoeff] = 0.0;
  }
  if(cpLsda==1&&numStateDnProc!=0){
    memcpy(&(cre_dn[1]),wfInReDn,numCoeffDnTotal);
    memcpy(&(cim_dn[1]),wfInImDn,numCoeffDnTotal);
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      fcre_dn[iCoeff] = 0.0;
      fcim_dn[iCoeff] = 0.0;
    }
  }
  
/*==========================================================================*/
/* 1) Calculate the H/sigma|phi> */
  //control_vps_atm_list will be done somewhere else (perhaps in density calculation?)

  control_cp_eext_recip(clatoms_info,clatoms_pos,cpcoeffs_info,
                       cpcoeffs_pos,cpewald,cpscr,cpopts,pseudo,
		       ewd_scr,atommaps,cell,ewald,ptens,&(stat_avg->vrecip),
		       &(stat_avg->cp_enl),communicate,for_scr,cp_dual_grid_opt_on,
                       cp_para_fft_pkg3d_lg);

  coef_force_control(cpopts,cpcoeffs_info,cpcoeffs_pos,cpscr,ewald,cpewald,
		    cell,stat_avg,pseudo->vxc_typ,ptens->pvten_tmp,pseudo->gga_cut,
		    pseudo->alpha_conv_dual,pseudo->n_interp_pme_dual,cp_min_on,
		    communicate,cp_comm_state_pkg_up,
                     cp_comm_state_pkg_dn,cp_para_fft_pkg3d_lg,cp_sclr_fft_pkg3d_lg,
                     cp_para_fft_pkg3d_dens_cp_box,cp_sclr_fft_pkg3d_dens_cp_box,
		     cp_para_fft_pkg3d_sm,cp_sclr_fft_pkg3d_sm,cp_dual_grid_opt_on);

  for(iState=0;iState<numStateUpProc;iState++){
    iCoeffStart = iState*numCoeff;
    for(iCoeff=0;iCoeff<numCoeffM1;iCoeff++){
      index1 = iCoeffStart+iCoeff;
      wfOutUp[index1] = fcre_up[index1+1]*scale1+fcim_up[index1+1]*scale1*I;
    }//endfor iCoeff
    index1 = iCoeffStart+numCoeffM1;
    wfOutUp[index1] = fcre_up[index1+1]*scale2;
  }//endfor iState
  if(cpLsda==1&&numStateDnProc!=0){
    for(iState=0;iState<numStateDnProc;iState++){
      iCoeffStart = iState*numCoeff;
      for(iCoeff=0;iCoeff<numCoeffM1;iCoeff++){
	index1 = iCoeffStart+iCoeff;
	wfOutDn[index1] = fcre_dn[index1+1]*scale1+fcim_dn[index1+1]*scale1*I;
      }//endfor iCoeff
      index1 = iCoeffStart+numCoeffM1;
      index2 = index1+1;
      wfOutDn[index1] = fcre_dn[index1+1]*scale2+fcim_dn[index1+1]*scale2*I;
    }//endfor iState
  }

/*==========================================================================*/
/* 2) Calculate P_(n+1)(H)|phi> */

  ZAXPY(&numCoeffUpTotal,&prefact,wfInReUp,&incx,wfOutReUp,&incy);
  ZAXPY(&numCoeffUpTotal,&prefact,wfInImUp,&incx,wfOutImUp,&incy);
  if(cpLsda==1&&numStateDnProc!=0){
    ZAXPY(&numCoeffDnTotal,&prefact,wfInReDn,&incx,wfOutReDn,&incy);
    ZAXPY(&numCoeffDnTotal,&prefact,wfInImDn,&incx,wfOutImDn,&incy);
  }
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void normHNewtonNoHerm(CP *cp,CPCOEFFS_POS *cpcoeffs_pos,CPCOEFFS_INFO *cpcoeffs_info,
	 CELL *cell,CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
	 EWALD *ewald,EWD_SCR *ewd_scr,ATOMMAPS *atommaps,FOR_SCR *for_scr,
	 STAT_AVG *stat_avg,PTENS *ptens,double complex zn)


/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Let's first build a simple version. We will first apply this on       */
/* Hermitian operator. So all the sample points will locate on the real  */
/* axis. zn is the interpolation point. For details please read	     */
/* Ashkenazi et.al. JChemPhys 103, 10005(1995).		     */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPEWALD *cpewald              = &(cp->cpewald);
  CPSCR *cpscr                  = &(cp->cpscr);
  CPOPTS *cpopts                = &(cp->cpopts);
  PSEUDO *pseudo                = &(cp->pseudo);
  COMMUNICATE *communicate      = &(cp->communicate);

  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm;            = &(cp->cp_sclr_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_sm;	       = &(cp->cp_para_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_dens_cp_box;   = &(cp->cp_sclr_fft_pkg3d_dens_cp_box);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box;   = &(cp->cp_para_fft_pkg3d_dens_cp_box);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_lg;	       = &(cp->cp_sclr_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg;	       = &(cp->cp_para_fft_pkg3d_lg);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_up          = &(cp_comm_state_pkg_up);
  CP_COMM_STATE_PKG *cp_comm_state_pkg_dn          = &(cp_comm_state_pkg_dn);

  NWETONINFO *newtonInfo = stodftInfo->newtonInfo;
  double Smin	    = newtonInfo->Smin;
  double Smax	    = newtonInfo->Smax;
  double scale	    = newtonInfo->scale;
  double energyMean = stodftInfo->energyMean;
  double energyDiff = stodftInfo->energyDiff;
  double prefact    = -scale*energyMean-zn;
  double scale1	    = -scale*0.5;
  double scale2	    = -scale;

  int cpLsda         = cpopts->cp_lsda;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int cp_dual_grid_opt_on  = cpopts->cp_dual_grid_opt;
  int numCoeffM1     = numCoeff-1;
  int incx = 1;
  int incy = 1;
  int iState,iCoeff,iCoeffStart,index1,index2;
  
  double complex *expanCoeff = (double complex*)stodftCoefPos->expanCoeff;
  double complex *wfInUp   = stodftCoefPos->wfInUp;
  double complex *wfInDn   = stodftCoefPos->wfInDn;
  double complex *wfOutUp  = stodftCoefPos->wfOutUp;
  double complex *wfOutDn  = stodftCoefPos->wfOutDn;
  
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;

/*==========================================================================*/
/* 0) Copy the input wave function to CP coeff and zero the force */
 
  memcpy(&(cre_up[1]),wfInReUp,numCoeffUpTotal);
  memcpy(&(cim_up[1]),wfInImUp,numCoeffUpTotal);
  for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
    fcre_up[iCoeff] = 0.0;
    fcim_up[iCoeff] = 0.0;
  }
  if(cpLsda==1&&numStateDnProc!=0){
    memcpy(&(cre_dn[1]),wfInReDn,numCoeffDnTotal);
    memcpy(&(cim_dn[1]),wfInImDn,numCoeffDnTotal);
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      fcre_dn[iCoeff] = 0.0;
      fcim_dn[iCoeff] = 0.0;
    }
  }
  
/*==========================================================================*/
/* 1) Calculate the H/sigma|phi> */
  //control_vps_atm_list will be done somewhere else (perhaps in density calculation?)

  control_cp_eext_recip(clatoms_info,clatoms_pos,cpcoeffs_info,
                       cpcoeffs_pos,cpewald,cpscr,cpopts,pseudo,
	       ewd_scr,atommaps,cell,ewald,ptens,&(stat_avg->vrecip),
	       &(stat_avg->cp_enl),communicate,for_scr,cp_dual_grid_opt_on,
                       cp_para_fft_pkg3d_lg);

  coef_force_control(cpopts,cpcoeffs_info,cpcoeffs_pos,cpscr,ewald,cpewald,
	    cell,stat_avg,pseudo->vxc_typ,ptens->pvten_tmp,pseudo->gga_cut,
	    pseudo->alpha_conv_dual,pseudo->n_interp_pme_dual,cp_min_on,
	    communicate,cp_comm_state_pkg_up,
                     cp_comm_state_pkg_dn,cp_para_fft_pkg3d_lg,cp_sclr_fft_pkg3d_lg,
                     cp_para_fft_pkg3d_dens_cp_box,cp_sclr_fft_pkg3d_dens_cp_box,
	     cp_para_fft_pkg3d_sm,cp_sclr_fft_pkg3d_sm,cp_dual_grid_opt_on);

  for(iState=0;iState<numStateUpProc;iState++){
    iCoeffStart = iState*numCoeff;
    for(iCoeff=0;iCoeff<numCoeffM1;iCoeff++){
      index1 = iCoeffStart+iCoeff;
      wfOutUp[index1] = fcre_up[index1+1]*scale1+fcim_up[index1+1]*scale1*I;
    }//endfor iCoeff
    index1 = iCoeffStart+numCoeffM1;
    wfOutUp[index1] = fcre_up[index1+1]*scale2;
  }//endfor iState
  if(cpLsda==1&&numStateDnProc!=0){
    for(iState=0;iState<numStateDnProc;iState++){
      iCoeffStart = iState*numCoeff;
      for(iCoeff=0;iCoeff<numCoeffM1;iCoeff++){
    index1 = iCoeffStart+iCoeff;
    wfOutDn[index1] = fcre_dn[index1+1]*scale1+fcim_dn[index1+1]*scale1*I;
      }//endfor iCoeff
      index1 = iCoeffStart+numCoeffM1;
      index2 = index1+1;
      wfOutDn[index1] = fcre_dn[index1+1]*scale2+fcim_dn[index1+1]*scale2*I;
    }//endfor iState
  }

/*==========================================================================*/
/* 2) Calculate P_(n+1)(H)|phi> */

  ZAXPY(&numCoeffUpTotal,&prefact,wfInReUp,&incx,wfOutReUp,&incy);
  ZAXPY(&numCoeffUpTotal,&prefact,wfInImUp,&incx,wfOutImUp,&incy);
  if(cpLsda==1&&numStateDnProc!=0){
    ZAXPY(&numCoeffDnTotal,&prefact,wfInReDn,&incx,wfOutReDn,&incy);
    ZAXPY(&numCoeffDnTotal,&prefact,wfInImDn,&incx,wfOutImDn,&incy);
  }
  
/*==========================================================================*/
}/*end Routine*/
/*=======================================================================*/

