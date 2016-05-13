/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: init.c                                         */
/*                                                                          */
/* This routine initialize the stochastic dft calculation.                  */
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
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_stodft_local.h"


#define TIME_CP_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initStodft(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp,
		int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate expanCoeff for different chemical potentials                */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS  *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  ATOMMAPS     *atommaps     = &(class->atommaps);
  CELL         *cell         = &(general_data->cell);
  FOR_SCR      *for_scr      = &(class->for_scr);
  EWD_SCR      *ewd_scr      = &(class->ewd_scr);



  CPOPTS *cpopts = &(cp->cpopts);
  PSEUDO *pseudo = &(cp->pseudo);
  CPCOEFFS_INFO *cpcoeffs_info = &(cp->cpcoeffs_info);
  STODFTINFO *stodftInfo       = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos = cp->stodftCoefPos;
  NEWTONINFO *newtonInfo;

  int cpLsda         = cpopts->cp_lsda;
  int cpGga         = cpopts->cp_gga;
  int expanType      = stodftInfo->expanType;
  int numOrbital     = stodftInfo->numOrbital;
  int polynormLength = stodftInfo->polynormLength;
  int numChemPot     = stodftInfo->numChemPot;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numStateUpTot  = numStateUpProc*numCoeff;
  int numStateDnTot  = numStateDnProc*numCoeff;
  int totalPoly	     = polynormLength*numChemPot;
  int filterFunType   = stodftInfo->filterFunType;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int iChem,iSamp;

  char *ggaxTyp     = pseudo->ggax_typ;
  char *ggacTyp     = pseudo->ggac_typ;
  
  double Smin = -2.0;
  double Smax = 2.0;
  double energyDiff;

/*==========================================================================*/
/* I) General parameters and malloc					    */
  
  stodftInfo->vpsAtomListFlag = 0;

  stodftCoefPos->chemPot = (double*)cmalloc(numChemPot*sizeof(double));

  stodftCoefPos->stoWfUpRe = (double**)cmalloc(numChemPot*sizeof(double*));
  stodftCoefPos->stoWfUpIm = (double**)cmalloc(numChemPot*sizeof(double*));

  for(iChem=0;iChem<numChemPot;iChem++){
    stodftCoefPos->stoWfUpRe[iChem] = (double*)cmalloc(numStateUpTot*sizeof(double))-1;
    stodftCoefPos->stoWfUpIm[iChem] = (double*)cmalloc(numStateUpTot*sizeof(double))-1;
  }
  if(cpLsda==1&&numStateDnProc!=0){
    stodftCoefPos->stoWfUpRe = (double**)cmalloc(numChemPot*sizeof(double*));
    stodftCoefPos->stoWfUpIm = (double**)cmalloc(numChemPot*sizeof(double*));
    for(iChem=0;iChem<numChemPot;iChem++){
      stodftCoefPos->stoWfDnRe[iChem] = (double*)cmalloc(numStateUpTot*sizeof(double))-1;
      stodftCoefPos->stoWfDnIm[iChem] = (double*)cmalloc(numStateUpTot*sizeof(double))-1;
    }//endfor iChem
  }//endif

  if(expanType==2&&filterFunType==1)stodftInfo->fermiFunctionReal = &fermiExpReal;
  if(expanType==2&&filterFunType==2)stodftInfo->fermiFunctionReal = &fermiErfcReal;
  if(expanType==2&&filterFunType==3)stodftInfo->fermiFunctionReal = &gaussianReal;
  if(expanType==3&&filterFunType==1)stodftInfo->fermiFunctionComplex = &fermiExpComplex;
  if(expanType==3&&filterFunType==2){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("We haven't implement erfc type of Fermi function \n");
    printf("for non-Hermitian KS Hamiltonian. Please use the \n");
    printf("exponential type in this case!\n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(0);
  }
  if(expanType==3&&filterFunType==3){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("We haven't implement gaussian type of Fermi function \n");
    printf("for non-Hermitian KS Hamiltonian. Please use the \n");
    printf("exponential type in this case!\n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(0);
  }



/*==========================================================================*/
/* II) Malloc by expension type						    */

  switch(expanType){
    case 2:
      stodftCoefPos->expanCoeff = (double *)cmalloc(totalPoly*sizeof(double));
      stodftInfo->newtonInfo = (NEWTONINFO *)cmalloc(sizeof(NEWTONINFO));
      newtonInfo = stodftInfo->newtonInfo;
      newtonInfo->sampPoint = (double *)cmalloc(polynormLength*sizeof(double));
      newtonInfo->sampPointUnscale = (double *)cmalloc(polynormLength*sizeof(double));
      newtonInfo->Smin = Smin;
      newtonInfo->Smax = Smax;
      //newtonInfo->scale = (Smax-Smin)/energyDiff;      
      
      break;
  }

  //For debug only
  //stodftCoefPos->chemPot[0] = -0.17435045;
  stodftCoefPos->chemPot[0] = 0.075726635;
  

/*==========================================================================*/
/* III) Initialize utility data						    */
  FILE *fileSampPoint = fopen("samp-point","r");
  double *sampLocal = (double*)(newtonInfo->sampPoint);

  switch(expanType){
    case 2:
      
      for(iSamp=0;iSamp<polynormLength;iSamp++){
	fscanf(fileSampPoint,"%lg",&(sampLocal[iSamp]));
	printf("samp %lg\n",sampLocal[iSamp]);
      }
      
      //genSampNewtonHermit(stodftInfo,stodftCoefPos);
      break;
    // I'll do chebyshev and non-Hermitain Newtonian later    
  }
  fclose(fileSampPoint);
/*==========================================================================*/
/* IV) Initialize Flags							    */

  general_data->stat_avg.count_diag_srot      = 0.0;
  general_data->stat_avg.fatm_mag = 10000.0;
  general_data->stat_avg.fatm_max = 10000.0;

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
  if(cpGga==1){
    if(strcasecmp(ggaxTyp,"becke"   )==0){cpopts->cp_becke=1;}
    if(strcasecmp(ggaxTyp,"pw91x"   )==0){cpopts->cp_pw91x=1;}
    if(strcasecmp(ggaxTyp,"fila_1x" )==0){cpopts->cp_fila_1x=1;}
    if(strcasecmp(ggaxTyp,"fila_2x" )==0){cpopts->cp_fila_2x=1;}
    if(strcasecmp(ggaxTyp,"pbe_x"   )==0){cpopts->cp_pbe_x=1;}
    if(strcasecmp(ggaxTyp,"revpbe_x")==0){cpopts->cp_revpbe_x=1;}
    if(strcasecmp(ggaxTyp,"rpbe_x"  )==0){cpopts->cp_rpbe_x=1;}
    if(strcasecmp(ggaxTyp,"xpbe_x"  )==0){cpopts->cp_xpbe_x=1;}
    if(strcasecmp(ggaxTyp,"brx89"   )==0){cpopts->cp_brx89=1;}
    if(strcasecmp(ggaxTyp,"brx2k"   )==0){cpopts->cp_brx2k=1;}
    if(strcasecmp(ggacTyp,"lyp"     )==0){cpopts->cp_lyp=1;  }
    if(strcasecmp(ggacTyp,"lypm1"   )==0){cpopts->cp_lypm1=1;  }
    if(strcasecmp(ggacTyp,"pw91c"   )==0){cpopts->cp_pw91c=1;}
    if(strcasecmp(ggacTyp,"pbe_c"   )==0){cpopts->cp_pbe_c=1;}
    if(strcasecmp(ggacTyp,"xpbe_c"  )==0){cpopts->cp_xpbe_c=1;}
    if(strcasecmp(ggacTyp,"tau1_c"  )==0){cpopts->cp_tau1_c=1;}
    if(strcasecmp(ggacTyp,"debug97x")==0){cpopts->cp_debug_xc=1;}
  }/*endif*/


/*==========================================================================*/
/* V) Calculate the non-local pseudopotential list                          */
  /*
  if(stodftInfo->vpsAtomListFlag==0||cpDualGridOptOn>= 1){
    control_vps_atm_list(pseudo,cell,clatoms_pos,clatoms_info,
                         atommaps,ewd_scr,for_scr,cpDualGridOptOn,
                         stodftInfo->vpsAtomListFlag);
    stodftInfo->vpsAtomListFlag = 1;
  }
  */
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcRhoInit(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                   CP *cp,int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the routine to calculate the initial density			 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  STODFTINFO *stodftInfo       = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos = cp->stodftCoefPos;
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]); 
 
  int reInitFlag = cp->reInitFlag;
  
  if(reInitFlag==0) calcRhoSto(class,bonded,general_data,cp,cpcoeffs_pos); 
  if(reInitFlag==1) calcRhoDet(class,bonded,general_data,cp,cpcoeffs_pos);
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void reInitWaveFunMin(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                   CP *cp,int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the routine to reset number of wave functions after           */
/* calculating the initial density. This is the sp/opt version.		 */
/* check files:#coords_cp/mall_properties.c				 */
/*	       #cp_ewald/control_set_cp_ewald.c				 */
/*	       #parse/parse.c						 */
/*	       #parse/zero_cp.c						 */
/*	       scratch/mall_scratch.c					 */
/*	       				 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos    = &(cp->cpcoeffs_pos[ip_now]);
  COMMUNICATE  *communicate     = &(cp->communicate);
  CPOPTS       *cpopts          = &(cp->cpopts);
  CPSCR        *cpscr           = &(cp->cpscr);

  int numStateStoUp = stodftInfo->numStateStoUp;
  int numStateStoDn = stodftInfo->numStateStoDn;
  int numProcstates = communicate->np_states;
  int cpLsda = cpopts->cp_lsda;
  int piBeadsProc  = cpcoeffs_info->pi_beads_proc;
  int iState;


/*==========================================================================*/
/* I) Initialize Check                                                      */
     
  cpcoeffs_info->nstate_up = numStateStoUp;
  cpcoeffs_info->nstate_dn = numStateStoDn;

  if(cpLsda==0){
    if(numStateStoUp<numProcstates){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("Number of states less than number of processors\n");
      printf("If possible, reduce number of processors to be\n");
      printf("less than the number of states or run a bigger system.\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }//endif
  }else{
    if(numStateStoUp+numStateStoDn<numProcstates){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("Number of states less than number of processors\n");
      printf("If possible, reduce number of processors to be\n");
      printf("less than the number of states or run a bigger system.\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }//endif
  }//endif

/*==========================================================================*/
/* II) Free old arrays, except scratch					    */
/*     (//N means don't remalloc these arrays)				    */

  free(&(cpcoeffs_pos->cre_up[1]));
  free(&(cpcoeffs_pos->cim_up[1]));
  free(&(cpcoeffs_pos->cre_dn[1]));
  free(&(cpcoeffs_pos->cim_dn[1]));
  free(&(cpcoeffs_pos->vcre_up[1]));//N
  free(&(cpcoeffs_pos->vcim_up[1]));//N
  free(&(cpcoeffs_pos->vcre_dn[1]));//N
  free(&(cpcoeffs_pos->vcim_dn[1]));//N
  free(&(cpcoeffs_pos->fcre_up[1]));
  free(&(cpcoeffs_pos->fcim_up[1]));
  free(&(cpcoeffs_pos->fcre_dn[1]));
  free(&(cpcoeffs_pos->fcim_dn[1]));
  free(&(cpcoeffs_pos->kfcre_up[1]));//N
  free(&(cpcoeffs_pos->kfcim_up[1]));//N
  free(&(cpcoeffs_pos->ksmat_up[1]));//N
  free(&(cpcoeffs_pos->ksmat_dn[1]));//N
  free(&(cpcoeffs_pos->ksmat_eig_up[1]));//N
  free(&(cpcoeffs_pos->ksmat_eig_dn[1]));//N
  free(&(cpcoeffs_pos->norbmat_up[1]));//N
  free(&(cpcoeffs_pos->norbmat_dn[1]));//N
  free(&(cpcoeffs_pos->norbmati_up[1]));//N
  free(&(cpcoeffs_pos->norbmati_dn[1]));//N
  free(&(cpcoeffs_pos->ovmat_eigv_up[1]));//N
  free(&(cpcoeffs_pos->ovmat_eigv_dn[1]));//N

  // cp_min_on>0 I need fake cp_min_on to cheat the code 
  // when I call coef_force_control
  free(&(cpcoeffs_pos->cp_hess_re_up[1]));//N
  free(&(cpcoeffs_pos->cp_hess_im_up[1]));//N
  free(&(cpcoeffs_pos->cp_hess_re_dn[1]));//N
  free(&(cpcoeffs_pos->cp_hess_im_dn[1]));//N


  free(&(cpcoeffs_info->ioff_up[1]));
  free(&(cpcoeffs_info->ioff_dn[1]));
  free(&(cpcoeffs_info->ioff_upt[1]));
  free(&(cpcoeffs_info->ioff_dnt[1]));

  free(&(cpopts->occ_up[1]));//N
  free(&(cpopts->occ_dn[1]));//N
  free(&(cpopts->rocc_sum_up[1]));//N
  free(&(cpopts->rocc_sum_dn[1]));//N

/*==========================================================================*/
/* III) Reinit arrays, except scratch					    */
  
  //need set up mpi first
  stoRealloc(cp,cpcoeffs_pos);

/*==========================================================================*/
/* IV) Free all scratch		                                            */


/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void stoRealloc(CP *cp,CPCOEFFS_POS *cpcoeffs_pos)
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*    Local Variables   */
#include "../typ_defs/typ_mask.h"
  int myid = cp->communicate.myid;
  int i,ip,nread,is;
  int par_size_up,par_size_dn,ncoef_up_tot;
  int nstate,nstate2,ncoef_dn_tot;
  double *cre,*cim,*vcre,*vcim,*fcre,*fcim;
  int *ioff_up,*ioff_upt,*ioff_dn,*ioff_dnt;
  double mem_test;

/*  Local Pointers */
  int pi_beads       = cp->cpcoeffs_info.pi_beads;
  int pi_beads_proc  = cp->cpcoeffs_info.pi_beads_proc;
  int nstate_up_proc = cp->cpcoeffs_info.nstate_up_proc;
  int nstate_dn_proc = cp->cpcoeffs_info.nstate_dn_proc;
  int nstate_up      = cp->cpcoeffs_info.nstate_up;
  int nstate_dn      = cp->cpcoeffs_info.nstate_dn;
  int cp_lda         = cp->cpopts.cp_lda;
  int cp_lsda        = cp->cpopts.cp_lsda;
  int ncoef_up       = cp->cpcoeffs_info.ncoef;
  int ncoef_dn       = cp->cpcoeffs_info.ncoef;
  int cp_norb        = cp->cpopts.cp_norb;

  int np_states               = cp->communicate.np_states;
  int ncoef_up_proc           = cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
  int ncoef_dn_proc           = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;
  int nstate_max_up           =cp->cp_comm_state_pkg_up.nstate_max;
  int nstate_ncoef_proc_max_up=cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
  int nstate_max_dn           =cp->cp_comm_state_pkg_dn.nstate_max;
  int nstate_ncoef_proc_max_dn=cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;

  if(cp_lda == 1){ncoef_dn = 0;}
/*==========================================================================*/
/* 0) Calculate the sizes */

  par_size_up = nstate_max_up*(nstate_ncoef_proc_max_up);
  ncoef_up_tot = nstate_up_proc*ncoef_up;
  ncoef_up_tot = MAX(ncoef_up_tot,par_size_up);

  par_size_dn = nstate_max_dn*(nstate_ncoef_proc_max_dn);
  if(cp_lda ==1){par_size_dn = 0;}
  ncoef_dn_tot = MAX(nstate_dn_proc*ncoef_dn,1);
  ncoef_dn_tot = MAX(ncoef_dn_tot,par_size_dn);

  nstate  = MAX(nstate_up,nstate_dn);
  nstate2 = nstate*nstate;

/*==========================================================================*/
/* I) Malloc the variables */

  for(i=1;i<=pi_beads_proc;i++){
    cpcoeffs_pos->cre_up =(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
    cpcoeffs_pos->cim_up =(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
    cpcoeffs_pos->cre_dn =(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
    cpcoeffs_pos->cim_dn =(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
    cpcoeffs_pos->fcre_up=(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
    cpcoeffs_pos->fcim_up=(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
    cpcoeffs_pos->fcre_dn=(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
    cpcoeffs_pos->fcim_dn=(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  }/*endfor*/

  cp->cpcoeffs_info.ioff_up  = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_dn  = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_upt = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_dnt = (int *)cmalloc(nstate*sizeof(int))-1;

/*==========================================================================*/
/* II) Assign the offsets */

  ioff_up  = cp->cpcoeffs_info.ioff_up;
  ioff_upt = cp->cpcoeffs_info.ioff_upt;
  ioff_dn  = cp->cpcoeffs_info.ioff_dn;
  ioff_dnt = cp->cpcoeffs_info.ioff_dnt;

  for(is=1;is<=nstate;is++){ioff_up[is]=(is-1)*ncoef_up;}
  for(is=1;is<=nstate;is++){ioff_dn[is]=(is-1)*ncoef_dn;}

  if(np_states==1){
    for(is=1;is<=nstate;is++){ioff_upt[is]=(is-1)*ncoef_up;}
    for(is=1;is<=nstate;is++){ioff_dnt[is]=(is-1)*ncoef_dn;}
  }else{
    for(is=1;is<=nstate;is++){ioff_upt[is]=(is-1)*ncoef_up_proc;}
    for(is=1;is<=nstate;is++){ioff_dnt[is]=(is-1)*ncoef_dn_proc;}
  }/*endif*/

/*========================================================================*/
/* III) Initialize coeficient arrays and form flags                       */


  cpcoeffs_pos->icoef_form_up  = 0;
  cpcoeffs_pos->ivcoef_form_up = 0;
  cpcoeffs_pos->ifcoef_form_up = 1;
  cpcoeffs_pos->icoef_orth_up  = 1;
  if(cp_norb>0){cpcoeffs_pos->icoef_orth_up = 0;}
  cpcoeffs_pos->ivcoef_orth_up  = 1;
  if(cp_norb>0){cpcoeffs_pos->ivcoef_orth_up = 0;}
  cre = cpcoeffs_pos->cre_up;
  cim = cpcoeffs_pos->cim_up;
  fcre = cpcoeffs_pos->fcre_up;
  fcim = cpcoeffs_pos->fcim_up;
  nread = ncoef_up_tot;
  for(i=1;i<=nread;i++){
    cre[i] = 0.0;
    cim[i] = 0.0;
    fcre[i] = 0.0;
    fcim[i] = 0.0;
  }/*endfor*/
  if(cp_lsda==1){
    cpcoeffs_pos->icoef_form_dn  = 0;
    cpcoeffs_pos->ivcoef_form_dn = 0;
    cpcoeffs_pos->ifcoef_form_dn = 1;
    cpcoeffs_pos->icoef_orth_dn  = 1;
    if(cp_norb>0)cpcoeffs_pos->icoef_orth_dn = 0;
    cpcoeffs_pos->ivcoef_orth_dn  = 1;
    if(cp_norb>0)cp->cpcoeffs_pos[ip].ivcoef_orth_dn = 0;

    cre = cpcoeffs_pos->cre_dn;
    cim = cpcoeffs_pos->cim_dn;
    fcre = cp->cpcoeffs_pos[ip].fcre_dn;
    fcim = cp->cpcoeffs_pos[ip].fcim_dn;
    nread = ncoef_dn_tot;
    for(i=1;i<=nread;i++){
      cre[i] = 0.0;
      cim[i] = 0.0;
      fcre[i] = 0.0;
      fcim[i] = 0.0;
    }/*endfor*/
  }/*endif:lsda*/


/*-----------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/

