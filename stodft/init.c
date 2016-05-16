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
/* II) Reinit communication group for new number of wave function           */
/*     (//N means don't remalloc these arrays)				    */

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
void reInitComm(CP *cp)
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*    Local Variables   */
  int irem,idiv,iii;

/*==========================================================================*/
/* I) Up states, state per process                                          */

  idiv = cp->cpcoeffs_info.nstate_up/cp->communicate.np_states;
  irem = (cp->cpcoeffs_info.nstate_up % cp->communicate.np_states);
  cp->cpcoeffs_info.nstate_up_proc = idiv;
  if(cp->communicate.myid_state < irem) {
     cp->cpcoeffs_info.nstate_up_proc = idiv+1;
  }/*endif*/
  if(cp->communicate.myid_state <= irem) {
    cp->cpcoeffs_info.istate_up_st = cp->communicate.myid_state*(idiv+1)+1;
  } else {
    cp->cpcoeffs_info.istate_up_st = irem*(idiv+1)
                                   + (cp->communicate.myid_state-irem)*idiv+1;
  }/*endif*/
  cp->cpcoeffs_info.istate_up_end = cp->cpcoeffs_info.istate_up_st +
                                    cp->cpcoeffs_info.nstate_up_proc-1;

  cp->cp_comm_state_pkg_up.nstate     = cp->cpcoeffs_info.nstate_up;
  cp->cp_comm_state_pkg_up.nstate_proc= cp->cpcoeffs_info.nstate_up_proc;


  irem             = (cp->cp_comm_state_pkg_up.nstate %
                      cp->cp_comm_state_pkg_up.num_proc);
  cp->cp_comm_state_pkg_up.nstate_proc_max  = (irem > 0 ? idiv+1 : idiv);
  cp->cp_comm_state_pkg_up.nstate_max = (irem > 0 ?
                          ((idiv+1)*cp->communicate.np_states) :
                          (idiv*cp->communicate.np_states)) ;

/*==========================================================================*/
/* I) Down states, state per process                                        */

  idiv = cp->cpcoeffs_info.nstate_dn/cp->communicate.np_states;
  irem = (cp->cpcoeffs_info.nstate_dn % cp->communicate.np_states);
  cp->cpcoeffs_info.nstate_dn_proc = idiv;
  if(cp->communicate.myid_state < irem) {
     cp->cpcoeffs_info.nstate_dn_proc = idiv+1;
  }/*endif*/
  if(cp->communicate.myid_state <= irem) {
    cp->cpcoeffs_info.istate_dn_st = cp->communicate.myid_state*(idiv+1)+1;
  } else {
    cp->cpcoeffs_info.istate_dn_st = irem*(idiv+1)
                                   + (cp->communicate.myid_state-irem)*idiv+1;
  }/*endif*/
  cp->cpcoeffs_info.istate_dn_end = cp->cpcoeffs_info.istate_dn_st +
                                    cp->cpcoeffs_info.nstate_dn_proc-1;

  cp->cp_comm_state_pkg_dn.nstate     = cp->cpcoeffs_info.nstate_dn;
  cp->cp_comm_state_pkg_dn.nstate_proc= cp->cpcoeffs_info.nstate_dn_proc;

  irem             = (cp->cp_comm_state_pkg_dn.nstate %
                      cp->cp_comm_state_pkg_dn.num_proc);
  cp->cp_comm_state_pkg_dn.nstate_proc_max  = (irem > 0 ? idiv+1 : idiv);
  cp->cp_comm_state_pkg_dn.nstate_max = (irem > 0 ?
                          ((idiv+1)*cp->communicate.np_states) :
                          (idiv*cp->communicate.np_states)) ;


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
/* I) Calculate the sizes */

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
/* II) Free everything */

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
/* III) Malloc the variables */

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
/* IV) Assign the offsets */

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
/* V) Initialize coeficient arrays and form flags                       */


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


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void reallocScratch(CP *cp)
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*    Local Variables   */
  CPSCR *cpscr = cp->cpscr;
  CPSCR_LOC    *cpscr_loc = cpscr->cpscr_loc;
  CPSCR_NONLOC *cpscr_nonloc = cpscr->cpscr_nonloc;
  CPSCR_RHO    *cpscr_rho = cpscr->cpscr_rho;
  CPSCR_GRHO   *cpscr_grho = cpscr->cpscr_grho;
  CPSCR_OVMAT  *cpscr_ovmat = cpscr->cpscr_ovmat;
  CPSCR_WAVE   *cpscr_wave = cpscr->cpscr_wave;
  CPSCR_THERM  *cpscr_therm = cpscr->cpscr_therm;
  CPSCR_DUAL_PME *cpscr_dual_pme = cpscr->cpscr_dual_pme;
  CPSCR_ATOM_PME *cpscr_atom_pme = cpscr->cpscr_atom_pme;
  CPSCR_WANNIER *cpscr_wannier = cpscr->cpscr_wannier;

  int nstate_up         = cpcoeffs_info->nstate_up_proc;
  int nstate_dn         = cpcoeffs_info->nstate_dn_proc;
  int nstate_up_tot     = cpcoeffs_info->nstate_up;
  int nstate_dn_tot     = cpcoeffs_info->nstate_dn;
  int ncoef             = cpcoeffs_info->ncoef;
  int ncoef_l           = cpcoeffs_info->ncoef_l;
  int pi_beads          = cpcoeffs_info->pi_beads;
  int cp_laplacian_on   = cpcoeffs_info->cp_laplacian_on;
  int cp_tau_functional = cpcoeffs_info->cp_tau_functional;
  int num_c_nhc_proc    = cptherm_info->num_c_nhc_proc;
  int massiv_flag       = cptherm_info->massiv_flag;
  int n_ang_max         = pseudo->n_ang_max;
  int n_rad_max         = pseudo->n_rad_max;
  int natm_nls_max      = cpscr->cpscr_nonloc.natm_nls_max;
  int cp_lsda           = cpopts->cp_lsda;
  int cp_ptens_calc     = cpopts->cp_ptens_calc;
  int cp_hess_calc      = cpopts->cp_hess_calc;
  int cp_gga            = cpopts->cp_gga;
  int cp_ke_dens_on     = cpcoeffs_info->cp_ke_dens_on;
  int cp_elf_calc_frq   = cpcoeffs_info->cp_elf_calc_frq;
  int cp_norb           = cpopts->cp_norb;
  int cp_para_opt       = cpopts->cp_para_opt;
  int np_states         = cp_comm_state_pkg_up->num_proc;
  int nfft_up_proc      = cp_para_fft_pkg3d_lg->nfft_proc;
  int nfft_up           = cp_para_fft_pkg3d_lg->nfft;
  int nstate_max_up     = cp_comm_state_pkg_up->nstate_max;
  int nstate_ncoef_proc_max_up = cp_comm_state_pkg_up->nstate_ncoef_proc_max;
  int nstate_max_dn     = cp_comm_state_pkg_up->nstate_max;
  int nstate_ncoef_proc_max_dn = cp_comm_state_pkg_up->nstate_ncoef_proc_max;
  int num_c_nhc1        = num_c_nhc_proc+1;

  int ncoef_l_pme_dual,ncoef_l_pme_dual_proc;
  int ncoef_l_proc_max_mall;
  int ncoef_l_proc_max_mall_ke;

  int ncoef_l_dens_cp_box;
  int ncoef_l_proc_max_dens_cp_box;
  int nfft_up_proc_dens_cp_box,nfft_up_dens_cp_box,nfft2_up_dens_cp_box;
  int nfft_dn_proc_dens_cp_box,nfft_dn_dens_cp_box ;
  int ncoef_l_proc_max_dn;

  int nkf1_cp_box = cp_para_fft_pkg3d_dens_cp_box->nkf1;
  int nkf2_cp_box = cp_para_fft_pkg3d_dens_cp_box->nkf2;
  int nkf3_cp_box = cp_para_fft_pkg3d_dens_cp_box->nkf3;
  int n_interp_pme_dual = pseudo->n_interp_pme_dual;

/*--------------------------------------------------------------------------*/
/*         Local variable declarations                                      */

  double now_memory;
  int i,iii;
  int nlscr_up,nlscr_dn,nlscr_up_pv,nlscr_dn_pv,ncoef_l_pv,ncoef_l_proc_max;
  int ncoef_l_proc_max_mall_cp_box,ncoef_l_proc_max_mall_cp_box_dn;
  int nfft2_mall_up,nfft2_mall_dn,nfft2_mall_up_proc,nfft2_mall_dn_proc;
  int nfft2_up,nfft2_dn,irem;
  int nfft2_up_proc,nfft2_dn_proc,nfft_dn;
  int nfft2_up_gga,nfft2_dn_gga,nlap_up,nlap_dn,nlap_up_ptens;
  int nfft2_up_ke_dens,nfft2_dn_ke_dens;
  int nlap_dn_ptens,ngga_up,ngga_dn,nlap_g_up,nlap_g_dn;
  int nlap_g_up_ptens, nlap_g_dn_ptens;
  int ncoef_up,ncoef_dn;
  int par_size_up,par_size_dn;
  int ncoef2_up_c,ncoef2_dn_c;
  int ncoef2_up,ncoef2_dn,ncoef2_up_spec,ncoef2_dn_spec;
  int ncoef2_up_par,nstate,nstate2,nstate_tot,nstate2_tot;
  int num=0;
  int zero=0;
  int ncoef_l_mall_proc_max_dual,ncoef_l_mall_proc_max_dual_dn;
  int map_count;
  int mtemp;
  int nlen_pme,pme_nkf3,ninterp_pme,nmall;
  int mall_size;

  int cp_wan_opt        = cpopts->cp_wan_opt;
  int cp_wan_min_opt    = cpopts->cp_wan_min_opt;
  int cp_wan_init_opt   = cpopts->cp_wan_init_opt;
  int ndim_wannier;
  int mm=5;

/*==========================================================================*/
/* Free everything*/

  free(&(cpscr_nonloc->vnlre_up[1]));
  free(&(cpscr_nonloc->vnlim_up[1]));
  free(&(cpscr_nonloc->vnlre_dn[1]));
  free(&(cpscr_nonloc->vnlim_dn[1]));

  free(&(cpscr_nonloc->dvnlre_x_up[1]));
  free(&(cpscr_nonloc->dvnlre_y_up[1]));
  free(&(cpscr_nonloc->dvnlre_z_up[1]));
  free(&(cpscr_nonloc->dvnlim_x_up[1]));
  free(&(cpscr_nonloc->dvnlim_y_up[1]));
  free(&(cpscr_nonloc->dvnlim_z_up[1]));

  free(&(cpscr_nonloc->dvnlre_x_dn[1]));
  free(&(cpscr_nonloc->dvnlre_y_dn[1]));
  free(&(cpscr_nonloc->dvnlre_z_dn[1]));
  free(&(cpscr_nonloc->dvnlim_x_dn[1]));
  free(&(cpscr_nonloc->dvnlim_y_dn[1]));
  free(&(cpscr_nonloc->dvnlim_z_dn[1]));

  free(&(cpscr_nonloc->dvnlre_gxgx_up[1]));
  free(&(cpscr_nonloc->dvnlre_gygy_up[1]));
  free(&(cpscr_nonloc->dvnlre_gzgz_up[1]));
  free(&(cpscr_nonloc->dvnlre_gxgy_up[1]));
  free(&(cpscr_nonloc->dvnlre_gygz_up[1]));
  free(&(cpscr_nonloc->dvnlre_gxgz_up[1]));
  free(&(cpscr_nonloc->dvnlim_gxgx_up[1]));
  free(&(cpscr_nonloc->dvnlim_gygy_up[1]));
  free(&(cpscr_nonloc->dvnlim_gzgz_up[1]));
  free(&(cpscr_nonloc->dvnlim_gxgy_up[1]));
  free(&(cpscr_nonloc->dvnlim_gygz_up[1]));
  free(&(cpscr_nonloc->dvnlim_gxgz_up[1]));

  free(&(cpscr_nonloc->dvnlre_gxgx_dn[1]));
  free(&(cpscr_nonloc->dvnlre_gygy_dn[1]));
  free(&(cpscr_nonloc->dvnlre_gzgz_dn[1]));
  free(&(cpscr_nonloc->dvnlre_gxgy_dn[1]));
  free(&(cpscr_nonloc->dvnlre_gygz_dn[1]));
  free(&(cpscr_nonloc->dvnlre_gxgz_dn[1]));
  free(&(cpscr_nonloc->dvnlim_gxgx_dn[1]));
  free(&(cpscr_nonloc->dvnlim_gygy_dn[1]));
  free(&(cpscr_nonloc->dvnlim_gzgz_dn[1]));
  free(&(cpscr_nonloc->dvnlim_gxgy_dn[1]));
  free(&(cpscr_nonloc->dvnlim_gygz_dn[1]));
  free(&(cpscr_nonloc->dvnlim_gxgz_dn[1]));

  free(&(cpscr_wave->cre_up[1]));
  free(&(cpscr_wave->cim_up[1]));
  free(&(cpscr_wave->cre_dn[1]));
  free(&(cpscr_wave->cim_dn[1]));

  free(&(cpscr_ovmat->ovlap1[1]));
  free(&(cpscr_ovmat->ovlap2[1]));
  free(&(cpscr_ovmat->ovlap3[1]));
  free(&(cpscr_ovmat->ovlap4[1]));
  free(&(cpscr_ovmat->ovlap5[1]));
  free(&(cpscr_ovmat->ovlap6[1]));
  free(&(cpscr_ovmat->ovlap7[1]));
  free(&(cpscr_ovmat->ovlap8[1]));

  free(&(cpscr_ovmat->state_vec1[1]));
  free(&(cpscr_ovmat->state_vec2[1]));
  free(&(cpscr_ovmat->state_vec3[1]));
  free(&(cpscr_ovmat->state_vec4[1]));
  free(&(cpscr_ovmat->state_vec5[1]));
  free(&(cpscr_ovmat->state_vec6[1]));
  free(&(cpscr_ovmat->state_vec7[1]));

  if(cpscr_wannier->cp_wannier_on==1){
    free(cpscr_wannier->A[1]);
    free(cpscr_wannier->R_real[1]);
    free(cpscr_wannier->R_imag[1]);
    free(cpscr_wannier->U_real[1]);
    free(cpscr_wannier->U_tmp1[1]);
    free(cpscr_wannier->U_tmp2[1]);
    free(cpscr_wannier->U_tmp3[1]);
    free(cpscr_wannier->U_tmp4[1]);
    free(cpscr_wannier->Bt_real[1]);
    free(cpscr_wannier->Bt_imag[1]);
    free(cpscr_wannier->M_real[1]);
    free(cpscr_wannier->phi[1]);
    free(&(cpscr_wannier->A[1]));
    free(&(cpscr_wannier->R_real[1]));
    free(&(cpscr_wannier->R_imag[1]));
    free(&(cpscr_wannier->U_real[1]));
    free(&(cpscr_wannier->U_tmp1[1]));
    free(&(cpscr_wannier->U_tmp2[1]));
    free(&(cpscr_wannier->U_tmp3[1]));
    free(&(cpscr_wannier->U_tmp4[1]));
    free(&(cpscr_wannier->Bt_real[1]));
    free(&(cpscr_wannier->Bt_imag[1]));
    free(&(cpscr_wannier->M_real[1]));
    free(&(cpscr_wannier->phi[1]));

    free(&(cpscr_wannier->real));
    free(&(cpscr_wannier->imag));
    free(&(cpscr_wannier->D));
    free(&(cpscr_wannier->norm));
  }

/*=========================================================================*/
/* I) Malloc size calculation */

 /*-------------------------------------------------------------------------*/
 /* i) Dual grid CP : Define the small dense grid sizes */

  if(cp_dual_grid_opt_on >= 1){

    nfft_up_proc_dens_cp_box   = cp_para_fft_pkg3d_dens_cp_box->nfft_proc;
    nfft_up_dens_cp_box        = cp_para_fft_pkg3d_dens_cp_box->nfft;
    nfft2_up_dens_cp_box       = nfft_up_dens_cp_box/2;

    nfft_dn_proc_dens_cp_box   = (cp_lsda == 1 ? nfft_up_proc_dens_cp_box : 0);
    nfft_dn_dens_cp_box        = (cp_lsda == 1 ? nfft_up_dens_cp_box : 0);

    ncoef_l_dens_cp_box          = cpcoeffs_info->ncoef_l_dens_cp_box;
    ncoef_l_proc_max_dens_cp_box = ncoef_l_dens_cp_box/np_states;
    irem                         = (ncoef_l_proc_max_dens_cp_box % np_states);
    if(irem>0){ncoef_l_proc_max_dens_cp_box++;}

  }/*endif cp_dual_grid_opt_on */

 if(cp_dual_grid_opt_on == 2){
   ncoef_l_pme_dual = cp_para_fft_pkg3d_lg->ncoef;
 }/*endif cp_dual_grid_opt_on*/

 /*-------------------------------------------------------------------------*/
 /* ii) Normal CP : Define the grid size              */
 /*     Dual   CP : Define the large sparse grid size */

  nfft2_up      = nfft_up/2;
  nfft2_up_proc = nfft_up_proc/2;

  nfft_dn       = (cp_lsda == 1 ? nfft_up : 0);
  nfft2_dn      = (cp_lsda == 1 ? nfft2_up : 0);
  nfft2_dn_proc = (cp_lsda == 1 ? nfft2_up_proc : 0);

 /*-------------------------------------------------------------------------*/
 /* iii) Choose the correct size for your application                       */
 /*      This is always the small dense grid                                */

  nfft2_mall_up      = (cp_dual_grid_opt_on >= 1 ?
                        nfft_up_dens_cp_box/2 : nfft2_up);
  nfft2_mall_up_proc = (cp_dual_grid_opt_on >= 1 ?
                        nfft_up_proc_dens_cp_box/2 : nfft2_up_proc);

  nfft2_mall_dn      = (cp_dual_grid_opt_on >= 1 ?
                        nfft_dn_dens_cp_box/2 : nfft2_dn);
  nfft2_mall_dn_proc = (cp_dual_grid_opt_on >= 1 ?
                        nfft_dn_proc_dens_cp_box/2 : nfft2_dn_proc);

 /*-------------------------------------------------------------------------*/
 /* iv) Wave function size (spherically cutoff small g-space for dense box) */

  ncoef_up  = ncoef;
  ncoef_dn  = (cp_lsda == 1 ? ncoef : 0);

 /*-------------------------------------------------------------------------*/
 /* v) Normal CP: The sphere cut large g-space for the dense box            */
 /*    Dual CP  : The sphere cut large g-space for the large sparse box     */

  ncoef_l_proc_max = ncoef_l/np_states;
  irem = (ncoef_l % np_states);
  if(irem>0){ncoef_l_proc_max++;}

 if(cp_dual_grid_opt_on == 2){
   ncoef_l_pme_dual_proc = ncoef_l_pme_dual/np_states;
   irem = (ncoef_l_pme_dual_proc % np_states);
   if(irem>0){ncoef_l_pme_dual_proc++;}
  }/*endif cp_dual_grid_opt_on */

  ncoef_l_proc_max_mall = (cp_dual_grid_opt_on == 2 ? ncoef_l_pme_dual_proc
                                                    : ncoef_l_proc_max);
  ncoef_l_proc_max_mall_ke = (cp_ke_dens_on == 1 ? ncoef_l_proc_max_mall:0);
  ncoef_l_proc_max_dn      = (cp_lsda == 1 ? ncoef_l_proc_max_mall : 0);

 /*-------------------------------------------------------------------------*/
 /* vi) Choose the large g-space malloc size based on the dual or normal opt*/
 /*     The malloc size is the always the dense grid.                       */
 /*     Set special dual malloc sizes to zero to avoid mallocing extra      */
 /*     memory during normal CP.                                            */

  if(cp_dual_grid_opt_on==0){

    ncoef_l_proc_max_mall_cp_box = ncoef_l_proc_max;

  }else{

    ncoef_l_proc_max_mall_cp_box = ncoef_l_dens_cp_box/np_states;
    irem = (ncoef_l_dens_cp_box % np_states);
    if(irem>0){ncoef_l_proc_max_mall_cp_box++;}
    ncoef_l_proc_max_mall_cp_box_dn =
                   (cp_lsda == 1 ? ncoef_l_proc_max_mall_cp_box : 0);

  }/*endif cp_dual_grid_opt_on*/

  ncoef_l_mall_proc_max_dual    = (cp_dual_grid_opt_on >= 1 ?
                                   ncoef_l_proc_max_mall_cp_box : 0);
  ncoef_l_mall_proc_max_dual_dn = (cp_dual_grid_opt_on >= 1 ?
                                   ncoef_l_proc_max_mall_cp_box_dn : 0);

 /*-------------------------------------------------------------------*/
 /* vii) Wavefunction scratch sizes : always on small dense grid     */

  par_size_up = nstate_max_up*nstate_ncoef_proc_max_up;
  par_size_dn = nstate_max_dn*nstate_ncoef_proc_max_dn;

  if(massiv_flag==0){
    ncoef2_up_c = MAX(ncoef_up*nstate_up,num_c_nhc_proc);
  }else{
    ncoef2_up_c = MAX(ncoef_up*nstate_up,2*ncoef_up);
  }/*endif*/

  if(np_states>1){ncoef2_up_c = MAX(ncoef2_up_c,par_size_up);}
  ncoef2_up      = ncoef_up*nstate_up;
  ncoef2_dn      = ncoef_dn*nstate_dn;
  ncoef2_up      = MAX(ncoef2_up,par_size_up);
  ncoef2_dn      = MAX(ncoef2_dn,par_size_dn);
  ncoef2_up_spec = MAX(ncoef2_up,ncoef2_dn);
  ncoef2_dn_spec = ncoef2_up_spec;
  if(np_states == 1 &&  cp_norb==0){ncoef2_up_spec = 1;}
  ncoef2_up_par = ncoef2_up_spec;
  if(np_states == 1){ncoef2_up_par = 1;}

  nstate        = MAX(nstate_up, nstate_dn);
  nstate2       = nstate*nstate;
  nstate_tot    = MAX(nstate_up_tot,nstate_dn_tot);
  nstate2_tot   = nstate_tot*nstate_tot;

 /*-------------------------------------------------------------------*/
 /* viii) Nonlocal sizes : always on small dense grid */

   nlscr_up  = nstate_up*(n_ang_max+1)*(n_ang_max+1)
                        *(n_rad_max)*natm_nls_max;
  nlscr_dn  = 0;
  if(cp_lsda==1){
     nlscr_dn = nstate_dn*(n_ang_max+1)*(n_ang_max+1)
                         *(n_rad_max)*(natm_nls_max);
  }/*endif*/
  nlscr_up_pv = 0;
  nlscr_dn_pv = 0;
  ncoef_l_pv  = 0;
  if(cp_ptens_calc == 1){
    nlscr_up_pv = nlscr_up;
    nlscr_dn_pv = nlscr_dn;
    ncoef_l_pv  = ncoef_l_proc_max_mall;
  }/* endif */
  if(atm_hess_calc == 3){
    nlscr_up_pv = nlscr_up;
    nlscr_dn_pv = nlscr_dn;
  }/* endif */

 /*-------------------------------------------------------------------*/
 /* ix) GGA sizes : Always the small dense grid                       */

  nfft2_up_gga  = ((cp_gga == 1 || cp_elf_calc_frq > 0) ? nfft2_mall_up_proc:0);
  nfft2_dn_gga  = (((cp_gga == 1 || cp_elf_calc_frq > 0) && cp_lsda == 1)
                ? nfft2_mall_dn_proc:0);
  nfft2_up_ke_dens = (cp_ke_dens_on == 1 ? nfft2_mall_up_proc:0);
  nfft2_dn_ke_dens = ((cp_ke_dens_on == 1 && cp_lsda == 1) ? nfft2_mall_dn_proc:0);
  nlap_up       = (cp_laplacian_on == 1 ? nfft2_up_gga : 1);
  nlap_dn       = (cp_laplacian_on == 1 ? nfft2_dn_gga : 1);
  nlap_up_ptens = (cp_laplacian_on == 1&&cp_ptens_calc == 1?nfft2_up_gga : 1);
  nlap_dn_ptens = (cp_laplacian_on == 1&&cp_ptens_calc == 1?nfft2_dn_gga : 1);

  nlap_up       = ( cp_laplacian_on == 1 ? nfft2_up_gga : 1);
  nlap_dn       = ( cp_laplacian_on == 1 ? nfft2_dn_gga : 1);
  nlap_up_ptens = ((cp_laplacian_on == 1 && cp_ptens_calc == 1) ?
                                           nfft2_up_gga : 1);
  nlap_dn_ptens = ((cp_laplacian_on == 1 && cp_ptens_calc == 1) ?
                                           nfft2_dn_gga : 1);

  ngga_up         = ((cp_gga == 1 || cp_elf_calc_frq > 0)
                  ? ncoef_l_proc_max_mall_cp_box:0);
  ngga_dn         = (cp_lsda == 1 ? ngga_up:0);
  nlap_g_up       = (cp_laplacian_on == 1 ? ngga_up : 1);
  nlap_g_dn       = (cp_laplacian_on == 1 ? ngga_dn : 1);
  nlap_g_up_ptens = (cp_laplacian_on == 1 && cp_ptens_calc == 1 ? ngga_up : 1);
  nlap_g_dn_ptens = (cp_laplacian_on == 1 && cp_ptens_calc == 1 ? ngga_dn : 1);

  /*------------------------------------------------------------------------*/
  /* x) Wannier scratch size   */

  if(cp_wan_opt ==1 || cp_wan_min_opt==1 || cp_wan_init_opt==1){
    cpscr->cpscr_wannier.cp_wannier_on=1;
    ndim_wannier=nstate_up_tot*nstate_up_tot;
  }else{
    cpscr->cpscr_wannier.cp_wannier_on=0;
    ndim_wannier=0;
  }

/*===========================================================================*/
/* II) Malloc the vectors  */

/*------------------------------------------------------------------*/
/* Non_local  : Always small dense grid                             */

  cpscr->cpscr_nonloc.vnlre_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr->cpscr_nonloc.vnlim_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;

  cpscr->cpscr_nonloc.vnlre_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr->cpscr_nonloc.vnlim_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;

  cpscr->cpscr_nonloc.dvnlre_x_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_y_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_z_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_x_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_y_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_z_up
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;

  num += 10*nlscr_up;

  cpscr->cpscr_nonloc.dvnlre_x_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_y_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_z_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_x_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_y_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_z_dn
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;

  num += 6*nlscr_dn;

  cpscr->cpscr_nonloc.dvnlre_gxgx_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gygy_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gzgz_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gxgy_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gygz_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gxgz_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gxgx_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gygy_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gzgz_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gxgy_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gygz_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gxgz_up
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;

  num += 12*nlscr_up_pv;

  cpscr->cpscr_nonloc.dvnlre_gxgx_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gygy_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gzgz_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gxgy_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gygz_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gxgz_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;

  cpscr->cpscr_nonloc.dvnlim_gxgx_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gygy_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gzgz_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gxgy_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gygz_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gxgz_dn
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;

  num += 12*nlscr_dn_pv;

/*------------------------------------------------------------------*/
/* wave */
  cpscr->cpscr_wave.cre_up
                    = (double *)cmalloc(ncoef2_up_c*sizeof(double))-1;
  cpscr->cpscr_wave.cim_up
                    = (double *)cmalloc(ncoef2_up*sizeof(double))-1;
  cpscr->cpscr_wave.cre_dn
                    = (double *)cmalloc(ncoef2_dn*sizeof(double))-1;
  cpscr->cpscr_wave.cim_dn
                    = (double *)cmalloc(ncoef2_dn*sizeof(double))-1;


/*------------------------------------------------------------------*/
/* ovmat */
  cpscr->cpscr_ovmat.ovlap1
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.ovlap2
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.ovlap3
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.ovlap4
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.ovlap5
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.ovlap6
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.ovlap7
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.ovlap8
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;

  num += 8*nstate2_tot;

  cpscr->cpscr_ovmat.state_vec1
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.state_vec2
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.state_vec3
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.state_vec4
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.state_vec5
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.state_vec6
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.state_vec7
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;

  num += 7*nstate_tot;

/*-----------------------------------------------------------------------*/
/* Wannier scratch                                                       */

  if(cpscr->cpscr_wannier.cp_wannier_on==1){
    cpscr->cpscr_wannier.U_final = (double *) cmalloc(ndim_wannier*sizeof(double))-1;
    cpscr->cpscr_wannier.g       = (double *) cmalloc(ndim_wannier*sizeof(double))-1;
    cpscr->cpscr_wannier.gg      = (double *) cmalloc(ndim_wannier*sizeof(double))-1;
    cpscr->cpscr_wannier.diag    = (double *) cmalloc(ndim_wannier*sizeof(double))-1;
    cpscr->cpscr_wannier.scr     =
                   (double *) cmalloc((ndim_wannier*(2*mm+1)+2*mm)*sizeof(double))-1;
    cpscr->cpscr_wannier.iprint  = (int *)cmalloc(2*sizeof(int))-1;


    num += 4*ndim_wannier+ndim_wannier*(2*mm+1)+2*mm;

    cpscr->cpscr_wannier.A       = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.R_real  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.R_imag  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.U_real  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.U_imag  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.U_tmp1  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.U_tmp2  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.U_tmp3  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.U_tmp4  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.Bt_real = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.Bt_imag = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.M_real  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);

    num += 12*ndim_wannier;

    cpscr->cpscr_wannier.real    = (double *)cmalloc(nstate_up_tot*sizeof(double))-1;
    cpscr->cpscr_wannier.imag    = (double *)cmalloc(nstate_up_tot*sizeof(double))-1;
    cpscr->cpscr_wannier.D       = (double *)cmalloc(nstate_up_tot*sizeof(double))-1;
    cpscr->cpscr_wannier.norm    = (double *)cmalloc(nstate_up_tot*sizeof(double))-1;
    cpscr->cpscr_wannier.phi     = cmall_mat(1,nstate_up_tot,1,3);
    cpscr->cpscr_wannier.HMatrix = cmall_mat(1,3,1,3);
    printf("CP_WANNIER memory allocation \n");
    
    num += 7*nstate_up_tot;
  } 

  
/*-----------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/

