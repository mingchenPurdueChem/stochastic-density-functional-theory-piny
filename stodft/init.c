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







