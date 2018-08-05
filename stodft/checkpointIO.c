/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: checkpointIO.c                                 */
/*                                                                          */
/*   This routine reads and writes checkpoint file			    */
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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void checkpointOutput(CP *cp, GENERAL_DATA *general_data)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*************************************************************************/
/* This function output the checkpoint file, including current SCF step, */
/* number of steps in diis bank, diis and error bank.                    */
/*************************************************************************/

  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  COMMUNICATE *communicate      = &(cp->communicate);
  CPSCR *cpscr			= &(cp->cpscr);
  CELL *cell			= &(general_data->cell);
  CPOPTS *cpopts		= &(cp->cpopts);
 
  int myidState		= communicate->myid_state;
  int numProcStates	= communicate->np_states;
  int numDiis		= stodftInfo->numDiis;
  int numDiisNow	= stodftInfo->numDiisNow;
  int numDiisOutput;
  int rhoRealGridNum	= stodftInfo->rhoRealGridNum;
  int rhoRealGridTot	= stodftInfo->rhoRealGridTot;
  int cpLsda            = cpopts->cp_lsda;
  int iScf		= stodftInfo->iScf;
  int iDiis,iGrid;
  MPI_Comm commStates = communicate->comm_states;
 
  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;
    
  FILE *fileCheckpoint;

  double vol;
  double *hmat_cp        = cell->hmat_cp;
  double *rhoTemp;
  double *rhoUp	     = cpscr->cpscr_rho.rho_up;
  double *rhoDn          = cpscr->cpscr_rho.rho_dn;  
  double **rhoUpBank = stodftCoefPos->rhoUpBank;
  double **rhoDnBank = stodftCoefPos->rhoDnBank;
  double **rhoUpErr  = stodftCoefPos->rhoUpErr;
  double **rhoDnErr  = stodftCoefPos->rhoDnErr;

/*======================================================================*/
/* I) First output SCF step and DIIS step                               */

  if(myidState==0){
    fileCheckpoint = fopen("density-checkpoint","w");
    fprintf(fileCheckpoint,"%i %i\n",iScf,numDiisNow);  
  }

/*======================================================================*/
/* II) Output current density in r space	                        */

  if(myidState==0)rhoTemp = (double*)cmalloc(rhoRealGridTot*sizeof(double));
  
  if(numProcStates>1){
    Barrier(commStates);
    Gatherv(&rhoUp[1],rhoRealGridNum,MPI_DOUBLE,rhoTemp,
	    rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,0,commStates);
    Barrier(commStates);
    printf("myidState %i rhoRealGridNum %i rhoRealSendCounts %i rhoRealDispls %i rhoUp[1] %lg\n",
	   myidState,rhoRealGridNum,rhoRealSendCounts[myidState],rhoRealDispls[myidState],rhoUp[1]);
    if(myidState==0)printf("rhoTemp %lg %lg %lg\n",rhoTemp[0],rhoTemp[200000],rhoTemp[400000]);
  } 
  else{
    memcpy(rhoTemp,&rhoUp[1],rhoRealGridTot*sizeof(double));
  } 
  if(myidState==0){
    // Scale it back to original volumn factor
    vol = getdeth(hmat_cp);
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      rhoTemp[iGrid] *= vol;
      fprintf(fileCheckpoint,"%.16lg\n",rhoTemp[iGrid]);
    }
  }
  if(cpLsda==1){
    if(numProcStates>1){
      Barrier(commStates);
      Gatherv(&rhoDn[1],rhoRealGridNum,MPI_DOUBLE,rhoTemp,
	      rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,0,commStates);
      Barrier(commStates);
    }
    else{
      memcpy(rhoTemp,&rhoDn[1],rhoRealGridTot*sizeof(double));
    }
    if(myidState==0){
      // Scale it back to original volumn factor
      vol = getdeth(hmat_cp);
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	rhoTemp[iGrid] *= vol;
	fprintf(fileCheckpoint,"%.16lg\n",rhoTemp[iGrid]);
      }
    }
  }
  

/*======================================================================*/
/* III) Output DIIS bank and error bank		                        */

  
  if(numDiisNow<numDiis)numDiisOutput = numDiisNow+1;
  else numDiisOutput = numDiisNow;

  for(iDiis=0;iDiis<numDiisOutput;iDiis++){
    if(numProcStates>1){
      Barrier(commStates);
      Gatherv(rhoUpBank[iDiis],rhoRealGridNum,MPI_DOUBLE,rhoTemp,
	      rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,0,commStates);
      Barrier(commStates);
    }
    else{
      memcpy(rhoTemp,rhoUpBank[iDiis],rhoRealGridTot*sizeof(double));
    }
    if(myidState==0){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	fprintf(fileCheckpoint,"%.16lg\n",rhoTemp[iGrid]);
      }
    }
  }
  for(iDiis=0;iDiis<numDiisNow;iDiis++){
    if(numProcStates>1){
      Barrier(commStates);
      Gatherv(rhoUpErr[iDiis],rhoRealGridNum,MPI_DOUBLE,rhoTemp,
              rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,0,commStates);
      Barrier(commStates);
    }
    else{
      memcpy(rhoTemp,rhoUpErr[iDiis],rhoRealGridTot*sizeof(double));
    }
    if(myidState==0){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
        fprintf(fileCheckpoint,"%.16lg\n",rhoTemp[iGrid]);
      }
    }
  }
  if(cpLsda==1){
    for(iDiis=0;iDiis<numDiisOutput;iDiis++){
      if(numProcStates>1){
        Barrier(commStates);
	Gatherv(rhoDnBank[iDiis],rhoRealGridNum,MPI_DOUBLE,rhoTemp,
		rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,0,commStates); 
        Barrier(commStates);
      }
      else{
	memcpy(rhoTemp,rhoDnBank[iDiis],rhoRealGridTot*sizeof(double));
      }//endif
      if(myidState==0){
	for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	  fprintf(fileCheckpoint,"%.16lg\n",rhoTemp[iGrid]);
	}//endfor
      }//endif
    }//endfor iDiis
    for(iDiis=0;iDiis<numDiisNow;iDiis++){
      if(numProcStates>1){
	Barrier(commStates);
	Gatherv(rhoDnErr[iDiis],rhoRealGridNum,MPI_DOUBLE,rhoTemp,
		rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,0,commStates);
	Barrier(commStates);
      }
      else{
	memcpy(rhoTemp,rhoDnErr[iDiis],rhoRealGridTot*sizeof(double));
      }//endif
      if(myidState==0){
	for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	  fprintf(fileCheckpoint,"%.16lg\n",rhoTemp[iGrid]);
	}//endfor iGrid
      }//endif myidState
    }//endfor iDiis
  }//endif cpLsda
  

  if(myidState==0){
    fclose(fileCheckpoint);
  }

  if(myidState==0)free(rhoTemp);

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void checkpointInput(CP *cp,GENERAL_DATA *general_data,CLASS *class)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*************************************************************************/
/* This function output the checkpoint file, including current SCF step, */
/* number of steps in diis bank, diis and error bank.                    */
/*************************************************************************/

  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  COMMUNICATE *communicate      = &(cp->communicate);
  CPEWALD *cpewald		= &(cp->cpewald);
  CPSCR *cpscr			= &(cp->cpscr);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS  *cpcoeffs_pos   = &(cp->cpcoeffs_pos[1]);
  CPOPTS *cpopts		= &(cp->cpopts);
  EWALD *ewald			= &(general_data->ewald);
  CELL *cell		        = &(general_data->cell);
  PSEUDO *pseudo		= &(cp->pseudo);
  CLATOMS_POS *clatoms_pos      = &(class->clatoms_pos[1]);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  ATOMMAPS *atommaps            = &(class->atommaps);
  EWD_SCR      *ewd_scr         = &(class->ewd_scr);
  FOR_SCR      *for_scr         = &(class->for_scr);  
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

  int myidState         = communicate->myid_state;
  int numProcStates     = communicate->np_states;
  int numDiisNow,numDiisOutput;
  int numDiis		= stodftInfo->numDiis;
  int rhoRealGridNum    = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot    = stodftInfo->rhoRealGridTot;
  int cpGGA		= cpopts->cp_gga;
  int cpLsda            = cpopts->cp_lsda;  
  int cpDualGridOptOn	= cpopts->cp_dual_grid_opt;
  int numInterpPmeDual	= pseudo->n_interp_pme_dual;
  int numCoeff		= cpcoeffs_info->ncoef;
  int iScf,iDiis,iGrid,iCoeff;

  MPI_Comm commStates = communicate->comm_states;

  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;

  FILE *fileCheckpoint;

  double *rhoTemp = (double*)cmalloc(rhoRealGridTot*sizeof(double));
  double *rhoCoeffReUp   = cpscr->cpscr_rho.rhocr_up;
  double *rhoCoeffImUp   = cpscr->cpscr_rho.rhoci_up;
  double *rhoCoeffReDn   = cpscr->cpscr_rho.rhocr_dn;
  double *rhoCoeffImDn   = cpscr->cpscr_rho.rhoci_dn;
  double *rhoUp		 = cpscr->cpscr_rho.rho_up;
  double *rhoDn          = cpscr->cpscr_rho.rho_dn;
  double *rhoCoeffReUpDensCpBox = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
  double *rhoCoeffImUpDensCpBox = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
  double *rhoCoeffReDnDensCpBox = cpscr->cpscr_rho.rhocr_dn_dens_cp_box;
  double *rhoCoeffImDnDensCpBox = cpscr->cpscr_rho.rhoci_dn_dens_cp_box;
  double *divRhoxUp       = cpscr->cpscr_grho.d_rhox_up;
  double *divRhoyUp       = cpscr->cpscr_grho.d_rhoy_up;
  double *divRhozUp       = cpscr->cpscr_grho.d_rhoz_up;
  double *d2RhoUp         = cpscr->cpscr_grho.d2_rho_up;
  double *divRhoxDn       = cpscr->cpscr_grho.d_rhox_dn;
  double *divRhoyDn       = cpscr->cpscr_grho.d_rhoy_dn;
  double *divRhozDn       = cpscr->cpscr_grho.d_rhoz_dn;
  double *d2RhoDn         = cpscr->cpscr_grho.d2_rho_dn;
  double *testWfMaxRe	  = stodftCoefPos->testWfMaxRe;
  double *testWfMaxIm	  = stodftCoefPos->testWfMaxIm;
  double *testWfMinRe     = stodftCoefPos->testWfMinRe;
  double *testWfMinIm     = stodftCoefPos->testWfMinIm;
  double *rhoUpOld = stodftCoefPos->rhoUpOld;
  double *rhoDnOld = stodftCoefPos->rhoDnOld;

  double **rhoUpBank = stodftCoefPos->rhoUpBank;
  double **rhoDnBank = stodftCoefPos->rhoDnBank;
  double **rhoUpErr  = stodftCoefPos->rhoUpErr;
  double **rhoDnErr  = stodftCoefPos->rhoDnErr;

/*======================================================================*/
/* I) Read SCF step and DIIS step                                       */

  if(myidState==0){
    fileCheckpoint = fopen("density-checkpoint","r");
    fscanf(fileCheckpoint,"%i",&stodftInfo->iScf);
    fscanf(fileCheckpoint,"%i",&stodftInfo->numDiisNow);   
  }
  Bcast(&stodftInfo->iScf,1,MPI_INT,0,commStates);
  Bcast(&stodftInfo->numDiisNow,1,MPI_INT,0,commStates);

  numDiisNow = stodftInfo->numDiisNow;
  iScf = stodftInfo->iScf;

/*======================================================================*/
/* II) Read real space density                                          */

  if(myidState==0){
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      fscanf(fileCheckpoint,"%lg",&rhoTemp[iGrid]);
    }
  }
  if(numProcStates>1){
    Barrier(commStates);
    Scatterv(rhoTemp,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	     &rhoUp[1],rhoRealGridNum,MPI_DOUBLE,0,commStates);
    Barrier(commStates);
  }
  else{
    memcpy(&rhoUp[1],rhoTemp,rhoRealGridTot*sizeof(double));
  }
  memcpy(rhoUpOld,&rhoUp[1],rhoRealGridNum*sizeof(double));
  calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
                     rhoCoeffReUp,rhoCoeffImUp,rhoUp,
		     rhoCoeffReUpDensCpBox,rhoCoeffImUpDensCpBox,
                     divRhoxUp,divRhoyUp,divRhozUp,d2RhoUp,cpGGA,
		     cpDualGridOptOn,numInterpPmeDual,
                     communicate,&(cp->cp_para_fft_pkg3d_lg),
		     &(cp->cp_para_fft_pkg3d_dens_cp_box));
  if(cpLsda==1){
    if(myidState==0){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	fscanf(fileCheckpoint,"%lg",&rhoTemp[iGrid]);
      }
    }
    if(numProcStates>1){
      Barrier(commStates);
      Scatterv(rhoTemp,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	       &rhoDn[1],rhoRealGridNum,MPI_DOUBLE,0,commStates);
      Barrier(commStates);
    }
    else{
      memcpy(&rhoDn[1],rhoTemp,rhoRealGridTot*sizeof(double));
    }
    memcpy(rhoDnOld,&rhoDn[1],rhoRealGridNum*sizeof(double));
    calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
			 rhoCoeffReDn,rhoCoeffImDn,rhoDn,
			 rhoCoeffReDnDensCpBox,rhoCoeffImDnDensCpBox,
			 divRhoxDn,divRhoyDn,divRhozDn,d2RhoDn,cpGGA,
			 cpDualGridOptOn,numInterpPmeDual,
			 communicate,&(cp->cp_para_fft_pkg3d_lg),
			 &(cp->cp_para_fft_pkg3d_dens_cp_box));
  }//endif cpLsda

/*======================================================================*/
/* II) Generating Pseudopotential List                                  */
/*     (This should be done in calcRhoInit. However, we skip this	*/
/*      function when we read checkpoint file. So we need do it here.)  */

  
  if(myidState==0){
    PRINT_LINE_STAR;
    printf("Start Generating Pseudopotential List\n");
    PRINT_LINE_DASH;
  }

  if(stodftInfo->vpsAtomListFlag==0||cpDualGridOptOn>= 1){
    control_vps_atm_list(pseudo,cell,clatoms_pos,clatoms_info,
                         atommaps,ewd_scr,for_scr,cpDualGridOptOn,
                         stodftInfo->vpsAtomListFlag);
    stodftInfo->vpsAtomListFlag = 1;
  }

  if(myidState==0){
    PRINT_LINE_DASH;
    printf("Finish Generating Pseudopotential List\n");
    PRINT_LINE_STAR;
  }

/*======================================================================*/
/* II) Generating Occupatation number                                   */
/*     (This should be done in calcRhoInit. However, we skip this       */
/*      function when we read checkpoint file. So we need do it here.)  */

  if(cpLsda==1)stodftInfo->occNumber = 1;
  else stodftInfo->occNumber = 2;

/*======================================================================*/
/* III) Read Diis                                                       */

  for(iDiis=0;iDiis<numDiis;iDiis++){
    rhoUpBank[iDiis] = (double*)cmalloc(rhoRealGridNum*sizeof(double));
    rhoUpErr[iDiis] = (double*)cmalloc(rhoRealGridNum*sizeof(double));
  }
  if(numDiisNow<numDiis)numDiisOutput = numDiisNow+1;
  else numDiisOutput = numDiisNow;

  for(iDiis=0;iDiis<numDiisOutput;iDiis++){
    if(myidState==0){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	fscanf(fileCheckpoint,"%lg",&rhoTemp[iGrid]);
      }//endfor iGrid
    }//endif myidState
    if(numProcStates>1){
      Barrier(commStates);
      Scatterv(rhoTemp,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	       rhoUpBank[iDiis],rhoRealGridNum,MPI_DOUBLE,0,commStates);
      Barrier(commStates);
    }
    else{
      memcpy(rhoUpBank[iDiis],rhoTemp,rhoRealGridTot*sizeof(double));
    }
  }
  for(iDiis=0;iDiis<numDiisNow;iDiis++){
    if(myidState==0){    
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
        fscanf(fileCheckpoint,"%lg",&rhoTemp[iGrid]);
      }//endfor iGrid
    }
    if(numProcStates>1){
      Barrier(commStates);
      Scatterv(rhoTemp,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	       rhoUpErr[iDiis],rhoRealGridNum,MPI_DOUBLE,0,commStates);
      Barrier(commStates);
    }
    else{
      memcpy(rhoUpErr[iDiis],rhoTemp,rhoRealGridTot*sizeof(double));
    }
  }
  if(cpLsda==1){
    for(iDiis=0;iDiis<numDiisOutput;iDiis++){
      if(myidState==0){
	for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	  fscanf(fileCheckpoint,"%lg",&rhoTemp[iGrid]);
	}//endfor iGrid
      }//endif myidState
      if(numProcStates>1){
        Barrier(commStates);
	Scatterv(rhoTemp,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
		 rhoDnErr[iDiis],rhoRealGridNum,MPI_DOUBLE,0,commStates);
	Barrier(commStates);
      }
      else{
        memcpy(rhoDnBank[iDiis],rhoTemp,rhoRealGridTot*sizeof(double));
      }
    }
    for(iDiis=0;iDiis<numDiisNow;iDiis++){
      if(myidState==0){ 
	for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	  fscanf(fileCheckpoint,"%lg",&rhoTemp[iGrid]);
	}//endfor iGrid
      }
      if(numProcStates>1){
	Barrier(commStates);
	Scatterv(rhoTemp,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
		 rhoDnErr[iDiis],rhoRealGridNum,MPI_DOUBLE,0,commStates);
	Barrier(commStates);
      }
      else{
	memcpy(rhoUpErr[iDiis],rhoTemp,rhoRealGridTot*sizeof(double));
      }
    }
  }
  
  if(myidState==0){
    fclose(fileCheckpoint);
    printf("Finish reading checkpoint file...\n");
  }
  
  free(rhoTemp);

/*======================================================================*/
/* IV) Generate a random guess for calculate H_KS spetral range         */

  if(myidState==0){
    for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
      testWfMaxRe[iCoeff] = 1.0;
      testWfMaxIm[iCoeff] = -1.0;
    }
    testWfMaxRe[numCoeff-1] = 1.0;
    testWfMaxIm[numCoeff-1] = 0.0;
    for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
      testWfMinRe[iCoeff] = 1.0;
      testWfMinIm[iCoeff] = -1.0;
    }
    testWfMinRe[numCoeff-1] = 1.0;
    testWfMinIm[numCoeff-1] = 0.0;
  }


/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

