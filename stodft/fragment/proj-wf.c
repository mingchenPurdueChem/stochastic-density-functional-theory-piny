/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: proj-wf.c                                    */
/*                                                                          */
/* This routine project the stochastic wave function on to fragments        */
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
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_frag_local.h"

//#define DMAT
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void projRhoMiniMol(CP *cp,GENERAL_DATA *general_data,CLASS *class,
		 CP *cpMini,GENERAL_DATA *generalDataMini,CLASS *classMini,
		 int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  STODFTINFO *stodftInfo = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  CPOPTS *cpOpts = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cp->cpcoeffs_info);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  COMMUNICATE *commCP = &(cp->communicate);
  CELL *cell	      = &(general_data->cell);
 
  int iFrag,jFrag,iGrid,iState,jState,iStateFrag,iProc;
  int iAtom;
  int gridIndex;
  int numGrid;
  int numStateUpMini;
  int numStateDn;
  int numStateDnMini;
  int cpLsda		    = cpOpts->cp_lsda;
  int numFragProc	    = fragInfo->numFragProc;
  int numFragTot	    = fragInfo->numFragTot;
  int occNumber		    = stodftInfo->occNumber;
  int myidState		    = commCP->myid_state;
  int numProcStates	    = commCP->np_states;
  int rhoRealGridNum	    = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot	    = stodftInfo->rhoRealGridTot;
  int numStateUpProc	    = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc	    = cpcoeffs_info->nstate_dn_proc;
  int numStateStoUp     = stodftInfo->numStateStoUp;
  int numStateStoDn     = stodftInfo->numStateStoDn;
  int numCoeff		= cpcoeffs_info->ncoef;
  int countWf;
  int fragOpt		    = stodftInfo->fragOpt;

  MPI_Comm commStates   =    commCP->comm_states;


  int *numGridFragProc	    = fragInfo->numGridFragProc;
  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;
  int *numStateUpAllProc;
  int *numStateDnAllProc;

  int **gridMapProc	    = fragInfo->gridMapProc;
 
  double proj;
  double numElecProj = 0.0;
  double numGridTotInv = 1.0/rhoRealGridTot;
  double vol,volInv;
  double volMini;
  double pre,preNe;
  double preDot;
  double *rhoTemp,*rhoTempReduce;
  double *wfTemp,*wfFragTemp,*rhoFragTemp;
  double *hmat_cp = cell->hmat_cp;
  double *hmatCpMini;
  double **rhoUpFragProc;
  double **coefUpFragProc;
  double **rhoDnFragProc;
  double **coefDnFragProc;
  double *rhoUpFragSum;
  double *rhoDnFragSum;
  double *noiseWfUpReal,*noiseWfDnReal;

  //debug
  //double *fragWfCpy = (double*)cmalloc(rhoRealGridNum*sizeof(double));
  
/*======================================================================*/
/* I) Allocated memories				                */

  vol = getdeth(hmat_cp);
  volInv = 1.0/vol;
  
  // prefactor for prjection part
  pre = (double)(occNumber)*vol*vol/
	((double)(rhoRealGridTot)*(double)(rhoRealGridTot)*(double)(numStateStoUp));
  // prefactor for number of e in proj part
  preNe = pre/(double)(rhoRealGridTot);
  preDot = 1.0/sqrt(vol);
  /*
  printf("vol %lg\n",vol);
  printf("rhoRealGridTot %i pre %lg preNe %lg\n",rhoRealGridTot,pre,preNe);
  */
  
  fragInfo->rhoUpFragProc = (double**)cmalloc(numFragProc*sizeof(double*));
  //fragInfo->rhoDnFragProc = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->coefUpFragProc = (double**)cmalloc(numFragProc*sizeof(double*));
  //fragInfo->coefDnFragProc = (double**)cmalloc(numFragProc*sizeof(double*));
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    numGrid = numGridFragProc[iFrag];
    numStateUpMini = cpMini[iFrag].cpcoeffs_info.nstate_up_proc;
    fragInfo->rhoUpFragProc[iFrag] = (double*)cmalloc(numGrid*sizeof(double));
    fragInfo->coefUpFragProc[iFrag] = (double*)cmalloc(numStateUpMini*numGrid*sizeof(double));
  }
  rhoUpFragProc = fragInfo->rhoUpFragProc;
  coefUpFragProc = fragInfo->coefUpFragProc;
  if(cpLsda==1&&numStateDn!=0){
    fragInfo->rhoDnFragProc = (double**)cmalloc(numFragProc*sizeof(double*));
    fragInfo->coefDnFragProc = (double**)cmalloc(numFragProc*sizeof(double*));
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      numGrid = numGridFragProc[iFrag];
      numStateDnMini = cpMini[iFrag].cpcoeffs_info.nstate_dn_proc;
      fragInfo->rhoUpFragProc[iFrag] = (double*)cmalloc(numGrid*sizeof(double));
      fragInfo->coefUpFragProc[iFrag] = (double*)cmalloc(numStateDnMini*numGrid*sizeof(double));
    }
    rhoDnFragProc = fragInfo->rhoDnFragProc;
    coefDnFragProc = fragInfo->coefDnFragProc;
  }

  fragInfo->rhoUpFragSum = (double*)cmalloc(rhoRealGridNum*sizeof(double));
  fragInfo->rhoUpFragSumCpy = (double*)cmalloc(rhoRealGridNum*sizeof(double));
  rhoUpFragSum = fragInfo->rhoUpFragSum;
  if(cpLsda==1&&numStateDn!=0){
    fragInfo->rhoDnFragSum = (double*)cmalloc(rhoRealGridNum*sizeof(double));
    fragInfo->rhoDnFragSumCpy = (double*)cmalloc(rhoRealGridNum*sizeof(double));
    rhoDnFragSum = fragInfo->rhoDnFragSum;
  }

/*======================================================================*/
/* II) Transform fragment wf to real spaces and store them              */

  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInfo->iFrag = iFrag;
    if(fragOpt==1){
      rhoRealCalcDriverFragMol(&generalDataMini[iFrag],&cpMini[iFrag],&classMini[iFrag],cp);
    }
    if(fragOpt==3){
      rhoRealCalcDriverFragUnitCell(&generalDataMini[iFrag],&cpMini[iFrag],&classMini[iFrag],cp);
    }
  }
  //debug
  /*
  char fileNameFragMOReal[100];
  char fileNameFragRhoReal[100];
  FILE *fileFragMOReal,*fileFragRhoReal;
  double preTest = 1.0/sqrt(2.0*vol);
  int junk;
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    sprintf(fileNameFragMOReal,"wf-frag-real-%i",iFrag);
    sprintf(fileNameFragRhoReal,"rho-frag-%i",iFrag);
    fileFragMOReal = fopen(fileNameFragMOReal,"r");
    fileFragRhoReal = fopen(fileNameFragRhoReal,"r");
    numStateUpMini = cpMini[iFrag].cpcoeffs_info.nstate_up_proc;
    numGrid = numGridFragProc[iFrag];
    for(iState=0;iState<numStateUpMini;iState++){
      for(iGrid=0;iGrid<numGrid;iGrid++){
	fscanf(fileFragMOReal,"%lg",&(fragInfo->coefUpFragProc[iFrag][iState*numGrid+iGrid]));
	fragInfo->coefUpFragProc[iFrag][iState*numGrid+iGrid] *= preTest;
      }
    }
    for(iGrid=0;iGrid<numGrid;iGrid++){
      fscanf(fileFragRhoReal,"%i",&junk);
      fscanf(fileFragRhoReal,"%i",&junk);
      fscanf(fileFragRhoReal,"%i",&junk);
      fscanf(fileFragRhoReal,"%lg",&(fragInfo->rhoUpFragProc[iFrag][iGrid]));
      fragInfo->rhoUpFragProc[iFrag][iGrid] *= volInv;
    }
  }
  */

/*======================================================================*/
/* III) Assemble fragments densities		                        */

  // Allocate rhoTemp for each proc
  rhoTemp = (double*)cmalloc(rhoRealGridTot*sizeof(double));
  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
    rhoTemp[iGrid] = 0.0;
  }//endfor iGrid
  // Reduce all fragment densities in each proc
  /*
  for(iGrid=0;iGrid<numGrid;iGrid++){
    printf("1111111 rhofrag1 %lg\n",rhoUpFragProc[0][iGrid]);
    //rhoTemp[gridMapProc[iFrag][iGrid]] += rhoUpFragProc[0][iGrid];
  }//endfor iGrid
  fflush(stdout);
  exit(0);
  */

  for(iFrag=0;iFrag<numFragProc;iFrag++){
    numGrid = numGridFragProc[iFrag];
    for(iGrid=0;iGrid<numGrid;iGrid++){
      //printf("1111111 rhofrag1 %lg\n",rhoUpFragProc[iFrag][iGrid]);
      rhoTemp[gridMapProc[iFrag][iGrid]] += rhoUpFragProc[iFrag][iGrid];
    }//endfor iGrid
  }//endfor iFrag
  
  // Allocate rhoTempReduce on the master proc
  if(myidState==0&&numProcStates>1){
    rhoTempReduce = (double*)cmalloc(rhoRealGridTot*sizeof(double));
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoTempReduce[iGrid] = 0.0;
  }
  // Reduce fragment densities on different procs and scatter them
  if(numProcStates>1){
    Barrier(commStates);
    Reduce(rhoTemp,rhoTempReduce,rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Scatterv(rhoTempReduce,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
             rhoUpFragSum,rhoRealGridNum,MPI_DOUBLE,0,commStates);
    // I need to multiply the density by volumn to match the |ksi(r)|^2
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoUpFragSum[iGrid] *= vol;
    //if(myidState==0)free(&rhoTempReduce[0]);
  }
  // Sequential, just copy and scale the density
  else{
    //memcpy(rhoUpFragSum,rhoTemp,rhoRealGridNum*sizeof(double));
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
      rhoUpFragSum[iGrid] = rhoTemp[iGrid]*vol;
      //printf("111111 rhofrag %.8lg\n",rhoTemp[iGrid]*vol);
    }
  }
  //fflush(stdout);
  //exit(0);
  memcpy(&(fragInfo->rhoUpFragSumCpy[0]),rhoUpFragSum,rhoRealGridNum*sizeof(double));
  // Do the same thing for spin down states
  if(cpLsda==1&&numStateDn!=0){
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      rhoTemp[iGrid] = 0.0;
    }
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      numGrid = numGridFragProc[iFrag];
      for(iGrid=0;iGrid<numGrid;iGrid++){
	rhoTemp[gridMapProc[iFrag][iGrid]] += rhoDnFragProc[iFrag][iGrid];
      }
    }

    if(myidState==0&&numProcStates>1){
      rhoTempReduce = (double*)cmalloc(rhoRealGridTot*sizeof(double));
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoTempReduce[iGrid] = 0.0;
    }
    if(numProcStates>1){
      Barrier(commStates);
      Reduce(rhoTemp,rhoTempReduce,rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
      Scatterv(rhoTempReduce,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	       rhoDnFragSum,rhoRealGridNum,MPI_DOUBLE,0,commStates);
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoUpFragSum[iGrid] *= vol;
      //if(myidState==0)free(&rhoTempReduce[0]);
    }
    else{
      //memcpy(rhoDnFragSum,rhoTemp,rhoRealGridNum*sizeof(double));
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
	rhoUpFragSum[iGrid] = rhoTemp[iGrid]*vol;
      }//endfor iGrid
    }//endif 
    memcpy(&(fragInfo->rhoDnFragSumCpy[0]),rhoDnFragSum,rhoRealGridNum*sizeof(double));
  }//endif

  //debug
  //memcpy(fragWfCpy,rhoUpFragSum,rhoRealGridNum*sizeof(double));

/*======================================================================*/
/* IV) Free/allocate some memory for next step                          */

  for(iFrag=0;iFrag<numFragProc;iFrag++){
    free(fragInfo->rhoUpFragProc[iFrag]);
  }
  free(fragInfo->rhoUpFragProc);
  if(cpLsda==1&&numStateDn!=0){
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      free(fragInfo->rhoDnFragProc[iFrag]);
    }
    free(fragInfo->rhoDnFragProc);   
  }
  //free(rhoTemp);
  fragInfo->noiseWfUpReal = (double*)cmalloc(numStateUpProc*rhoRealGridTot*sizeof(double));
  noiseWfUpReal = fragInfo->noiseWfUpReal;
  if(cpLsda==1&&numStateDn!=0){
    fragInfo->noiseWfDnReal = (double*)cmalloc(numStateDnProc*rhoRealGridTot*sizeof(double));
    noiseWfDnReal = fragInfo->noiseWfDnReal;
  }
  
  fragInfo->wfProjUp = (double**)cmalloc(numFragProc*sizeof(double*));
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    numStateUpMini = cpMini[iFrag].cpcoeffs_info.nstate_up_proc;
    fragInfo->wfProjUp[iFrag] = (double*)cmalloc(numStateStoUp*numStateUpMini*sizeof(double));
  }
  if(cpLsda==1&&numStateDn!=0){
    fragInfo->wfProjDn = (double**)cmalloc(numFragProc*sizeof(double*));
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      numStateDnMini = cpMini[iFrag].cpcoeffs_info.nstate_dn_proc;
      fragInfo->wfProjDn[iFrag] = (double*)cmalloc(numStateStoDn*numStateDnMini*sizeof(double));
    }
  }

  
/*======================================================================*/
/* IV) Calculate the real space noise wave function                     */

  //debug copy the deterministic wf here and backup the stochastic ones
  /*
  int numCoeffTotal = numStateStoUp*numCoeff;
  int iCoeff;
  double *noiseReUpBackup = (double*)cmalloc(numCoeffTotal*sizeof(double));
  double *noiseImUpBackup = (double*)cmalloc(numCoeffTotal*sizeof(double));
  
  for(iCoeff=0;iCoeff<numCoeffTotal;iCoeff++){
    noiseReUpBackup[iCoeff] = cp->cpcoeffs_pos[1].cre_up[iCoeff+1];
    noiseImUpBackup[iCoeff] = cp->cpcoeffs_pos[1].cim_up[iCoeff+1];
  }
  for(iState=0;iState<4;iState++){
    FILE *fwfread = fopen("wf-det","r");
    for(jState=0;jState<4;jState++){
      for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
        fscanf(fwfread,"%lg",&(cp->cpcoeffs_pos[1].cre_up[iState*4*numCoeff+jState*numCoeff+iCoeff]));
        fscanf(fwfread,"%lg",&(cp->cpcoeffs_pos[1].cim_up[iState*4*numCoeff+jState*numCoeff+iCoeff]));
      }      
    }
    fclose(fwfread);
  }
  double presq = sqrt(2.0);
  for(iCoeff=1;iCoeff<=numCoeffTotal;iCoeff++){
    cp->cpcoeffs_pos[1].cre_up[iCoeff] *= presq;
    cp->cpcoeffs_pos[1].cim_up[iCoeff] *= presq;
  }
  */

  
  /*
  char wfname[100];
  //sprintf(wfname,"/scratch/mingchen/tmp/sto-wf-save-%i",myidState);
  printf("Read in stochastic orbitals...\n");
  sprintf(wfname,"sto-wf-save-%i",myidState);

  FILE *filePrintWF = fopen(wfname,"r");
  int iCoeff;
  for(iState=0;iState<numStateUpProc;iState++){
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      fscanf(filePrintWF,"%lg",&(cp->cpcoeffs_pos[1].cre_up[iState*numCoeff+iCoeff]));
      fscanf(filePrintWF,"%lg",&(cp->cpcoeffs_pos[1].cim_up[iState*numCoeff+iCoeff]));
    }//endfor iCoeff
  }//endfor iState
  fclose(filePrintWF);
  printf("myid %i finish reading in WF.\n",myidState);
  */
  // Copy the deterministic wf to the stochastic wf
  //memcpy(&(cp->cpcoeffs_pos[1].cre_up[1]),&(stodftCoefPos->wfDetBackupUpRe[0]),numStateUpProc*numCoeff*sizeof(double));
  //memcpy(&(cp->cpcoeffs_pos[1].cim_up[1]),&(stodftCoefPos->wfDetBackupUpIm[0]),numStateUpProc*numCoeff*sizeof(double));
  /*
  int iCoeff;
  for(iCoeff=0;iCoeff<numStateUpProc*numCoeff;iCoeff++){
    cp->cpcoeffs_pos[1].cre_up[iCoeff+1] = stodftCoefPos->wfDetBackupUpRe[iCoeff];
    cp->cpcoeffs_pos[1].cim_up[iCoeff+1] = stodftCoefPos->wfDetBackupUpIm[iCoeff];
  }

  rhoRealCalcDriverNoise(general_data,cp,class,ip_now);
  */
  
  //fflush(stdout);
  //exit(0);
  

  noiseRealReGen(general_data,cp,class,ip_now);

/*======================================================================*/
/* IV) Project the real space noise wave function                       */

  numStateUpAllProc = (int*)cmalloc(numProcStates*sizeof(int));
  if(numProcStates>1){
    Allgather(&numStateUpProc,1,MPI_INT,numStateUpAllProc,1,MPI_INT,0,commStates);
  }
  else numStateUpAllProc[0] = numStateUpProc;
  if(cpLsda==1&&numStateDn!=0){
    numStateDnAllProc = (int*)cmalloc(numProcStates*sizeof(int));
    if(numProcStates>1){
      Allgather(&numStateDnProc,1,MPI_INT,numStateDnAllProc,1,MPI_INT,0,commStates);
    }
    else numStateDnAllProc[0] = numStateDnProc;
  }
  wfTemp = (double*)cmalloc(rhoRealGridTot*sizeof(double));
  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
    rhoTemp[iGrid] = 0.0;
  }
  countWf = 0;
  //debug
  /*
  double **fix_frag = (double**)cmalloc(numFragProc*sizeof(double*));
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fix_frag[iFrag] = (double*)calloc(rhoRealGridTot,sizeof(double));
  }
  */
  //debug
  /*
  double *projMatrixTest = (double*)calloc(numStateUpProc*numStateUpProc,sizeof(double));
  double *wfProjjTemp = (double*)calloc(rhoRealGridTot,sizeof(double));
  for(iState=0;iState<numStateUpProc;iState++){
    memcpy(wfTemp,&noiseWfUpReal[iState*rhoRealGridTot],rhoRealGridTot*sizeof(double));
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)wfProjjTemp[iGrid] = 0.0;
    if(numProcStates>1)Bcast(wfTemp,rhoRealGridTot,MPI_DOUBLE,iProc,commStates);
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      double *coefUpFragExt = (double*)calloc(numStateUpMini*rhoRealGridTot,sizeof(double));
      numGrid = numGridFragProc[iFrag];
      numStateUpMini = cpMini[iFrag].cpcoeffs_info.nstate_up_proc;
      hmatCpMini = generalDataMini[iFrag].cell.hmat_cp;
      volMini = getdeth(hmatCpMini);
      volMini /= numGrid;
      wfFragTemp = (double*)cmalloc(numGrid*sizeof(double));
      rhoFragTemp = (double*)cmalloc(numGrid*sizeof(double));
      for(iGrid=0;iGrid<numGrid;iGrid++){
	gridIndex = gridMapProc[iFrag][iGrid];
	wfFragTemp[iGrid] = wfTemp[gridIndex];
	rhoFragTemp[iGrid] = 0.0;
      }
      for(iStateFrag=0;iStateFrag<numStateUpMini;iStateFrag++){
	//double *tempExt = (double*)calloc(rhoRealGridTot,sizeof(double));
	double preCopyFrag = 1.0/sqrt(numFragProc);
	for(jFrag=0;jFrag<numFragProc;jFrag++){   
	  for(iGrid=0;iGrid<numGrid;iGrid++){
	    gridIndex = gridMapProc[jFrag][iGrid];
	    coefUpFragExt[iStateFrag*rhoRealGridTot+gridIndex] = coefUpFragProc[iFrag][iStateFrag*numGrid+iGrid]*preCopyFrag;
	  }
	}
        proj = ddotBlasWrapper(rhoRealGridTot,&wfTemp[0],1,&coefUpFragExt[iStateFrag*rhoRealGridTot],1);
        //free(tempExt); 
	//proj = ddotBlasWrapper(numGrid,&wfFragTemp[0],1,&coefUpFragProc[iFrag][iStateFrag*numGrid],1);
	fragInfo->wfProjUp[iFrag][countWf*numStateUpMini+iStateFrag] = proj*preDot*volMini;
	//daxpyBlasWrapper(numGrid,proj,&coefUpFragProc[iFrag][iStateFrag*numGrid],1,&rhoFragTemp[0],1);
	//daxpyBlasWrapper(rhoRealGridTot,proj,&coefUpFragProc[iFrag][iStateFrag*rhoRealGridTot],1,&wfProjjTemp[0],1);
	for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	  wfProjjTemp[iGrid] += proj*coefUpFragExt[iStateFrag*rhoRealGridTot+iGrid];
	}
      }//endfor iStateFrag
      for(jState=iState;jState<numStateUpProc;jState++){
	//for(iGrid=0;iGrid<numGrid;iGrid++){
	  //gridIndex = gridMapProc[iFrag][iGrid];
	  ////rhoTemp[gridIndex] += rhoFragTemp[iGrid]*rhoFragTemp[iGrid];
	  //projMatrixTest[iState*numStateUpProc+jState] += noiseWfUpReal[jState*rhoRealGridTot+gridIndex]*rhoFragTemp[iGrid]*volMini*volMini;
	  ////fix_frag[iFrag][gridIndex] += rhoFragTemp[iGrid]*rhoFragTemp[iGrid]*pre;
	//}
	for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	  projMatrixTest[iState*numStateUpProc+jState] += noiseWfUpReal[jState*rhoRealGridTot+iGrid]*wfProjjTemp[iGrid]*volMini*volMini;
	}
      }//endfor jState
      free(wfFragTemp);
      free(rhoFragTemp);
      free(coefUpFragExt);
    }//endfor iFrag
    countWf += 1;
    for(jState=iState;jState<numStateUpProc;jState++){
      projMatrixTest[iState*numStateUpProc+jState] *= 0.5*volInv;
      projMatrixTest[jState*numStateUpProc+iState] = projMatrixTest[iState*numStateUpProc+jState];
    }
  }//endfor iState
  FILE *testMatrix = fopen("test-matrix","w");
  for(iState=0;iState<numStateUpProc;iState++){
    for(jState=0;jState<numStateUpProc;jState++){
      if(iState==jState)projMatrixTest[iState*numStateUpProc+jState] -= 1.0;
      fprintf(testMatrix,"%.8lg ",projMatrixTest[iState*numStateUpProc+jState]);
    }
    fprintf(testMatrix,"\n");
    printf("111111 matrix diag %.8lg\n",projMatrixTest[iState*numStateUpProc+iState]+1.0);
  }
  fflush(stdout);
  exit(0);   
  */
  //end debug  

  for(iProc=0;iProc<numProcStates;iProc++){
    for(iState=0;iState<numStateUpAllProc[iProc];iState++){
      if(myidState==iProc){
	memcpy(wfTemp,&noiseWfUpReal[iState*rhoRealGridTot],rhoRealGridTot*sizeof(double));
      }
      if(numProcStates>1)Bcast(wfTemp,rhoRealGridTot,MPI_DOUBLE,iProc,commStates);
      //double *rho_frag_ext = (double*)calloc(rhoRealGridTot,sizeof(double));
      for(iFrag=0;iFrag<numFragProc;iFrag++){
	numGrid = numGridFragProc[iFrag];
	numStateUpMini = cpMini[iFrag].cpcoeffs_info.nstate_up_proc;
	hmatCpMini = generalDataMini[iFrag].cell.hmat_cp;
	volMini = getdeth(hmatCpMini);
	volMini /= numGrid;
	wfFragTemp = (double*)cmalloc(numGrid*sizeof(double));
	rhoFragTemp = (double*)cmalloc(numGrid*sizeof(double));
	for(iGrid=0;iGrid<numGrid;iGrid++){
	  gridIndex = gridMapProc[iFrag][iGrid];
	  wfFragTemp[iGrid] = wfTemp[gridIndex];
	  rhoFragTemp[iGrid] = 0.0;
	}
	for(iStateFrag=0;iStateFrag<numStateUpMini;iStateFrag++){
	  // debug Let's try dot product in k space. But I can do this in real space just for test
	  /*
	  double *tempExt = (double*)calloc(rhoRealGridTot,sizeof(double));
	  double preCopyFrag = 1.0/sqrt(numFragProc);
	  for(jFrag=0;jFrag<numFragProc;jFrag++){   
	    for(iGrid=0;iGrid<numGrid;iGrid++){
	      gridIndex = gridMapProc[jFrag][iGrid];
	      tempExt[gridIndex] = coefUpFragProc[iFrag][iStateFrag*numGrid+iGrid]*preCopyFrag;
	    }
	  }
	  */
	  /*
	  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	    rho_frag_ext[iGrid] += tempExt[iGrid]*tempExt[iGrid];
	  }
	  */
	  proj = ddotBlasWrapper(numGrid,&wfFragTemp[0],1,&coefUpFragProc[iFrag][iStateFrag*numGrid],1);
	  //proj = ddotBlasWrapper(rhoRealGridTot,&wfTemp[0],1,&tempExt[0],1);
	  //free(tempExt);
	  fragInfo->wfProjUp[iFrag][countWf*numStateUpMini+iStateFrag] = proj*preDot*volMini;
	  daxpyBlasWrapper(numGrid,proj,&coefUpFragProc[iFrag][iStateFrag*numGrid],1,&rhoFragTemp[0],1);
	}//endfor iStateFrag
	/*
	char fname[100];
	sprintf(fname,"wf-%i",iState);
	FILE *testwf = fopen(fname,"w");
	double *temp = (double*)calloc(numGrid,sizeof(double));
	for(iGrid=0;iGrid<numGrid;iGrid++)temp[gridMapProc[iFrag][iGrid]] = rhoFragTemp[iGrid];
	for(iGrid=0;iGrid<numGrid;iGrid++)fprintf(testwf,"%i %.16lg\n",iGrid,-temp[iGrid]*volMini);
	fclose(testwf);
	free(temp);
	*/
	for(iGrid=0;iGrid<numGrid;iGrid++){
	  gridIndex = gridMapProc[iFrag][iGrid];
	  rhoTemp[gridIndex] += rhoFragTemp[iGrid]*rhoFragTemp[iGrid];
	  //fix_frag[iFrag][gridIndex] += rhoFragTemp[iGrid]*rhoFragTemp[iGrid]*pre;
	}
	free(wfFragTemp);
	free(rhoFragTemp);
      }//endfor iFrag
      /*
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	printf("2222222 rhoFragExt %.8lg\n",rho_frag_ext[iGrid]*vol*2.0);
      }
      free(rho_frag_ext);
      exit(0);
      */
      countWf += 1;
    }//endfor iState   
  }//endfor iProc
  //fflush(stdout);
  //exit(0);
  //debug
  /*
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    char fname[100];
    sprintf(fname,"proj-%i",iFrag);
    FILE *test_proj = fopen(fname,"w");
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      fprintf(test_proj,"%.8lg\n",fix_frag[iFrag][iGrid]);
    }
    fclose(test_proj);
  }
  */
  /*
  numStateUpMini = cpMini[0].cpcoeffs_info.nstate_up_proc;
  int iCoeff;
  int numCoeffMini = cpMini->cpcoeffs_info.ncoef;
  double dott;
  double presqrt2 = 1.0/sqrt(2.0);
  //double presqrt2 = 1.0;
  double *cree = cpMini->cpcoeffs_pos[1].cre_up;
  double *cimm = cpMini->cpcoeffs_pos[1].cim_up;
  double *creComb = (double *)cmalloc(16*numCoeffMini*sizeof(double));
  double *cimComb = (double *)cmalloc(16*numCoeffMini*sizeof(double));
  double *creback = (double *)cmalloc(4*numCoeffMini*sizeof(double));
  double *cimback = (double *)cmalloc(4*numCoeffMini*sizeof(double));

  memcpy(&creback[0],&(cpMini->cpcoeffs_pos[1].cre_up[1]),4*numCoeffMini*sizeof(double));
  memcpy(&cimback[0],&(cpMini->cpcoeffs_pos[1].cim_up[1]),4*numCoeffMini*sizeof(double));

  for(iGrid=0;iGrid<16*numCoeffMini;iGrid++){
    creComb[iGrid] = 0.0;
    cimComb[iGrid] = 0.0;
    
  }
  for(iState=0;iState<16;iState++){
    for(jState=0;jState<numStateUpMini;jState++){
      for(iCoeff=0;iCoeff<numCoeffMini;iCoeff++){
	creComb[iState*numCoeffMini+iCoeff] += 
	    fragInfo->wfProjUp[0][iState*numStateUpMini+jState]*cree[jState*numCoeffMini+iCoeff+1]*presqrt2;
        cimComb[iState*numCoeffMini+iCoeff] +=             
            fragInfo->wfProjUp[0][iState*numStateUpMini+jState]*cimm[jState*numCoeffMini+iCoeff+1]*presqrt2;
      }
    }
    dott = 0.0;
    for(iCoeff=0;iCoeff<numCoeffMini-1;iCoeff++){
      dott += creComb[iState*numCoeffMini+iCoeff]*cree[iCoeff+1]*presqrt2
	      +cimComb[iState*numCoeffMini+iCoeff]*cimm[iCoeff+1]*presqrt2;
    }
    dott *= 2.0;
    dott += creComb[iState*numCoeffMini+numCoeffMini-1]*cree[numCoeffMini]*presqrt2;
    printf("iState %i dott %lg\n",iState,dott);
  }
  for(iAtom=1;iAtom<=3;iAtom++){
    classMini->clatoms_pos[1].fx[iAtom] = 0.0;
    classMini->clatoms_pos[1].fy[iAtom] = 0.0;
    classMini->clatoms_pos[1].fz[iAtom] = 0.0;
  }
  cpMini->cpcoeffs_info.nstate_up = 16;
  cpMini->cpcoeffs_info.nstate_dn = 16;
  memcpy(&(cree[1]),&creComb[0],16*numCoeffMini*sizeof(double));
  memcpy(&(cimm[1]),&cimComb[0],16*numCoeffMini*sizeof(double));
  //for(iGrid=1;iGrid<=16*numCoeffMini;iGrid++){
  //  cpMini->cpcoeffs_pos[1].cre_up[iGrid] = 0.0;
  //  cpMini->cpcoeffs_pos[1].cim_up[iGrid] = 0.0;
  //}
  //memcpy(&(cpMini->cpcoeffs_pos[1].cre_up[1]),&creback[0],4*numCoeffMini*sizeof(double));
  //memcpy(&(cpMini->cpcoeffs_pos[1].cim_up[1]),&cimback[0],4*numCoeffMini*sizeof(double));

  
  cp_ks_energy_ctrl(cpMini,1,&(generalDataMini->ewald),&(classMini->ewd_scr),
                       &(generalDataMini->cell),&(classMini->clatoms_info),
                       &(classMini->clatoms_pos[1]),&(classMini->atommaps),
                       &(generalDataMini->stat_avg),&(generalDataMini->ptens),&(generalDataMini->simopts),
                       &(classMini->for_scr));

  cpMini->cpcoeffs_info.nstate_up = 4;
  cpMini->cpcoeffs_info.nstate_dn = 4;
  
  memcpy(&(cree[1]),&creback[0],4*numCoeffMini*sizeof(double));
  memcpy(&(cimm[1]),&cimback[0],4*numCoeffMini*sizeof(double));
  */
    
  //enddebug
  
  if(numProcStates>1){
    Barrier(commStates);
    Reduce(rhoTemp,rhoTempReduce,rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Scatterv(rhoTempReduce,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	     rhoTemp,rhoRealGridNum,MPI_DOUBLE,0,commStates);
  } 
  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)numElecProj += rhoTemp[iGrid];
 
  //printf("numElecProj %lg\n",numElecProj);
  //debug
  
  double sumElecFrag = 0.0;
  double sumElecProj = 0.0;
  for(iProc=0;iProc<numProcStates;iProc++){
    if(myidState==iProc){
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
	sumElecFrag += rhoUpFragSum[iGrid];
	sumElecProj += pre*rhoTemp[iGrid];
        //printf("rhofraggggg %lg rhoproj %lg\n",rhoUpFragSum[iGrid],pre*rhoTemp[iGrid]);
      }
      sumElecFrag /= rhoRealGridTot;
      sumElecProj /= rhoRealGridTot;
    }
    if(numProcStates>1)Barrier(commStates);
  }
  //printf("sumElecFrag %lg sumElecProj %lg rhoRealGridNum %i\n",sumElecFrag,sumElecProj,rhoRealGridNum);
  if(numProcStates>1)Barrier(commStates);
  //exit(0);
    
  //debug
  /*
  FILE *file_test = fopen("rho-test","w");
  double sum1 = 0.0;
  double sum2 = 0.0;
  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
    fprintf(file_test,"%.8lg %.8lg\n",rhoUpFragSum[iGrid],rhoTemp[iGrid]*pre);
    sum1 += rhoUpFragSum[iGrid];
    sum2 += rhoTemp[iGrid]*pre;
  }
  printf("sum1 %.8lg sum2 %.8lg\n",sum1,sum2);
  fclose(file_test);
  */
  
  daxpyBlasWrapper(rhoRealGridNum,-pre,&rhoTemp[0],1,&rhoUpFragSum[0],1);

  if(cpLsda==1&&numStateDn!=0){
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      rhoTemp[iGrid] = 0.0;
    }
    countWf = 0;
    for(iProc=0;iProc<numProcStates;iProc++){
      for(iState=0;iState<numStateDnAllProc[iProc];iState++){
	if(myidState==iProc){
	  memcpy(wfTemp,&noiseWfDnReal[iState*rhoRealGridTot],rhoRealGridTot*sizeof(double));
	}
	if(numProcStates>1)Bcast(wfTemp,rhoRealGridTot,MPI_DOUBLE,iProc,commStates);
	for(iFrag=0;iFrag<numFragProc;iFrag++){
	  numGrid = numGridFragProc[iFrag];
	  numStateDnMini = cpMini[iFrag].cpcoeffs_info.nstate_dn_proc;
	  hmatCpMini = generalDataMini[iFrag].cell.hmat_cp;
	  volMini = getdeth(hmatCpMini);
	  volMini /= numGrid;
	  wfFragTemp = (double*)cmalloc(numGrid*sizeof(double));
	  rhoFragTemp = (double*)cmalloc(numGrid*sizeof(double));
	  for(iGrid=0;iGrid<numGrid;iGrid++){
	    gridIndex = gridMapProc[iFrag][iGrid];
	    wfFragTemp[iGrid] = wfTemp[gridIndex];
	    rhoFragTemp[iGrid] = 0.0;
	  }
	  for(iStateFrag=0;iStateFrag<numStateDnMini;iStateFrag++){
	    proj = ddotBlasWrapper(numGrid,&wfFragTemp[0],1,&coefDnFragProc[iFrag][iStateFrag*numGrid],1);
	    fragInfo->wfProjUp[iFrag][countWf*numStateDnMini+iStateFrag] = proj*preDot*volMini;
	    daxpyBlasWrapper(numGrid,proj,&coefDnFragProc[iFrag][iStateFrag*numGrid],1,&rhoFragTemp[0],1);
	  }//endfor iGrid
	  for(iGrid=0;iGrid<numGrid;iGrid++){
	    gridIndex = gridMapProc[iFrag][iGrid];
	    rhoTemp[gridIndex] += rhoFragTemp[iGrid]*rhoFragTemp[iGrid];
	  }
	  free(wfFragTemp);
	  free(rhoFragTemp);
	}//endfor iFrag
	countWf += 1;
      }//endfor iState   
    }//endfor iProc
    if(numProcStates>1){
      Barrier(commStates);
      Reduce(rhoTemp,rhoTempReduce,rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
      Scatterv(rhoTempReduce,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	       rhoTemp,rhoRealGridNum,MPI_DOUBLE,0,commStates);
    }
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)numElecProj += rhoTemp[iGrid];
    daxpyBlasWrapper(rhoRealGridNum,-pre,&rhoTemp[0],1,&rhoDnFragSum[0],1);
  }

  
  if(numProcStates>1){
    //printf("numElecProj %lg\n",numElecProj*preNe);
    Allreduce(&numElecProj,&(stodftInfo->numElecTrueFrag),1,MPI_DOUBLE,MPI_SUM,0,commStates);
    stodftInfo->numElecTrueFrag *= preNe;
    if(myidState==0)printf("Number of Electron for StoWf is %.16lg\n",stodftInfo->numElecTrueFrag);
  }
  else{
    //printf("numElecTrue %.16lg\n",stodftInfo->numElecTrueFrag);
    stodftInfo->numElecTrueFrag = numElecProj*preNe;
    printf("Number of Electron for StoWf is %.16lg\n",stodftInfo->numElecTrueFrag);
  }
  
  // I would like to seperate them, but this will make life easier
  stodftInfo->numElecTrue = stodftInfo->numElecTrueFrag;
  //debug
  /*
  for(iProc=0;iProc<numProcStates;iProc++){
    if(myidState==iProc){
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
	printf("rhodiffff %lg %lg\n",rhoUpFragSum[iGrid],fragWfCpy[iGrid]);
      }
    }
    if(numProcStates>1)Barrier(commStates);
  }
  if(numProcStates>1)Barrier(commStates);
  exit(0);
  */

/*======================================================================*/
/* VI) Free memories                                                    */

  //debug
  /*
  for(iCoeff=0;iCoeff<numCoeffTotal;iCoeff++){
    cp->cpcoeffs_pos[1].cre_up[iCoeff+1] = noiseReUpBackup[iCoeff];
    cp->cpcoeffs_pos[1].cim_up[iCoeff+1] = noiseImUpBackup[iCoeff];
  } 
  free(noiseReUpBackup);
  free(noiseImUpBackup);
  */


  if(numProcStates>1&&myidState==0)free(rhoTempReduce);
  free(rhoTemp);
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    //free(coefUpFragProc[iFrag]);   
  }
  //free(coefUpFragProc);
  //free(noiseWfUpReal);
  free(numStateUpAllProc);
  if(cpLsda==1&&numStateDn!=0){
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      //free(coefDnFragProc[iFrag]);
    }
    //free(coefDnFragProc);
    //free(noiseWfDnReal);
    free(numStateDnAllProc);
  }
  //fflush(stdout);
  //exit(0);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void projRhoMiniUnitCell(CP *cp,GENERAL_DATA *general_data,CLASS *class,
		 CP *cpMini,GENERAL_DATA *generalDataMini,CLASS *classMini,
		 int ip_now)
/*========================================================================*/
{/*begin routine*/
/*************************************************************************/
/* This is the driver for calculating the fragmentation correction to    */
/* sDFT. This part used to be a whole function. Now it is been splited   */
/* into two functions so that energy window */
/*************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  
  STODFTINFO *stodftInfo = cp->stodftInfo;

  int energyWindowOn = stodftInfo->energyWindowOn;
  
  combineRhoUC(cp,general_data,class,
                     cpMini,generalDataMini,classMini,
                     ip_now);
  if(energyWindowOn==0){
    combineStoUC(cp,general_data,class,
                 cpMini,generalDataMini,classMini,
                 ip_now);
  }
  // The combineStoUCEnergyWindow is moved to genStoOrbitalEnergyWindow
  /*
  else{
    combineStoUCEnergyWindow(cp,general_data,class,cpMini,
                 generalDataMini,classMini,ip_now);

  }
  */

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void combineRhoUC(CP *cp,GENERAL_DATA *general_data,CLASS *class,
		 CP *cpMini,GENERAL_DATA *generalDataMini,CLASS *classMini,
		 int ip_now)
/*========================================================================*/
{/*begin routine*/
/*************************************************************************/
/* This first part calculate the density from deterministic calculation  */
/*************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  STODFTINFO *stodftInfo = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  CPOPTS *cpOpts = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cp->cpcoeffs_info);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  COMMUNICATE *commCP = &(cp->communicate);
  CELL *cell	      = &(general_data->cell);
 
  int iFrag,jFrag,iGrid,jGrid,iState,jState,iStateFrag,iProc;
  int iAtom;
  int gridIndex;
  int numGrid,numGridSmall;
  int numStateUpMini;
  int numStateDn;
  int numStateDnMini;
  int cpLsda		    = cpOpts->cp_lsda;
  int numFragProc	    = fragInfo->numFragProc;
  int numFragTot	    = fragInfo->numFragTot;
  int occNumber		    = stodftInfo->occNumber;
  int myidState		    = commCP->myid_state;
  int numProcStates	    = commCP->np_states;
  int rhoRealGridNum	    = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot	    = stodftInfo->rhoRealGridTot;
  int numStateUpProc	    = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc	    = cpcoeffs_info->nstate_dn_proc;
  int numStateStoUp     = stodftInfo->numStateStoUp;
  int numStateStoDn     = stodftInfo->numStateStoDn;
  int numCoeff		= cpcoeffs_info->ncoef;
  int countWf;
  int fragOpt		    = stodftInfo->fragOpt;

  MPI_Comm commStates   =    commCP->comm_states;

  int *numGridFragProc	    = fragInfo->numGridFragProc;
  int *numGridFragProcSmall = fragInfo->numGridFragProcSmall;
  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;
  int *numStateUpAllProc;
  int *numStateDnAllProc;

  int **gridMapProc	    = fragInfo->gridMapProc;
  int **gridMapProcSmall    = fragInfo->gridMapProcSmall;
 
  double proj;
  double numElecProj = 0.0;
  double numElecFrag = 0.0;
  double numElecProjTot = 0.0;
  double numElecFragTot = 0.0;
  double numElecSys = stodftInfo->numElecTrue;
  double numGridTotInv = 1.0/rhoRealGridTot;
  double vol,volInv;
  double volMini;
  double pre,preNe;
  double preDot;
  double *rhoTemp,*rhoTempReduce;
  double *wfTemp,*wfFragTemp,*rhoFragTemp;
  double *hmat_cp = cell->hmat_cp;
  double *hmatCpMini;
  double **rhoUpFragProc;
  double **coefUpFragProc;
  double **rhoDnFragProc;
  double **coefDnFragProc;
  double **coefUpFragCoreProc;
  double **coefDnFragCoreProc;
  double *rhoUpFragSum;
  double *rhoDnFragSum;
  double *noiseWfUpReal,*noiseWfDnReal;

  
/*======================================================================*/
/* I) Allocated memories				                */

  vol = getdeth(hmat_cp);
  volInv = 1.0/vol;
  
  // prefactor for prjection part
  pre = (double)(occNumber)*vol*vol/
	((double)(rhoRealGridTot)*(double)(rhoRealGridTot)*(double)(numStateStoUp));
  // prefactor for number of e in proj part
  preNe = pre/(double)(rhoRealGridTot);
  preDot = 1.0/sqrt(vol);
  
  fragInfo->rhoUpFragProc = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->coefUpFragProc = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->coefUpFragCoreProc = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->noiseWfUpReal = (double*)cmalloc(numStateUpProc*rhoRealGridTot*sizeof(double));
  /*
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    numGrid = numGridFragProc[iFrag];
    numStateUpMini = cpMini[iFrag].cpcoeffs_info.nstate_up_proc;
    fragInfo->rhoUpFragProc[iFrag] = (double*)cmalloc(numGrid*sizeof(double));
    fragInfo->coefUpFragProc[iFrag] = (double*)cmalloc(numStateUpMini*numGrid*sizeof(double));
  }
  */
  rhoUpFragProc = fragInfo->rhoUpFragProc;
  coefUpFragProc = fragInfo->coefUpFragProc;
  coefUpFragCoreProc = fragInfo->coefUpFragCoreProc;
  if(cpLsda==1&&numStateDn!=0){
    fragInfo->rhoDnFragProc = (double**)cmalloc(numFragProc*sizeof(double*));
    fragInfo->coefDnFragProc = (double**)cmalloc(numFragProc*sizeof(double*));
    fragInfo->coefDnFragCoreProc = (double**)cmalloc(numFragProc*sizeof(double*));
    fragInfo->noiseWfDnReal = (double*)cmalloc(numStateDnProc*rhoRealGridTot*sizeof(double));
    /*
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      numGrid = numGridFragProc[iFrag];
      numStateDnMini = cpMini[iFrag].cpcoeffs_info.nstate_dn_proc;
      fragInfo->rhoDnFragProc[iFrag] = (double*)cmalloc(numGrid*sizeof(double));
      fragInfo->coefDnFragProc[iFrag] = (double*)cmalloc(numStateDnMini*numGrid*sizeof(double));
      fragInfo->coefDnFragCoreProc = (double**)cmalloc(numFragProc*sizeof(double*));

    }
    */
    rhoDnFragProc = fragInfo->rhoDnFragProc;
    coefDnFragProc = fragInfo->coefDnFragProc;
    coefDnFragCoreProc = fragInfo->coefDnFragCoreProc;
  }//endif

  fragInfo->rhoUpFragSum = (double*)cmalloc(rhoRealGridNum*sizeof(double));
  fragInfo->rhoUpFragSumCpy = (double*)cmalloc(rhoRealGridNum*sizeof(double));
  rhoUpFragSum = fragInfo->rhoUpFragSum;
  if(cpLsda==1&&numStateDn!=0){
    fragInfo->rhoDnFragSum = (double*)cmalloc(rhoRealGridNum*sizeof(double));
    fragInfo->rhoDnFragSumCpy = (double*)cmalloc(rhoRealGridNum*sizeof(double));
    rhoDnFragSum = fragInfo->rhoDnFragSum;
  }

/*======================================================================*/
/* II) Transform fragment wf to real spaces and store them              */

  // Allocate rhoTemp for each proc
  rhoTemp = (double*)cmalloc(rhoRealGridTot*sizeof(double));
  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
    rhoTemp[iGrid] = 0.0;
  }//endfor iGrid
  // Reduce all fragment densities in each proc


  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInfo->iFrag = iFrag;
    numGrid = numGridFragProc[iFrag];
    numGridSmall = numGridFragProcSmall[iFrag];
    numStateUpMini = cpMini[iFrag].cpcoeffs_info.nstate_up_proc;
    numStateDnMini = cpMini[iFrag].cpcoeffs_info.nstate_dn_proc;
    fragInfo->rhoUpFragProc[iFrag] = (double*)cmalloc(numGrid*sizeof(double));
    fragInfo->coefUpFragProc[iFrag] = (double*)cmalloc(numStateUpMini*numGrid*sizeof(double));
    coefUpFragCoreProc[iFrag] = (double*)cmalloc(numStateUpMini*numGridSmall*sizeof(double));
    if(cpLsda==1&&numStateDn!=0){
      numStateDnMini = cpMini[iFrag].cpcoeffs_info.nstate_dn_proc;
      fragInfo->rhoDnFragProc[iFrag] = (double*)cmalloc(numGrid*sizeof(double));
      fragInfo->coefDnFragProc[iFrag] = (double*)cmalloc(numStateDnMini*numGrid*sizeof(double));
      coefDnFragCoreProc[iFrag] = (double*)cmalloc(numStateDnMini*numGridSmall*sizeof(double));
    }

    rhoRealCalcDriverFragMol(&generalDataMini[iFrag],&cpMini[iFrag],&classMini[iFrag],cp);
    for(iGrid=0;iGrid<numGridSmall;iGrid++){
      rhoTemp[gridMapProc[iFrag][gridMapProcSmall[iFrag][iGrid]]] +=
                rhoUpFragProc[iFrag][gridMapProcSmall[iFrag][iGrid]];
    }//endfor iGrid
    for(iState=0;iState<numStateUpMini;iState++){
      for(iGrid=0;iGrid<numGridSmall;iGrid++){
	coefUpFragCoreProc[iFrag][iState*numGridSmall+iGrid] = coefUpFragProc[iFrag][iState*numGrid+gridMapProcSmall[iFrag][iGrid]];
      }
    }
    if(cpLsda==1&&numStateDn!=0){
      for(iState=0;iState<numStateDnMini;iState++){
	for(iGrid=0;iGrid<numGridSmall;iGrid++){
	  coefDnFragCoreProc[iFrag][iState*numGridSmall+iGrid] = coefDnFragProc[iFrag][iState*numGrid+gridMapProcSmall[iFrag][iGrid]];
	}
      }
    }
#ifndef DMAT
    free(fragInfo->rhoUpFragProc[iFrag]);
    free(fragInfo->coefUpFragProc[iFrag]);
#endif
  }//endfor iFrag
  Barrier(commStates);
  //fflush(stdout);
  //exit(0);
  //debug

/*======================================================================*/
/* III) Assemble fragments densities		                        */

  /*
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    //numGrid = numGridFragProc[iFrag];
    numGrid = numGridFragProcSmall[iFrag];
    for(iGrid=0;iGrid<numGrid;iGrid++){
      rhoTemp[gridMapProc[iFrag][gridMapProcSmall[iFrag][iGrid]]] += 
		rhoUpFragProc[iFrag][gridMapProcSmall[iFrag][iGrid]];
    }//endfor iGrid
  }//endfor iFrag
  */
  /*
  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
    printf("111111 rhofrag %.8lg\n",rhoTemp[iGrid]*vol);
  }
  */
  //fflush(stdout);
  //exit(0);
  
  
  // Allocate rhoTempReduce on the master proc
  if(myidState==0&&numProcStates>1){
    rhoTempReduce = (double*)cmalloc(rhoRealGridTot*sizeof(double));
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoTempReduce[iGrid] = 0.0;
  }
  // Reduce fragment densities on different procs and scatter them
  if(numProcStates>1){
    Barrier(commStates);
    Reduce(rhoTemp,rhoTempReduce,rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Scatterv(rhoTempReduce,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
             rhoUpFragSum,rhoRealGridNum,MPI_DOUBLE,0,commStates);
    // I need to multiply the density by volumn to match the |ksi(r)|^2
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
      rhoUpFragSum[iGrid] *= vol;
    }
    //if(myidState==0)free(&rhoTempReduce[0]);
  }
  // Sequential, just copy and scale the density
  else{
    //memcpy(rhoUpFragSum,rhoTemp,rhoRealGridNum*sizeof(double));
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
      rhoUpFragSum[iGrid] = rhoTemp[iGrid]*vol;
    }
  }
  memcpy(&(fragInfo->rhoUpFragSumCpy[0]),rhoUpFragSum,rhoRealGridNum*sizeof(double));
  // Do the same thing for spin down states
  if(cpLsda==1&&numStateDn!=0){
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      rhoTemp[iGrid] = 0.0;
    }
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      //numGrid = numGridFragProc[iFrag];
      numGrid = numGridFragProcSmall[iFrag];
      for(iGrid=0;iGrid<numGrid;iGrid++){
	//rhoTemp[gridMapProc[iFrag][iGrid]] += rhoDnFragProc[iFrag][iGrid];
        rhoTemp[gridMapProc[iFrag][gridMapProcSmall[iFrag][iGrid]]] +=
                  rhoDnFragProc[iFrag][gridMapProcSmall[iFrag][iGrid]];

      }
    }

    if(myidState==0&&numProcStates>1){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoTempReduce[iGrid] = 0.0;
    }
    if(numProcStates>1){
      Barrier(commStates);
      Reduce(rhoTemp,rhoTempReduce,rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
      Scatterv(rhoTempReduce,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	       rhoDnFragSum,rhoRealGridNum,MPI_DOUBLE,0,commStates);
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++)rhoUpFragSum[iGrid] *= vol;
    }
    else{
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
	rhoDnFragSum[iGrid] = rhoTemp[iGrid]*vol;
      }//endfor iGrid
    }//endif 
    memcpy(&(fragInfo->rhoDnFragSumCpy[0]),rhoDnFragSum,rhoRealGridNum*sizeof(double));
  }//endif cpLsda

/*======================================================================*/
/* IV) Free/allocate some memory for next step                          */

  /*
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    free(fragInfo->rhoUpFragProc[iFrag]);
  }
  
  free(fragInfo->rhoUpFragProc);
  if(cpLsda==1&&numStateDn!=0){
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      free(fragInfo->rhoDnFragProc[iFrag]);
    }
    free(fragInfo->rhoDnFragProc);   
  }
  */
  free(rhoTemp);
  if(myidState==0&&numProcStates>1)free(rhoTempReduce);
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void combineStoUC(CP *cp,GENERAL_DATA *general_data,CLASS *class,
                 CP *cpMini,GENERAL_DATA *generalDataMini,CLASS *classMini,
                 int ip_now)
/*========================================================================*/
{/*begin routine*/
/*************************************************************************/
/* The second part calculate the density from stochastic calculation     */
/*************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  STODFTINFO *stodftInfo = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  CPOPTS *cpOpts = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cp->cpcoeffs_info);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  COMMUNICATE *commCP = &(cp->communicate);
  CELL *cell	      = &(general_data->cell);
 
  int iFrag,jFrag,iGrid,jGrid,iState,jState,iStateFrag,iProc;
  int iAtom;
  int gridIndex;
  int numGrid,numGridSmall;
  int numStateUpMini;
  int numStateDn;
  int numStateDnMini;
  int cpLsda		    = cpOpts->cp_lsda;
  int numFragProc	    = fragInfo->numFragProc;
  int numFragTot	    = fragInfo->numFragTot;
  int occNumber		    = stodftInfo->occNumber;
  int myidState		    = commCP->myid_state;
  int numProcStates	    = commCP->np_states;
  int rhoRealGridNum	    = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot	    = stodftInfo->rhoRealGridTot;
  int numStateUpProc	    = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc	    = cpcoeffs_info->nstate_dn_proc;
  int numStateStoUp     = stodftInfo->numStateStoUp;
  int numStateStoDn     = stodftInfo->numStateStoDn;
  int numCoeff		= cpcoeffs_info->ncoef;
  int countWf;
  int fragOpt		    = stodftInfo->fragOpt;

  MPI_Comm commStates   =    commCP->comm_states;

  int *numGridFragProc	    = fragInfo->numGridFragProc;
  int *numGridFragProcSmall = fragInfo->numGridFragProcSmall;
  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;
  int *numStateUpAllProc;
  int *numStateDnAllProc;

  int **gridMapProc	    = fragInfo->gridMapProc;
  int **gridMapProcSmall    = fragInfo->gridMapProcSmall;
 
  double proj;
  double numElecProj = 0.0;
  double numElecFrag = 0.0;
  double numElecProjTot = 0.0;
  double numElecFragTot = 0.0;
  double numElecSys = stodftInfo->numElecTrue;
  double numGridTotInv = 1.0/rhoRealGridTot;
  double vol,volInv;
  double volMini;
  double pre,preNe;
  double preDot;
  double *rhoTemp,*rhoTempReduce;
  double *wfTemp,*wfFragTemp,*rhoFragTemp;
  double *hmat_cp = cell->hmat_cp;
  double *hmatCpMini;
  double **rhoUpFragProc;
  double **coefUpFragProc;
  double **rhoDnFragProc;
  double **coefDnFragProc;
  double **coefUpFragCoreProc;
  double **coefDnFragCoreProc;
  double *rhoUpFragSum;
  double *rhoDnFragSum;
  double *noiseWfUpReal,*noiseWfDnReal;

  
/*======================================================================*/
/* I) Asign local variables				                */

  vol = getdeth(hmat_cp);
  volInv = 1.0/vol;

  // prefactor for prjection part
  pre = (double)(occNumber)*vol*vol/
        ((double)(rhoRealGridTot)*(double)(rhoRealGridTot)*(double)(numStateStoUp));
  // prefactor for number of e in proj part
  preNe = pre/(double)(rhoRealGridTot);
  preDot = 1.0/sqrt(vol);

  rhoUpFragProc = fragInfo->rhoUpFragProc;
  coefUpFragProc = fragInfo->coefUpFragProc;
  coefUpFragCoreProc = fragInfo->coefUpFragCoreProc;
  rhoUpFragSum = fragInfo->rhoUpFragSum;
  noiseWfUpReal = fragInfo->noiseWfUpReal;

  if(cpLsda==1&&numStateDn!=0){
    rhoDnFragProc = fragInfo->rhoDnFragProc;
    coefDnFragProc = fragInfo->coefDnFragProc;
    coefDnFragCoreProc = fragInfo->coefDnFragCoreProc;    
    rhoDnFragSum = fragInfo->rhoDnFragSum;
    noiseWfDnReal = fragInfo->noiseWfDnReal;
  }

  rhoTemp = (double*)cmalloc(rhoRealGridTot*sizeof(double));
  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
    rhoTemp[iGrid] = 0.0;
  }//endfor iGrid


  /*
  fragInfo->wfProjUp = (double**)cmalloc(numFragProc*sizeof(double*));

  for(iFrag=0;iFrag<numFragProc;iFrag++){
    numStateUpMini = cpMini[iFrag].cpcoeffs_info.nstate_up_proc;
    fragInfo->wfProjUp[iFrag] = (double*)cmalloc(numStateStoUp*numStateUpMini*sizeof(double));
  }

  if(cpLsda==1&&numStateDn!=0){
    fragInfo->wfProjDn = (double**)cmalloc(numFragProc*sizeof(double*));
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      numStateDnMini = cpMini[iFrag].cpcoeffs_info.nstate_dn_proc;
      fragInfo->wfProjDn[iFrag] = (double*)cmalloc(numStateStoDn*numStateDnMini*sizeof(double));
    }
  }
  */

/*======================================================================*/
/* IV) Calculate the real space noise wave function                     */

  // Copy the deterministic wf to the stochastic wf
  //memcpy(&(cp->cpcoeffs_pos[1].cre_up[1]),&(stodftCoefPos->wfDetBackupUpRe[0]),numStateUpProc*numCoeff*sizeof(double));
  //memcpy(&(cp->cpcoeffs_pos[1].cim_up[1]),&(stodftCoefPos->wfDetBackupUpIm[0]),numStateUpProc*numCoeff*sizeof(double));
  
#ifdef DMAT
  printf("Test Density Matrix\n");
  int iCoeff;
  for(iCoeff=0;iCoeff<numStateUpProc*numCoeff;iCoeff++){
    cp->cpcoeffs_pos[1].cre_up[iCoeff+1] = stodftCoefPos->wfDetBackupUpRe[iCoeff];
    cp->cpcoeffs_pos[1].cim_up[iCoeff+1] = stodftCoefPos->wfDetBackupUpIm[iCoeff];
  }
  
  rhoRealCalcDriverNoise(general_data,cp,class,ip_now);
  free(stodftCoefPos->wfDetBackupUpRe);
  free(stodftCoefPos->wfDetBackupUpIm);
#else
  noiseRealReGen(general_data,cp,class,ip_now);
#endif

/*======================================================================*/
/* IV) Project the real space noise wave function                       */


  numStateUpAllProc = (int*)cmalloc(numProcStates*sizeof(int));
  if(numProcStates>1){
    Allgather(&numStateUpProc,1,MPI_INT,numStateUpAllProc,1,MPI_INT,0,commStates);
  }
  else numStateUpAllProc[0] = numStateUpProc;
  if(cpLsda==1&&numStateDn!=0){
    numStateDnAllProc = (int*)cmalloc(numProcStates*sizeof(int));
    if(numProcStates>1){
      Allgather(&numStateDnProc,1,MPI_INT,numStateDnAllProc,1,MPI_INT,0,commStates);
    }
    else numStateDnAllProc[0] = numStateDnProc;
  }
  wfTemp = (double*)cmalloc(rhoRealGridTot*sizeof(double));
  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
    rhoTemp[iGrid] = 0.0;
  }
  //countWf = 0;
  //end debug  
  //debug
  /*
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInfo->rhoUpFragProc[iFrag] = (double*)calloc(numGrid,sizeof(double));
    fragInfo->coefUpFragProc[iFrag] = (double*)calloc(numStateUpMini*numGrid,sizeof(double));
    rhoRealCalcDriverFragMol(&generalDataMini[iFrag],&cpMini[iFrag],&classMini[iFrag],cp);
  }
  */
  
#ifdef DMAT
  /*
  // Test1, map the fragment density matrix on the system occupied orbitals
  double *projMatrixTest = (double*)calloc(numStateUpProc*numStateUpProc,sizeof(double));
  double *wfProjjTemp;
  double *coefUpFragExt;
  double preCopyFrag = sqrt(volInv);
  for(iState=0;iState<numStateUpProc;iState++){
    printf("iState %i\n",iState);
    memcpy(wfTemp,&noiseWfUpReal[iState*rhoRealGridTot],rhoRealGridTot*sizeof(double));
    if(numProcStates>1)Bcast(wfTemp,rhoRealGridTot,MPI_DOUBLE,iProc,commStates);
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      numGrid = numGridFragProc[iFrag];
      numGridSmall = numGridFragProcSmall[iFrag];

      wfProjjTemp = (double*)calloc(numGridSmall,sizeof(double));
      //coefUpFragExt = (double*)calloc(numStateUpMini*numGridSmall,sizeof(double));

      numStateUpMini = cpMini[iFrag].cpcoeffs_info.nstate_up_proc;
      hmatCpMini = generalDataMini[iFrag].cell.hmat_cp;
      volMini = getdeth(hmatCpMini);
      volMini /= numGrid;
      wfFragTemp = (double*)cmalloc(numGrid*sizeof(double));
      rhoFragTemp = (double*)cmalloc(numGrid*sizeof(double));
      //rhoRealCalcDriverFragMol(&generalDataMini[iFrag],&cpMini[iFrag],&classMini[iFrag],cp);
      for(iGrid=0;iGrid<numGrid;iGrid++){
        gridIndex = gridMapProc[iFrag][iGrid];
        wfFragTemp[iGrid] = wfTemp[gridIndex];
        rhoFragTemp[iGrid] = 0.0;
      }
      for(iStateFrag=0;iStateFrag<numStateUpMini;iStateFrag++){
	proj = 0.0;
	for(iGrid=0;iGrid<numGrid;iGrid++){
	  proj += wfFragTemp[iGrid]*coefUpFragProc[iFrag][iStateFrag*numGrid+iGrid];
	}
	proj *= volMini;
	for(iGrid=0;iGrid<numGridSmall;iGrid++){
	  gridIndex = gridMapProcSmall[iFrag][iGrid];
	  wfProjjTemp[iGrid] += proj*coefUpFragProc[iFrag][iStateFrag*numGrid+gridIndex];
	}
      }//endfor iStateFrag
      #pragma omp parallel for private(jState,iGrid,gridIndex)
      for(jState=0;jState<numStateUpProc;jState++){
        for(iGrid=0;iGrid<numGridSmall;iGrid++){
	  gridIndex = gridMapProc[iFrag][gridMapProcSmall[iFrag][iGrid]];
          projMatrixTest[iState*numStateUpProc+jState] += noiseWfUpReal[jState*rhoRealGridTot+gridIndex]*wfProjjTemp[iGrid]*volMini;
        }
      }//endfor jState
      free(wfFragTemp);
      free(rhoFragTemp);
      free(wfProjjTemp);
    }//endfor iFrag
    countWf += 1;
    for(jState=0;jState<numStateUpProc;jState++){
      projMatrixTest[iState*numStateUpProc+jState] *= 0.5*volInv;
      //projMatrixTest[jState*numStateUpProc+iState] = projMatrixTest[iState*numStateUpProc+jState];
    }
  }//endfor iState
  FILE *testMatrix = fopen("test-matrix","w");
  for(iState=0;iState<numStateUpProc;iState++){
    for(jState=0;jState<numStateUpProc;jState++){
      if(iState==jState)projMatrixTest[iState*numStateUpProc+jState] -= 1.0;
      fprintf(testMatrix,"%.8lg ",projMatrixTest[iState*numStateUpProc+jState]);
    }
    fprintf(testMatrix,"\n");
    printf("111111 matrix diag %.8lg\n",projMatrixTest[iState*numStateUpProc+iState]+1.0);
  }
  */
  // Test2, test the bound for y=z=0
  int numGridTest = cp_para_fft_pkg3d_lg->nkf1;
  int gridTest,fragInd;
  int *fragIndTest = (int*)cmalloc(numGridTest*sizeof(int));
  int *gridIndDress = (int*)cmalloc(numGridTest*sizeof(int));
  // searching fragment containing testing grid point in core region
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    numGridSmall = numGridFragProcSmall[iFrag];
    for(iGrid=0;iGrid<numGridSmall;iGrid++){
      gridTest = gridMapProc[iFrag][gridMapProcSmall[iFrag][iGrid]];
      if(gridTest<numGridTest){
        fragIndTest[gridTest] = iFrag;
        gridIndDress[gridTest] = gridMapProcSmall[iFrag][iGrid];
      }
    }
  }
  double *dmDet = (double*)calloc(numGridTest*rhoRealGridTot,sizeof(double));
  double *dmFrag = (double*)calloc(numGridTest*rhoRealGridTot,sizeof(double));
  //deterministic part
  for(iState=0;iState<numStateUpProc;iState++){
    #pragma omp parallel for private(iGrid,jGrid)
    for(iGrid=0;iGrid<numGridTest;iGrid++){
      for(jGrid=0;jGrid<rhoRealGridTot;jGrid++){
        dmDet[iGrid*rhoRealGridTot+jGrid] += noiseWfUpReal[iState*rhoRealGridTot+iGrid]*
                                        noiseWfUpReal[iState*rhoRealGridTot+jGrid];
      }
    }
  }
  #pragma omp parallel for private(iGrid)
  for(iGrid=0;iGrid<numGridTest*rhoRealGridTot;iGrid++){
    dmDet[iGrid] *= 0.5*volInv;
  }
  //next stochastic part
  for(iGrid=0;iGrid<numGridTest;iGrid++){
    fragInd = fragIndTest[iGrid];
    numStateUpMini = cpMini[fragInd].cpcoeffs_info.nstate_up_proc;
    hmatCpMini = generalDataMini[fragInd].cell.hmat_cp;
    volMini = getdeth(hmatCpMini);
    volMini /= numGrid;
    numGrid = numGridFragProc[fragInd];
    for(iStateFrag=0;iStateFrag<numStateUpMini;iStateFrag++){
      #pragma omp parallel for private(jGrid)
      for(jGrid=0;jGrid<numGrid;jGrid++){
        dmFrag[iGrid*rhoRealGridTot+gridMapProc[fragInd][jGrid]] += coefUpFragProc[fragInd][iStateFrag*numGrid+gridIndDress[iGrid]]*coefUpFragProc[fragInd][iStateFrag*numGrid+jGrid];
      }
    }
     
  }
  //next lets worry about max(<r|rho|r'>^2)
  double *dmDetMax = (double*)calloc(numGridTest,sizeof(double));
  double dmMax,dmSq;
  #pragma omp parallel for private(iGrid,jGrid,dmMax,dmSq)
  for(iGrid=0;iGrid<numGridTest;iGrid++){
    dmMax = -10000.0;
    for(jGrid=0;jGrid<rhoRealGridTot;jGrid++){
      dmSq = dmDet[iGrid*rhoRealGridTot+jGrid]*dmDet[iGrid*rhoRealGridTot+jGrid];
      if(dmSq>dmMax)dmMax = dmSq;
    }
    dmDetMax[iGrid] = dmMax;
  }
  //next lets work on <r|delta rho|r'>
  double *dmDiffSqInt = (double*)calloc(numGridTest,sizeof(double));
  double diff;
  #pragma omp parallel for private(iGrid,jGrid,diff)
  for(iGrid=0;iGrid<numGridTest;iGrid++){
    for(jGrid=0;jGrid<rhoRealGridTot;jGrid++){
      diff = dmDet[iGrid*rhoRealGridTot+jGrid]-dmFrag[iGrid*rhoRealGridTot+jGrid];
      dmDiffSqInt[iGrid] += diff*diff;
    }
    dmDiffSqInt[iGrid] *= volMini;
  }  
  for(iGrid=0;iGrid<numGridTest;iGrid++){
    printf("1111111 dmxxxxxx %.16lg %.16lg\n",dmDetMax[iGrid],dmDiffSqInt[iGrid]);
  }

  fflush(stdout);
  exit(0);   
#endif  
  //end debug  

  int res = numFragTot%numProcStates;
  int div = numFragTot/numProcStates;
  int numFragProcMax;
  if(res==0)numFragProcMax = div;
  else numFragProcMax = div+1;

  //DEBUG
  FILE *frho;
  FILE *frho_frag;

  //debug, output all wf filtered by frag wf
  //double *stowffrag = (double*)cmalloc(numStateUpAllProc[0]*rhoRealGridTot*sizeof(double));
  for(iFrag=0;iFrag<numFragProcMax;iFrag++){
    if(iFrag<numFragProc){
      fragInfo->iFrag = iFrag;
      numGrid = numGridFragProc[iFrag];
      numStateUpMini = cpMini[iFrag].cpcoeffs_info.nstate_up_proc;
      fragInfo->rhoUpFragProc[iFrag] = (double*)cmalloc(numGrid*sizeof(double));
      fragInfo->coefUpFragProc[iFrag] = (double*)cmalloc(numStateUpMini*numGrid*sizeof(double));
      rhoRealCalcDriverFragMol(&generalDataMini[iFrag],&cpMini[iFrag],&classMini[iFrag],cp);
    }//endif iFrag
    countWf = 0;
    for(iProc=0;iProc<numProcStates;iProc++){
      for(iState=0;iState<numStateUpAllProc[iProc];iState++){
	if(myidState==iProc){
	  memcpy(wfTemp,&noiseWfUpReal[iState*rhoRealGridTot],rhoRealGridTot*sizeof(double));
	}
	if(numProcStates>1)Bcast(wfTemp,rhoRealGridTot,MPI_DOUBLE,iProc,commStates);
	//double *rho_frag_ext = (double*)calloc(rhoRealGridTot,sizeof(double));
	//for(iFrag=0;iFrag<numFragProc;iFrag++){
	if(iFrag<numFragProc){
	  numGrid = numGridFragProc[iFrag];
	  numGridSmall = numGridFragProcSmall[iFrag];
	  numStateUpMini = cpMini[iFrag].cpcoeffs_info.nstate_up_proc;
	  hmatCpMini = generalDataMini[iFrag].cell.hmat_cp;
	  volMini = getdeth(hmatCpMini);
	  volMini /= numGrid;
	  wfFragTemp = (double*)calloc(numGrid,sizeof(double));
	  rhoFragTemp = (double*)calloc(numGridSmall,sizeof(double));
	  for(iGrid=0;iGrid<numGrid;iGrid++){
	    gridIndex = gridMapProc[iFrag][iGrid];
	    wfFragTemp[iGrid] = wfTemp[gridIndex];
	    //wfFragTemp[iGrid] = coefUpFragProc[iFrag][iGrid];
	  }
	  for(iGrid=0;iGrid<numGridSmall;iGrid++)rhoFragTemp[iGrid] = 0.0;
	  //double testsum = 0.0;
	  for(iStateFrag=0;iStateFrag<numStateUpMini;iStateFrag++){
	    //proj = ddotBlasWrapper(numGrid,&wfFragTemp[0],1,&coefUpFragProc[iFrag][iStateFrag*numGrid],1);
            
	    proj = 0.0;
	    for(iGrid=0;iGrid<numGrid;iGrid++){
	      proj += wfFragTemp[iGrid]*coefUpFragProc[iFrag][iStateFrag*numGrid+iGrid];
	    }
	    fragInfo->wfProjUp[iFrag][countWf*numStateUpMini+iStateFrag] = proj*preDot*volMini;
	    //daxpyBlasWrapper(numGrid,proj,&coefUpFragProc[iFrag][iStateFrag*numGrid],1,&rhoFragTemp[0],1);
	    for(iGrid=0;iGrid<numGridSmall;iGrid++){
	      rhoFragTemp[iGrid] += proj*coefUpFragProc[iFrag][iStateFrag*numGrid+gridMapProcSmall[iFrag][iGrid]];
	    }
	    //testsum += proj*coefUpFragProc[iFrag][iStateFrag*numGrid+gridMapProcSmall[iFrag][0]];
	  }//endfor iStateFrag
	  //debug
	  /*
	  for(iGrid=0;iGrid<numGridSmall;iGrid++){
	    gridIndex = gridMapProc[iFrag][gridMapProcSmall[iFrag][iGrid]];
	    stowffrag[iState*rhoRealGridTot+gridIndex] = rhoFragTemp[iGrid];
	  }
	  */
	  //printf("iFrag %i testsum %lg\n",iFrag,testsum*testsum);
	  for(iGrid=0;iGrid<numGridSmall;iGrid++){
	    gridIndex = gridMapProc[iFrag][gridMapProcSmall[iFrag][iGrid]];
	    rhoTemp[gridIndex] += rhoFragTemp[iGrid]*rhoFragTemp[iGrid];
	    //fix_frag[iFrag][gridIndex] += rhoFragTemp[iGrid]*rhoFragTemp[iGrid]*pre;
	  }
	  free(wfFragTemp);
	  free(rhoFragTemp);
	}//endif iFrag<numFragProc
	//}//endfor iFrag
	countWf += 1;
      }//endfor iState   
    }//endfor iProc
    if(iFrag<numFragProc){
      free(fragInfo->rhoUpFragProc[iFrag]);
      free(fragInfo->coefUpFragProc[iFrag]);
    }
  }//endfor iFrag

  /*
  char name[100];
  FILE *fstowf;
  for(iState=0;iState<numStateUpAllProc[0];iState++){
    sprintf(name,"stowf-frag-%i",iState);
    fstowf = fopen(name,"w");
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      fprintf(fstowf,"%.16lg\n",stowffrag[iState*rhoRealGridTot+iGrid]);
    }
    fclose(fstowf);
  }
  */
  /*--------------------------------------------------*/
  /* We are going to calculate number of electrons.   */
  /* The Formula should be N_e^{frag}-N_e^{proj}+     */
  /* N_e^{sto} = N_e.                                 */
  /*--------------------------------------------------*/
  
  if(numProcStates>1){
    if(myidState==0){
      rhoTempReduce = (double*)cmalloc(rhoRealGridTot*sizeof(double));
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoTempReduce[iGrid] = 0.0;
    }
    Barrier(commStates);
    Reduce(rhoTemp,rhoTempReduce,rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    //DEBUG
    /*
    if(myidState==0){
      frho = fopen("density-cor","w");
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
        fprintf(frho,"%.16lg\n",rhoTempReduce[iGrid]*pre);
      }
      fclose(frho);
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoTempReduce[iGrid] = 0.0;
    }
    Barrier(commStates);
    Gatherv(rhoUpFragSum,rhoRealGridNum,MPI_DOUBLE,rhoTempReduce,rhoRealSendCounts,
            rhoRealDispls,MPI_DOUBLE,0,commStates);
    if(myidState==0){
      frho_frag = fopen("density-frag","w");
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
        fprintf(frho_frag,"%.16lg\n",rhoTempReduce[iGrid]);
      }
      fclose(frho_frag);
    }
    */
    Barrier(commStates);
    //fflush(stdout);
    //exit(0);
    Scatterv(rhoTempReduce,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	     rhoTemp,rhoRealGridNum,MPI_DOUBLE,0,commStates);
  } 
  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
    numElecProj += rhoTemp[iGrid];
    numElecFrag += rhoUpFragSum[iGrid];
  }
  
  double sumElecFrag = 0.0;
  double sumElecProj = 0.0;
  for(iProc=0;iProc<numProcStates;iProc++){
    if(myidState==iProc){
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
	sumElecFrag += rhoUpFragSum[iGrid];
	sumElecProj += pre*rhoTemp[iGrid];
      }
      sumElecFrag /= rhoRealGridTot;
      sumElecProj /= rhoRealGridTot;
    }
    if(numProcStates>1)Barrier(commStates);
  }
  //printf("11111111111111 sumElecFrag %.16lg sumElecProj %.16lg\n",sumElecFrag,sumElecProj);
  if(numProcStates>1)Barrier(commStates);
  
  daxpyBlasWrapper(rhoRealGridNum,-pre,&rhoTemp[0],1,&rhoUpFragSum[0],1);
#ifdef TUNE_NOISE
  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
    rhoUpFragSum[iGrid] *= NOISE_FACTOR;
  }
  sumElecFrag *= NOISE_FACTOR;
  sumElecProj *= NOISE_FACTOR;
#endif
  if(numProcStates>1)Barrier(commStates);
  
  /*
  for(iProc=0;iProc<numProcStates;iProc++){
    if(myidState==iProc){
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
        printf("111111111111 rhofix %i %lg %lg\n",iGrid,rhoUpFragSum[iGrid],rhoTemp[iGrid]*pre);
      }
    }
    if(numProcStates>1)Barrier(commStates);
  }
  */
  
  //fflush(stdout);
  //exit(0);

  if(cpLsda==1&&numStateDn!=0){
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      rhoTemp[iGrid] = 0.0;
    }
    countWf = 0;
    for(iProc=0;iProc<numProcStates;iProc++){
      for(iState=0;iState<numStateDnAllProc[iProc];iState++){
	if(myidState==iProc){
	  memcpy(wfTemp,&noiseWfDnReal[iState*rhoRealGridTot],rhoRealGridTot*sizeof(double));
	}
	if(numProcStates>1)Bcast(wfTemp,rhoRealGridTot,MPI_DOUBLE,iProc,commStates);
	for(iFrag=0;iFrag<numFragProc;iFrag++){
	  numGrid = numGridFragProc[iFrag];
	  numGridSmall = numGridFragProcSmall[iFrag];
	  numStateDnMini = cpMini[iFrag].cpcoeffs_info.nstate_dn_proc;
	  hmatCpMini = generalDataMini[iFrag].cell.hmat_cp;
	  volMini = getdeth(hmatCpMini);
	  volMini /= numGrid;
	  wfFragTemp = (double*)cmalloc(numGrid*sizeof(double));
	  rhoFragTemp = (double*)cmalloc(numGridSmall*sizeof(double));
	  for(iGrid=0;iGrid<numGrid;iGrid++){
	    gridIndex = gridMapProc[iFrag][iGrid];
	    wfFragTemp[iGrid] = wfTemp[gridIndex];
	  }
	  for(iGrid=0;iGrid<numGridSmall;iGrid++)rhoFragTemp[iGrid] = 0.0;
	  for(iStateFrag=0;iStateFrag<numStateDnMini;iStateFrag++){
	    proj = ddotBlasWrapper(numGrid,&wfFragTemp[0],1,&coefDnFragProc[iFrag][iStateFrag*numGrid],1);
	    fragInfo->wfProjUp[iFrag][countWf*numStateDnMini+iStateFrag] = proj*preDot*volMini;
	    //daxpyBlasWrapper(numGrid,proj,&coefDnFragProc[iStateFrag*numGrid],1,&rhoFragTemp[0],1);
	    for(iGrid=0;iGrid<numGridSmall;iGrid++){
	      rhoFragTemp[iGrid] += proj*coefDnFragProc[iFrag][iStateFrag*numGrid+gridMapProcSmall[iFrag][iGrid]];
	    } 
	  }//endfor iGrid
	  for(iGrid=0;iGrid<numGridSmall;iGrid++){
	    gridIndex = gridMapProc[iFrag][gridMapProcSmall[iFrag][iGrid]];
	    rhoTemp[gridIndex] += rhoFragTemp[iGrid]*rhoFragTemp[iGrid];
	  }
	  free(wfFragTemp);
	  free(rhoFragTemp);
	}//endfor iFrag
	countWf += 1;
      }//endfor iState   
    }//endfor iProc
    if(numProcStates>1){
      Barrier(commStates);
      Reduce(rhoTemp,rhoTempReduce,rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
      Scatterv(rhoTempReduce,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	       rhoTemp,rhoRealGridNum,MPI_DOUBLE,0,commStates);
    }
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
      numElecProj += rhoTemp[iGrid];
      numElecFrag += rhoDnFragSum[iGrid];
    }
    daxpyBlasWrapper(rhoRealGridNum,-pre,&rhoTemp[0],1,&rhoDnFragSum[0],1);
  }//endif cpLsda

  
  if(numProcStates>1){
    //printf("numElecProj %lg\n",numElecProj*preNe);
    //Allreduce(&numElecProj,&(stodftInfo->numElecTrueFrag),1,MPI_DOUBLE,MPI_SUM,0,commStates);
    //stodftInfo->numElecTrueFrag *= preNe;
   
    Allreduce(&numElecProj,&numElecProjTot,1,MPI_DOUBLE,MPI_SUM,0,commStates);
    Allreduce(&numElecFrag,&numElecFragTot,1,MPI_DOUBLE,MPI_SUM,0,commStates);
    stodftInfo->numElecTrueFrag = numElecProjTot*preNe-numElecFragTot/rhoRealGridTot+
                                  numElecSys;
    if(myidState==0)printf("Number of Electron for StoWf is %.16lg %.16lg %.16lg %lg\n",
                           numElecProjTot*preNe,numElecFragTot/rhoRealGridTot,
                           stodftInfo->numElecTrueFrag,numElecSys);
  }
  else{
    //printf("numElecTrue %.16lg\n",stodftInfo->numElecTrueFrag);
    //stodftInfo->numElecTrueFrag = numElecProj*preNe;

    stodftInfo->numElecTrueFrag = numElecProj*preNe-numElecFrag/rhoRealGridTot+
                                  numElecSys;
    printf("Number of Electron for StoWf is %.16lg %.16lg %.16lg\n",
           numElecProj*preNe,numElecFrag/rhoRealGridTot,
           stodftInfo->numElecTrueFrag);
  }
  
  // I would like to seperate them, but this will make life easier
  stodftInfo->numElecTrue = stodftInfo->numElecTrueFrag;
  //debug

/*======================================================================*/
/* VI) Free memories                                                    */

  if(numProcStates>1&&myidState==0)free(rhoTempReduce);
  free(rhoTemp);
  free(numStateUpAllProc);
  if(cpLsda==1&&numStateDn!=0){
    free(numStateDnAllProc);
  }
  if(myidState==0){
    printf("Finish Calculating Density Correction.\n");
    fflush(stdout);
  }
  if(numProcStates>1)Barrier(commStates);
  //fflush(stdout);
  //exit(0);



/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void combineStoUCEnergyWindow(CP *cp,GENERAL_DATA *general_data,CLASS *class,
                 CP *cpMini,GENERAL_DATA *generalDataMini,CLASS *classMini,
                 int ip_now)
/*========================================================================*/
{/*begin routine*/
/*************************************************************************/
/* The second part calculate the density from stochastic calculation     */
/* with energy window.                                                   */
/*************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  STODFTINFO *stodftInfo = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  CPOPTS *cpOpts = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cp->cpcoeffs_info);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  COMMUNICATE *commCP = &(cp->communicate);
  CELL *cell	      = &(general_data->cell);
 
  int iFrag,jFrag,iGrid,jGrid,iState,jState,iStateFrag,iProc;
  int iChem;
  int iAtom;
  int gridIndex;
  int numGrid,numGridSmall;
  int numStateUpMini;
  int numStateDn;
  int numStateDnMini;
  int cpLsda		    = cpOpts->cp_lsda;
  int numFragProc	    = fragInfo->numFragProc;
  int numFragTot	    = fragInfo->numFragTot;
  int occNumber		    = stodftInfo->occNumber;
  int myidState		    = commCP->myid_state;
  int numProcStates	    = commCP->np_states;
  int rhoRealGridNum	    = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot	    = stodftInfo->rhoRealGridTot;
  int numStateUpProc	    = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc	    = cpcoeffs_info->nstate_dn_proc;
  int numStateStoUp     = stodftInfo->numStateStoUp;
  int numStateStoDn     = stodftInfo->numStateStoDn;
  int numCoeff		= cpcoeffs_info->ncoef;
  int countWf;
  int fragOpt		    = stodftInfo->fragOpt;
  int energyWindowOn        = stodftInfo->energyWindowOn;
  int numChemPot                = stodftInfo->numChemPot; // I already +1

  MPI_Comm commStates   =    commCP->comm_states;

  int *numGridFragProc	    = fragInfo->numGridFragProc;
  int *numGridFragProcSmall = fragInfo->numGridFragProcSmall;
  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;
  int *numStateUpAllProc;
  int *numStateDnAllProc;

  int **gridMapProc	    = fragInfo->gridMapProc;
  int **gridMapProcSmall    = fragInfo->gridMapProcSmall;
 
  double proj;
  double numElecProj = 0.0;
  double numElecFrag = 0.0;
  double numElecProjTot = 0.0;
  double numElecFragTot = 0.0;
  double numElecSys = stodftInfo->numElecSys;
  double numGridTotInv = 1.0/rhoRealGridTot;
  double vol,volInv;
  double volMini;
  double pre,preNe;
  double preDot;
  double *rhoTemp,*rhoTempReduce;
  double *wfTemp,*wfFragTemp,*rhoFragTemp;
  double *hmat_cp = cell->hmat_cp;
  double *hmatCpMini;
  double **rhoUpFragProc;
  double **coefUpFragProc;
  double **rhoDnFragProc;
  double **coefDnFragProc;
  double **coefUpFragCoreProc;
  double **coefDnFragCoreProc;
  double *rhoUpFragSum;
  double *rhoDnFragSum;
  double *rhoUpFragSumCpy;
  double *rhoDnFragSumCpy;
  double *noiseWfUpReal,*noiseWfDnReal;
  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;
  
/*======================================================================*/
/* I) Asign local variables				                */

  vol = getdeth(hmat_cp);
  volInv = 1.0/vol;

  // prefactor for prjection part
  pre = (double)(occNumber)*vol*vol/
        ((double)(rhoRealGridTot)*(double)(rhoRealGridTot)*(double)(numStateStoUp));
  // prefactor for number of e in proj part
  preNe = pre/(double)(rhoRealGridTot);
  preDot = 1.0/sqrt(vol);

  rhoUpFragProc = fragInfo->rhoUpFragProc;
  coefUpFragProc = fragInfo->coefUpFragProc;
  coefUpFragCoreProc = fragInfo->coefUpFragCoreProc;
  rhoUpFragSum = fragInfo->rhoUpFragSum;
  rhoUpFragSumCpy = fragInfo->rhoUpFragSumCpy;
  memcpy(rhoUpFragSum,rhoUpFragSumCpy,rhoRealGridNum*sizeof(double));

  if(cpLsda==1&&numStateDn!=0){
    rhoDnFragProc = fragInfo->rhoDnFragProc;
    coefDnFragProc = fragInfo->coefDnFragProc;
    coefDnFragCoreProc = fragInfo->coefDnFragCoreProc;    
    rhoDnFragSum = fragInfo->rhoDnFragSum;
    memcpy(rhoDnFragSum,rhoDnFragSumCpy,rhoRealGridNum*sizeof(double));
  }

  rhoTemp = (double*)cmalloc(rhoRealGridTot*sizeof(double));
  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
    rhoTemp[iGrid] = 0.0;
  }//endfor iGrid

  noiseWfUpReal = fragInfo->noiseWfUpReal;
  if(cpLsda==1&&numStateDn!=0){
    noiseWfDnReal = fragInfo->noiseWfDnReal;
  }

/*======================================================================*/
/* IV) Calculate the real space noise wave function                     */

  // Copy the deterministic wf to the stochastic wf
  //memcpy(&(cp->cpcoeffs_pos[1].cre_up[1]),&(stodftCoefPos->wfDetBackupUpRe[0]),numStateUpProc*numCoeff*sizeof(double));
  //memcpy(&(cp->cpcoeffs_pos[1].cim_up[1]),&(stodftCoefPos->wfDetBackupUpIm[0]),numStateUpProc*numCoeff*sizeof(double));
  
  //noiseRealReGen(general_data,cp,class,ip_now);
  //stodftInfo->iScf = 1;
  
  //noiseFilterGen(general_data,cp,class,ip_now);
  //noiseFilterGenFake(general_data,cp,class,ip_now);

/*======================================================================*/
/* IV) Project the real space noise wave function                       */


  numStateUpAllProc = (int*)cmalloc(numProcStates*sizeof(int));
  if(numProcStates>1){
    Allgather(&numStateUpProc,1,MPI_INT,numStateUpAllProc,1,MPI_INT,0,commStates);
  }
  else numStateUpAllProc[0] = numStateUpProc;
  if(cpLsda==1&&numStateDn!=0){
    numStateDnAllProc = (int*)cmalloc(numProcStates*sizeof(int));
    if(numProcStates>1){
      Allgather(&numStateDnProc,1,MPI_INT,numStateDnAllProc,1,MPI_INT,0,commStates);
    }
    else numStateDnAllProc[0] = numStateDnProc;
  }
  wfTemp = (double*)cmalloc(rhoRealGridTot*sizeof(double));
  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
    rhoTemp[iGrid] = 0.0;
  }

  int res = numFragTot%numProcStates;
  int div = numFragTot/numProcStates;
  int numFragProcMax;
  if(res==0)numFragProcMax = div;
  else numFragProcMax = div+1;

  //debug, output all wf filtered by frag wf
  //double *stowffrag = (double*)cmalloc(numStateUpAllProc[0]*rhoRealGridTot*sizeof(double));
  //double testgrid;
  //DEBUG 
  double **rhoCorrectWindow;
  double *rhoTempTemp = (double*)cmalloc(rhoRealGridTot*sizeof(double));
  double *rhoTempTempReduce = (double*)cmalloc(rhoRealGridTot*sizeof(double));
  FILE *fileout;
  char fileName[100];

  if(myidState==0){
    rhoCorrectWindow = (double**)cmalloc(numChemPot*sizeof(double*));
    for(iChem=0;iChem<numChemPot;iChem++){
      rhoCorrectWindow[iChem] = (double*)calloc(rhoRealGridTot,sizeof(double));
    }
  }

  for(iFrag=0;iFrag<numFragProcMax;iFrag++){
    if(iFrag<numFragProc){
      fragInfo->iFrag = iFrag;
      numGrid = numGridFragProc[iFrag];
      numStateUpMini = cpMini[iFrag].cpcoeffs_info.nstate_up_proc;
      fragInfo->rhoUpFragProc[iFrag] = (double*)cmalloc(numGrid*sizeof(double));
      fragInfo->coefUpFragProc[iFrag] = (double*)cmalloc(numStateUpMini*numGrid*sizeof(double));
      rhoRealCalcDriverFragMol(&generalDataMini[iFrag],&cpMini[iFrag],&classMini[iFrag],cp);
    }//endif iFrag
    //for(iChem=1;iChem<numChemPot;iChem++){ // debug
    for(iChem=0;iChem<numChemPot;iChem++){
      //DEBUG
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
        rhoTempTemp[iGrid] = 0.0;
        rhoTempTempReduce[iGrid] = 0.0;
      }
      noiseFilterRealReGen(class,general_data,cp,stoWfUpRe[iChem],stoWfUpIm[iChem],
                           noiseWfUpReal,numStateUpProc);
      countWf = 0;
      //testgrid = 0.0;
      for(iProc=0;iProc<numProcStates;iProc++){
        for(iState=0;iState<numStateUpAllProc[iProc];iState++){
          if(myidState==iProc){
            memcpy(wfTemp,&noiseWfUpReal[iState*rhoRealGridTot],rhoRealGridTot*sizeof(double));
          }
          if(numProcStates>1)Bcast(wfTemp,rhoRealGridTot,MPI_DOUBLE,iProc,commStates);
          //double *rho_frag_ext = (double*)calloc(rhoRealGridTot,sizeof(double));
          //for(iFrag=0;iFrag<numFragProc;iFrag++){
          if(iFrag<numFragProc){
            numGrid = numGridFragProc[iFrag];
            numGridSmall = numGridFragProcSmall[iFrag];
            numStateUpMini = cpMini[iFrag].cpcoeffs_info.nstate_up_proc;
            hmatCpMini = generalDataMini[iFrag].cell.hmat_cp;
            volMini = getdeth(hmatCpMini);
            volMini /= numGrid;
            wfFragTemp = (double*)calloc(numGrid,sizeof(double));
            rhoFragTemp = (double*)calloc(numGridSmall,sizeof(double));
            for(iGrid=0;iGrid<numGrid;iGrid++){
              gridIndex = gridMapProc[iFrag][iGrid];
              wfFragTemp[iGrid] = wfTemp[gridIndex];
              //wfFragTemp[iGrid] = coefUpFragProc[iFrag][iGrid];
            }
            for(iGrid=0;iGrid<numGridSmall;iGrid++)rhoFragTemp[iGrid] = 0.0;
            //double testsum = 0.0;
            for(iStateFrag=0;iStateFrag<numStateUpMini;iStateFrag++){
              //proj = ddotBlasWrapper(numGrid,&wfFragTemp[0],1,&coefUpFragProc[iFrag][iStateFrag*numGrid],1);
              
              proj = 0.0;
              for(iGrid=0;iGrid<numGrid;iGrid++){
                proj += wfFragTemp[iGrid]*coefUpFragProc[iFrag][iStateFrag*numGrid+iGrid];
              }
              fragInfo->wfProjUp[iFrag][iChem*numStateStoUp*numStateUpMini+countWf*numStateUpMini+iStateFrag] 
                                      = proj*preDot*volMini;
              //printf("ppppproj iChem %i countwf %i iStateFrag %i %lg %lg %lg\n",
              //       iChem,countWf,iStateFrag,
              //       fragInfo->wfProjUp[iFrag][iChem*numStateStoUp*numStateUpMini+countWf*numStateUpMini+iStateFrag],
              //       preDot,volMini);
              //daxpyBlasWrapper(numGrid,proj,&coefUpFragProc[iFrag][iStateFrag*numGrid],1,&rhoFragTemp[0],1);
              for(iGrid=0;iGrid<numGridSmall;iGrid++){
                rhoFragTemp[iGrid] += proj*coefUpFragProc[iFrag][iStateFrag*numGrid+gridMapProcSmall[iFrag][iGrid]];
              }
              //testsum += proj*coefUpFragProc[iFrag][iStateFrag*numGrid+gridMapProcSmall[iFrag][0]];
            }//endfor iStateFrag
            //debug
            /*
            for(iGrid=0;iGrid<numGridSmall;iGrid++){
              gridIndex = gridMapProc[iFrag][gridMapProcSmall[iFrag][iGrid]];
              stowffrag[iState*rhoRealGridTot+gridIndex] = rhoFragTemp[iGrid];
            }
            */
            //printf("iFrag %i testsum %lg\n",iFrag,testsum*testsum);
            for(iGrid=0;iGrid<numGridSmall;iGrid++){
              gridIndex = gridMapProc[iFrag][gridMapProcSmall[iFrag][iGrid]];
              rhoTemp[gridIndex] += rhoFragTemp[iGrid]*rhoFragTemp[iGrid];
              //rhoTemp[gridIndex] = rhoFragTemp[iGrid]*rhoFragTemp[iGrid]; //debug     
              //fix_frag[iFrag][gridIndex] += rhoFragTemp[iGrid]*rhoFragTemp[iGrid]*pre; 
              //DEBUG
              rhoTempTemp[gridIndex] += rhoFragTemp[iGrid]*rhoFragTemp[iGrid];
            }
            //testgrid += rhoFragTemp[0]*rhoFragTemp[0];
            free(wfFragTemp);
            free(rhoFragTemp);
          }//endif iFrag<numFragProc
          //}//endfor iFrag
          countWf += 1;
        }//endfor iState   
      }//endfor iProc
      //printf("1111111111 testgrid %lg\n",testgrid);
      //DEBUG
      Reduce(rhoTempTemp,rhoTempTempReduce,rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
      
      if(myidState==0){
        for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoCorrectWindow[iChem][iGrid] += rhoTempTempReduce[iGrid];
      }
    }//endfor iChem
    if(iFrag<numFragProc){
      free(fragInfo->rhoUpFragProc[iFrag]);
      free(fragInfo->coefUpFragProc[iFrag]);
    }
  }//endfor iFrag

  //DEBUG
  /*
  if(myidState==0){
    for(iChem=0;iChem<numChemPot;iChem++){
      sprintf(fileName,"density-window-correction-%i",iChem);
      fileout = fopen(fileName,"w");
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
        fprintf(fileout,"%.16lg\n",rhoCorrectWindow[iChem][iGrid]*pre);
      }
      fclose(fileout);
    }
  } 
  */
  
  /*
  for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
    printf("gridddddddddd %lg\n",rhoTemp[iGrid]);
  }
  
  fflush(stdout);
  exit(0);
  */
  

  /*--------------------------------------------------*/
  /* We are going to calculate number of electrons.   */
  /* The Formula should be N_e^{frag}-N_e^{proj}+     */
  /* N_e^{sto} = N_e.                                 */
  /*--------------------------------------------------*/
  
  if(numProcStates>1){
    if(myidState==0){
      rhoTempReduce = (double*)cmalloc(rhoRealGridTot*sizeof(double));
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++)rhoTempReduce[iGrid] = 0.0;
    }
    Barrier(commStates);
    Reduce(rhoTemp,rhoTempReduce,rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Scatterv(rhoTempReduce,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	     rhoTemp,rhoRealGridNum,MPI_DOUBLE,0,commStates);
  } 
  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
    numElecProj += rhoTemp[iGrid];
    numElecFrag += rhoUpFragSum[iGrid];
  }
  
  double sumElecFrag = 0.0;
  double sumElecProj = 0.0;
  double sumElecFragTest;
  double sumElecProjTest;
  for(iProc=0;iProc<numProcStates;iProc++){
    if(myidState==iProc){
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
	sumElecFrag += rhoUpFragSum[iGrid];
	sumElecProj += pre*rhoTemp[iGrid];
      }
      sumElecFragTest = sumElecFrag;
      sumElecProjTest = sumElecProj;
      sumElecFrag /= rhoRealGridTot;
      sumElecProj /= rhoRealGridTot;
    }
    if(numProcStates>1)Barrier(commStates);
  }
  //printf("11111111111111 sumElecFrag %.16lg %.16lg sumElecProj %.16lg %.16lg %.16lg\n",sumElecFragTest,sumElecFrag,sumElecProjTest,sumElecProj,pre);
  if(numProcStates>1)Barrier(commStates);
  
  daxpyBlasWrapper(rhoRealGridNum,-pre,&rhoTemp[0],1,&rhoUpFragSum[0],1);
#ifdef TUNE_NOISE
  for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
    rhoUpFragSum[iGrid] *= NOISE_FACTOR;
  }
  sumElecFrag *= NOISE_FACTOR;
  sumElecProj *= NOISE_FACTOR;
#endif
  if(numProcStates>1)Barrier(commStates);
  
  /*
  for(iProc=0;iProc<numProcStates;iProc++){
    if(myidState==iProc){
      for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
        printf("111111111111 rhofix %i %lg %lg\n",iGrid,rhoUpFragSum[iGrid],rhoTemp[iGrid]*pre);
      }
    }
    if(numProcStates>1)Barrier(commStates);
  }
  */
  
  //fflush(stdout);
  //exit(0);

  if(cpLsda==1&&numStateDn!=0){
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      rhoTemp[iGrid] = 0.0;
    }
    countWf = 0;
    for(iChem=0;iChem<numChemPot;iChem++){
      noiseFilterRealReGen(class,general_data,cp,stoWfDnRe[iChem],stoWfDnIm[iChem],
                           noiseWfDnReal,numStateDnProc);
      for(iProc=0;iProc<numProcStates;iProc++){
        for(iState=0;iState<numStateDnAllProc[iProc];iState++){
          if(myidState==iProc){
            memcpy(wfTemp,&noiseWfDnReal[iState*rhoRealGridTot],rhoRealGridTot*sizeof(double));
          }
          if(numProcStates>1)Bcast(wfTemp,rhoRealGridTot,MPI_DOUBLE,iProc,commStates);
          for(iFrag=0;iFrag<numFragProc;iFrag++){
            numGrid = numGridFragProc[iFrag];
            numGridSmall = numGridFragProcSmall[iFrag];
            numStateDnMini = cpMini[iFrag].cpcoeffs_info.nstate_dn_proc;
            hmatCpMini = generalDataMini[iFrag].cell.hmat_cp;
            volMini = getdeth(hmatCpMini);
            volMini /= numGrid;
            wfFragTemp = (double*)cmalloc(numGrid*sizeof(double));
            rhoFragTemp = (double*)cmalloc(numGridSmall*sizeof(double));
            for(iGrid=0;iGrid<numGrid;iGrid++){
              gridIndex = gridMapProc[iFrag][iGrid];
              wfFragTemp[iGrid] = wfTemp[gridIndex];
            }
            for(iGrid=0;iGrid<numGridSmall;iGrid++)rhoFragTemp[iGrid] = 0.0;
            for(iStateFrag=0;iStateFrag<numStateDnMini;iStateFrag++){
              proj = ddotBlasWrapper(numGrid,&wfFragTemp[0],1,&coefDnFragProc[iFrag][iStateFrag*numGrid],1);
              fragInfo->wfProjDn[iFrag][iChem*numStateStoDn*numStateDnMini+countWf*numStateDnMini+iStateFrag] 
                                            = proj*preDot*volMini;
              //daxpyBlasWrapper(numGrid,proj,&coefDnFragProc[iStateFrag*numGrid],1,&rhoFragTemp[0],1);
              for(iGrid=0;iGrid<numGridSmall;iGrid++){
                rhoFragTemp[iGrid] += proj*coefDnFragProc[iFrag][iStateFrag*numGrid+gridMapProcSmall[iFrag][iGrid]];
              } 
            }//endfor iGrid
            for(iGrid=0;iGrid<numGridSmall;iGrid++){
              gridIndex = gridMapProc[iFrag][gridMapProcSmall[iFrag][iGrid]];
              rhoTemp[gridIndex] += rhoFragTemp[iGrid]*rhoFragTemp[iGrid];
            }
            free(wfFragTemp);
            free(rhoFragTemp);
          }//endfor iFrag
          countWf += 1;
        }//endfor iState   
      }//endfor iProc
    }//endfor iChem
    if(numProcStates>1){
      Barrier(commStates);
      Reduce(rhoTemp,rhoTempReduce,rhoRealGridTot,MPI_DOUBLE,MPI_SUM,0,commStates);
      Scatterv(rhoTempReduce,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	       rhoTemp,rhoRealGridNum,MPI_DOUBLE,0,commStates);
    }
    for(iGrid=0;iGrid<rhoRealGridNum;iGrid++){
      numElecProj += rhoTemp[iGrid];
      numElecFrag += rhoDnFragSum[iGrid];
    }
    daxpyBlasWrapper(rhoRealGridNum,-pre,&rhoTemp[0],1,&rhoDnFragSum[0],1);
  }//endif cpLsda

  
  if(numProcStates>1){
    //printf("numElecProj %lg\n",numElecProj*preNe);
    //Allreduce(&numElecProj,&(stodftInfo->numElecTrueFrag),1,MPI_DOUBLE,MPI_SUM,0,commStates);
    //stodftInfo->numElecTrueFrag *= preNe;
   
    Allreduce(&numElecProj,&numElecProjTot,1,MPI_DOUBLE,MPI_SUM,0,commStates);
    Allreduce(&numElecFrag,&numElecFragTot,1,MPI_DOUBLE,MPI_SUM,0,commStates);
    stodftInfo->numElecTrueFrag = numElecProjTot*preNe-numElecFragTot/rhoRealGridTot+
                                  numElecSys;
    if(myidState==0)printf("Number of Electron for StoWf is %.16lg %.16lg %.16lg %lg\n",
                           numElecProjTot*preNe,numElecFragTot/rhoRealGridTot,
                           stodftInfo->numElecTrueFrag,numElecSys);
  }
  else{
    //printf("numElecTrue %.16lg\n",stodftInfo->numElecTrueFrag);
    //stodftInfo->numElecTrueFrag = numElecProj*preNe;

    stodftInfo->numElecTrueFrag = numElecProj*preNe-numElecFrag/rhoRealGridTot+
                                  numElecSys;
    printf("Number of Electron for StoWf is %.16lg %.16lg %.16lg %.16lg\n",
           numElecProj*preNe,numElecFrag/rhoRealGridTot,
           stodftInfo->numElecTrueFrag,numElecSys);
  }
  
  // I would like to seperate them, but this will make life easier
  stodftInfo->numElecTrue = stodftInfo->numElecTrueFrag;
  //debug

/*======================================================================*/
/* VI) Free memories                                                    */

  if(numProcStates>1&&myidState==0)free(rhoTempReduce);
  free(rhoTemp);
  free(numStateUpAllProc);
  if(cpLsda==1&&numStateDn!=0){
    free(numStateDnAllProc);
  }
  if(myidState==0){
    printf("Finish Calculating Density Correction.\n");
    fflush(stdout);
  }
  if(numProcStates>1)Barrier(commStates);
  //fflush(stdout);
  //exit(0);



/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/



