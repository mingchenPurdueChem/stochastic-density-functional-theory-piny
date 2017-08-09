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
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void projRhoMini(CP *cp,GENERAL_DATA *general_data,CLASS *class,
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
 
  int iFrag,iGrid,iState,jState,iStateFrag,iProc;
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
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    numGrid = numGridFragProc[iFrag];
    for(iGrid=0;iGrid<numGrid;iGrid++){
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
      printf("111111 rhofrag %.8lg\n",rhoUpFragSum[iGrid]);
    }
  }
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

  rhoRealCalcDriverNoise(general_data,cp,class,ip_now);

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
  for(iProc=0;iProc<numProcStates;iProc++){
    for(iState=0;iState<numStateUpAllProc[iProc];iState++){
      if(myidState==iProc){
	memcpy(wfTemp,&noiseWfUpReal[iState*rhoRealGridTot],rhoRealGridTot*sizeof(double));
      }
      if(numProcStates>1)Bcast(wfTemp,rhoRealGridTot,MPI_DOUBLE,iProc,commStates);
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
	  proj = ddotBlasWrapper(numGrid,&wfFragTemp[0],1,&coefUpFragProc[iFrag][iStateFrag*numGrid],1);
	  fragInfo->wfProjUp[iFrag][countWf*numStateUpMini+iStateFrag] = proj*preDot*volMini;
	  daxpyBlasWrapper(numGrid,proj,&coefUpFragProc[iFrag][iStateFrag*numGrid],1,&rhoFragTemp[0],1);
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

  //debug
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
	    daxpyBlasWrapper(numGrid,proj,&coefDnFragProc[iStateFrag*numGrid],1,&rhoFragTemp[0],1);
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

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/



