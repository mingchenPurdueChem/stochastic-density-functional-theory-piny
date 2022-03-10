/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: filter-diag.c                                  */
/*                                                                          */
/* This routine constructs filters of Fermi function.                       */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_stodft_local.h"

#include "../typ_defs/typ_mask.h"
#define MKL_THREADS
// scatterv has problem on bell. This is a replaced version
#define TEST_SCATTERV 

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void orthDiagDriver(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                          int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */ 
  CLATOMS_POS *clatoms_pos	  = &(class->clatoms_pos[ip_now]);
  CPCOEFFS_POS  *cpcoeffs_pos     = &(cp->cpcoeffs_pos[ip_now]);
  STODFTINFO    *stodftInfo       = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos    = cp->stodftCoefPos;

  orthNormStoWf(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
  buildKSMatrix(cp,class,general_data,cpcoeffs_pos,clatoms_pos);
  diagKSMatrix(cp,class,general_data,cpcoeffs_pos,clatoms_pos);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void orthNormStoWf(CP *cp,CLASS *class,GENERAL_DATA *general_data,
		    CPCOEFFS_POS  *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  CPCOEFFS_INFO *cpcoeffs_info    = &(cp->cpcoeffs_info);
  STODFTINFO    *stodftInfo       = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos    = cp->stodftCoefPos;
  COMMUNICATE *communicate      = &(cp->communicate);

  int iChem,iCoeff,iState,jState,indCoeff1,indCoeff2;
  int iProc;
  int indexCut;
  int iCoeff1,iCoeff2;
  int numChemPot     = stodftInfo->numChemPot;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int numStateUpTotal = stodftInfo->numStateStoUp;
  int numStateUpAllProc = numChemPot*numStateUpTotal;
  int numCoeffUpAllProc = numChemPot*numStateUpTotal*numCoeff;
  int numStateStoUp = stodftInfo->numStateStoUp;
  int numCeoffUpAllProcOneChem = numStateStoUp*numCoeff;
  int numProcStates = communicate->np_states;
  int myidState = communicate->myid_state;
  int numThreadsMKL = communicate->numThreadsMKL;
  MPI_Comm comm_states   =    communicate->comm_states;

  int *stowfRecvCounts = stodftInfo->stowfRecvCounts;
  int *stowfDispls = stodftInfo->stowfDispls;
  int *stowfRecvCountsComplex = stodftInfo->stowfRecvCountsComplex;
  int *stowfDisplsComplex = stodftInfo->stowfDisplsComplex;
  int *stowfRecvCountsComplex2 = stodftInfo->stowfRecvCountsComplex2;
  int *stowfDisplsComplex2 = stodftInfo->stowfDisplsComplex2;
  int *numStates2 = stodftInfo->numStates2;
  int *dsplStates2 = stodftInfo->dsplStates2;

  
  double pre = sqrt(2.0);
  double sum;
  //double *wfBfOrthUp = stodftCoefPos->wfBfOrthUp;  
  double *wfBfOrthUp;
  double *wfOrthUpRe = stodftCoefPos->wfOrthUpRe;
  double *wfOrthUpIm = stodftCoefPos->wfOrthUpIm;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  double *allWF = (double *)cmalloc(2*numCoeffUpTotal*
						    sizeof(double));

  //for debug
  int *numStates = stodftInfo->numStates;
  int *dsplStates = stodftInfo->dsplStates;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  //double *dot1 = (double*)cmalloc(numChemPot*numStateUpProc*sizeof(double));
  double *dotAll;

  int numStatesDet = stodftInfo->numStatesDet;
  double *wfUpReDet = stodftCoefPos->wfUpReDet;
  double *wfUpImDet = stodftCoefPos->wfUpImDet;

  double time_st,time_end;


/*==========================================================================*/
/* I) Pack all stochastic wave functions */

  /*
  if(myidState==0){
    for(iProc=0;iProc<numProcStates;iProc++){
      printf("%i recvCounts %i displs %i\n",iProc,stowfRecvCountsComplex[iProc],stowfDisplsComplex[iProc]);
    }    
  }
  */
  if(myidState==0)printf("Start Gathering WF\n");
  //DEBUG
  for(iChem=numChemPot-1;iChem>0;iChem--){
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      stoWfUpRe[iChem][iCoeff] -= stoWfUpRe[iChem-1][iCoeff];
      stoWfUpIm[iChem][iCoeff] -= stoWfUpIm[iChem-1][iCoeff];
    }//endfor iCoeff
  }//endfor iChem
  

  /*
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iState=0;iState<numStateUpProc;iState++){
      sum = 0.0;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	iCoeff1 = iState*numCoeff+iCoeff;
	sum += stoWfUpRe[iChem][iCoeff1]*stoWfUpRe[iChem][iCoeff1]+
		stoWfUpIm[iChem][iCoeff1]*stoWfUpIm[iChem][iCoeff1];
      }
      sum *= 2.0;
      iCoeff1 = iState*numCoeff+numCoeff;
      sum += stoWfUpRe[iChem][iCoeff1]*stoWfUpRe[iChem][iCoeff1];
      sum = sqrt(sum);
      for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
	iCoeff1 = iState*numCoeff+iCoeff;
	stoWfUpRe[iChem][iCoeff1] /= sum;
	stoWfUpIm[iChem][iCoeff1] /= sum;
      }
      printf("1111111111 iChem %i iState %i wflnth %lg\n",iChem,iState,sum);
    }
    
    calcForceWrapper(cp,class,general_data,cpcoeffs_pos,clatoms_pos,
                    stoWfUpRe[iChem],stoWfUpIm[iChem]);   
    for(iState=0;iState<numStateUpProc;iState++){
      sum = 0.0;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	iCoeff1 = iState*numCoeff+iCoeff;
	sum += fcre_up[iCoeff1]*fcre_up[iCoeff1]+
		fcim_up[iCoeff1]*fcim_up[iCoeff1];
      }
      sum *= 2.0;
      iCoeff1 = iState*numCoeff+numCoeff;
      sum += fcre_up[iCoeff1]*fcre_up[iCoeff1];
      dot1[iChem*numStateUpProc+iState] = sum;
      printf("iState %i dotttt %lg\n",iState,sum);
    }
  }
  if(myidState==0){
    dotAll = (double*)cmalloc(numStateUpAllProc*sizeof(double));
    for(iProc=0;iProc<numProcStates;iProc++){
      printf("iProc %i numStates %i dsplStates %i\n",iProc,numStates[iProc],dsplStates[iProc]);
    }
  }

  if(numProcStates>1){
    Barrier(comm_states);
    Gatherv(dot1,numChemPot*numStateUpProc,MPI_DOUBLE,dotAll,numStates,dsplStates,
	    MPI_DOUBLE,0,comm_states);
    Barrier(comm_states);
  }
  else{
    memcpy(&dotAll[0],&dot1[0],numStateUpAllProc*sizeof(double));
  }
  if(myidState==0){
    printf("before orth\n");
    for(iState=0;iState<numStateUpAllProc;iState++){
      printf("bo iState %i Eave %.8lg\n",iState,dotAll[iState]);
    }
  }
  */
  
  
  int numCoeffUpAllProc2 = numCoeffUpAllProc*2;
  if(myidState==0){
    wfBfOrthUp = (double*)calloc(numCoeffUpAllProc*2,sizeof(double));
  }
  //DEBUG
  
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iState=0;iState<numStateUpProc;iState++){
      for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
	indCoeff2 = iState*numCoeff+iCoeff;
	//indCoeff1 = indCoeff2+iChem*numCoeffUpTotal;
	allWF[2*indCoeff2] = pre*stoWfUpRe[iChem][indCoeff2+1];
	allWF[2*indCoeff2+1] = pre*stoWfUpIm[iChem][indCoeff2+1];
	//allWF[indCoeff1] = pre*stoWfUpRe[iChem][indCoeff2+1]+
	//		    I*pre*stoWfUpIm[iChem][indCoeff2+1];
      }//endfor iCoeff
      indCoeff2 = iState*numCoeff+numCoeff-1;
      //indCoeff1 = indCoeff2+iChem*numCoeffUpTotal-1;
      //allWF[indCoeff1] = stoWfUpRe[iChem][indCoeff2]+I*stoWfUpIm[iChem][indCoeff2];
      allWF[2*indCoeff2] = stoWfUpRe[iChem][indCoeff2+1];
      allWF[2*indCoeff2+1] = 0.0;
    }//endfor iState

    // Gather to master proc
    if(numProcStates>1){
      Gatherv(allWF,2*numCoeffUpTotal,MPI_DOUBLE,
              &wfBfOrthUp[iChem*numCeoffUpAllProcOneChem*2],
              stowfRecvCountsComplex,stowfDisplsComplex,
              MPI_DOUBLE,0,comm_states);
      //printf("iChem %i allWF[0] %lg\n",iChem,allWF[0]);
    }
    else{
      memcpy(&wfBfOrthUp[iChem*numCeoffUpAllProcOneChem*2],&allWF[0],
             2*numCeoffUpAllProcOneChem*sizeof(double));
    }
    free(&stoWfUpRe[iChem][1]);
    free(&stoWfUpIm[iChem][1]);
  }//endfor iChem
  

  /*
  if(myidState==0){
    for(iChem=0;iChem<numChemPot;iChem++){
      printf("iiiiiiChem %i wfBfOrthUp %lg\n",iChem,wfBfOrthUp[23760*iChem]);
    }
  }
  */
  
/*==========================================================================*/
/* II) Gather to master proc */

  /*
  if(numProcStates>1){
    Gatherv(allWF,2*numChemPot*numCoeffUpTotal,MPI_DOUBLE,
  	    wfBfOrthUp,stowfRecvCountsComplex,stowfDisplsComplex,MPI_DOUBLE,
	    0,comm_states);
  }
  else{
    memcpy(&wfBfOrthUp[0],&allWF[0],2*numCoeffUpAllProc*sizeof(double));
  }
  */
 
/*==========================================================================*/
/* III) SVD */
  
  if(myidState==0){
    time_st = omp_get_wtime();
    printf("Start SVD\n");
    //DEBUG    
    indexCut = orthogSVD(numStateUpAllProc,2*numCoeff,wfBfOrthUp,numThreadsMKL);
    //indexCut = 1250;
    //indexCut = orthogSVD(numCoeff,numStateUpAllProc,wfBfOrthUp);
    printf("indexCut %i\n",indexCut);
    if(indexCut<stodftInfo->numStatePrintUp){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("You require %i of orbitals. However, I can only find %i\n",
              stodftInfo->numStatePrintUp,indexCut);
      printf("independent orbitals. Please increase the number of windows.\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }
    time_end = omp_get_wtime();
    printf("svd time %lg\n",time_end-time_st);
  }
  if(numProcStates>1)Barrier(comm_states);
  if(numProcStates>1)Bcast(&indexCut,1,MPI_INT,0,comm_states);
  // I need do something for indexCut
  // Save wave functios on master proc
  stodftInfo->numStateUpIdp = indexCut;
  if(numProcStates>1)Barrier(comm_states);
  //exit(0);
/*==========================================================================*/
/* IV) Scatter back to all proc */

  stodftInfo->numStateUpProcBackup = numStateUpProc;
  stodftInfo->numStateDnProcBackup = numStateDnProc;

  int div,res;
  div = (int)(indexCut/numProcStates);
  res = indexCut%numProcStates;
  if(myidState<res)cpcoeffs_info->nstate_up_proc = div+1;
  else cpcoeffs_info->nstate_up_proc = div;
  numStateUpProc = cpcoeffs_info->nstate_up_proc;
  numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  numCoeffUpTotal = numStateUpProc*numCoeff;
  numCoeffDnTotal = numStateDnProc*numCoeff;
  numStateUpTotal = stodftInfo->numStateStoUp;
  numStateUpAllProc = numChemPot*numStateUpTotal;
  numCoeffUpAllProc = numChemPot*numStateUpTotal*numCoeff;

  if(numProcStates>1){
    Allgather(&numCoeffUpTotal,1,MPI_INT,stowfRecvCountsComplex2,1,MPI_INT,0,comm_states);
    Allgather(&numCoeffUpTotal,1,MPI_INT,stowfRecvCounts,1,MPI_INT,0,comm_states);
    Allgather(&numStateUpProc,1,MPI_INT,numStates2,1,MPI_INT,0,comm_states);    
    //Allgather(&numCoeffUpAllChemPot,1,MPI_INT,stowfRecvCounts,1,MPI_INT,0,comm_states);
    for(iProc=0;iProc<numProcStates;iProc++){
      stowfRecvCountsComplex2[iProc] *= 2;
    }
    stowfDisplsComplex2[0] = 0;
    dsplStates2[0] = 0;
    stowfDispls[0] = 0;
    for(iProc=1;iProc<numProcStates;iProc++){
      stowfDisplsComplex2[iProc] = stowfRecvCountsComplex2[iProc-1]+stowfDisplsComplex2[iProc-1];
      dsplStates2[iProc] = numStates2[iProc-1]+dsplStates2[iProc-1];
      stowfDispls[iProc] = stowfDispls[iProc-1]+stowfRecvCounts[iProc-1];
    }
  }

  /*
  for(iProc=0;iProc<numProcStates;iProc++){
    printf("count %i %i %i\n",stowfRecvCountsComplex2[iProc],stowfRecvCounts[iProc],numStates2[iProc]);
  }
  for(iProc=0;iProc<numProcStates;iProc++){
    printf("dspl %i %i %i\n",stowfDisplsComplex2[iProc],stowfDispls[iProc],dsplStates2[iProc]);
  }
  */


  //Gather(&numCoeffUpAllProc,1,MPI_INT,stowfRecvCounts,numProcStates,MPI_INT,0,comm_states); 
  //Bcast(stowfRecvCounts,numProcStates,MPI_INT,0,comm_states);
  Barrier(comm_states);
  free(&(cpcoeffs_pos->cre_up[1]));
  free(&(cpcoeffs_pos->cim_up[1]));
  free(&(cpcoeffs_pos->fcre_up[1]));
  free(&(cpcoeffs_pos->fcim_up[1]));

  cpcoeffs_pos->cre_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  cpcoeffs_pos->cim_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  cpcoeffs_pos->fcre_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  cpcoeffs_pos->fcim_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;

  free(allWF);
  allWF = (double*)calloc(2*numCoeffUpTotal,sizeof(double));

  if(myidState==0)printf("Start Scatter Data\n");

  if(numProcStates>1){
    for(iProc=0;iProc<numProcStates;iProc++){
      printf("stowfRecvCountsComplex2 %i stowfDisplsComplex2 %i %i\n",stowfRecvCountsComplex2[iProc],stowfDisplsComplex2[iProc],2*numCoeffUpTotal);
    }
    printf("2*numCoeffUpTotal %i numCoeffUpAllProc*2 %i %i\n",2*numCoeffUpTotal,numCoeffUpAllProc2,numCoeffUpAllProc*2);
#ifdef TEST_SCATTERV
    if(myidState!=0){
      Recv(allWF,2*numCoeffUpTotal,MPI_DOUBLE,0,myidState,comm_states);
    }
    else{
      for(iProc=1;iProc<numProcStates;iProc++){
        Send(&wfBfOrthUp[stowfDisplsComplex2[iProc]],stowfRecvCountsComplex2[iProc],
             MPI_DOUBLE,iProc,iProc,comm_states);
      }
      memcpy(allWF,wfBfOrthUp,2*numCoeffUpTotal*sizeof(double));
    }
#else
    Scatterv(wfBfOrthUp,stowfRecvCountsComplex2,stowfDisplsComplex2,MPI_DOUBLE,
  	      allWF,2*numCoeffUpTotal,MPI_DOUBLE,0,comm_states);
#endif
  }
  else{
    memcpy(&allWF[0],&wfBfOrthUp[0],2*numCoeffUpTotal*sizeof(double));
  }

  if(myidState==0)free(wfBfOrthUp);

/*==========================================================================*/
/* IV) Unpack all stochastic wave functions */

  /*
  if(numProcStates>1)Barrier(comm_states);
  //printf("111111111111\n");
  for(iChem=0;iChem<numChemPot;iChem++){
    //printf("myid %i iChem %i\n",myidState,iChem);
    for(iState=0;iState<numStateUpProc;iState++){
      for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
        indCoeff2 = iState*numCoeff+iCoeff;
        indCoeff1 = indCoeff2+iChem*numCoeffUpTotal;
	stoWfUpRe[iChem][indCoeff2+1] = allWF[2*indCoeff1];
        stoWfUpIm[iChem][indCoeff2+1] = allWF[2*indCoeff1+1];
	//stoWfUpRe[iChem][indCoeff2+1] = creal(allWF[indCoeff1])/pre;
	//stoWfUpIm[iChem][indCoeff2+1] = cimag(allWF[indCoeff1])/pre;
      }//endfor iCoeff
      indCoeff2 = iState*numCoeff+numCoeff-1;
      indCoeff1 = indCoeff2+iChem*numCoeffUpTotal;
      stoWfUpRe[iChem][indCoeff2+1] = allWF[2*indCoeff1];
      stoWfUpIm[iChem][indCoeff2+1] = 0.0;
    }//endfor iState
  }//endfor iChem
  */
  if(numProcStates>1)Barrier(comm_states);
  for(iState=0;iState<numStateUpProc;iState++){
    for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
      indCoeff2 = iState*numCoeff+iCoeff;
      cre_up[indCoeff2+1] = allWF[2*indCoeff2]/pre;
      cim_up[indCoeff2+1] = allWF[2*indCoeff2+1]/pre;
    }
    indCoeff2 = iState*numCoeff+numCoeff-1;
    cre_up[indCoeff2+1] = allWF[2*indCoeff2];
    cim_up[indCoeff2+1] = 0.0;
  }

/*==========================================================================*/
/* V) Normalization all stochastic wave functions */
  
  if(myidState==0)printf("Start Normalization\n");  

  /*
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iState=0;iState<numStateUpProc;iState++){
      sum = 0.0;
      iCoeff2 = iState*numCoeff;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	sum += stoWfUpRe[iChem][iCoeff2+iCoeff]*stoWfUpRe[iChem][iCoeff2+iCoeff]+
		stoWfUpIm[iChem][iCoeff2+iCoeff]*stoWfUpIm[iChem][iCoeff2+iCoeff];
      }//endfor iCoeff
      sum *= 2.0;
      sum += stoWfUpRe[iChem][iCoeff2+numCoeff]*stoWfUpRe[iChem][iCoeff2+numCoeff];
      sum = sqrt(sum);
      for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
	stoWfUpRe[iChem][iCoeff2+iCoeff] /= sum;
	stoWfUpIm[iChem][iCoeff2+iCoeff] /= sum;
      }//endfor iCoeff
    }//endfor iState
  }//endfor iChem
  */
  
  for(iState=0;iState<numStateUpProc;iState++){
    sum = 0.0;
    iCoeff2 = iState*numCoeff;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      sum += cre_up[iCoeff2+iCoeff]*cre_up[iCoeff2+iCoeff]+
             cim_up[iCoeff2+iCoeff]*cim_up[iCoeff2+iCoeff];
    }//endfor iCoeff
    sum *= 2.0;
    sum += cre_up[iCoeff2+numCoeff]*cre_up[iCoeff2+numCoeff];
    sum = sqrt(sum);
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      cre_up[iCoeff2+iCoeff] /= sum;
      cim_up[iCoeff2+iCoeff] /= sum;
    }//endfor iCoeff
  }//endfor iState

  
  /*
  for(iState=0;iState<numStateUpProc;iState++){
    indCoeff1 = iState*numCoeff;
    for(jState=0;jState<numStateUpProc;jState++){
      indCoeff2 = jState*numCoeff;
      sum = 0.0;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
        sum += cre_up[indCoeff1+iCoeff]*cre_up[indCoeff2+iCoeff]+
               cim_up[indCoeff1+iCoeff]*cim_up[indCoeff2+iCoeff];
      }
      sum *= 2.0;
      sum += cre_up[indCoeff1+numCoeff]*cre_up[indCoeff2+numCoeff];
      printf("%lg ",sum);
    }
    printf("\n");
  }
  */
  

  //debug 
  /*
  for(iChem=0;iChem<numChemPot;iChem++){
    calcForceWrapper(cp,class,general_data,cpcoeffs_pos,clatoms_pos,
                    stoWfUpRe[iChem],stoWfUpIm[iChem]);
    for(iState=0;iState<numStateUpProc;iState++){
      sum = 0.0;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
        iCoeff1 = iState*numCoeff+iCoeff;
        sum += stoWfUpRe[iChem][iCoeff1]*fcre_up[iCoeff1]+
                stoWfUpIm[iChem][iCoeff1]*fcim_up[iCoeff1];
      }
      sum *= 2.0;
      iCoeff1 = iState*numCoeff+numCoeff;
      sum += stoWfUpRe[iChem][iCoeff1]*fcre_up[iCoeff1];
      dot1[iChem*numStateUpProc+iState] = sum;
    }
  }
  if(myidState==0){
    dotAll = (double*)cmalloc(numStateUpAllProc*sizeof(double));
    for(iProc=0;iProc<numProcStates;iProc++){
      printf("iProc %i numStates %i dsplStates %i\n",iProc,numStates[iProc],dsplStates[iProc]);
    }
  }

  Barrier(comm_states);
  Gatherv(dot1,numChemPot*numStateUpProc,MPI_DOUBLE,dotAll,numStates,dsplStates,
          MPI_DOUBLE,0,comm_states);
  Barrier(comm_states);
  if(myidState==0){
    printf("before orth\n");
    for(iState=0;iState<numStateUpAllProc;iState++){
      printf("ao iState %i Eave %.8lg\n",iState,dotAll[iState]);
    }
  }
  */

  free(allWF);
  //exit(0);
   
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
int orthogSVD(int m,int n,double *matrix,int nthreads)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  // m is wf number, n is PW number
  int i,j,k;
  char jobu = 'O';
  char jobvt = 'N';
  int lda = n;
  double *S = (double*)cmalloc(m*sizeof(double));
  int ldu = 1;
  int ldvt = 1;
  int lwork = 5*(m+n);
  //int lwork = 3*(m+n);
  //MKL_Complex16 *work = (MKL_Complex16 *)cmalloc(lwork*sizeof(MKL_Complex16));
  double *work = (double*)cmalloc(lwork*sizeof(double));
  double *rwork = (double*)cmalloc(5*m*sizeof(double));
  int info;
  int indexCut;
  
/*==========================================================================*/
/* I) SVD */

  /*
  double complex *matrixComplex = (double complex*)cmalloc(m*n*sizeof(double complex));
  for(i=0;i<m*n;i++){
    matrixComplex[i] = matrix[2*i]+I*matrix[2*i+1];
  }
  */

  //zgesvd_(&jobu,&jobvt,&n,&m,&matrixComplex[0],&lda,&S[0],NULL,&ldu,NULL,&ldvt,
  //	    &work[0],&lwork,&rwork[0],&info);
#ifdef MKL_THREADS
  printf("nthreads %i\n",nthreads);
  mkl_set_dynamic(0);
  mkl_set_num_threads(nthreads);
  //omp_set_nested(1);
  omp_set_max_active_levels(1);
#endif
  {
  dgesvd_(&jobu,&jobvt,&n,&m,&matrix[0],&lda,&S[0],NULL,&ldu,NULL,&ldvt,
	  &work[0],&lwork,&info);
  }
  if(info!=0){
    printf("Bad SVD results %i\n",info);
    fflush(stdout);
    exit(0);
  }
/*==========================================================================*/
/* I) Get Cutoff index if singular value is too small */

  indexCut = m;
  for(i=0;i<m;i++){
    //printf("i %i SV %lg\n",i,S[i]);
    if(S[i]/S[0]<1.0e-10){indexCut = i;break;}//endif
  }//endfor i

  /*
  for(i=0;i<m*n;i++){
    //matrix[2*i] = creal(matrixComplex[i]);
    //matrix[2*i+1] = cimag(matrixComplex[i]);
    printf("SVDDDDD %.8lg %.8lg\n",matrix[2*i],matrix[2*i+1]);
  }
  */

  free(S);
  free(work);
  free(rwork);
  //free(matrixComplex);
  
  return indexCut;

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void buildKSMatrix(CP *cp,CLASS *class,GENERAL_DATA *general_data,
		   CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
 
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  COMMUNICATE *communicate      = &(cp->communicate);

  int iChem,iCoeff,iState,jState,iOff;
  int iProc;
  int index1,index2;
  int indexState1,indexState2;
  int numChemPot     = stodftInfo->numChemPot;  
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numStateUpTotal = stodftInfo->numStateStoUp;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int numCoeffUpAllProc = numChemPot*numStateUpTotal*numCoeff;
  int numStateUpAllProc = numChemPot*numStateUpTotal;
  int numStateUpIdp = stodftInfo->numStateUpIdp;
  int numProcStates = communicate->np_states;
  int myidState = communicate->myid_state;
  MPI_Comm comm_states   =    communicate->comm_states;

  int *stowfRecvCounts = stodftInfo->stowfRecvCounts;
  int *stowfDispls = stodftInfo->stowfDispls;
  int *numStates2 = stodftInfo->numStates2;
  int *dsplStates2 = stodftInfo->dsplStates2;

  double sum;

  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  double *forceUpRe = NULL;
  double *forceUpIm = NULL;
  double *wfOrthUpRe = stodftCoefPos->wfOrthUpRe;
  double *wfOrthUpIm = stodftCoefPos->wfOrthUpIm;
  double *KSMatrix = stodftCoefPos->KSMatrix;
  double *KSMatrixReduce,*KSMatrixReduceAll;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

/*==========================================================================*/
/* I)  Calculate force	*/
 
  /*   
  if(myidState==0){
    allForceRe = (double*)malloc(numCoeffUpAllProc*sizeof(double))-1;
    allForceIm = (double*)malloc(numCoeffUpAllProc*sizeof(double))-1;
  }
  //allForceRe = (double*)malloc(numCoeffUpAllProc*sizeof(double))-1;
  //allForceIm = (double*)malloc(numCoeffUpAllProc*sizeof(double))-1;
  //Barrier(comm_states);
  forceUpRe = (double*)cmalloc(numChemPot*numCoeffUpTotal*sizeof(double))-1;
  forceUpIm = (double*)cmalloc(numChemPot*numCoeffUpTotal*sizeof(double))-1;
  */

  /*
  if(numProcStates>1)Barrier(comm_states);
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
      cre_up[iCoeff] = stoWfUpRe[iChem][iCoeff];
      cim_up[iCoeff] = stoWfUpIm[iChem][iCoeff];
      // In case I forget to zero them
      fcre_up[iCoeff] = 0.0;
      fcim_up[iCoeff] = 0.0;
    }
    if(numProcStates>1)Barrier(comm_states);
    calcCoefForceWrapSCF(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    //calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
    //calcCoefForceWrapReduce(class,general_data,cp,cpcoeffs_pos,clatoms_pos);    
    for(iState=0;iState<numStateUpProc;iState++){
      iOff = iState*numCoeff;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	fcre_up[iOff+iCoeff] *= -0.25;
	fcim_up[iOff+iCoeff] *= -0.25;
      }//endfor iCoeff
      fcre_up[iOff+numCoeff] *= -0.5;
      fcim_up[iOff+numCoeff] *= -0.5;
    }//endfor iState
    memcpy(&forceUpRe[iChem*numCoeffUpTotal+1],&fcre_up[1],numCoeffUpTotal*sizeof(double));
    memcpy(&forceUpIm[iChem*numCoeffUpTotal+1],&fcim_up[1],numCoeffUpTotal*sizeof(double));
  }//endfor iChem
  */
  for(iState=0;iState<numStateUpProc;iState++){
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      fcre_up[iState*numCoeff+iCoeff] = 0.0;
      fcim_up[iState*numCoeff+iCoeff] = 0.0;
    }
  }  

  calcCoefForceWrapSCF(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  for(iState=0;iState<numStateUpProc;iState++){
    iOff = iState*numCoeff;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      fcre_up[iOff+iCoeff] *= -0.25;
      fcim_up[iOff+iCoeff] *= -0.25;
    }
    fcre_up[iOff+numCoeff] *= -0.5;
    fcim_up[iOff+numCoeff] *= -0.5;
  }

/*==========================================================================*/
/* II)  Gather all forces to the master proc  */

  if(myidState==0)printf("Start gather force\n");

  /*  
  if(numProcStates>1){
    Gatherv(&forceUpRe[1],numChemPot*numCoeffUpTotal,MPI_DOUBLE,
            &allForceRe[1],stowfRecvCounts,stowfDispls,MPI_DOUBLE,
	    0,comm_states);
    Gatherv(&forceUpIm[1],numChemPot*numCoeffUpTotal,MPI_DOUBLE,
            &allForceIm[1],stowfRecvCounts,stowfDispls,MPI_DOUBLE,
	    0,comm_states);
  }
  else{
    memcpy(&allForceRe[1],&forceUpRe[1],numCoeffUpAllProc*sizeof(double));
    memcpy(&allForceIm[1],&forceUpIm[1],numCoeffUpAllProc*sizeof(double));
  }
  */

/*==========================================================================*/
/* III)  Construct KS matrix  */

  stodftCoefPos->KSMatrix = (double*)crealloc(stodftCoefPos->KSMatrix,
                                numStateUpIdp*numStateUpIdp*sizeof(double));
  KSMatrix = stodftCoefPos->KSMatrix;
  for(iState=0;iState<numStateUpIdp*numStateUpIdp;iState++){
    KSMatrix[iState] = 0.0;
  }
  if(myidState==0&&numProcStates>1){
    KSMatrixReduce = (double*)calloc(numStateUpIdp*numStateUpIdp,sizeof(double));
    KSMatrixReduceAll = (double*)calloc(numStateUpIdp*numStateUpIdp,sizeof(double));
  }

  printf("aaaaa numCoeffUpTotal %i numStates2[iProc] %i\n",numCoeffUpTotal,numStates2[myidState]*numCoeff);
  for(iProc=0;iProc<numProcStates;iProc++){
    if(numProcStates>1)Barrier(comm_states);
    forceUpRe = (double*)crealloc(forceUpRe,numStates2[iProc]*numCoeff*sizeof(double)); 
    forceUpIm = (double*)crealloc(forceUpIm,numStates2[iProc]*numCoeff*sizeof(double));
    for(iCoeff=0;iCoeff<numStates2[iProc]*numCoeff;iCoeff++){
      forceUpRe[iCoeff] = 0.0;
      forceUpIm[iCoeff] = 0.0;
    }
    for(iState=0;iState<numStateUpIdp*numStateUpIdp;iState++)KSMatrix[iState] = 0.0;
    if(myidState==0){
      for(iState=0;iState<numStateUpIdp*numStateUpIdp;iState++)KSMatrixReduce[iState] = 0.0;
    }
    if(numProcStates>1)Barrier(comm_states);
    if(myidState==iProc){
      memcpy(forceUpRe,&fcre_up[1],numCoeffUpTotal*sizeof(double));
      memcpy(forceUpIm,&fcim_up[1],numCoeffUpTotal*sizeof(double));
    }
    if(numProcStates>1){
      Barrier(comm_states);
      Bcast(forceUpRe,numStates2[iProc]*numCoeff,MPI_DOUBLE,iProc,comm_states);
      Bcast(forceUpIm,numStates2[iProc]*numCoeff,MPI_DOUBLE,iProc,comm_states);
      Barrier(comm_states);
    }
    for(iState=0;iState<numStateUpProc;iState++){
      indexState1 = dsplStates2[myidState]+iState;
      for(jState=0;jState<numStates2[iProc];jState++){
        indexState2 = dsplStates2[iProc]+jState;
        sum = 0.0;
        for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
          index1 = iState*numCoeff+iCoeff;
          index2 = jState*numCoeff+iCoeff;
          sum += cre_up[index1+1]*forceUpRe[index2]
                 +cim_up[index1+1]*forceUpIm[index2];
        }//endfor iCoeff
        sum *= 2.0;
        index1 = iState*numCoeff+numCoeff-1;
        index2 = jState*numCoeff+numCoeff-1;
        sum += cre_up[index1+1]*forceUpRe[index2];
        KSMatrix[indexState1*numStateUpIdp+indexState2] = sum;
        //KSMatrix[indexState2*numStateUpIdp+indexState1] = sum;
      }//endfor jState
    }//endfor iState
    if(numProcStates>1){
      Barrier(comm_states);
      Reduce(KSMatrix,KSMatrixReduce,numStateUpIdp*numStateUpIdp,
             MPI_DOUBLE,MPI_SUM,0,comm_states);
      if(myidState==0){
        for(iState=0;iState<numStateUpIdp*numStateUpIdp;iState++){
          KSMatrixReduceAll[iState] += KSMatrixReduce[iState];
        }
      }
      Barrier(comm_states);
    }
  }//endfor iProc

  // debug
  /*
  numCoeffUpAllProc = numStateUpIdp*numCoeff;
  double *forceUpReAll = (double*)cmalloc(numCoeffUpAllProc*sizeof(double));
  double *forceUpImAll = (double*)cmalloc(numCoeffUpAllProc*sizeof(double));
  double *creAll = (double*)cmalloc(numCoeffUpAllProc*sizeof(double));
  double *cimAll = (double*)cmalloc(numCoeffUpAllProc*sizeof(double));
  double *KSMatrixTest = (double*)cmalloc(numStateUpIdp*numStateUpIdp*sizeof(double));
  Gatherv(&cre_up[1],numCoeffUpTotal,MPI_DOUBLE,
            creAll,stowfRecvCounts,stowfDispls,MPI_DOUBLE,
            0,comm_states);
  Gatherv(&cim_up[1],numCoeffUpTotal,MPI_DOUBLE,
            cimAll,stowfRecvCounts,stowfDispls,MPI_DOUBLE,
            0,comm_states);
  Gatherv(&fcre_up[1],numCoeffUpTotal,MPI_DOUBLE,
            forceUpReAll,stowfRecvCounts,stowfDispls,MPI_DOUBLE,
            0,comm_states);
  Gatherv(&fcim_up[1],numCoeffUpTotal,MPI_DOUBLE,
            forceUpImAll,stowfRecvCounts,stowfDispls,MPI_DOUBLE,
            0,comm_states);

  if(myidState==0){
    for(iState=0;iState<numStateUpIdp;iState++){
      for(jState=iState;jState<numStateUpIdp;jState++){
        sum = 0.0;
        for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
          index1 = iState*numCoeff+iCoeff;
          index2 = jState*numCoeff+iCoeff;
          sum += creAll[index1]*forceUpReAll[index2]
                  +cimAll[index1]*forceUpImAll[index2];
        }//endfor iCoeff
        sum *= 2.0;
        index1 = iState*numCoeff+numCoeff-1;
        index2 = jState*numCoeff+numCoeff-1;
        sum += creAll[index1]*forceUpReAll[index2];
        //sum *= 0.5; I dont need to do this since our wave functions are normalized to 1
        KSMatrixTest[iState*numStateUpIdp+jState] = sum;
        KSMatrixTest[jState*numStateUpIdp+iState] = sum;
      }//endfor jState
    }//endfor iState
    for(iState=0;iState<numStateUpIdp;iState++){
      for(jState=iState;jState<numStateUpIdp;jState++){
        printf("iState %i jState %i KSMatrix Diff %lg\n",iState,jState,
                KSMatrixTest[iState*numStateUpIdp+jState]-KSMatrixReduce[iState*numStateUpIdp+jState]);
      }
    }
  }
  */
   

  /*
  for(iState=0;iState<numStateUpIdp;iState++){
    for(jState=0;jState<numStateUpIdp;jState++){
      printf("%lg ",KSMatrix[iState*numStateUpIdp+jState]);
    }
    printf("\n");
  }
  */

  
  if(numProcStates>1){
    Barrier(comm_states);
    if(myidState==0){
      //memcpy(KSMatrix,KSMatrixReduce,numStateUpIdp*numStateUpIdp*sizeof(double));
      memcpy(KSMatrix,KSMatrixReduceAll,numStateUpIdp*numStateUpIdp*sizeof(double));
      /*
      for(iState=0;iState<numStateUpIdp;iState++){
        for(jState=0;jState<numStateUpIdp;jState++){
          printf("%lg ",KSMatrix[iState*numStateUpIdp+jState]);
        }
        printf("\n");
      }
      */
    }//endif myidStatea
    Barrier(comm_states);
  }//endif numProcStates
  

  /*
  if(myidState==0){
    printf("Start construct KS matrix\n");
    for(iState=0;iState<numStateUpIdp;iState++){
      for(jState=iState;jState<numStateUpIdp;jState++){
	sum = 0.0;
	for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	  index1 = iState*numCoeff+iCoeff;
	  index2 = jState*numCoeff+iCoeff;
	  sum += wfOrthUpRe[index1]*allForceRe[index2]
		  +wfOrthUpIm[index1]*allForceIm[index2];
	}//endfor iCoeff
	sum *= 2.0;
	index1 = iState*numCoeff+numCoeff;
	index2 = jState*numCoeff+numCoeff;
	sum += wfOrthUpRe[index1]*allForceRe[index2];
	//sum *= 0.5; I dont need to do this since our wave functions are normalized to 1
	KSMatrix[iState*numStateUpIdp+jState] = sum;
	KSMatrix[jState*numStateUpIdp+iState] = sum;
      }//endfor jState
    }//endfor iState
  }
  */

/*==========================================================================*/
/* IV)  Free local memories  */

  free(forceUpRe);
  free(forceUpIm);
  if(myidState==0&&numProcStates>1){
    //free(&allForceRe[1]);
    //free(&allForceIm[1]);
    free(KSMatrixReduce);
    free(KSMatrixReduceAll);
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void diagKSMatrix(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                   CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  COMMUNICATE *communicate      = &(cp->communicate);
  
  int iChem,iCoeff,iState,jState,iComb,iProc;
  int index1,index2;
  int missMatchFlag = 0;
  int numChemPot     = stodftInfo->numChemPot;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numStateUpTotal = stodftInfo->numStateStoUp;
  int numStateUpIdp = stodftInfo->numStateUpIdp;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int numCoeffUpAllProc = numStateUpIdp*numCoeff;
  int numStatesDet = stodftInfo->numStatesDet;
  int myidState = communicate->myid_state;
  int numProcStates         = communicate->np_states;
  int numStatePrintUp = stodftInfo->numStatePrintUp;
  int numStatePrintDn = stodftInfo->numStatePrintDn;
  int numElecTrueUp = stodftInfo->numElecTrueUp;
  int smearOpt = stodftInfo->smearOpt;
  int numThreadsMKL = communicate->numThreadsMKL;
  MPI_Comm comm_states = communicate->comm_states;
  
  int *stowfRecvCounts = stodftInfo->stowfRecvCounts;
  int *stowfDispls = stodftInfo->stowfDispls;
  int *numStates2 = stodftInfo->numStates2;
  int *dsplStates2 = stodftInfo->dsplStates2;

  double energyMOSq;

  double *KSMatrix = stodftCoefPos->KSMatrix;
  //double *wfOrthUpRe = stodftCoefPos->wfOrthUpRe;
  //double *wfOrthUpIm = stodftCoefPos->wfOrthUpIm;
  double *wfOrthUpRe,*wfOrthUpIm;
  double *moUpRe = stodftCoefPos->moUpRe;
  double *moUpIm = stodftCoefPos->moUpIm;
  double *moUpReAll,*moUpImAll;
  double *coeffUpReBackup,*coeffUpImBackup;
  double *energyTest,*energyTestAll;
  double *energyLevel = stodftCoefPos->energyLevel;
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  double *numOccDetProc,*numOccDetAll;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;
  double time_st,time_end;

  if(myidState==0){
    wfOrthUpRe = (double*)cmalloc(numCoeffUpAllProc*sizeof(double));
    wfOrthUpIm = (double*)cmalloc(numCoeffUpAllProc*sizeof(double));
  }
  stodftCoefPos->energyLevel = (double*)crealloc(stodftCoefPos->energyLevel,
                                           numStateUpIdp*sizeof(double));
  energyLevel = stodftCoefPos->energyLevel;

  
  if(numProcStates>1){
    Gatherv(&cre_up[1],numCoeffUpTotal,MPI_DOUBLE,
            wfOrthUpRe,stowfRecvCounts,stowfDispls,MPI_DOUBLE,
            0,comm_states);
    Gatherv(&cim_up[1],numCoeffUpTotal,MPI_DOUBLE,
            wfOrthUpIm,stowfRecvCounts,stowfDispls,MPI_DOUBLE,
            0,comm_states);    
  }
  else{
    memcpy(wfOrthUpRe,&cre_up[1],numCoeffUpTotal*sizeof(double));
    memcpy(wfOrthUpIm,&cim_up[1],numCoeffUpTotal*sizeof(double));
  }
  cfree(&cre_up[1]);
  cfree(&cim_up[1]);
  cfree(&fcre_up[1]);
  cfree(&fcim_up[1]);

  if(myidState==0){
/*==========================================================================*/
/* I)  Diagonal KS Matrix  */
    printf("Start diagonal KS Matrix\n");
#ifdef FAST_FILTER
    stodftCoefPos->moUpRe = (double*)crealloc(stodftCoefPos->moUpRe,
                                 numCoeffUpAllProc*sizeof(double));
    stodftCoefPos->moUpIm = (double*)crealloc(stodftCoefPos->moUpIm,
                                 numCoeffUpAllProc*sizeof(double));
#else
    stodftCoefPos->moUpRe = (double*)cmalloc(numCoeffUpAllProc*sizeof(double));
    stodftCoefPos->moUpIm = (double*)cmalloc(numCoeffUpAllProc*sizeof(double));
#endif
    moUpRe = stodftCoefPos->moUpRe;
    moUpIm = stodftCoefPos->moUpIm;

    time_st = omp_get_wtime();
    diagSymMatWrapper(numStateUpIdp,KSMatrix,energyLevel);
    time_end = omp_get_wtime();
    printf("diag time %lg\n",time_end-time_st);
    for(iCoeff=0;iCoeff<numCoeffUpAllProc;iCoeff++){
      moUpRe[iCoeff] = 0.0;
      moUpIm[iCoeff] = 0.0;
    }//endfor iCoeff

/*==========================================================================*/
/* II)  Calculate MO in the PW basis  wfOrthUpRe is numCoeff*		    */
/*     numStateUpAllProc (in Fortran) so we need wfOrthUpRe times KSMatrix( */
/*     after diag). Call blas to do it.					    */
    
    printf("Calculate MO\n");
    fflush(stdout);

    time_st = omp_get_wtime();
    genMatrixMulWrapper(numCoeff,numStateUpIdp,wfOrthUpRe,
			KSMatrix,moUpRe,numThreadsMKL);
    genMatrixMulWrapper(numCoeff,numStateUpIdp,wfOrthUpIm,
			KSMatrix,moUpIm,numThreadsMKL);
    time_end = omp_get_wtime();
    printf("rotate orbital time %lg\n",time_end-time_st);
    fflush(stdout);

    cfree(wfOrthUpRe);
    cfree(wfOrthUpIm);
    printf("Start scattering wf\n");
    fflush(stdout);
  }//endif myid

/*==========================================================================*/
/* III) Scatter the wave functions  */

  coeffUpReBackup = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  coeffUpImBackup = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;

  if(numProcStates>1){
    Bcast(energyLevel,numStateUpIdp,MPI_DOUBLE,0,comm_states);

    Scatterv(moUpRe,stowfRecvCounts,stowfDispls,MPI_DOUBLE,
             &coeffUpReBackup[1],numCoeffUpTotal,MPI_DOUBLE,0,comm_states);

    Scatterv(moUpIm,stowfRecvCounts,stowfDispls,MPI_DOUBLE,
             &coeffUpImBackup[1],numCoeffUpTotal,MPI_DOUBLE,0,comm_states);
  }
  else{
    memcpy(&coeffUpReBackup[1],moUpRe,numCoeffUpAllProc*sizeof(double));
    memcpy(&coeffUpImBackup[1],moUpIm,numCoeffUpAllProc*sizeof(double));
  }

  if(myidState==0){
    cfree(moUpRe);
    cfree(moUpIm);
  } 

/*==========================================================================*/
/* IV) Double check the MO energy */

/*--------------------------------------------------------------------------*/
/* i)  Calculate H^2|phi>  */
  if(myidState==0){
    printf("Double check orbitals...\n");
    fflush(stdout);
  }
  cpcoeffs_pos->cre_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  cpcoeffs_pos->cim_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  cpcoeffs_pos->fcre_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  cpcoeffs_pos->fcim_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  cre_up = cpcoeffs_pos->cre_up;
  cim_up = cpcoeffs_pos->cim_up;
  fcre_up = cpcoeffs_pos->fcre_up;
  fcim_up = cpcoeffs_pos->fcim_up;


  calcForceWrapper(cp,class,general_data,cpcoeffs_pos,clatoms_pos,
                  coeffUpReBackup,coeffUpImBackup);
  for(iState=0;iState<numStateUpProc;iState++){
    fcim_up[iState*numCoeff+numCoeff] = 0.0; //exact 0
  }
  calcForceWrapper(cp,class,general_data,cpcoeffs_pos,clatoms_pos,
                  fcre_up,fcim_up);

/*--------------------------------------------------------------------------*/
/* iii)  Calculate <phi|H^2|phi> and compare with KS   */

  energyTest = (double*)cmalloc(numStateUpProc*sizeof(double));

  for(iState=0;iState<numStateUpProc;iState++){
  //printf("iState %i numStateUpIdp %i\n",iState,numStateUpIdp);
    energyTest[iState] = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      index1 = iState*numCoeff+iCoeff;
      energyTest[iState] += coeffUpReBackup[index1]*fcre_up[index1]+
                            coeffUpImBackup[index1]*fcim_up[index1];
    }//endfor iCoeff
    energyTest[iState] *= 2.0;
    index1 = iState*numCoeff+numCoeff;
    energyTest[iState] += coeffUpReBackup[index1]*fcre_up[index1];
    //printf("11111111 %lg %lg %lg\n",energyTest[iState],coeffUpReBackup[index1],fcre_up[index1]);
  }
  
  if(myidState==0){
    energyTestAll = (double*)cmalloc(numStateUpIdp*sizeof(double));
  }
  if(numProcStates>1){
    Gatherv(energyTest,numStateUpProc,MPI_DOUBLE,energyTestAll,
            numStates2,dsplStates2,MPI_DOUBLE,0,comm_states);
  }
  
  if(myidState==0){
    printf("Double Check: Calculate and Compare energy level\n");
    for(iState=0;iState<numStatePrintUp;iState++){
      energyMOSq = energyLevel[iState]*energyLevel[iState];
      if(fabs(energyMOSq-energyTestAll[iState])>1.0e-4){
        printf("Bad State %i\n",iState);
        missMatchFlag += 1;
      }//endif      
    }//endfor iState
    if(missMatchFlag>0){
      printf("Energy missmatch! I'll print them all.\n");
      for(iState=0;iState<numStatePrintUp;iState++){
        energyMOSq = energyLevel[iState]*energyLevel[iState];
        printf("eigvvv %i %.8lg %.8lg %.8lg\n",iState,energyLevel[iState],
                energyMOSq,fabs(energyTestAll[iState]-energyMOSq));
      }//endfor iState
      fflush(stdout);
      exit(0);
    }//endif missMatch
    else{
      for(iState=0;iState<numStatePrintUp;iState++){
        energyMOSq = energyLevel[iState]*energyLevel[iState];
        printf("eigvvv %i %.8lg %.8lg %.8lg\n",iState,energyLevel[iState],
                energyMOSq,fabs(energyTestAll[iState]-energyMOSq));
      }//endfor iState
      fflush(stdout);
    }//endif missMatchFlag
  }//endif myidState

  //copy the LUMO 
  if(numStateUpIdp>numStatesDet)stodftInfo->eigValLUMO = energyLevel[numStatesDet];
  else stodftInfo->eigValLUMO = energyLevel[numStatesDet-1]+0.1;
    
  if(numProcStates>1)Barrier(comm_states);


/*--------------------------------------------------------------------------*/
/* iv) Determine occupatation number */

  numOccDetProc = (double*)cmalloc(numStateUpProc*sizeof(double));

  if(smearOpt>0){
    if(myidState==0)printf("**Determine Chemical Potential...\n");
    calcChemPotMetal(cp,numOccDetProc);
#ifndef FAST_FILTER
    // We do this for metallic system to calculate entropy
    for(iState=0;iState<numStateUpProc;iState++){
      stodftInfo->numOccDetProc[iState] = numOccDetProc[iState];
    }
#endif 
    if(myidState==0)printf("**Finish Determining Chemical Potential...\n");
  }
  else{
    if(myidState==0){
      numOccDetAll = (double*)calloc(numStateUpIdp,sizeof(double));
      for(iState=0;iState<numElecTrueUp;iState++)numOccDetAll[iState] = sqrt(2.0);
    }
    if(numProcStates>1){
      Barrier(comm_states);
      Scatterv(numOccDetAll,numStates2,dsplStates2,MPI_DOUBLE,
               numOccDetProc,numStateUpProc,MPI_DOUBLE,0,comm_states);      
      Barrier(comm_states);
    }
    else{
      memcpy(numOccDetProc,numOccDetAll,numStateUpProc*sizeof(double));
    }
    cfree(numOccDetAll);
  }
  if(numProcStates>1)Barrier(comm_states);
/*--------------------------------------------------------------------------*/
/* v) Scale by occupied number */
  
  for(iState=0;iState<numStateUpProc;iState++){
    for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
      index1 = iState*numCoeff+iCoeff;
      cre_up[index1] = coeffUpReBackup[index1]*numOccDetProc[iState];
      cim_up[index1] = coeffUpImBackup[index1]*numOccDetProc[iState];
    }
  }

/*--------------------------------------------------------------------------*/
/* vi) Prepare for calculating electron friction if needed                  */

  METALLIC *metallic = stodftInfo->metallic;
  int electronFricFlag = metallic->electronFricFlag;
  int numStateFric;
  int countState;
  int div,res;
  int procInd,stateInd;
  int *stateFricIndex;
  double Emin,Emax;
  double sigma = metallic->sigma;
  double chemPotUpMetallic = stodftInfo->chemPotUpMetallic;
  double *ksStateChemPotRe;
  double *ksStateChemPotIm;
  double *ksEnergyFric;

  /* We believe only a small number of orbitals are needed */
  /* Therefore we will broadcast ksStateChemPot to all process */
  if(electronFricFlag==1){
    if(myidState==0){
      Emin = chemPotUpMetallic-4.0*sigma;
      Emax = chemPotUpMetallic+4.0*sigma;
      countState = 0;
      for(iState=0;iState<numStatePrintUp;iState++){
        if(energyLevel[iState]>Emin&&energyLevel[iState]<Emax)countState += 1;
      }//endfor iState
    }//endif myidState
    Bcast(&(countState),1,MPI_INT,0,comm_states);
    metallic->numStateFric = countState;
    numStateFric = countState;
    printf("numStateFric %i\n",metallic->numStateFric);

    metallic->stateFricIndex = (int*)realloc(metallic->stateFricIndex,countState*sizeof(int));
    metallic->ksEnergyFric = (double*)realloc(metallic->ksEnergyFric,countState*sizeof(double));
    metallic->ksStateChemPotRe = (double*)realloc(metallic->ksStateChemPotRe,countState*numCoeff*sizeof(double));
    metallic->ksStateChemPotIm = (double*)realloc(metallic->ksStateChemPotIm,countState*numCoeff*sizeof(double));
    stateFricIndex = metallic->stateFricIndex;
    ksEnergyFric = metallic->ksEnergyFric;
    ksStateChemPotRe = metallic->ksStateChemPotRe;
    ksStateChemPotIm = metallic->ksStateChemPotIm;

    if(myidState==0){
      countState = 0;
      for(iState=0;iState<numStatePrintUp;iState++){
        if(energyLevel[iState]>Emin&&energyLevel[iState]<Emax){
          //printf("iState %i energyLevel %lg Emin %lg Emax %lg\n",iState,energyLevel[iState],Emin,Emax);
          stateFricIndex[countState] = iState;
          ksEnergyFric[countState] = energyLevel[iState];
          countState += 1;
        }//endif
      }//endfor iState
    }//endif myidState
    Bcast(stateFricIndex,countState,MPI_INT,0,comm_states);
    Bcast(ksEnergyFric,countState,MPI_DOUBLE,0,comm_states);

    div = numStateUpIdp/numProcStates;
    res = numStateUpIdp%numProcStates;
    for(iState=0;iState<numStateFric;iState++){
      if(stateFricIndex[iState]<res*(div+1))procInd = stateFricIndex[iState]/(div+1);
      else procInd = (stateFricIndex[iState]-res*(div+1))/div+res;
      stateInd = stateFricIndex[iState]-dsplStates2[procInd];
      printf("myidState %i iState %i stateIndex %i procInd %i stateInd %i\n",myidState,iState,stateFricIndex[iState],procInd,stateInd);
      if(myidState==procInd){
        memcpy(&ksStateChemPotRe[iState*numCoeff],&coeffUpReBackup[stateInd*numCoeff],
               numCoeff*sizeof(double));
        memcpy(&ksStateChemPotIm[iState*numCoeff],&coeffUpImBackup[stateInd*numCoeff],
               numCoeff*sizeof(double));
      }
      if(numProcStates>1){
        Bcast(&ksStateChemPotRe[iState*numCoeff],numCoeff,MPI_DOUBLE,procInd,comm_states);
        Bcast(&ksStateChemPotIm[iState*numCoeff],numCoeff,MPI_DOUBLE,procInd,comm_states);
        Barrier(comm_states);
      }//endif numProcStates
    }//endfor iState
  }

#ifdef FAST_FILTER
  cfree(&fcre_up[1]);
  cfree(&fcim_up[1]);
  cfree(&cre_up[1]);
  cfree(&cim_up[1]);

  if(myidState==0){
    stodftCoefPos->moUpRe = (double*)cmalloc(numCoeffUpAllProc*sizeof(double));
    stodftCoefPos->moUpIm = (double*)cmalloc(numCoeffUpAllProc*sizeof(double));
    moUpRe = stodftCoefPos->moUpRe;
    moUpIm = stodftCoefPos->moUpIm;
  }

  if(numProcStates>1){
    Gatherv(&coeffUpReBackup[1],numCoeffUpTotal,MPI_DOUBLE,
            moUpRe,stowfRecvCounts,stowfDispls,MPI_DOUBLE,
            0,comm_states);
    Gatherv(&coeffUpImBackup[1],numCoeffUpTotal,MPI_DOUBLE,
            moUpIm,stowfRecvCounts,stowfDispls,MPI_DOUBLE,
            0,comm_states);
  }
  else{
    memcpy(moUpRe,&coeffUpReBackup[1],numCoeffUpTotal*sizeof(double));
    memcpy(moUpIm,&coeffUpImBackup[1],numCoeffUpTotal*sizeof(double));
  }
#endif 
  if(numProcStates>1)Barrier(comm_states);
  //debug
  /*
  for(iProc=0;iProc<numProcStates;iProc++){
    if(myidState==iProc){
      for(iState=0;iState<numChemPot*numStateUpProc;iState++){
	printf("myid %i occ %lg\n",myidState,numOccDetProc[iState]);
      }
    }
    Barrier(comm_states);
  }
  */

  cfree(energyTest);
  cfree(numOccDetProc);
  cfree(&coeffUpReBackup[1]);
  cfree(&coeffUpImBackup[1]);
  if(myidState==0){
    cfree(energyTestAll);
  }
  
  //exit(0);
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void diagSymMatWrapper(int n,double *A,double *w)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
  char jobz = 'V';
  char uplo = 'U';
  int lda = n;
  int lwork = 10*n;
  int info;
  double *work = (double*)cmalloc(lwork*sizeof(double));
  
  dsyev_(&jobz,&uplo,&n,&A[0],&lda,&w[0],&work[0],&lwork,&info);  

  //dsyev_("V","U",&ist.ms,&H[0],&ist.ms,&eval[0],&work[0],&lwk,&info);


  if(info!=0){
    printf("Something is wrong with the diagnolization!\n");
    printf("Check the error code %i\n",info);
    fflush(stdout);
    exit(0);
  }
  free(work);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genMatrixMulWrapper(int m,int n,double *A,double *B,double *C,int nthreads)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  char transa = 'N';
  char transb = 'N';
  int k = n;
  int lda = m;
  int ldb = n;
  int ldc = m;
  double alpha = 1.0;
  double beta = 0.0;

#ifdef MKL_THREADS
  mkl_set_dynamic(0);
  mkl_set_num_threads(nthreads);
  //omp_set_nested(1);
  omp_set_max_active_levels(1);
#endif
  {
    //DGEMM(&transa,&transb,&m,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc); 
    dgemm_(&transa,&transb,&m,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcForceWrapper(CP *cp,CLASS *class,GENERAL_DATA *general_data,
                   CPCOEFFS_POS *cpcoeffs_pos,CLATOMS_POS *clatoms_pos,
		    double *cre,double *cim)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  COMMUNICATE *communicate      = &(cp->communicate);

  int iChem,iCoeff,iState,jState,iComb,iOff;
  int index1,index2;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int myidState = communicate->myid_state;
  MPI_Comm comm_states = communicate->comm_states;

  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;

  for(iCoeff=1;iCoeff<=numCoeffUpTotal;iCoeff++){
    cre_up[iCoeff] = cre[iCoeff];
    cim_up[iCoeff] = cim[iCoeff];
    // In case I forget to zero them
    fcre_up[iCoeff] = 0.0;
    fcim_up[iCoeff] = 0.0;
  }
  //debug
  /*
  double sum;
  for(iState=0;iState<numStateUpProc;iState++){
    sum = 0.0;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      sum += cre_up[iState*numCoeff+iCoeff]*cre_up[iState*numCoeff+iCoeff]+
	     cim_up[iState*numCoeff+iCoeff]*cim_up[iState*numCoeff+iCoeff];
    }
    sum *= 2.0;
    sum += cre_up[iState*numCoeff+numCoeff]*cre_up[iState*numCoeff+numCoeff];
    printf("sum sum %lg\n",sum);
  }
  */
  //calcKSPotExtRecipWrap(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  //calcCoefForceWrapReduce(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  calcCoefForceWrapSCF(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
  for(iState=0;iState<numStateUpProc;iState++){
    iOff = iState*numCoeff;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      fcre_up[iOff+iCoeff] *= -0.25;
      fcim_up[iOff+iCoeff] *= -0.25;
    }//endfor iCoeff
    fcre_up[iOff+numCoeff] *= -0.5;
    fcim_up[iOff+numCoeff] *= -0.5;
  }//endfor iState

  //printf("fcre %lg\n",fcre_up[numCoeff]);
  /*
  double sum;
  int ioff;
  for(iState=0;iState<numStateUpProc;iState++){
    sum = 0.0;
    ioff = iState*numCoeff;
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      sum += fcre_up[ioff+iCoeff]*fcre_up[ioff+iCoeff]+fcim_up[ioff+iCoeff]*fcim_up[ioff+iCoeff];
      //printf("force test 111111111111 %lg %lg\n",fcre_up[iCoeff],fcim_up[iCoeff]);
    }
    sum *= 2.0;
    sum += fcre_up[ioff+numCoeff]*fcre_up[ioff+numCoeff];
    printf("force dot force %lg\n",sum);
  }
  */
  //exit(0);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


#ifdef FAST_FILTER   

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void broadcastWfDet(CP *cp2,CLASS *class2,GENERAL_DATA *general_data2,CP *cp)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp2->cpcoeffs_info);
  CPCOEFFS_POS  *cpcoeffs_pos     = &(cp2->cpcoeffs_pos[1]);  
  STODFTINFO *stodftInfo        = cp2->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp2->stodftCoefPos;
  CPOPTS *cpopts                = &(cp2->cpopts);
  COMMUNICATE *communicate      = &(cp2->communicate);

  int numChemPot     = stodftInfo->numChemPot;
  int numStateUpPrintProc = cpcoeffs_info->nstate_up_proc;
  int numStateUpTotal = stodftInfo->numStateStoUp;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal;
  int numCoeffDnTotal;
  int numCoeffUpAllProc = numChemPot*numStateUpTotal*numCoeff;
  int numStateUpAllProc = numChemPot*numStateUpTotal;
  int numStateUpIdp = stodftInfo->numStateUpIdp;
  int numStatesDet = stodftInfo->numStatesDet;
  int myidState = communicate->myid_state;
  int numProcStates         = communicate->np_states;
  int numStatePrintUp = stodftInfo->numStatePrintUp;
  int numStatePrintDn = stodftInfo->numStatePrintDn;
  int numStateUpProc,numStateDnProc;
  int numStatesNow;
  int numStatePrintUpProc;
  int iChem,iProc;
  int div,res;

  MPI_Comm comm_states = communicate->comm_states;

  int *stowfRecvCountsLocal = (int*)cmalloc(numProcStates*sizeof(int));
  int *stowfDisplsLocal = (int*)cmalloc(numProcStates*sizeof(int));
  int *numStates22 = cp->stodftInfo->numStates2;
  int *dsplStates22 = cp->stodftInfo->dsplStates2;
  int *stowfRecvCounts = stodftInfo->stowfRecvCounts;
  int *stowfDispls = stodftInfo->stowfDispls;

  double *moUpRe = stodftCoefPos->moUpRe;
  double *moUpIm = stodftCoefPos->moUpIm;
  double *moUpRePrint = cp->stodftCoefPos->moUpRePrint;
  double *moUpImPrint = cp->stodftCoefPos->moUpImPrint;
  double *energyLevel2 = cp2->stodftCoefPos->energyLevel;
  double *energyLevel = cp->stodftCoefPos->energyLevel;
  double *moUpReAll,*moUpImAll;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

/*--------------------------------------------------------------------------*/
/* i) Gather all wavefunctions to the master process */

  div = (int)(numStatePrintUp/numProcStates);
  res = numStatePrintUp%numProcStates;
  if(myidState<res) numStatePrintUpProc = div+1;
  else numStatePrintUpProc = div;
  if(numProcStates>1){
    Allgather(&numStatePrintUpProc,1,MPI_INT,numStates22,1,MPI_INT,0,comm_states);
  }
  else{
    numStates22[0] = numStatePrintUp;
  }
  
  dsplStates22[0] = 0;
  for(iProc=1;iProc<numProcStates;iProc++){
    dsplStates22[iProc] = dsplStates22[iProc-1]+numStates22[iProc-1];
  }  
  
  numStateUpPrintProc = numStates22[myidState];
  numCoeffUpTotal = numStateUpPrintProc*numCoeff;
  
  if(numProcStates>1)Barrier(comm_states);
  if(myidState==0){
    memcpy(energyLevel,energyLevel2,numStatePrintUp*sizeof(double));
  }

  if(numProcStates>1){
    Allgather(&numCoeffUpTotal,1,MPI_INT,stowfRecvCounts,1,MPI_INT,0,comm_states);
  }
  stowfDispls[0] = 0;
  for(iProc=1;iProc<numProcStates;iProc++){
    stowfDispls[iProc] = stowfDispls[iProc-1]+stowfRecvCounts[iProc-1];
  }
  
/*--------------------------------------------------------------------------*/
/* ii) Scatter numStatePrint orbitals to all processors */

  //cp->stodftInfo.numStatePrintUp = numStatePrintUp;

  if(numProcStates>1){
    Scatterv(moUpRe,stowfRecvCounts,stowfDispls,MPI_DOUBLE,
             moUpRePrint,numCoeffUpTotal,MPI_DOUBLE,0,comm_states);
    Scatterv(moUpIm,stowfRecvCounts,stowfDispls,MPI_DOUBLE,
             moUpImPrint,numCoeffUpTotal,MPI_DOUBLE,0,comm_states);
    Bcast(energyLevel,numStatePrintUp,MPI_DOUBLE,0,comm_states);
  }
  else{
    memcpy(moUpRePrint,moUpRe,numCoeffUpTotal*sizeof(double));
    memcpy(moUpImPrint,moUpIm,numCoeffUpTotal*sizeof(double));
  }

/*--------------------------------------------------------------------------*/
/* iii) Free local memory and allocate others */

  //cfree(stowfRecvCounts);
  //cfree(stowfDispls);
  if(myidState==0){
    cfree(moUpRe);
    cfree(moUpIm);
    stodftCoefPos->moUpRe = NULL;
    stodftCoefPos->moUpIm = NULL;
  }
  cpcoeffs_info->nstate_up_proc = stodftInfo->numStateUpProcBackup;
  cpcoeffs_info->nstate_dn_proc = stodftInfo->numStateDnProcBackup;

  numStateUpProc = cpcoeffs_info->nstate_up_proc;
  numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  numCoeffUpTotal = numStateUpProc*numCoeff;
  numCoeffDnTotal = numStateDnProc*numCoeff;
  
  //cpcoeffs_pos->cre_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  //cpcoeffs_pos->cim_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  //cpcoeffs_pos->fcre_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  //cpcoeffs_pos->fcim_up = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;

  /*   
  for(iChem=0;iChem<numChemPot;iChem++){
    stoWfUpRe[iChem] = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
    stoWfUpIm[iChem] = (double*)cmalloc(numCoeffUpTotal*sizeof(double))-1;
  }
  */

  if(numProcStates>1){
    Allgather(&numCoeffUpTotal,1,MPI_INT,stowfRecvCounts,1,MPI_INT,0,comm_states);
    stowfDispls[0] = 0;
    for(iProc=1;iProc<numProcStates;iProc++){
      stowfDispls[iProc] = stowfDispls[iProc-1]+stowfRecvCounts[iProc-1];
    }
  }

  //free(stowfRecvCounts);
  //free(stowfDispls);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

#endif
