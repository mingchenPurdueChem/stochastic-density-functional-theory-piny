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
  int numProcStates = communicate->np_states;
  int myidState = communicate->myid_state;
  MPI_Comm comm_states   =    communicate->comm_states;

  int *stowfRecvCounts = stodftInfo->stowfRecvCounts;
  int *stowfDispls = stodftInfo->stowfDispls;
  int *stowfRecvCountsComplex = stodftInfo->stowfRecvCountsComplex;
  int *stowfDisplsComplex = stodftInfo->stowfDisplsComplex;
  
  double pre = sqrt(2.0);
  double sum;
  double *wfBfOrthUp = stodftCoefPos->wfBfOrthUp;  
  double *wfOrthUpRe = stodftCoefPos->wfOrthUpRe;
  double *wfOrthUpIm = stodftCoefPos->wfOrthUpIm;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  double *allWF = (double *)cmalloc(2*numChemPot*numCoeffUpTotal*
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
  
  
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iState=0;iState<numStateUpProc;iState++){
      for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
	indCoeff2 = iState*numCoeff+iCoeff;
	indCoeff1 = indCoeff2+iChem*numCoeffUpTotal;
	allWF[2*indCoeff1] = pre*stoWfUpRe[iChem][indCoeff2+1];
	allWF[2*indCoeff1+1] = pre*stoWfUpIm[iChem][indCoeff2+1];
	//allWF[indCoeff1] = pre*stoWfUpRe[iChem][indCoeff2+1]+
	//		    I*pre*stoWfUpIm[iChem][indCoeff2+1];
      }//endfor iCoeff
      indCoeff2 = iState*numCoeff+numCoeff;
      indCoeff1 = indCoeff2+iChem*numCoeffUpTotal-1;
      //allWF[indCoeff1] = stoWfUpRe[iChem][indCoeff2]+I*stoWfUpIm[iChem][indCoeff2];
      allWF[2*indCoeff1] = stoWfUpRe[iChem][indCoeff2];
      allWF[2*indCoeff1+1] = 0.0;
    }//endfor iState
  }//endfor iChem
  
/*==========================================================================*/
/* II) Gather to master proc */

  if(numProcStates>1){
    Gatherv(allWF,2*numChemPot*numCoeffUpTotal,MPI_DOUBLE,
  	    wfBfOrthUp,stowfRecvCountsComplex,stowfDisplsComplex,MPI_DOUBLE,
	    0,comm_states);
  }
  else{
    memcpy(&wfBfOrthUp[0],&allWF[0],2*numCoeffUpAllProc*sizeof(double));
  }
 
/*==========================================================================*/
/* III) SVD */

  if(myidState==0){
    printf("Start SVD\n");
    /*
    for(iState=0;iState<numStateUpAllProc;iState++){
      printf("%lg %lg\n",creal(wfBfOrthUp[iState*numCoeff]),cimag(wfBfOrthUp[iState*numCoeff]));
    }
    
    for(iProc=0;iProc<numProcStates;iProc++){
      printf("%i recvCounts %i displs %i\n",iProc,stowfRecvCountsComplex[iProc],stowfDisplsComplex[iProc]);
    }
    */
    // test determ wf orthog, make sure we keep the correct wf
    /*
    for(iState=0;iState<numStatesDet;iState++){
      for(jState=iState;jState<numStatesDet;jState++){
	sum = 0.0;
	for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
	  iCoeff1 = iState*numCoeff+iCoeff;
	  iCoeff2 = jState*numCoeff+iCoeff;
	  sum += wfUpReDet[iCoeff1]*wfUpReDet[iCoeff2]+
		wfUpImDet[iCoeff1]*wfUpImDet[iCoeff2];
	}//endfor iCoeff
	sum *=2.0;
	iCoeff1 = iState*numCoeff+numCoeff-1;
	iCoeff2 = jState*numCoeff+numCoeff-1;
	sum += wfUpReDet[iCoeff1]*wfUpReDet[iCoeff2];
	printf("iState %i jState %i dot %lg\n",iState,jState,sum);
      }//endfor jState
    }//endfor iState
    */
    //Normalize them!
    /*
    pre = 1.0/sqrt(2.0);
    for(iState=0;iState<numStatesDet;iState++){
      for(iCoeff=0;iCoeff<numCoeff;iCoeff++){
	iCoeff1 = iState*numCoeff+iCoeff;
	wfUpReDet[iCoeff1] *= pre;
	wfUpImDet[iCoeff1] *= pre;
      }//endfor iCoeff
    }
    //correct k!=0 terms
    pre = sqrt(2.0);
    for(iState=0;iState<numStatesDet;iState++){
      for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
        iCoeff1 = iState*numCoeff+iCoeff;
        wfUpReDet[iCoeff1] *= pre;
        wfUpImDet[iCoeff1] *= pre;
      }//endfor iCoeff
    }//endfor iState
    */
    // Done, next let's just copy those WF to SVD to check
    /*
    int div = numStateUpAllProc/numStatesDet;
    int res = numStateUpAllProc%numStatesDet;
    int iCpy;
    for(iCpy=0;iCpy<div;iCpy++){
      for(iState=0;iState<numStatesDet;iState++){
	for(iCoeff=0;iCoeff<numCoeff;iCoeff++){
	  iCoeff1 = iCpy*numStatesDet*numCoeff+iState*numCoeff+iCoeff;
	  iCoeff2 = iState*numCoeff+iCoeff;
	  //printf("iCoeff1 %i iCoeff2 %i\n",iCoeff1,iCoeff2);
	  wfBfOrthUp[2*iCoeff1] = wfUpReDet[iCoeff2];
          wfBfOrthUp[2*iCoeff1+1] = wfUpImDet[iCoeff2];
	}//endfor iCoeff
      }//endfor iState
    }//endfor iCpy
    for(iState=0;iState<res;iState++){
      for(iCoeff=0;iCoeff<numCoeff;iCoeff++){
	iCoeff1 = div*numStatesDet*numCoeff+iState*numCoeff+iCoeff;
	iCoeff2 = iState*numCoeff+iCoeff;
	//printf("iCoeff1 %i iCoeff2 %i\n",iCoeff1,iCoeff2);
	wfBfOrthUp[2*iCoeff1] = wfUpReDet[iCoeff2];
	wfBfOrthUp[2*iCoeff1+1] = wfUpImDet[iCoeff2];
      }//endfor iCoeff
    }//endfor iState
    */
    indexCut = orthogSVD(numStateUpAllProc,2*numCoeff,wfBfOrthUp);
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

    // I'll keep 
    for(iState=0;iState<numStateUpAllProc;iState++){
      sum = 0.0;
      for(iCoeff=0;iCoeff<numCoeff-1;iCoeff++){
	iCoeff1 = iState*numCoeff+iCoeff;
	wfBfOrthUp[2*iCoeff1] /= pre;
	wfBfOrthUp[2*iCoeff1+1] /= pre;
	sum += wfBfOrthUp[2*iCoeff1]*wfBfOrthUp[2*iCoeff1]+
		wfBfOrthUp[2*iCoeff1+1]*wfBfOrthUp[2*iCoeff1+1];
	//wfOrthUpRe[iCoeff1+1] = wfBfOrthUp[2*iCoeff1];
	//wfOrthUpIm[iCoeff1+1] = wfBfOrthUp[2*iCoeff1+1];
      }
      sum *= 2.0;
      iCoeff1 = iState*numCoeff+numCoeff-1; 
      sum += wfBfOrthUp[2*iCoeff1]*wfBfOrthUp[2*iCoeff1];
      sum = sqrt(sum);
      for(iCoeff=0;iCoeff<numCoeff;iCoeff++){
	iCoeff1 = iState*numCoeff+iCoeff;
	wfBfOrthUp[2*iCoeff1] /= sum;
	wfBfOrthUp[2*iCoeff1+1] /= sum;
	wfOrthUpRe[iCoeff1+1] = wfBfOrthUp[2*iCoeff1];
	wfOrthUpIm[iCoeff1+1] = wfBfOrthUp[2*iCoeff1+1];
      }
      wfOrthUpIm[iState*numCoeff+numCoeff] = 0.0;
      //wfOrthUpRe[iCoeff1+1] = wfBfOrthUp[2*iCoeff1];
      //wfOrthUpIm[iCoeff1+1] = 0.0;  
    }//endfor iState
    //test orthogal
    /*
    for(iState=0;iState<indexCut;iState++){
      for(jState=iState;jState<indexCut;jState++){
	sum = 0.0;
	for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	  iCoeff1 = iState*numCoeff+iCoeff;
	  iCoeff2 = jState*numCoeff+iCoeff;
	  sum += wfOrthUpRe[iCoeff1]*wfOrthUpRe[iCoeff2]+
		 wfOrthUpIm[iCoeff1]*wfOrthUpIm[iCoeff2];
	}
	sum *= 2.0;
	iCoeff1 = iState*numCoeff+numCoeff;
	iCoeff2 = jState*numCoeff+numCoeff;
	sum += wfOrthUpRe[iCoeff1]* wfOrthUpRe[iCoeff2];
	printf("iState %i jState %i %lg\n",iState,jState,sum);
      }//enfor jState
    }//endfor iState
    */
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

  if(myidState==0)printf("Start Scatter Data\n");

  if(numProcStates>1){
    Scatterv(wfBfOrthUp,stowfRecvCountsComplex,stowfDisplsComplex,MPI_DOUBLE,
  	      allWF,2*numChemPot*numCoeffUpTotal,MPI_DOUBLE,0,comm_states);
  }
  else{
    memcpy(&allWF[0],&wfBfOrthUp[0],2*numCoeffUpAllProc*sizeof(double));
  }

/*==========================================================================*/
/* IV) Unpack all stochastic wave functions */

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

/*==========================================================================*/
/* V) Normalization all stochastic wave functions */
  
  if(myidState==0)printf("Start Normalization\n");  

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
int orthogSVD(int m,int n,double *matrix)
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
  dgesvd_(&jobu,&jobvt,&n,&m,&matrix[0],&lda,&S[0],NULL,&ldu,NULL,&ldvt,
	  &work[0],&lwork,&info);
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
  int index1,index2;
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

  double sum;

  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;
  double *forceUpRe,*forceUpIm;
  double *allForceRe,*allForceIm;
  double *wfOrthUpRe = stodftCoefPos->wfOrthUpRe;
  double *wfOrthUpIm = stodftCoefPos->wfOrthUpIm;
  double *KSMatrix = stodftCoefPos->KSMatrix;

  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

/*==========================================================================*/
/* I)  Calculate force	*/
    
  if(myidState==0){
    allForceRe = (double*)malloc(numCoeffUpAllProc*sizeof(double))-1;
    allForceIm = (double*)malloc(numCoeffUpAllProc*sizeof(double))-1;
  }
  allForceRe = (double*)malloc(numCoeffUpAllProc*sizeof(double))-1;
  allForceIm = (double*)malloc(numCoeffUpAllProc*sizeof(double))-1;
  //Barrier(comm_states);
  forceUpRe = (double*)cmalloc(numChemPot*numCoeffUpTotal*sizeof(double))-1;
  forceUpIm = (double*)cmalloc(numChemPot*numCoeffUpTotal*sizeof(double))-1;

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

/*==========================================================================*/
/* II)  Gather all forces to the master proc  */

  if(myidState==0)printf("Start gather force\n");
  
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

/*==========================================================================*/
/* III)  Construct KS matrix  */

  if(myidState==0){
    /*
    for(iState=0;iState<numStateUpAllProc;iState++){
      printf("iState %i fre %lg fim %lg\n",iState,allForceRe[iState*numCoeff+numCoeff],
	     allForceIm[iState*numCoeff+numCoeff]);
    }
    */
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
    /*
    for(iState=0;iState<numStateUpIdp*numStateUpIdp;iState++){
      printf("KS Mat %i %lg\n",iState,KSMatrix[iState]);
    }
    */
  }

/*==========================================================================*/
/* IV)  Free local memories  */

  free(&forceUpRe[1]);
  free(&forceUpIm[1]);
  if(myidState==0){
    free(&allForceRe[1]);
    free(&allForceIm[1]);
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
  int missMatchFlag;
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
  int numStatesDet = stodftInfo->numStatesDet;
  int myidState = communicate->myid_state;
  int numProcStates         = communicate->np_states;
  int numStatePrintUp = stodftInfo->numStatePrintUp;
  int numStatePrintDn = stodftInfo->numStatePrintDn;
  int smearOpt = stodftInfo->smearOpt;
  MPI_Comm comm_states = communicate->comm_states;
  
  int *stowfRecvCounts = stodftInfo->stowfRecvCounts;
  int *stowfDispls = stodftInfo->stowfDispls;

  double energyMOSq;

  double *KSMatrix = stodftCoefPos->KSMatrix;
  double *wfOrthUpRe = stodftCoefPos->wfOrthUpRe;
  double *wfOrthUpIm = stodftCoefPos->wfOrthUpIm;
  double *moUpRe = stodftCoefPos->moUpRe;
  double *moUpIm = stodftCoefPos->moUpIm;
  double *moUpReAll,*moUpImAll;
  double *allForceRe,*allForceIm,*forceUpRe,*forceUpIm;
  double *energyTest;
  double *energyLevel = stodftCoefPos->energyLevel;
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *fcre_up = cpcoeffs_pos->fcre_up;
  double *fcim_up = cpcoeffs_pos->fcim_up;
  double *fcre_dn = cpcoeffs_pos->fcre_dn;
  double *fcim_dn = cpcoeffs_pos->fcim_dn;


  double **stoWfUpRe = stodftCoefPos->stoWfUpRe;
  double **stoWfUpIm = stodftCoefPos->stoWfUpIm;
  double **stoWfDnRe = stodftCoefPos->stoWfDnRe;
  double **stoWfDnIm = stodftCoefPos->stoWfDnIm;

  if(myidState==0){
/*==========================================================================*/
/* I)  Diagonal KS Matrix  */
    printf("Start diagonal KS Matrix\n");

    moUpReAll = (double*)cmalloc(numCoeffUpAllProc*sizeof(double))-1;
    moUpImAll = (double*)cmalloc(numCoeffUpAllProc*sizeof(double))-1;

    diagSymMatWrapper(numStateUpIdp,KSMatrix,energyLevel);
    for(iCoeff=1;iCoeff<=numCoeffUpAllProc;iCoeff++){
      moUpReAll[iCoeff] = 0.0;
      moUpImAll[iCoeff] = 0.0;
    }//endfor iCoeff

/*==========================================================================*/
/* II)  Calculate MO in the PW basis  wfOrthUpRe is numCoeff*		    */
/*     numStateUpAllProc (in Fortran) so we need wfOrthUpRe times KSMatrix( */
/*     after diag). Call blas to do it.					    */
    
    printf("Calculate MO\n");

    genMatrixMulWrapper(numCoeff,numStateUpIdp,&wfOrthUpRe[1],
			KSMatrix,&moUpReAll[1]);
    genMatrixMulWrapper(numCoeff,numStateUpIdp,&wfOrthUpIm[1],
			KSMatrix,&moUpImAll[1]);

    printf("Start scattering wf\n");
  }//endif myid

/*==========================================================================*/
/* III) Scatter the wave functions  */

  if(numProcStates>1){
    Bcast(energyLevel,numStateUpAllProc,MPI_DOUBLE,0,comm_states);

    Scatterv(&moUpReAll[1],stowfRecvCounts,stowfDispls,MPI_DOUBLE,
             &moUpRe[1],numChemPot*numCoeffUpTotal,MPI_DOUBLE,0,comm_states);

    Scatterv(&moUpImAll[1],stowfRecvCounts,stowfDispls,MPI_DOUBLE,
             &moUpIm[1],numChemPot*numCoeffUpTotal,MPI_DOUBLE,0,comm_states);
  }
  else{
    memcpy(&moUpRe[1],&moUpReAll[1],numCoeffUpAllProc*sizeof(double));
    memcpy(&moUpIm[1],&moUpImAll[1],numCoeffUpAllProc*sizeof(double));
  }

/*==========================================================================*/
/* IV) Double check the MO energy */

/*--------------------------------------------------------------------------*/
/* i)  Calculate force  */

  if(myidState==0){
    printf("Double Check: Calculate force\n");
    allForceRe = (double*)cmalloc(numCoeffUpAllProc*sizeof(double))-1;
    allForceIm = (double*)cmalloc(numCoeffUpAllProc*sizeof(double))-1;
  }
  forceUpRe = (double*)cmalloc(numChemPot*numCoeffUpTotal*sizeof(double))-1;
  forceUpIm = (double*)cmalloc(numChemPot*numCoeffUpTotal*sizeof(double))-1;

  for(iChem=0;iChem<numChemPot;iChem++){
    calcForceWrapper(cp,class,general_data,cpcoeffs_pos,clatoms_pos,
                    moUpRe+iChem*numCoeffUpTotal,moUpIm+iChem*numCoeffUpTotal);
    for(iState=0;iState<numStateUpProc;iState++){
      fcim_up[iState*numCoeff+numCoeff] = 0.0; //exact 0
    }
    calcForceWrapper(cp,class,general_data,cpcoeffs_pos,clatoms_pos,
		    fcre_up,fcim_up);
    memcpy(&forceUpRe[iChem*numCoeffUpTotal+1],&fcre_up[1],numCoeffUpTotal*sizeof(double));
    memcpy(&forceUpIm[iChem*numCoeffUpTotal+1],&fcim_up[1],numCoeffUpTotal*sizeof(double));
  }//endfor iChem

/*--------------------------------------------------------------------------*/
/* ii)  Gather all forces to the master proc  */

  if(numProcStates>1)Barrier(comm_states);
  if(myidState==0)printf("Double Check: gather force\n");

  if(numProcStates>1){
    Gatherv(&forceUpRe[1],numChemPot*numCoeffUpTotal,MPI_DOUBLE,
              &allForceRe[1],stowfRecvCounts,stowfDispls,MPI_DOUBLE,0,comm_states);
    Gatherv(&forceUpIm[1],numChemPot*numCoeffUpTotal,MPI_DOUBLE,
              &allForceIm[1],stowfRecvCounts,stowfDispls,MPI_DOUBLE,0,comm_states);
  }
  else{
    memcpy(&allForceRe[1],&forceUpRe[1],numCoeffUpAllProc*sizeof(double));
    memcpy(&allForceIm[1],&forceUpIm[1],numCoeffUpAllProc*sizeof(double));
  }
  //printf("11111111111111111\n");
  if(numProcStates>1)Barrier(comm_states);

/*--------------------------------------------------------------------------*/
/* iii)  Calculate <phi|H|phi> and compare with KS   */
    
  if(myidState==0){
    printf("Double Check: Calculate and Compare energy level\n");
    energyTest = (double*)cmalloc(numStateUpIdp*sizeof(double));
    missMatchFlag = 0;
    for(iState=0;iState<numStatePrintUp;iState++){
      //printf("iState %i numStateUpIdp %i\n",iState,numStateUpIdp);
      energyTest[iState] = 0.0;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	index1 = iState*numCoeff+iCoeff;
	energyTest[iState] += moUpReAll[index1]*allForceRe[index1]+
			      moUpImAll[index1]*allForceIm[index1];
      }//endfor iCoeff
      energyTest[iState] *= 2.0;
      index1 = iState*numCoeff+numCoeff;
      energyTest[iState] += moUpReAll[index1]*allForceRe[index1];   
      energyMOSq = energyLevel[iState]*energyLevel[iState];
      if(fabs(energyMOSq-energyTest[iState])>1.0e-4){
	printf("Bad State %i\n",iState);
	missMatchFlag += 1;
      }
    }//endfor iState
    if(missMatchFlag>0){
      printf("Energy missmatch! I'll print them all.\n");
      for(iState=0;iState<numStatePrintUp;iState++){
	energyMOSq = energyLevel[iState]*energyLevel[iState];
	printf("eigvvv %i %.8lg %.8lg %.8lg\n",iState,energyLevel[iState],
		energyMOSq,energyTest[iState]);
      }//endfor iState
      fflush(stdout);
      exit(0);
    }//endif missMatch
    //print something
    for(iState=0;iState<numStatePrintUp;iState++){
      energyMOSq = energyLevel[iState]*energyLevel[iState];
      printf("eigvvv %i %.8lg %.8lg %.8lg\n",iState,energyLevel[iState],
              energyMOSq,fabs(energyTest[iState]-energyMOSq));
    }//endfor iState
    free(&allForceRe[1]);
    free(&allForceIm[1]);
    free(energyTest);
    free(&moUpReAll[1]);
    free(&moUpImAll[1]);
  }//endif myidState
  free(&forceUpRe[1]);
  free(&forceUpIm[1]);

  //copy the LUMO 
  if(numStateUpIdp>numStatesDet)stodftInfo->eigValLUMO = energyLevel[numStatesDet];
  else stodftInfo->eigValLUMO = energyLevel[numStatesDet-1]+0.1;
    
  if(numProcStates>1)Barrier(comm_states);


/*--------------------------------------------------------------------------*/
/* iv) Determine occupatation number for metallic systems */

  if(smearOpt>0){
    if(myidState==0)printf("**Determine Chemical Potential...\n");
    calcChemPotMetal(cp);
    if(myidState==0)printf("**Finish Determining Chemical Potential...\n");
  }
  if(numProcStates>1)Barrier(comm_states);
/*--------------------------------------------------------------------------*/
/* v) Scale by occupied number */
  
  double *numOccDetProc = stodftInfo->numOccDetProc;
  for(iChem=0;iChem<numChemPot;iChem++){
    for(iState=0;iState<numStateUpProc;iState++){
      for(iCoeff=1;iCoeff<=numCoeff;iCoeff++){
	index1 = iState*numCoeff+iCoeff;
	index2 = iChem*numCoeffUpTotal+index1;
	stoWfUpRe[iChem][index1] = moUpRe[index2]*numOccDetProc[iChem*numStateUpProc+iState];
        stoWfUpIm[iChem][index1] = moUpIm[index2]*numOccDetProc[iChem*numStateUpProc+iState];
      }
    }
  }
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
void genMatrixMulWrapper(int m,int n,double *A,double *B,double *C)
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

  DGEMM(&transa,&transb,&m,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc); 

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
  STODFTINFO *stodftInfo        = cp2->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp2->stodftCoefPos;
  CPOPTS *cpopts                = &(cp2->cpopts);
  COMMUNICATE *communicate      = &(cp2->communicate);

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
  int numStatesDet = stodftInfo->numStatesDet;
  int myidState = communicate->myid_state;
  int numProcStates         = communicate->np_states;
  int numStatePrintUp = stodftInfo->numStatePrintUp;
  int numStatePrintDn = stodftInfo->numStatePrintDn;

  MPI_Comm comm_states = communicate->comm_states;

  int *stowfRecvCounts = stodftInfo->stowfRecvCounts;
  int *stowfDispls = stodftInfo->stowfDispls;

  double *moUpRe = stodftCoefPos->moUpRe;
  double *moUpIm = stodftCoefPos->moUpIm;
  double *moUpRePrint = cp->stodftCoefPos->moUpRePrint;
  double *moUpImPrint = cp->stodftCoefPos->moUpImPrint;
  double *energyLevel2 = cp2->stodftCoefPos->energyLevel;
  double *energyLevel = cp->stodftCoefPos->energyLevel;
  double *moUpReAll,*moUpImAll;

/*--------------------------------------------------------------------------*/
/* i) Gather all wavefunctions to the master process */

  
  
  if(myidState==0){
    moUpReAll = (double*)cmalloc(numCoeffUpAllProc*sizeof(double))-1;
    moUpImAll = (double*)cmalloc(numCoeffUpAllProc*sizeof(double))-1;
  }
  
  if(numProcStates>1){
    Gatherv(&moUpRe[1],numChemPot*numCoeffUpTotal,MPI_DOUBLE,
            &moUpReAll[1],stowfRecvCounts,stowfDispls,MPI_DOUBLE,
            0,comm_states);
    Gatherv(&moUpIm[1],numChemPot*numCoeffUpTotal,MPI_DOUBLE,
            &moUpImAll[1],stowfRecvCounts,stowfDispls,MPI_DOUBLE,
            0,comm_states);
  }
  else{
    memcpy(&moUpReAll[1],&moUpRe[1],numCoeffUpAllProc*sizeof(double));
    memcpy(&moUpImAll[1],&moUpIm[1],numCoeffUpAllProc*sizeof(double));
  }
  if(numProcStates>1)Barrier(comm_states);
  if(myidState==0){
    memcpy(&moUpRePrint[1],&moUpReAll[1],numStatePrintUp*numCoeff*sizeof(double));
    memcpy(&moUpImPrint[1],&moUpImAll[1],numStatePrintUp*numCoeff*sizeof(double));
    memcpy(energyLevel,energyLevel2,numStatePrintUp*sizeof(double));
  }

/*--------------------------------------------------------------------------*/
/* ii) Broadcast numStatePrint orbitals to all processors */

  //cp->stodftInfo.numStatePrintUp = numStatePrintUp;

  if(numProcStates>1){
    Bcast(&moUpRePrint[1],numStatePrintUp*numCoeff,MPI_DOUBLE,0,comm_states);
    Bcast(&moUpImPrint[1],numStatePrintUp*numCoeff,MPI_DOUBLE,0,comm_states);
    Bcast(energyLevel,numStatePrintUp,MPI_DOUBLE,0,comm_states);
  }

/*--------------------------------------------------------------------------*/
/* iii) Free local memory */

  if(myidState==0){
    free(&moUpReAll[1]);
    free(&moUpImAll[1]);
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

#endif
