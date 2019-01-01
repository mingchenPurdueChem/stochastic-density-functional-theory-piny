/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: frag-scf.c                                   */
/*                                                                          */
/* Driver Routine for fragment SCF calculation                              */
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
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
//#include "../proto_defs/proto_stodft_local.h"
#include "../proto_defs/proto_integrate_cpmin_entry.h"
#include "../proto_defs/proto_interface_frag_entry.h"
#include "../proto_defs/proto_frag_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void fragScf(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
            CP *cp,ANALYSIS *analysis,GENERAL_DATA *generalDataMini,
	    CP *cpMini,CLASS *classMini,ANALYSIS *analysisMini,BONDED *bondedMini,
	    int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo	 = stodftInfo->fragInfo;
  COMMUNICATE   *communicate = &(cp->communicate);
  
  int numFragProc = fragInfo->numFragProc;
  int iFrag;
  //debug
  char fileNameFragMO[100];
  FILE *fileFragMO;
  int iState,iGrid;
  int ncoef,icoef;
  int myidState		= communicate->myid_state;
  int numProcStates     = communicate->np_states;
  int fragOpt		= stodftInfo->fragOpt;
  MPI_Comm commStates   = communicate->comm_states;

  int *numGridFragProc = fragInfo->numGridFragProc;
  int *numElecUpFragProc = fragInfo->numElecUpFragProc;

  //debug
  //ncoef = cpMini[0].cpcoeffs_info.ncoef;
  //int nstates = numElecUpFragProc[0];
  //int iAtom;
  //double *cre_temp = (double*)calloc(nstates*ncoef,sizeof(double));
  //double *cim_temp = (double*)calloc(nstates*ncoef,sizeof(double));

  /*
  if(myidState==0){
    for(iAtom=0;iAtom<classMini[0].clatoms_info.natm_tot;iAtom++){
      printf("ccccccccccccoord %.16lg %.16lg %.16lg\n",
	     classMini[0].clatoms_pos[1].x[iAtom+1],
             classMini[0].clatoms_pos[1].y[iAtom+1],
             classMini[0].clatoms_pos[1].z[iAtom+1]);
    }
  }
  Barrier(commStates);
  fflush(stdout);
  exit(0);
  */
  // This is to be used to read in a deterministic fragment and broadcast to all
  // fragments. Only use it for crystal/test
  /*
  printf("numFragProc %i nstates %i ncoef %i\n",numFragProc,nstates,ncoef);
  printf("before read coef %lg %lg\n",
         cpMini[0].cpcoeffs_pos[1].cre_up[1000],
         cpMini[0].cpcoeffs_pos[1].cim_up[1000]);
  if(myidState==0){
    fileFragMO = fopen("frag-MO-2","r");
    for(iState=0;iState<nstates;iState++){
      for(icoef=0;icoef<ncoef;icoef++){
        fscanf(fileFragMO,"%lg",&(cre_temp[iState*ncoef+icoef]));
        fscanf(fileFragMO,"%lg",&(cim_temp[iState*ncoef+icoef]));
      }
    }
    fclose(fileFragMO);
  }
  Barrier(commStates);
  MPI_Bcast(cre_temp,nstates*ncoef,MPI_DOUBLE,0,commStates);
  MPI_Bcast(cim_temp,nstates*ncoef,MPI_DOUBLE,0,commStates);
  Barrier(commStates);
  for(iFrag=0;iFrag<numFragProc;iFrag++){

    memcpy(&(cpMini[iFrag].cpcoeffs_pos[1].cre_up[1]),&cre_temp[0],nstates*ncoef*sizeof(double));
    memcpy(&(cpMini[iFrag].cpcoeffs_pos[1].cim_up[1]),&cim_temp[0],nstates*ncoef*sizeof(double));
  }
  free(cre_temp);
  free(cim_temp); 

  printf("after read coef %lg %lg\n",
         cpMini[0].cpcoeffs_pos[1].cre_up[1000],
         cpMini[0].cpcoeffs_pos[1].cim_up[1000]);
  */
 
  /* 
  sprintf(fileNameFragMO,"coef-si333-%i",myidState);
  if(numFragProc>0){
    //fileFragMO = fopen(fileNameFragMO,"r");
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      fileFragMO = fopen(fileNameFragMO,"r"); 	
      //printf("%p\n",fileFragMO);
      ncoef = cpMini[iFrag].cpcoeffs_info.ncoef;
      //printf("ncoef %i\n",ncoef);
      printf("before re %.16lg im %.16lg\n",cpMini[iFrag].cpcoeffs_pos[1].cre_up[1],cpMini[iFrag].cpcoeffs_pos[1].cim_up[1]);
      for(iState=0;iState<numElecUpFragProc[iFrag];iState++){
        for(icoef=1;icoef<=ncoef;icoef++){
          fscanf(fileFragMO,"%lg",&(cpMini[iFrag].cpcoeffs_pos[1].cre_up[iState*ncoef+icoef]));
          fscanf(fileFragMO,"%lg",&(cpMini[iFrag].cpcoeffs_pos[1].cim_up[iState*ncoef+icoef]));
        }
      }
      printf("after re %.16lg im %.16lg\n",cpMini[iFrag].cpcoeffs_pos[1].cre_up[1],cpMini[iFrag].cpcoeffs_pos[1].cim_up[1]);
      fclose(fileFragMO);
    }
    //fclose(fileFragMO);
  }
  */

    
  sprintf(fileNameFragMO,"frag-MO-%i",myidState);
  if(numFragProc>0){
    fileFragMO = fopen(fileNameFragMO,"r");
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      //printf("%p\n",fileFragMO);
      ncoef = cpMini[iFrag].cpcoeffs_info.ncoef;
      for(iState=0;iState<numElecUpFragProc[iFrag];iState++){
        for(icoef=1;icoef<=ncoef;icoef++){
          fscanf(fileFragMO,"%lg",&(cpMini[iFrag].cpcoeffs_pos[1].cre_up[iState*ncoef+icoef]));
          fscanf(fileFragMO,"%lg",&(cpMini[iFrag].cpcoeffs_pos[1].cim_up[iState*ncoef+icoef]));
        }
      }
    }
    fclose(fileFragMO);
  }
    
  for(iFrag=0;iFrag<numFragProc;iFrag++){
/*======================================================================*/
/* I) Allocate Mini Structures                                          */

/*======================================================================*/
/* I) Initialize Fragment SCF					        */
    
    /*
    parseFrag(class,bonded,general_data,cp,analysis,&classMini[iFrag],&bondedMini[iFrag],
           &generalDataMini[iFrag],&cpMini[iFrag],&analysisMini[iFrag]);
    */

/*======================================================================*/
/* II) SCF LOOP					                        */
  
    controlCpMinFrag(&classMini[iFrag],&bondedMini[iFrag],&generalDataMini[iFrag],
                     &cpMini[iFrag],&analysisMini[iFrag]);      
    
    /*
    if(myidState==0){
      fileFragMO = fopen("frag-MO-2","w");
      for(iState=0;iState<numElecUpFragProc[0];iState++){
        for(icoef=1;icoef<=ncoef;icoef++){
	  fprintf(fileFragMO,"%.16lg %.16lg\n",
	  	cpMini[0].cpcoeffs_pos[1].cre_up[iState*ncoef+icoef],
	  	cpMini[0].cpcoeffs_pos[1].cim_up[iState*ncoef+icoef]);
        }
      }
    }
    Barrier(commStates);
    fflush(stdout);
    exit(0);
    */
  }//endfor iFrag
  
  
  if(numProcStates>1)Barrier(commStates);
  
  /*  
  sprintf(fileNameFragMO,"frag-MO-%i",myidState);
  fileFragMO = fopen(fileNameFragMO,"w");
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    ncoef = cpMini[iFrag].cpcoeffs_info.ncoef;
    for(iState=0;iState<numElecUpFragProc[iFrag];iState++){
      for(icoef=1;icoef<=ncoef;icoef++){
	fprintf(fileFragMO,"%.16lg %.16lg\n",cpMini[iFrag].cpcoeffs_pos[1].cre_up[iState*ncoef+icoef],
	       cpMini[iFrag].cpcoeffs_pos[1].cim_up[iState*ncoef+icoef]);
      }
    }
  }
  fclose(fileFragMO);    
  */
  
  if(numProcStates>1)Barrier(commStates);
  //exit(0);

/*======================================================================*/
/* II) Transfer Data and Free Memory                                    */

  /*  
  sprintf(fileNameFragMO,"frag-MO-%i",myidState);
  if(numFragProc>0){
    fileFragMO = fopen(fileNameFragMO,"r");
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      //printf("%p\n",fileFragMO);
      ncoef = cpMini[iFrag].cpcoeffs_info.ncoef;
      for(iState=0;iState<numElecUpFragProc[iFrag];iState++){
	for(icoef=1;icoef<=ncoef;icoef++){
	  fscanf(fileFragMO,"%lg",&(cpMini[iFrag].cpcoeffs_pos[1].cre_up[iState*ncoef+icoef]));
	  fscanf(fileFragMO,"%lg",&(cpMini[iFrag].cpcoeffs_pos[1].cim_up[iState*ncoef+icoef]));
	}
      }
    }
    fclose(fileFragMO);
  }
  */

  if(fragOpt==1){
    projRhoMiniMol(cp,general_data,class,cpMini,generalDataMini,classMini,ip_now);
  }
  if(fragOpt==3){
    projRhoMiniUnitCell(cp,general_data,class,cpMini,generalDataMini,classMini,ip_now);
  }
  //fflush(stdout);
  //exit(0);
  energyCorrect(cpMini,generalDataMini,classMini,cp,class,ip_now);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void checkpointFragOutput(CP *cp,CLASS *class)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  STODFTINFO *stodftInfo = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CPOPTS *cpOpts	     = &(cp->cpopts);
  COMMUNICATE *commCP	     = &(cp->communicate);

  int iGrid,iAtom;
  int cpLsda                = cpOpts->cp_lsda;
  int myidState             = commCP->myid_state;
  int numProcStates         = commCP->np_states;
  int rhoRealGridNum        = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot        = stodftInfo->rhoRealGridTot;
  int numAtomTot            = clatoms_info->natm_tot;
  MPI_Comm commStates	    = commCP->comm_states;
    
  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;

  double *rhoUpFragSum = fragInfo->rhoUpFragSum;
  double *rhoDnFragSum = fragInfo->rhoDnFragSum;
  double *rhoUpFragSumCpy = fragInfo->rhoUpFragSumCpy;
  double *rhoDnFragSumCpy = fragInfo->rhoDnFragSumCpy;

  double *rhoTemp = (double*)cmalloc(rhoRealGridTot*sizeof(double));

  FILE *fileCheckpoint;

/*======================================================================*/
/* II) Density part		                                        */

  if(myidState==0){
    // check the aviability of checkpoint file
    fileCheckpoint = fopen("frag-checkpoint","r");
    if(fileCheckpoint!=NULL){
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      printf("Checkpoint file already exists. I will cover the\n");
      printf("old checkpoint file and write new data into it.\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      fflush(stdout);
      fclose(fileCheckpoint);
    }
    printf("Start writing new fragmentation checkpoint file...\n");
    fileCheckpoint = fopen("frag-checkpoint","w");
    fprintf(fileCheckpoint,"%.16lg\n",stodftInfo->numElecTrueFrag);
  }    

  if(numProcStates>1){
    Barrier(commStates);
    Gatherv(rhoUpFragSum,rhoRealGridNum,MPI_DOUBLE,rhoTemp,
	    rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,0,commStates);    
    Barrier(commStates);
  }
  else{
    memcpy(rhoTemp,rhoUpFragSum,rhoRealGridTot*sizeof(double));
  }
  if(myidState==0){
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      fprintf(fileCheckpoint,"%.16lg\n",rhoTemp[iGrid]);
    }
  }
  if(cpLsda==1){
    if(numProcStates>1){
      Barrier(commStates);
      Gatherv(rhoDnFragSum,rhoRealGridNum,MPI_DOUBLE,rhoTemp,
	      rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,0,commStates);
      Barrier(commStates);
    }
    else{
      memcpy(rhoTemp,rhoDnFragSum,rhoRealGridTot*sizeof(double));
    }
    if(myidState==0){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	fprintf(fileCheckpoint,"%.16lg\n",rhoTemp[iGrid]);
      }
    }    
  }

/*======================================================================*/
/* II) Initial density                                                  */

  if(numProcStates>1){
    Barrier(commStates);
    Gatherv(rhoUpFragSumCpy,rhoRealGridNum,MPI_DOUBLE,rhoTemp,
            rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,0,commStates);
    Barrier(commStates);
  }
  else{
    memcpy(rhoTemp,rhoUpFragSumCpy,rhoRealGridTot*sizeof(double));
  }
  if(myidState==0){
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      fprintf(fileCheckpoint,"%.16lg\n",rhoTemp[iGrid]);
    }
  }
  if(cpLsda==1){
    if(numProcStates>1){
      Barrier(commStates);
      Gatherv(rhoDnFragSumCpy,rhoRealGridNum,MPI_DOUBLE,rhoTemp,
              rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,0,commStates);
      Barrier(commStates);
    }
    else{
      memcpy(rhoTemp,rhoDnFragSumCpy,rhoRealGridTot*sizeof(double));
    }
    if(myidState==0){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
        fprintf(fileCheckpoint,"%.16lg\n",rhoTemp[iGrid]);
      }
    }
  }
  
/*======================================================================*/
/* II) Energy and Force part                                            */

  if(myidState==0){
    fprintf(fileCheckpoint,"%.16lg %.16lg\n",fragInfo->keCor,fragInfo->vnlCor);
    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      fprintf(fileCheckpoint,"%.16lg %.16lg %.16lg\n",fragInfo->vnlFxCor[iAtom],
	      fragInfo->vnlFyCor[iAtom],fragInfo->vnlFzCor[iAtom]);
    }
  }

  
  if(myidState==0){
    printf("Finish writing fragmentation checkpoint file...\n");
    fclose(fileCheckpoint);
  } 
  free(rhoTemp);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void checkpointFragInput(CP *cp,CLASS *class)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  STODFTINFO *stodftInfo = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  FRAGINFO *fragInfo;
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CPOPTS *cpOpts             = &(cp->cpopts);
  COMMUNICATE *commCP        = &(cp->communicate);

  int iGrid,iAtom;
  int cpLsda                = cpOpts->cp_lsda;
  int myidState             = commCP->myid_state;
  int numProcStates         = commCP->np_states;
  int rhoRealGridNum        = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot        = stodftInfo->rhoRealGridTot;
  int checkpointParFlag     = stodftInfo->checkpointParFlag;
  int checkpointFlag = 0;
  int checkpointFlagAll = 0;
  int numAtomTot            = clatoms_info->natm_tot;
  int readCoeffFlag;
  MPI_Comm commStates       = commCP->comm_states;

  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;
  char fileName[100];

  double junk;
  double *rhoUpFragSum;
  double *rhoDnFragSum;
  double *rhoUpFragSumCpy;
  double *rhoDnFragSumCpy;

  double *rhoTemp = (double*)cmalloc(rhoRealGridTot*sizeof(double));

  FILE *fileCheckpoint;
  FILE *fileCheckpointSys = NULL; //checkpoint file for the whole system

/*======================================================================*/
/* I) Allocation may needed                                             */
  if(myidState==0)fragInfo = stodftInfo->fragInfo;
  else{
    stodftInfo->fragInfo = (FRAGINFO*)cmalloc(sizeof(FRAGINFO));
    fragInfo = stodftInfo->fragInfo;
  }

  fragInfo->rhoUpFragSum = (double*)cmalloc(rhoRealGridNum*sizeof(double));
  fragInfo->rhoUpFragSumCpy = (double*)cmalloc(rhoRealGridNum*sizeof(double));
  rhoUpFragSum = fragInfo->rhoUpFragSum;
  rhoUpFragSumCpy = fragInfo->rhoUpFragSumCpy;
  if(cpLsda==1){
    fragInfo->rhoDnFragSum = (double*)cmalloc(rhoRealGridNum*sizeof(double));
    fragInfo->rhoDnFragSumCpy = (double*)cmalloc(rhoRealGridNum*sizeof(double));
    rhoDnFragSum = fragInfo->rhoDnFragSum;
    rhoDnFragSumCpy = fragInfo->rhoDnFragSumCpy;
  }
  fragInfo->vnlFxCor = (double*)cmalloc(numAtomTot*sizeof(double));
  fragInfo->vnlFyCor = (double*)cmalloc(numAtomTot*sizeof(double));
  fragInfo->vnlFzCor = (double*)cmalloc(numAtomTot*sizeof(double));

/*======================================================================*/
/* II) Read density part                                                */

  if(myidState==0){
    // check the aviability of checkpoint file
    printf("Start reading fragmentation checkpoint file...\n");
    fileCheckpoint = fopen("frag-checkpoint","r");
    if(fileCheckpoint==NULL){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("No checkpoint file. Please double check!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }
    fscanf(fileCheckpoint,"%lg",&stodftInfo->numElecTrueFrag);
    //stodftInfo->numElecTrue = stodftInfo->numElecTrueFrag;
    for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
      fscanf(fileCheckpoint,"%lg",&rhoTemp[iGrid]);
    }
  }//endif myidState
  Barrier(commStates);
  if(numProcStates>1){
    Barrier(commStates);
    Bcast(&stodftInfo->numElecTrueFrag,1,MPI_DOUBLE,0,commStates);
    Scatterv(rhoTemp,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	     rhoUpFragSum,rhoRealGridNum,MPI_DOUBLE,0,commStates);
    Barrier(commStates);
  }
  else{
    memcpy(rhoUpFragSum,rhoTemp,rhoRealGridTot*sizeof(double));
  }
  if(cpLsda==1){
    if(myidState==0){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	fscanf(fileCheckpoint,"%lg",&rhoTemp[iGrid]);
      }
    }
    if(numProcStates>1){
      Barrier(commStates);
      Scatterv(rhoTemp,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
	       rhoDnFragSum,rhoRealGridNum,MPI_DOUBLE,0,commStates);
      Barrier(commStates);
    }
    else{
      memcpy(rhoUpFragSum,rhoTemp,rhoRealGridTot*sizeof(double));
    }
  }
  stodftInfo->numElecTrue = stodftInfo->numElecTrueFrag;

/*======================================================================*/
/* II) Read initial density if necessary                                */

  if(checkpointParFlag==0){//sequential
    if(myidState==0){
      fileCheckpointSys = fopen("density-checkpoint","r");
      // I only have fragment checkpoint file. I need to remove the checkpoint 
      // file reading flag so that the system stochastic dft calculation 
      // is normal.
      if(fileCheckpointSys==NULL){
        printf("Can not find density-checkpoint. I'll read the initial density.\n");
        stodftInfo->readCoeffFlag = -1;
      }//endif fileCheckpointSys
      else{
        printf("Find system density-checkpoint file. I'll read that checkpoint file.\n");
      }
    }//endif myidState
  }//endif checkpointParFlag
  else{//parallel, make sure all processes can open its own checkpoint file
    sprintf(fileName,"density-checkpoint-%i",myidState);
    fileCheckpointSys = fopen(fileName,"r");
    if(fileCheckpointSys!=NULL){
      checkpointFlag = 1;
    }
    if(numProcStates>1){
      Barrier(commStates);
      Reduce(&checkpointFlag,&checkpointFlagAll,1,MPI_INT,MPI_SUM,0,commStates);
      Barrier(commStates);
    }
    else checkpointFlagAll = checkpointFlag;
    if(myidState==0){
      if(checkpointFlagAll==numProcStates){
        printf("Find system density-checkpoint files. I'll read those checkpoint files.\n");
      }
      else{
        printf("Can not file density-checkpoint files. There are total %i processes\n",numProcStates);
        printf("and total %i checkpoint files. I'll read the initial density instead!\n",checkpointFlagAll);
      }
    }
  }
  Bcast(&(stodftInfo->readCoeffFlag),1,MPI_INT,0,commStates);
  readCoeffFlag = stodftInfo->readCoeffFlag;
  if(readCoeffFlag==-1){
    if(myidState==0){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
        fscanf(fileCheckpoint,"%lg",&rhoTemp[iGrid]);
      }
    }
    if(numProcStates>1){
      Barrier(commStates);
      Scatterv(rhoTemp,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
  	       rhoUpFragSumCpy,rhoRealGridNum,MPI_DOUBLE,0,commStates);
      Barrier(commStates);
    }
    else{
      memcpy(rhoUpFragSumCpy,rhoTemp,rhoRealGridTot*sizeof(double));
    }
    if(cpLsda==1){
      if(myidState==0){
	for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	  fscanf(fileCheckpoint,"%lg",&rhoTemp[iGrid]);
	}
      }
      if(numProcStates>1){
        Barrier(commStates);
	Scatterv(rhoTemp,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
		 rhoDnFragSumCpy,rhoRealGridNum,MPI_DOUBLE,0,commStates);      
        Barrier(commStates);
      }
      else{
	memcpy(rhoDnFragSumCpy,rhoTemp,rhoRealGridTot*sizeof(double));
      }
    } 
    if(cpLsda==1)stodftInfo->occNumber = 1;
    else stodftInfo->occNumber = 2;
  }
  else{
    if(myidState==0){
      for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	fscanf(fileCheckpoint,"%lg",&junk);
      }
      if(cpLsda==1){
	for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
	  fscanf(fileCheckpoint,"%lg",&junk);
	}
      }
    }
  }


 
/*======================================================================*/
/* II) Read energy and force part                                       */

  if(myidState==0){
    fscanf(fileCheckpoint,"%lg",&fragInfo->keCor);
    fscanf(fileCheckpoint,"%lg",&fragInfo->vnlCor);
    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      fscanf(fileCheckpoint,"%lg",&fragInfo->vnlFxCor[iAtom]);
      fscanf(fileCheckpoint,"%lg",&fragInfo->vnlFyCor[iAtom]);
      fscanf(fileCheckpoint,"%lg",&fragInfo->vnlFzCor[iAtom]);
    }
    
    fclose(fileCheckpoint);
  }
    
  free(rhoTemp);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

