/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: energy-cor.c                                 */
/*                                                                          */
/* This routine calculate the kinectic energy, non-local pseudo-potential   */
/* energy and nuclei force correction.					    */
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
void energyCorrect(CP *cpMini,GENERAL_DATA *generalDataMini,CLASS *classMini,
		   CP *cp,CLASS *class,int ip_now)
/*========================================================================*/
{/*begin routine*/
/**************************************************************************/
/* This function calculate the energy correction, including kinetic	  */
/* energy, non-local pseudo potential energy and nuclei forces comes from */
/* non-local pseudo potential.						  */
/**************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  CPOPTS *cpOpts = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cp->cpcoeffs_info);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  COMMUNICATE *commCP = &(cp->communicate);
  PSEUDO *pseudo = &(cpMini->pseudo);

  int myidState             = commCP->myid_state;
  int numProcStates         = commCP->np_states;
  int numFragProc	    = fragInfo->numFragProc;
  int numAtomTot	    = clatoms_info->natm_tot;
  int iFrag,iAtom;
  int vnl_kb_flag = pseudo->vnl_kb_flag;
  int vnl_gh_flag = pseudo->vnl_gh_flag;
  double keCorProc = 0.0;
  double vnlCorProc = 0.0;
  double *vnlFxCorProc,*vnlFyCorProc,*vnlFzCorProc;
  MPI_Comm commStates = commCP->comm_states;

  vnlFxCorProc = (double*)cmalloc(numAtomTot*sizeof(double));
  vnlFyCorProc = (double*)cmalloc(numAtomTot*sizeof(double));
  vnlFzCorProc = (double*)cmalloc(numAtomTot*sizeof(double));
 

  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInfo->iFrag = iFrag;

/*======================================================================*/
/* I) Kinetic energy	                                                */
    
    calcKECor(cpMini,generalDataMini,cp,&keCorProc);

/*======================================================================*/
/* II) Non-local pseudo potential energy and force                      */

    printf("vnl_kb_flag %i\n",vnl_kb_flag);
    
    if(vnl_kb_flag==1){
      calcVnlCor(classMini,cpMini,generalDataMini,
		 cp,class,&vnlCorProc,vnlFxCorProc,
		 vnlFyCorProc,vnlFzCorProc);
    }
    

  }//endfor

/*======================================================================*/
/* I) Reduce everything                                                 */

  if(numProcStates>1){
    Reduce(&keCorProc,&(fragInfo->keCor),1,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&vnlCorProc,&(fragInfo->vnlCor),1,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&vnlFxCorProc[0],&(fragInfo->vnlFxCor[0]),numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&vnlFyCorProc[0],&(fragInfo->vnlFyCor[0]),numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
    Reduce(&vnlFzCorProc[0],&(fragInfo->vnlFzCor[0]),numAtomTot,MPI_DOUBLE,MPI_SUM,0,commStates);
  }
  else{
    fragInfo->keCor = keCorProc;
    fragInfo->vnlCor = vnlCorProc;
    memcpy(&(fragInfo->vnlFxCor[0]),&(vnlFxCorProc[0]),numAtomTot*sizeof(double));
    memcpy(&(fragInfo->vnlFyCor[0]),&(vnlFyCorProc[0]),numAtomTot*sizeof(double));
    memcpy(&(fragInfo->vnlFzCor[0]),&(vnlFzCorProc[0]),numAtomTot*sizeof(double));
  }
  
  if(numProcStates>1)Barrier(commStates);
  free(vnlFxCorProc);
  free(vnlFyCorProc);
  free(vnlFzCorProc);
  fflush(stdout);
  //exit(0);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcKECor(CP *cpMini,GENERAL_DATA *generalDataMini,CP *cp,double *keCorProc)
/*========================================================================*/
{/*begin routine*/
/**************************************************************************/
/* This function calculates the kinetic energy correction. \sum_f ke_f    */
/* -\sum_f a_f^T B_f a_f where ke_f = \sum_i <psi_i^f|K|psi_i^f> , a_f(i) */
/* = <kai|psi_f^i> , and B_f(i,j)=<psi_f^i|K|psi_f^j>.			  */
/**************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  CPOPTS *cpOpts = &(cpMini->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cpMini->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos = &(cpMini->cpcoeffs_pos[1]);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cpMini->cp_para_fft_pkg3d_lg);
  COMMUNICATE *commCP = &(cpMini->communicate);
  STAT_AVG *statAvg = &(generalDataMini->stat_avg);
  

  int iState,jState,iCoeff,iStoc;
  int iFrag = fragInfo->iFrag;
  int cpLsda = cpOpts->cp_lsda;
  int numFragProc           = fragInfo->numFragProc;
  int numFragTot            = fragInfo->numFragTot;
  int numStateUp = cpcoeffs_info->nstate_up_proc;
  int numStateDn = cpcoeffs_info->nstate_dn_proc;
  int numStateStoUp = stodftInfo->numStateStoUp;
  int numStateStoDn = stodftInfo->numStateStoDn;
  int occNumber = stodftInfo->occNumber;
  double *keMatrixUp,*keMatrixDn;
  double *wfProjUp,*wfProjDn;
  double *temp;
  double keCor,keCorUp,keCorDn;
  double ke = statAvg->cp_eke;

/*======================================================================*/
/* I) Calculate the matrix                                              */

  calcKEMatrix(cpMini,cp);

/*======================================================================*/
/* I) Allocate Local Memory                                             */

  //printf("ke %lg\n",ke);
  keMatrixUp = fragInfo->keMatrixUp[iFrag];
  wfProjUp = fragInfo->wfProjUp[iFrag];
  //debug
  /*
  for(iState=0;iState<numStateUp;iState++){
    for(jState=iState;jState<numStateUp;jState++){
      printf("iState %i jState %i Matrix %lg\n",iState,jState,keMatrixUp[iState*numStateUp+jState]);
    }
  }
  for(iStoc=0;iStoc<numStateStoUp;iStoc++){
    for(iState=0;iState<numStateUp;iState++){
      printf("iStoch %i iState %i wfProj %lg\n",iStoc,iState,wfProjUp[iStoc*numStateUp+iState]);
    }
  }
  */
  temp = (double*)cmalloc(numStateUp*sizeof(double));
  keCorUp = 0.0;
  for(iStoc=0;iStoc<numStateStoUp;iStoc++){
    dsymvWrapper('U',numStateUp,1.0,keMatrixUp,numStateUp,&wfProjUp[iStoc*numStateUp],1,0.0,temp,1);
    keCorUp += ddotBlasWrapper(numStateUp,temp,1,&wfProjUp[iStoc*numStateUp],1);
  }
  keCorUp /= numStateStoUp;
  keCor = keCorUp;
  free(temp);
  if(cpLsda==1&&numStateDn!=0){
    keMatrixDn = fragInfo->keMatrixDn[iFrag];
    wfProjDn = fragInfo->wfProjDn[iFrag];
    temp = (double*)cmalloc(numStateDn*sizeof(double));
    for(iStoc=0;iStoc<numStateStoDn;iStoc++){
      dsymvWrapper('U',numStateDn,1.0,keMatrixDn,numStateDn,&wfProjDn[iStoc*numStateDn],1,0.0,temp,1);
      keCorDn += ddotBlasWrapper(numStateDn,temp,1,&wfProjDn[iStoc*numStateDn],1);
    }
    keCorDn /= numStateStoDn;
    keCor += keCorDn;
    free(temp);
  }
  printf("ke %lg keCor %lg\n",ke,keCor);
  *keCorProc += ke-occNumber*keCor;

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcKEMatrix(CP *cpMini,CP *cp)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;  
  CPOPTS *cpOpts = &(cpMini->cpopts);
  CPEWALD *cpEwald = &(cpMini->cpewald);
  CPCOEFFS_INFO *cpcoeffs_info = &(cpMini->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos = &(cpMini->cpcoeffs_pos[1]);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cpMini->cp_para_fft_pkg3d_lg);
  COMMUNICATE *commCP = &(cpMini->communicate);

  int iState,jState,iCoeff;
  int numStateUp = cpcoeffs_info->nstate_up_proc;
  int numStateDn = cpcoeffs_info->nstate_dn_proc;
  int numCoeff = cpcoeffs_info->ncoef;
  int cpLsda = cpOpts->cp_lsda;
  int numCoeffUpTot = numStateUp*numCoeff;
  int numCoeffDnTot = numStateDn*numCoeff;
  int numAlloc = MAX(numCoeffUpTot,numCoeffDnTot);
  int index,index1,index2,index3;
  int iFrag = fragInfo->iFrag;
  
  double *cre_up = cpcoeffs_pos->cre_up;
  double *cim_up = cpcoeffs_pos->cim_up;
  double *cre_dn = cpcoeffs_pos->cre_dn;
  double *cim_dn = cpcoeffs_pos->cim_dn;
  double *ak2Small = cpEwald->ak2_sm;
  double *coefForceRe,*coefForceIm;
  double *keMatrixUp;
  double *keMatrixDn;
  double tpi = 2.0*M_PI;

/*======================================================================*/
/* I) Allocate Local Memory                                             */


  coefForceRe = (double*)cmalloc((numAlloc+1)*sizeof(double));
  coefForceIm = (double*)cmalloc((numAlloc+1)*sizeof(double));


/*======================================================================*/
/* II) Calculate Spin up matrix                                         */

  keMatrixUp = fragInfo->keMatrixUp[iFrag];

  for(iState=0;iState<numStateUp;iState++){
    for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
      index = iState*numCoeff+iCoeff;
      coefForceRe[index] = 0.5*ak2Small[iCoeff]*cre_up[index];
      coefForceIm[index] = 0.5*ak2Small[iCoeff]*cim_up[index];
    }//endfor iCoeff
    coefForceRe[iState*numCoeff+numCoeff] = 0.0;
    coefForceIm[iState*numCoeff+numCoeff] = 0.0;
  }//endfor iState

  for(iState=0;iState<numStateUp;iState++){
    for(jState=iState;jState<numStateUp;jState++){
      index = iState*numStateUp+jState;
      index1 = jState*numStateUp+iState;
      keMatrixUp[index] = 0.0;
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	index2 = iState*numCoeff+iCoeff;
	index3 = jState*numCoeff+iCoeff;
        // We should have *2 here since we have CC* and C*C but the wave functions
	// are normalized to 2, so there is another 0.5 to bring it back to normal 
	// But don't forget to scale everything by occupied number at the end of the day
	keMatrixUp[index] += coefForceRe[index2]*cre_up[index3]+coefForceIm[index2]*cim_up[index3];
      }//endfor iCoeff
      keMatrixUp[index1] = keMatrixUp[index];
    }//endfor jState
  }//endfor iState
  
/*======================================================================*/
/* III) Calculate Spin down matrix					*/


  if(cpLsda==1&&numStateDn!=0){// spin down
    keMatrixDn = fragInfo->keMatrixDn[iFrag];
    for(iState=0;iState<numStateDn;iState++){
      for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	index = iState*numCoeff+iCoeff;
	coefForceRe[index] = 0.5*ak2Small[iCoeff]*cre_dn[index];
	coefForceIm[index] = 0.5*ak2Small[iCoeff]*cim_dn[index];
      }
      coefForceRe[iState*numCoeff+numCoeff] = 0.0;
      coefForceIm[iState*numCoeff+numCoeff] = 0.0;
    }

    for(iState=0;iState<numStateDn;iState++){
      for(jState=iState;jState<numStateDn;jState++){
	index = iState*numStateDn+jState;
	index1 = jState*numStateDn+iState;
	keMatrixDn[index] = 0.0;
	for(iCoeff=1;iCoeff<numCoeff;iCoeff++){
	  index2 = iState*numCoeff+iCoeff;
	  index3 = jState*numCoeff+iCoeff;
	  // We should have *2 here since we have CC* and C*C but the wave functions
	  // are normalized to 2, so there is another 0.5 to bring it back to normal 
	  // But don't forget to scale everything by occupied number at the end of the day
	  keMatrixDn[index] += coefForceRe[index2]*cre_dn[index3]+coefForceIm[index2]*cim_dn[index3];
	}//endfor iCoeff
	keMatrixDn[index1] = keMatrixDn[index];
      }//endfor jState
    }//endfor iState   
  }//endif 


  /*
  //eke = 0.0;
  for(is=1 ; is<= nstate ; is++){
    ioff = (is-1)*ncoef;
    for(i=1; i<= ncoef1 ; i++){
      iis = ioff + i;
      fccreal[iis] -= 2.0*ak2_sm[i]*ccreal[iis];
      fccimag[iis] -= 2.0*ak2_sm[i]*ccimag[iis];
      //eke += (2.0*ak2_sm[i]*(ccreal[iis]*ccreal[iis] + ccimag[iis]*ccimag[iis]));
    }//endfor i
   nis = is*ncoef;
   fccimag[nis] = 0.0;
  }//endfor
  */

/*======================================================================*/
/* IV) free local memories						*/

  free(coefForceRe);
  free(coefForceIm);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcVnlCor(CLASS *classMini, CP *cpMini,GENERAL_DATA *generalDataMini,
		CP *cp,CLASS *class,double *vnlCorProc,double *vnlFxCorProc,
		double *vnlFyCorProc,double *vnlFzCorProc)
/*========================================================================*/
{/*begin routine*/
/**************************************************************************/
/* This function use the same way to correct the non-local pseudo-	  */
/* potential term							  */
/**************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  CPOPTS *cpOpts = &(cpMini->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cpMini->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos = &(cpMini->cpcoeffs_pos[1]);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cpMini->cp_para_fft_pkg3d_lg);
  COMMUNICATE *commCP = &(cpMini->communicate);
  STAT_AVG *statAvg = &(generalDataMini->stat_avg);
  CLATOMS_INFO *clatomsInfoMini = &(classMini->clatoms_info);
  CLATOMS_INFO *clatomsInfo = &(class->clatoms_info);
  CLATOMS_POS *clatomsPosMini = &(classMini->clatoms_pos[1]);
  

  int iState,jState,iCoeff,iStoc,iAtom;
  int iFrag = fragInfo->iFrag;
  int cpLsda = cpOpts->cp_lsda;
  int numFragProc           = fragInfo->numFragProc;
  int numFragTot            = fragInfo->numFragTot;
  int numStateUp = cpcoeffs_info->nstate_up_proc;
  int numStateDn = cpcoeffs_info->nstate_dn_proc;
  int numStateStoUp = stodftInfo->numStateStoUp;
  int numStateStoDn = stodftInfo->numStateStoDn;
  int occNumber = stodftInfo->occNumber;
  int numAtomFrag = clatomsInfoMini->natm_tot;
  int *atomFragMapProc = fragInfo->atomFragMapProc[iFrag];

  double vnlCor,vnlCorUp,vnlCorDn;
  double vnlFxCorTemp,vnlFyCorTemp,vnlFzCorTemp;
  double vnl = statAvg->cp_enl;

  double *vnlMatrixUp,*vnlMatrixDn;
  double *vnlFxMatrixUp,*vnlFxMatrixDn;
  double *vnlFyMatrixUp,*vnlFyMatrixDn;
  double *vnlFzMatrixUp,*vnlFzMatrixDn;
  double *wfProjUp,*wfProjDn;
  double *temp;
  double *vnlFxCorFragLoc,*vnlFyCorFragLoc,*vnlFzCorFragLoc; //local force correction
  double *Fx = fragInfo->Fx[iFrag];
  double *Fy = fragInfo->Fy[iFrag];
  double *Fz = fragInfo->Fz[iFrag];
  double *fx = clatomsPosMini->fx;
  double *fy = clatomsPosMini->fy;
  double *fz = clatomsPosMini->fz;

  

/*======================================================================*/
/* I) Calculate the matrix                                              */

  calcNonLocalMatrix(cp,cpMini,classMini,generalDataMini);

  for(iAtom=0;iAtom<numAtomFrag;iAtom++){
    Fx[iAtom] = fx[iAtom+1];
<<<<<<< HEAD
    Fy[iAtom] = fx[iAtom+1];
    Fz[iAtom] = fx[iAtom+1];
=======
    Fy[iAtom] = fy[iAtom+1];
    Fz[iAtom] = fz[iAtom+1];
>>>>>>> fragment-onebody-new
  }

/*======================================================================*/
/* I) Calculate vnl/force correction for spin up elections              */
/*    Pay attention to the 0.25 factor					*/
  
  //printf("ke %lg\n",ke);
  vnlMatrixUp = fragInfo->vnlMatrixUp[iFrag];
  wfProjUp = fragInfo->wfProjUp[iFrag];
  //debug
  /*
  for(iState=0;iState<numStateUp;iState++){
    for(jState=iState;jState<numStateUp;jState++){
      printf("iState %i jState %i Matrix %lg\n",iState,jState,keMatrixUp[iState*numStateUp+jState]);
    }
  }
  for(iStoc=0;iStoc<numStateStoUp;iStoc++){
    for(iState=0;iState<numStateUp;iState++){
      printf("iStoch %i iState %i wfProj %lg\n",iStoc,iState,wfProjUp[iStoc*numStateUp+iState]);
    }
  }
  */
  vnlFxCorFragLoc = (double*)cmalloc(numAtomFrag*sizeof(double));
  vnlFyCorFragLoc = (double*)cmalloc(numAtomFrag*sizeof(double));
  vnlFzCorFragLoc = (double*)cmalloc(numAtomFrag*sizeof(double));

  temp = (double*)cmalloc(numStateUp*sizeof(double));
  vnlCorUp = 0.0;

  printf("vnl matrix\n");
  for(iState=0;iState<numStateUp;iState++){
    for(jState=0;jState<numStateUp;jState++){
      printf("%lg ",vnlMatrixUp[iState*numStateUp+jState]);
    }
    printf("\n");
  }
  printf("end vnl matrix\n");
  for(iStoc=0;iStoc<numStateStoUp;iStoc++){
    dsymvWrapper('U',numStateUp,1.0,vnlMatrixUp,numStateUp,&wfProjUp[iStoc*numStateUp],1,0.0,temp,1);
    vnlCorUp += ddotBlasWrapper(numStateUp,temp,1,&wfProjUp[iStoc*numStateUp],1);
  }
  vnlCorUp /= numStateStoUp;
  vnlCor = vnlCorUp;
  // Force from spin up electrons
  for(iAtom=0;iAtom<numAtomFrag;iAtom++){
    vnlFxMatrixUp = &(fragInfo->vnlFxMatrixUp[iFrag][iAtom*numStateUp*numStateUp]);
    vnlFyMatrixUp = &(fragInfo->vnlFyMatrixUp[iFrag][iAtom*numStateUp*numStateUp]);
    vnlFzMatrixUp = &(fragInfo->vnlFzMatrixUp[iFrag][iAtom*numStateUp*numStateUp]);
    //debug
    for(iState=0;iState<numStateUp;iState++){
      for(jState=0;jState<numStateUp;jState++){
	printf("atom %i istate %i jstate %i vnlFxMatrixUp %lg vnlFyMatrixUp %lg vnlFzMatrixUp %lg\n",
		iAtom,iState,jState,vnlFxMatrixUp[iState*numStateUp+jState],
<<<<<<< HEAD
		vnlFzMatrixUp[iState*numStateUp+jState],vnlFzMatrixUp[iState*numStateUp+jState]);
=======
		vnlFyMatrixUp[iState*numStateUp+jState],vnlFzMatrixUp[iState*numStateUp+jState]);
>>>>>>> fragment-onebody-new
      }
    }
    
    vnlFxCorTemp = 0.0;
    vnlFyCorTemp = 0.0;
    vnlFzCorTemp = 0.0;
    for(iStoc=0;iStoc<numStateStoUp;iStoc++){
      dsymvWrapper('U',numStateUp,1.0,vnlFxMatrixUp,numStateUp,&wfProjUp[iStoc*numStateUp],1,0.0,temp,1);
      vnlFxCorTemp += ddotBlasWrapper(numStateUp,temp,1,&wfProjUp[iStoc*numStateUp],1);
      dsymvWrapper('U',numStateUp,1.0,vnlFyMatrixUp,numStateUp,&wfProjUp[iStoc*numStateUp],1,0.0,temp,1);
      vnlFyCorTemp += ddotBlasWrapper(numStateUp,temp,1,&wfProjUp[iStoc*numStateUp],1);
      dsymvWrapper('U',numStateUp,1.0,vnlFzMatrixUp,numStateUp,&wfProjUp[iStoc*numStateUp],1,0.0,temp,1);
      vnlFzCorTemp += ddotBlasWrapper(numStateUp,temp,1,&wfProjUp[iStoc*numStateUp],1);      
    }
    vnlFxCorFragLoc[iAtom] = vnlFxCorTemp/numStateUp;
    vnlFyCorFragLoc[iAtom] = vnlFyCorTemp/numStateUp;
    vnlFzCorFragLoc[iAtom] = vnlFzCorTemp/numStateUp;
  }

  free(temp);
  
/*======================================================================*/
/* I) Calculate vnl/force correction for spin down elections            */

  if(cpLsda==1&&numStateDn!=0){
    vnlMatrixDn = fragInfo->vnlMatrixDn[iFrag];
    wfProjDn = fragInfo->wfProjDn[iFrag];
    temp = (double*)cmalloc(numStateDn*sizeof(double));
    for(iStoc=0;iStoc<numStateStoDn;iStoc++){
      dsymvWrapper('U',numStateDn,1.0,vnlMatrixDn,numStateDn,&wfProjDn[iStoc*numStateDn],1,0.0,temp,1);
      vnlCorDn += ddotBlasWrapper(numStateDn,temp,1,&wfProjDn[iStoc*numStateDn],1);
    }
    vnlCorDn /= numStateStoDn;
    vnlCor += vnlCorDn;
    for(iAtom=0;iAtom<numAtomFrag;iAtom++){
      vnlFxMatrixDn = &(fragInfo->vnlFxMatrixDn[iFrag][iAtom*numStateDn*numStateDn]);
      vnlFyMatrixDn = &(fragInfo->vnlFyMatrixDn[iFrag][iAtom*numStateDn*numStateDn]);
      vnlFzMatrixDn = &(fragInfo->vnlFzMatrixDn[iFrag][iAtom*numStateDn*numStateDn]);
      vnlFxCorTemp = 0.0;
      vnlFyCorTemp = 0.0;
      vnlFzCorTemp = 0.0;
      for(iStoc=0;iStoc<numStateStoDn;iStoc++){
	dsymvWrapper('U',numStateDn,1.0,vnlFxMatrixDn,numStateDn,&wfProjDn[iStoc*numStateUp],1,0.0,temp,1);
	vnlFxCorTemp += ddotBlasWrapper(numStateDn,temp,1,&wfProjDn[iStoc*numStateUp],1);
	dsymvWrapper('U',numStateDn,1.0,vnlFyMatrixDn,numStateDn,&wfProjDn[iStoc*numStateUp],1,0.0,temp,1);
	vnlFyCorTemp += ddotBlasWrapper(numStateDn,temp,1,&wfProjDn[iStoc*numStateUp],1);
	dsymvWrapper('U',numStateDn,1.0,vnlFzMatrixDn,numStateDn,&wfProjDn[iStoc*numStateUp],1,0.0,temp,1);
	vnlFzCorTemp += ddotBlasWrapper(numStateDn,temp,1,&wfProjDn[iStoc*numStateUp],1);
      }
      vnlFxCorFragLoc[iAtom] += vnlFxCorTemp/numStateDn;
      vnlFyCorFragLoc[iAtom] += vnlFyCorTemp/numStateDn;
      vnlFzCorFragLoc[iAtom] += vnlFzCorTemp/numStateDn;
    }
    free(temp);
  }
  printf("vnl %lg vnlCor %lg\n",vnl,vnlCor);
  *vnlCorProc += vnl-vnlCor;
  for(iAtom=0;iAtom<numAtomFrag;iAtom++){
    printf("iAtom %i atomFragMapProc[iAtom] %i Fx[iAtom] %lg vnlFxCorFragLoc[iAtom] %lg\n",
    	    iAtom,atomFragMapProc[iAtom],Fx[iAtom],vnlFxCorFragLoc[iAtom]);
    printf("iAtom %i atomFragMapProc[iAtom] %i Fy[iAtom] %lg vnlFyCorFragLoc[iAtom] %lg\n",
            iAtom,atomFragMapProc[iAtom],Fy[iAtom],vnlFyCorFragLoc[iAtom]);
    printf("iAtom %i atomFragMapProc[iAtom] %i Fz[iAtom] %lg vnlFzCorFragLoc[iAtom] %lg\n",
            iAtom,atomFragMapProc[iAtom],Fz[iAtom],vnlFzCorFragLoc[iAtom]);
    vnlFxCorProc[atomFragMapProc[iAtom]-1] += Fx[iAtom]-vnlFxCorFragLoc[iAtom];
    vnlFyCorProc[atomFragMapProc[iAtom]-1] += Fy[iAtom]-vnlFyCorFragLoc[iAtom];
    vnlFzCorProc[atomFragMapProc[iAtom]-1] += Fz[iAtom]-vnlFzCorFragLoc[iAtom];
  }
  free(vnlFxCorFragLoc);
  free(vnlFyCorFragLoc);
  free(vnlFzCorFragLoc);
  //exit(0);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


