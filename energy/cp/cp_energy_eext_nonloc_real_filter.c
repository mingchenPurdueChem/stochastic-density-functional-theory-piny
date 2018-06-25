/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: cp-energy-nlpp-real.c                          */
/*                                                                          */
/* This routine wrapps all functions used within SCF. Nuclei forces are not */
/* calculated.                                                              */
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
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_stodft_local.h"

#include "complex.h"

//#define REAL_PP_DEBUG

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void controlPotentialEnergyFilter(CP *cp,CLASS *class,GENERAL_DATA *generalData,
                           double *v_ks,double *zfft,int flag, 
			   PARA_FFT_PKG3D *cpParaFftPkg3d)
/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
  PSEUDO *pseudo = &(cp->pseudo);
  CPOPTS *cpopts = &(cp->cpopts);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CPEWALD *cpewald = &(cp->cpewald);

  //PARA_FFT_PKG3D *cpParaFftPkg3dLgBigBox = &(cp->cp_para_fft_pkg3d_lg);
  //PARA_FFT_PKG3D *cpParaFftPkg3dSparse = &(cp->cp_sclr_fft_pkg3d_sparse);
  int realSparseOpt = cpewald->realSparseOpt;
  //int nfft,numGrid;
  int nfft = cpParaFftPkg3d->nfft;
  int numGrid = nfft/2;
  int iGrid;
  int forceCalcFlag = pseudoReal->forceCalcFlag;
  int threadFlag = cpopts->threadFlag;
  double *forceRealNlpp = pseudoReal->forceRealNlpp;
  //timing
  STODFTINFO *stodftInfo = cp->stodftInfo;
  double time_st,time_end;

  //printf("forceCalcFlag %i\n",forceCalcFlag);

  time_st = omp_get_wtime();
  #pragma omp parallel for private(iGrid)
  for(iGrid=0;iGrid<numGrid;iGrid++){
    forceRealNlpp[iGrid] = 0.0;
  }
  time_end = omp_get_wtime();
  stodftInfo->cputime_new[1] += time_end-time_st;  

  nlppKBRealEnergyFilter(cp,class,generalData,zfft,cpParaFftPkg3d,0);

  time_st = omp_get_wtime();  
  #pragma omp parallel for private(iGrid)
  for(iGrid=0;iGrid<numGrid;iGrid++){
    zfft[iGrid*2+1] = zfft[iGrid*2+1]*v_ks[iGrid+1]+forceRealNlpp[iGrid];
  }
  time_end = omp_get_wtime();
  stodftInfo->cputime_new[2] += time_end-time_st;

  if(flag==1){//double
    time_st = omp_get_wtime();
    #pragma omp parallel for private(iGrid)
    for(iGrid=0;iGrid<numGrid;iGrid++){
      forceRealNlpp[iGrid] = 0.0;
    }
    time_end = omp_get_wtime();
    stodftInfo->cputime_new[3] += time_end-time_st;

    nlppKBRealEnergyFilter(cp,class,generalData,zfft,cpParaFftPkg3d,1);
  
    time_st = omp_get_wtime();
    #pragma omp parallel for private(iGrid)
    for(iGrid=0;iGrid<numGrid;iGrid++){
      zfft[iGrid*2+2] = zfft[iGrid*2+2]*v_ks[iGrid+1]+forceRealNlpp[iGrid];
    }
    time_end = omp_get_wtime();
    stodftInfo->cputime_new[4] += time_end-time_st;
  }//endif

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void nlppKBRealEnergyFilter(CP *cp,CLASS *class,GENERAL_DATA *generalData,
	      double *zfft,PARA_FFT_PKG3D *cpParaFftPkg3d,int flag)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Real space nlpp, only used in filtering		 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  PSEUDO *pseudo = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CELL *cell = &(generalData->cell);
  CLATOMS_POS *clatoms_pos = &(class->clatoms_pos[1]);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  ATOMMAPS *atommaps = &(class->atommaps);
  STAT_AVG *stat_avg = &(generalData->stat_avg);
  COMMUNICATE *communicate = &(cp->communicate);
  CPEWALD *cpewald = &(cp->cpewald);

  int numAtomType = atommaps->natm_typ;
  int numAtom = clatoms_info->natm_tot;
  int numGrid;
  int numGridMax = pseudoReal->numGridMax;
  int countRad,countNlppRe,countNlppIm;
  int iPart,iRad,l,m,iType,iAtom,iGrid;
  int radIndex,gridIndex;
  int aIndex,bIndex,cIndex;
  int atomType;
  int realSparseOpt = cpewald->realSparseOpt;
  //int nkc,nkb,nka,numGridTot;
  int nkc = cpParaFftPkg3d->nkf3;
  int nkb = cpParaFftPkg3d->nkf2;
  int nka = cpParaFftPkg3d->nkf1;
  int numGridTot = nka*nkb*nkc;
  int div,res;
  int gridShiftNowRe,gridShiftNowIm;
  int energyCalcFlag = pseudoReal->energyCalcFlag;
  int numThreads = communicate->numThreads;
  int iThread;

  int *gridMap;
  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *gridStIndRe = pseudoReal->gridStIndRe;
  int *gridStIndIm = pseudoReal->gridStIndIm;
  int *atomLMax = pseudoReal->numLMax; //max L for each atom
  int **atomLRadNum = pseudoReal->atomLRadNum; //num of radical functions for each atom and each L
  int **atomRadMap = pseudoReal->atomRadMap;  //map of radical functions for each atom, starting from l=0 to l=lmax
  int *numGridNlppMap = pseudoReal->numGridNlppMap;
  int **gridNlppMap = pseudoReal->gridNlppMap;

  double vpsNorm;
  double vol,volInv;
  double volElem;
  double dotRe,dotIm;
  double energy = 0.0;
  double energyl;

  double *trig; // sin(theta),cos(theta),sin(phi),cos(phi)
  double *vpsNormList = pseudoReal->vpsNormList; // 1/vpsNorm 
  double *gridAtomNbhd;
  double *ylm;
  double a[3],b[3],c[3];
  double aGrid[3],bGrid[3],cGrid[3];
  double nucleiCoord[3];
  double *hmat = cell->hmat;
  double *x = clatoms_pos->x;
  double *y = clatoms_pos->y;
  double *z = clatoms_pos->z;
  double *wfNbhd = pseudoReal->wfNbhd;
  double *radFun;
  double *forceTemp = pseudoReal->forceTemp;
  double *forceRealNlpp = pseudoReal->forceRealNlpp;
  double *pseudoFunTemp;
  double *vnlPhiAtomGridRe = pseudoReal->vnlPhiAtomGridRe;
  double *vnlPhiAtomGridIm = pseudoReal->vnlPhiAtomGridIm;

  double **dotReAll = pseudoReal->dotReAll;
  double **dotImAll = pseudoReal->dotImAll;

  //Test multithread
  STODFTINFO *stodftInfo = cp->stodftInfo;
  double time_st,time_end;

/*======================================================================*/
/* I) Allocate local memory                                             */

  time_st = omp_get_wtime();
  /*
  numGridMax = numGridNlppMap[0];
  for(iPart=0;iPart<numAtom;iPart++){
    if(numGridNlppMap[iPart]>numGridMax)numGridMax = numGridNlppMap[iPart];
  }
  //printf("numThreads %i iThread %i\n",numThreads,iThread);

  forceTemp = (double*)cmalloc(numAtom*numGridMax*sizeof(double));
  wfNbhd = (double*)cmalloc(numThreads*numGridMax*sizeof(double));
  */
  
  vol = getdeth(hmat);
  volInv = 1.0/vol;
  volElem = vol/numGridTot;
  time_end = omp_get_wtime();
  stodftInfo->cputime_new[5] += time_end-time_st;
  

/*======================================================================*/
/* II) Loop over iPart/irad/ipart/m                                 */

  time_st = omp_get_wtime();
  #pragma omp parallel for private(iGrid)
  for(iGrid=0;iGrid<numAtom*numGridMax;iGrid++){
    forceTemp[iGrid] = 0.0;
  }
  time_end = omp_get_wtime();
  stodftInfo->cputime_new[6] += time_end-time_st;

  //cputime(&time_st);
  time_st = omp_get_wtime();
  omp_set_num_threads(numThreads);
  #pragma omp parallel private(iThread,iAtom,atomType,numGrid,countRad,countNlppRe,countNlppIm,gridShiftNowRe,gridShiftNowIm,iGrid,gridIndex,l,iRad,radIndex,energyl,m,dotRe,dotIm)
  {
    iThread = omp_get_thread_num();
    #pragma omp for
    for(iAtom=0;iAtom<numAtom;iAtom++){
      atomType = iAtomAtomType[iAtom+1]-1;
      numGrid = numGridNlppMap[iAtom];
      countRad = 0;
      countNlppRe = 0;
      countNlppIm = 0;
      gridShiftNowRe = gridStIndRe[iAtom];
      gridShiftNowIm = gridStIndIm[iAtom];
      /* cpy the wave function */
      if(numGrid>0){ //if numGrid=0, only local pp will be calculated
        if(flag==0){
	  for(iGrid=0;iGrid<numGrid;iGrid++){
	    gridIndex = gridNlppMap[iAtom][iGrid];
	    wfNbhd[iThread*numGridMax+iGrid] = zfft[2*gridIndex+1];
	  }
        }
        else{
          for(iGrid=0;iGrid<numGrid;iGrid++){
            gridIndex = gridNlppMap[iAtom][iGrid];
            wfNbhd[iThread*numGridMax+iGrid] = zfft[2*gridIndex+2];
          }
	}
	for(l=0;l<atomLMax[atomType];l++){
	  for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
	    radIndex = atomRadMap[atomType][countRad+iRad];
	    energyl = 0.0;
	    for(m=0;m<=l;m++){
	      if(m!=0){
		dotRe = ddotBlasWrapper(numGrid,&wfNbhd[iThread*numGridMax],1,
					&vnlPhiAtomGridRe[gridShiftNowRe],1)*volElem;
		dotIm = ddotBlasWrapper(numGrid,&wfNbhd[iThread*numGridMax],1,
					&vnlPhiAtomGridIm[gridShiftNowIm],1)*volElem;
		dotRe *= 2.0*vpsNormList[radIndex];
		dotIm *= 2.0*vpsNormList[radIndex];
		daxpyBlasWrapper(numGrid,dotRe,
				&vnlPhiAtomGridRe[gridShiftNowRe],1,
				&forceTemp[iAtom*numGridMax],1);
		daxpyBlasWrapper(numGrid,dotIm,
				 &vnlPhiAtomGridIm[gridShiftNowIm],1,
				 &forceTemp[iAtom*numGridMax],1);
		gridShiftNowRe += numGrid;
		gridShiftNowIm += numGrid;
	      }
	      else{
		dotRe = ddotBlasWrapper(numGrid,&wfNbhd[iThread*numGridMax],1,
					&vnlPhiAtomGridRe[gridShiftNowRe],1)*volElem;
		dotRe *= vpsNormList[radIndex];
		daxpyBlasWrapper(numGrid,dotRe,
				&vnlPhiAtomGridRe[gridShiftNowRe],1,
				&forceTemp[iAtom*numGridMax],1);
		gridShiftNowRe += numGrid;
	      }//endif m
	    }//endfor m
	    countNlppRe += l+1;
	    countNlppIm += l;
	  }//endfor iRad
	  countRad += atomLRadNum[atomType][l];
	}//endfor l
      }//endif numGrid
    }//endfor iAtom
  }//end omp
  time_end = omp_get_wtime();
  stodftInfo->cputime0 += time_end-time_st;

  int numUpdate,atomInd,gridInd2;
  int numGridNlppAll = pseudoReal->numGridNlppAll;
  int *gridNlppInd = pseudoReal->gridNlppInd;
  int **gridNlppMapInv = pseudoReal->gridNlppMapInv;

  /*
  for(iAtom=0;iAtom<numAtom;iAtom++){
    numGrid = numGridNlppMap[iAtom];
    if(numGrid>0){
      for(iGrid=0;iGrid<numGrid;iGrid++){
    gridIndex = gridNlppMap[iAtom][iGrid];
    forceRealNlpp[gridIndex] += forceTemp[iAtom*numGridMax+iGrid];
      }
    }
  }
  */

  
  //cputime(&time_st);
  time_st = omp_get_wtime();
  omp_set_num_threads(numThreads);
  #pragma omp parallel for private(iGrid,gridIndex,numUpdate,iAtom,atomInd,gridInd2)
  for(iGrid=0;iGrid<numGridNlppAll;iGrid++){
    gridIndex = gridNlppInd[iGrid];
    numUpdate = gridNlppMapInv[iGrid][0];
    for(iAtom=0;iAtom<numUpdate;iAtom++){
      atomInd = gridNlppMapInv[iGrid][iAtom*2+1];
      gridInd2 = gridNlppMapInv[iGrid][iAtom*2+2];
      forceRealNlpp[gridIndex] += forceTemp[atomInd*numGridMax+gridInd2];
    }
  }
  //cputime(&time_end);
  time_end = omp_get_wtime();
  stodftInfo->cputime1 += time_end-time_st;

/*======================================================================*/
/* III) free local memory                                               */

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/

