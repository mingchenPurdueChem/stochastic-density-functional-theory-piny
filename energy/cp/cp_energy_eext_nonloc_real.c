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
void controlEnergyNlppReal(CP *cp,CLASS *class,GENERAL_DATA *generalData,
			   double *zfft_tmp,double *zfft,int flag,double *energyNl,
			   double *fx,double *fy,double *fz,
			   PARA_FFT_PKG3D *cpParaFftPkg3d)
/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
  PSEUDO *pseudo = &(cp->pseudo);
  CPOPTS *cpopts = &(cp->cpopts);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  ATOMMAPS *atommaps = &(class->atommaps);
  CPEWALD *cpewald = &(cp->cpewald);

  int realSparseOpt = cpewald->realSparseOpt;
  //int nfft,numGrid;
  int nfft = cpParaFftPkg3d->nfft;
  int numGrid = nfft/2;
  int iGrid,iAtom,atomType;
  int forceCalcFlag = pseudoReal->forceCalcFlag;
  int threadFlag = cpopts->threadFlag;
  int numAtomTot = clatoms_info->natm_tot;
  int *numNlppAtom = pseudoReal->numNlppAtom;
  int *iAtomAtomType = atommaps->iatm_atm_typ;

  double *wfReal,*wfForceReal;
  double **dotReAll,**dotImAll;

  if(forceCalcFlag==1){
    dotReAll = (double**)cmalloc(numAtomTot*sizeof(double*));
    dotImAll = (double**)cmalloc(numAtomTot*sizeof(double*));
    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      atomType = iAtomAtomType[iAtom+1]-1;
      dotReAll[iAtom] = NULL;
      dotImAll[iAtom] = NULL;
      if(numNlppAtom[atomType]>=1){
	//printf("atomType %i %i\n",atomType,numNlppAtom[atomType]);
	dotReAll[iAtom] = (double*)cmalloc(numNlppAtom[atomType]*sizeof(double));
	dotImAll[iAtom] = (double*)cmalloc((numNlppAtom[atomType]-1)*sizeof(double));
      }
    }//endfor iAtom 
  }
  
  wfReal = (double*)cmalloc(numGrid*sizeof(double));
  wfForceReal = (double*)cmalloc(numGrid*sizeof(double));
  //printf("forceCalcFlag %i\n",forceCalcFlag);
  //printf("zzzzzzfft %lg\n",zfft_tmp[5481*2+1]);
  for(iGrid=0;iGrid<numGrid;iGrid++){
    wfReal[iGrid] = zfft_tmp[iGrid*2+1];
    wfForceReal[iGrid] = 0.0;
  }
  nlppKBRealEnergy(cp,class,generalData,wfReal,wfForceReal,energyNl,dotReAll,dotImAll,cpParaFftPkg3d);
  if(forceCalcFlag==1)nlppKBRealEnergyForce(cp,class,generalData,wfReal,fx,fy,fz,dotReAll,dotImAll,cpParaFftPkg3d);
  for(iGrid=0;iGrid<numGrid;iGrid++){
#ifdef REAL_PP_DEBUG  
    zfft[iGrid*2+1] = wfForceReal[iGrid];
    //printf("forceeeeeeeee grid 1 %.16lg\n",wfForceReal[iGrid]);
#else
    zfft[iGrid*2+1] += wfForceReal[iGrid];
#endif
  }
  if(flag==1){//double
    for(iGrid=0;iGrid<numGrid;iGrid++){
      wfReal[iGrid] = zfft_tmp[iGrid*2+2];
      wfForceReal[iGrid] = 0.0;
    }
    nlppKBRealEnergy(cp,class,generalData,wfReal,wfForceReal,energyNl,dotReAll,dotImAll,cpParaFftPkg3d);
    if(forceCalcFlag==1)nlppKBRealEnergyForce(cp,class,generalData,wfReal,fx,fy,fz,dotReAll,dotImAll,cpParaFftPkg3d);
    for(iGrid=0;iGrid<numGrid;iGrid++){
#ifdef REAL_PP_DEBUG  
      zfft[iGrid*2+2] = wfForceReal[iGrid];
      //printf("forceeeeeeeee grid 2 %.16lg\n",wfForceReal[iGrid]);
#else
      zfft[iGrid*2+2] += wfForceReal[iGrid];
#endif
    }
  }//endif
  free(wfReal);
  free(wfForceReal);
  
  if(forceCalcFlag==1){
    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      atomType = iAtomAtomType[iAtom+1]-1;
      if(numNlppAtom[atomType]>=1){
	cfree(dotReAll[iAtom]);
	cfree(dotImAll[iAtom]);
      }//endif
    }//endfor iAtom
    cfree(dotReAll);
    cfree(dotImAll);
  }//endif

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void nlppKBRealEnergy(CP *cp,CLASS *class,GENERAL_DATA *generalData,
		      double *wfReal,double *forceRealNlpp,double *energyNl,
		      double **dotReAll,double **dotImAll,
		      PARA_FFT_PKG3D *cpParaFftPkg3d)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Real space nlpp, only used in filtering				 */
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

  int realSparseOpt = cpewald->realSparseOpt;
  int numAtomType = atommaps->natm_typ;
  int numAtom = clatoms_info->natm_tot;
  int numGrid;
  int numGridMax;
  int countRad,countNlppRe,countNlppIm;
  int iPart,iRad,l,m,iType,iAtom,iGrid;
  int radIndex,gridIndex;
  int aIndex,bIndex,cIndex;
  int atomType;
  //int nkc,nkb,nka,numGridTot;
  int nkc = cpParaFftPkg3d->nkf3;
  int nkb = cpParaFftPkg3d->nkf2;
  int nka = cpParaFftPkg3d->nkf1;
  int numGridTot = nka*nkb*nkc;
  int div,res;
  int gridShiftNowRe,gridShiftNowIm;
  int energyCalcFlag = pseudoReal->energyCalcFlag;
  int numThreads = communicate->numThreads;
  int forceCalcFlag = pseudoReal->forceCalcFlag;

  int *gridMap;
  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *gridStIndRe = pseudoReal->gridStIndRe;
  int *gridStIndIm = pseudoReal->gridStIndIm;
  int *atomLMax = pseudoReal->numLMax; //max L for each atom
  int **atomLRadNum = pseudoReal->atomLRadNum; //num of radical functions for each atom and each L
  int **atomRadMap = pseudoReal->atomRadMap;  //map of radical functions for each atom, starting from l=0 to l=lmax
  int *numGridNlppMap = pseudoReal->numGridNlppMap;
  int *locOpt = pseudo->loc_opt;
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
  double *wfNbhd;
  double *radFun;
  double *forceTemp;
  double *pseudoFunTemp;
  double *vnlPhiAtomGridRe = pseudoReal->vnlPhiAtomGridRe;
  double *vnlPhiAtomGridIm = pseudoReal->vnlPhiAtomGridIm;

/*======================================================================*/
/* I) Allocate local memory                                             */

  numGridMax = numGridNlppMap[0];
  for(iPart=0;iPart<numAtom;iPart++){
    if(numGridNlppMap[iPart]>numGridMax)numGridMax = numGridNlppMap[iPart];
  }
  //printf("numThreads %i iThread %i\n",numThreads,iThread);
  forceTemp = (double*)calloc(numAtom*numGridMax,sizeof(double));
  wfNbhd = (double*)cmalloc(numGridMax*sizeof(double));
  
  vol = getdeth(hmat);
  volInv = 1.0/vol;
  volElem = vol/numGridTot;

/*======================================================================*/
/* II) Loop over iPart/irad/ipart/m                                 */

  for(iGrid=0;iGrid<numAtom*numGridMax;iGrid++){
    forceTemp[iGrid] = 0.0;
  }

  for(iAtom=0;iAtom<numAtom;iAtom++){
    atomType = iAtomAtomType[iAtom+1]-1;
    numGrid = numGridNlppMap[iAtom];
    countRad = 0;
    countNlppRe = 0;
    countNlppIm = 0;
    gridShiftNowRe = gridStIndRe[iAtom];
    gridShiftNowIm = gridStIndIm[iAtom];
    /* cpy the wave function */
    //printf("numGrid %i locOpt[atomType+1] %i\n",numGrid,locOpt[atomType+1]);
    if(numGrid>0){ //if numGrid=0, only local pp will be calculated
      for(iGrid=0;iGrid<numGrid;iGrid++){
	gridIndex = gridNlppMap[iAtom][iGrid];
	wfNbhd[iGrid] = wfReal[gridIndex];
      }
      for(l=0;l<=atomLMax[atomType];l++){
        if(locOpt[atomType+1]!=l){
          for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
            radIndex = atomRadMap[atomType][countRad+iRad];
            energyl = 0.0;
            for(m=0;m<=l;m++){
              if(m!=0){
                dotRe = ddotBlasWrapper(numGrid,&wfNbhd[0],1,
                        &vnlPhiAtomGridRe[gridShiftNowRe],1)*volElem;
                dotIm = ddotBlasWrapper(numGrid,&wfNbhd[0],1,
                        &vnlPhiAtomGridIm[gridShiftNowIm],1)*volElem;
                if(forceCalcFlag==1){
                  dotReAll[iAtom][countNlppRe+m] = dotRe;
                  dotImAll[iAtom][countNlppIm+m-1] = dotIm;
                }
                energyl += 2.0*(dotRe*dotRe+dotIm*dotIm)*vpsNormList[radIndex]*volInv;
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
                //dotRe = ddot1(numGrid,&wfNbhd[0],1,
                //	    &vnlPhiAtomGridRe[gridShiftNowRe],1)*volElem;
                dotRe = ddotBlasWrapper(numGrid,&wfNbhd[0],1,
                              &vnlPhiAtomGridRe[gridShiftNowRe],1)*volElem;
                //printf("volElem %lg\n",volElem);
                //printf("11111111 dotRe %.8lg\n",dotRe);
                
                //if(l==1&&m==0){
                  //printf("%lg\n",wfNbhd[4138]);
                  //for(iGrid=0;iGrid<numGrid;iGrid++){
                  //  printf("333333333 %i %i %i %lg %lg\n",iAtom,iGrid,gridShiftNowRe+iGrid,wfNbhd[iGrid],vnlPhiAtomGridRe[gridShiftNowRe+iGrid]);
                  //}
                  //fflush(stdout);
                  //exit(0);
                //}
                
                if(forceCalcFlag==1){
                  dotReAll[iAtom][countNlppRe] = dotRe;
                }
                energyl += dotRe*dotRe*vpsNormList[radIndex]*volInv;
                dotRe *= vpsNormList[radIndex];
                daxpyBlasWrapper(numGrid,dotRe,
                          &vnlPhiAtomGridRe[gridShiftNowRe],1,
                          &forceTemp[iAtom*numGridMax],1);
                gridShiftNowRe += numGrid;
              }//endif m
            }//endfor m
            energy += energyl;
            countNlppRe += l+1;
            countNlppIm += l;
          }//endfor iRad
          countRad += atomLRadNum[atomType][l];
        }
        else{
          /* We don't need to any calculation here, but we need to shift indecies. */
          for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
            radIndex = atomRadMap[atomType][countRad+iRad];
            for(m=0;m<=l;m++){
              if(m!=0){
                gridShiftNowRe += numGrid;
                gridShiftNowIm += numGrid;
              }
              else{
                gridShiftNowRe += numGrid;
              }//endif m
            }//endfor m
            countNlppRe += l+1;
            countNlppIm += l;
          }//endfor iRad
          countRad += atomLRadNum[atomType][l];
        }//endif locOpt
      }//endfor l
    }//endif numGrid
  }//endfor iAtom
  //fflush(stdout);
  //exit(0);
  //printf("energy %lg\n",energy);

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

  for(iGrid=0;iGrid<numGridNlppAll;iGrid++){
    gridIndex = gridNlppInd[iGrid];
    numUpdate = gridNlppMapInv[iGrid][0];
    for(iAtom=0;iAtom<numUpdate;iAtom++){
      atomInd = gridNlppMapInv[iGrid][iAtom*2+1];
      gridInd2 = gridNlppMapInv[iGrid][iAtom*2+2];
      forceRealNlpp[gridIndex] += forceTemp[atomInd*numGridMax+gridInd2];
    }
  }

  if(energyCalcFlag==1){
    *energyNl += energy;
    //stat_avg->cp_enl += energy;
    //printf("eeeeeeeeeeeeenergy %lg\n",stat_avg->cp_enl);
    //printf("energyNl %lg\n",*energyNl);
  }

/*======================================================================*/
/* III) free local memory                                               */

  free(wfNbhd);
  free(forceTemp);

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcPseudoWf(CP *cp,CLASS *class,GENERAL_DATA *generalData)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate pseudo wave function. Could be done 1. before SCF 2. in     */
/* filter. I need to test which one is faster.				 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  PSEUDO *pseudo = &(cp->pseudo);
  CPOPTS *cpopts = &(cp->cpopts);
  PARA_FFT_PKG3D *cpParaFftPkg3dLgBigBox = &(cp->cp_para_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cpParaFftPkg3dSparse = &(cp->cp_sclr_fft_pkg3d_sparse);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CELL *cell = &(generalData->cell);
  CLATOMS_POS *clatoms_pos = &(class->clatoms_pos[1]);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  ATOMMAPS *atommaps = &(class->atommaps);
  COMMUNICATE *communicate = &(cp->communicate);
  CPEWALD *cpewald = &(cp->cpewald);

  int numAtomType = atommaps->natm_typ;
  int numAtomTot = clatoms_info->natm_tot;
  int iAtom,iGrid;
  int numGrid;
  int numGridMax;
  int countRad;
  int iPart,iRad,l,m,iType;
  int radIndex,gridIndex;
  int aIndex,bIndex,cIndex;
  int atomType;
  int countWfTot;
  int gridShiftRe,gridShiftIm;
  int realSparseOpt = cpewald->realSparseOpt;
  int nkc,nkb,nka;
  //int nkc = cpParaFftPkg3dLgBigBox->nkf3;
  //int nkb = cpParaFftPkg3dLgBigBox->nkf2;
  //int nka = cpParaFftPkg3dLgBigBox->nkf1;
  int div,res;
  int ylmShift;
  int fftw3dFlag = cpopts->fftw3dFlag;
  int numThreads = communicate->numThreads;
  int iThread;
  
  int *gridMap;
  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *gridStIndRe = pseudoReal->gridStIndRe;
  int *gridStIndIm = pseudoReal->gridStIndIm;
  int *atomLMax = pseudoReal->numLMax; //max L for each atom
  int *locOpt = pseudo->loc_opt;
  int **atomLRadNum = pseudoReal->atomLRadNum; //num of radical functions for each atom and each L
  int **atomRadMap = pseudoReal->atomRadMap;  //map of radical functions for each atom, starting from l=0 to l=lmax
  int *numGridNlppMap = pseudoReal->numGridNlppMap;
  int **gridNlppMap = pseudoReal->gridNlppMap;

  double vpsNorm;
  double dotRe,dotIm;
  double xTemp,yTemp,zTemp;
  double xTemp2,yTemp2,zTemp2;

  double *trig; // sin(theta),cos(theta),sin(phi),cos(phi)
  double *vpsNormList = pseudoReal->vpsNormList; // 1/vpsNorm 
  double *gridAtomNbhd;
  //double *ylm;
  double a[3],b[3],c[3];
  double aGrid[3],bGrid[3],cGrid[3];
  double reduceCoord[3];
  double *nucleiCoord = (double*)cmalloc(numThreads*3*sizeof(double));
  double *hmat = cell->hmat;
  double *hmati = cell->hmati;
  double *x = clatoms_pos->x;
  double *y = clatoms_pos->y;
  double *z = clatoms_pos->z;
  double *radFun;
  double *pseudoFunTemp;

  double *vnlPhiAtomGridRe = pseudoReal->vnlPhiAtomGridRe;
  double *vnlPhiAtomGridIm = pseudoReal->vnlPhiAtomGridIm;
  //debug
  int nfft,numGridTot;
  double *testwfReal;
  
  if(realSparseOpt==0){
    nkc = cpParaFftPkg3dLgBigBox->nkf3;
    nkb = cpParaFftPkg3dLgBigBox->nkf2;
    nka = cpParaFftPkg3dLgBigBox->nkf1;
    nfft = cpParaFftPkg3dLgBigBox->nfft;
    numGridTot = nfft/2;
  }
  else{
    nkc = cpParaFftPkg3dSparse->nkf3;
    nkb = cpParaFftPkg3dSparse->nkf2;
    nka = cpParaFftPkg3dSparse->nkf1;
    nfft = cpParaFftPkg3dSparse->nfft;
    numGridTot = nfft/2;
  }
  testwfReal = (double*)calloc(numGridTot,sizeof(double));

  numGridMax = numGridNlppMap[0];
  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    if(numGridNlppMap[iAtom]>numGridMax)numGridMax = numGridNlppMap[iAtom];
  }
  //printf("numGridMax %i\n",numGridMax);
  gridAtomNbhd = (double*)cmalloc(numThreads*numGridMax*3*sizeof(double));
  radFun = (double*)cmalloc(numThreads*numGridMax*sizeof(double));
  trig = (double*)cmalloc(4*numGridMax*numThreads*sizeof(double));
  a[0] = hmat[1];a[1] = hmat[2];a[2] = hmat[3];
  b[0] = hmat[4];b[1] = hmat[5];b[2] = hmat[6];
  c[0] = hmat[7];c[1] = hmat[8];c[2] = hmat[9];
  aGrid[0] = a[0]/nka;aGrid[1] = a[1]/nka;aGrid[2] = a[2]/nka;
  bGrid[0] = b[0]/nkb;bGrid[1] = b[1]/nkb;bGrid[2] = b[2]/nkb;
  cGrid[0] = c[0]/nkc;cGrid[1] = c[1]/nkc;cGrid[2] = c[2]/nkc;

  countWfTot = 0;
  omp_set_num_threads(numThreads);
  #pragma omp parallel private(iThread,iAtom,atomType,numGrid,gridShiftRe,gridShiftIm,iGrid,gridIndex,cIndex,res,bIndex,aIndex,xTemp,yTemp,zTemp,xTemp2,yTemp2,zTemp2,countRad,l,iRad,radIndex,m,ylmShift)
  {
    iThread = omp_get_thread_num();
    #pragma omp for
    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      atomType = iAtomAtomType[iAtom+1]-1;
      numGrid = numGridNlppMap[iAtom];
      nucleiCoord[iThread*3] = x[iAtom+1];
      nucleiCoord[iThread*3+1] = y[iAtom+1];
      nucleiCoord[iThread*3+2] = z[iAtom+1];
      gridShiftRe = gridStIndRe[iAtom];
      gridShiftIm = gridStIndIm[iAtom];
      // calculate difference between nuclei and grid
      for(iGrid=0;iGrid<numGrid;iGrid++){
	gridIndex = gridNlppMap[iAtom][iGrid];
	if(fftw3dFlag==0){
	  cIndex = gridIndex/(nka*nkb);
	  res = gridIndex%(nka*nkb);
	  bIndex = res/nka;
	  aIndex = res%nka;
	}
	else{
	  aIndex = gridIndex/(nkb*nkc);
	  res = gridIndex%(nkb*nkc);
	  bIndex = res/nkc;
	  cIndex = res%nkc;
	}
	xTemp = aGrid[0]*aIndex+bGrid[0]*bIndex+cGrid[0]*cIndex-nucleiCoord[iThread*3];
	yTemp = aGrid[1]*aIndex+bGrid[1]*bIndex+cGrid[1]*cIndex-nucleiCoord[iThread*3+1];
	zTemp = aGrid[2]*aIndex+bGrid[2]*bIndex+cGrid[2]*cIndex-nucleiCoord[iThread*3+2];

	//fold the differences into the box
	xTemp2 = xTemp*hmati[1]+yTemp*hmati[4]+zTemp*hmati[7];
	yTemp2 = xTemp*hmati[2]+yTemp*hmati[5]+zTemp*hmati[8];
	zTemp2 = xTemp*hmati[3]+yTemp*hmati[6]+zTemp*hmati[9];
	xTemp = xTemp2-NINT(xTemp2);
	yTemp = yTemp2-NINT(yTemp2);
	zTemp = zTemp2-NINT(zTemp2);
	gridAtomNbhd[iThread*numGridMax*3+iGrid*3]   = xTemp*hmat[1]+yTemp*hmat[4]+zTemp*hmat[7];
	gridAtomNbhd[iThread*numGridMax*3+iGrid*3+1] = xTemp*hmat[2]+yTemp*hmat[5]+zTemp*hmat[8];
	gridAtomNbhd[iThread*numGridMax*3+iGrid*3+2] = xTemp*hmat[3]+yTemp*hmat[6]+zTemp*hmat[9];
      }
      calcTrig(&gridAtomNbhd[iThread*numGridMax*3],numGrid,&trig[4*numGridMax*iThread]);
      //printf("theta %.16lg phi %.16lg\n",acos(trig[5]),atan2(trig[6],trig[7]));
      countRad = 0;
      for(l=0;l<=atomLMax[atomType];l++){
        if(locOpt[atomType+1]!=l){
          /* calculate sherical harmonic */
          double *ylm = (double*)cmalloc((2*l+1)*numGridMax*numThreads*sizeof(double));
          calcSpHarm(&ylm[iThread*(2*l+1)*numGridMax],l,
                     &gridAtomNbhd[iThread*numGridMax*3],numGrid,
                     &trig[4*numGridMax*iThread]);
          for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
            radIndex = atomRadMap[atomType][countRad+iRad];
            /* calculate radial fun u*phi */
            calcRadFun(&gridAtomNbhd[iThread*numGridMax*3],radIndex,pseudoReal,
                        &radFun[iThread*numGridMax],numGrid);
            for(m=0;m<=l;m++){
              if(m!=0){
                ylmShift = (m*2-1)*numGrid;
                for(iGrid=0;iGrid<numGrid;iGrid++){
                  vnlPhiAtomGridRe[gridShiftRe+iGrid] = radFun[iThread*numGridMax+iGrid]*ylm[iThread*(2*l+1)*numGridMax+ylmShift+iGrid*2];
                  vnlPhiAtomGridIm[gridShiftIm+iGrid] = radFun[iThread*numGridMax+iGrid]*ylm[iThread*(2*l+1)*numGridMax+ylmShift+iGrid*2+1];
                }//endfor iGrid
                gridShiftRe += numGrid;
                gridShiftIm += numGrid;
              }else{
                for(iGrid=0;iGrid<numGrid;iGrid++){
                  vnlPhiAtomGridRe[gridShiftRe+iGrid] = radFun[iThread*numGridMax+iGrid]*ylm[iThread*(2*l+1)*numGridMax+iGrid];
                  //printf("radFun %lg ylm %lg vnl %lg\n",radFun[iThread*numGridMax+iGrid],ylm[iThread*(2*l+1)*numGridMax+iGrid],vnlPhiAtomGridRe[gridShiftRe+iGrid]);

                  //gridIndex = gridNlppMap[iAtom][iGrid];
                  //testwfReal[gridIndex] = vnlPhiAtomGridRe[gridShiftRe+iGrid];
                }
                gridShiftRe += numGrid;
              }//endif m
              /*
              if(l==1&&m==0){
                fflush(stdout);
                exit(0);
              }
              */
            }//endfor m
          }//endfor iRad
          countRad += atomLRadNum[atomType][l];
          free(ylm);
        }
        else{
          for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
            for(m=0;m<=l;m++){
              if(m!=0){
                gridShiftRe += numGrid;
                gridShiftIm += numGrid;
              }
              else{
                gridShiftRe += numGrid;
              }//endif m
            }//endfor m
          }//endfor iRad
          countRad += atomLRadNum[atomType][l];
        }//endif locOpt
      }//endfor l
    }//endfor iAtom
  }//endfor omp
  //debug
  numGrid = numGridNlppMap[0];
  //printf("vnllllllll m=0 %.16lg\n",vnlPhiAtomGridRe[1]);
  //printf("vnllllllll m=1 %.16lg %.16lg\n",vnlPhiAtomGridRe[numGrid+1],vnlPhiAtomGridIm[1]);

  /*
  for(iGrid=0;iGrid<numGridTot;iGrid++){
    printf("55555 %.16lg\n",testwfReal[iGrid]);
  }
  */
  //fflush(stdout);
  //exit(0);
  free(radFun);
  free(trig);
  free(gridAtomNbhd);
  free(nucleiCoord);
  free(testwfReal);
  

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcTrig(double *gridAtomNbhd,int numGrid,double *trig)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* KS potential from local pseudo pp can be calculated before SCF        */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int iGrid;
  int j,debugFlag;
  double x,y,z,rInv,rProj2,rProj,rProjInv,r2;
  
  for(iGrid=0;iGrid<numGrid;iGrid++){
    x = gridAtomNbhd[iGrid*3];
    y = gridAtomNbhd[iGrid*3+1];
    z = gridAtomNbhd[iGrid*3+2];
    rProj2 = x*x+y*y;
    r2 = rProj2+z*z;
    if(rProj2<1.0e-30){
      x += 1.0e-14;
      rProj2 = x*x+y*y;
      r2 = rProj2+z*z;
    }
    rInv = 1.0/sqrt(r2);
    rProj = sqrt(rProj2);
    rProjInv = 1.0/rProj;
    trig[iGrid*4] = rProj*rInv; //sin(theta)
    trig[iGrid*4+1] = z*rInv; //cos(theta)
    trig[iGrid*4+2] = y*rProjInv;
    trig[iGrid*4+3] = x*rProjInv;
    
    //debugFlag = 0;
    /*
    for(j=0;j<4;j++){
      //if(isinf(trig[iGrid*4+j])==1||isnan(trig[iGrid*4+j])==1){
	printf("iiiiiGrid %i j %i x %lg y %lg z %lg trig %lg\n",iGrid,j,x,y,z,trig[iGrid*4+j]);
      //}
    }
    */
  }
  //fflush(stdout);
  //exit(0);
  
/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcSpHarm(double *ylm,int l,double *gridAtomNbhd,int numGrid,
		double *trig)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Spheircal harmonic values at each grid point.		         */
/* For each l, array ylm stores yl0,yl1(Re),yl1(Im)...			 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int iGrid;
  int ind;
  double temp1;
  double piInv = 1.0/M_PI;
  double pre00 = 0.5*sqrt(piInv);
  double pre10 = 0.5*sqrt(3.0*piInv);
  double pre11 = -0.5*sqrt(1.5*piInv);
  double pre20 = 0.25*sqrt(5.0*piInv);
  double pre21 = -0.5*sqrt(7.5*piInv);
  double pre22 = 0.25*sqrt(7.5*piInv);
  double pre30 = 0.25*sqrt(7.0*piInv);
  double pre31 = -0.125*sqrt(21*piInv);
  double pre32 = 0.25*sqrt(52.5*piInv);
  double pre33 = -0.125*sqrt(35*piInv);


  switch(l){
    case 0:
      for(iGrid=0;iGrid<numGrid;iGrid++)ylm[iGrid] = pre00;
      break;
    case 1:
      for(iGrid=0;iGrid<numGrid;iGrid++)ylm[iGrid] = pre10*trig[iGrid*4+1];
      for(iGrid=0;iGrid<numGrid;iGrid++){
	temp1 = pre11*trig[iGrid*4];
	ind = numGrid+iGrid*2;
	ylm[ind] = temp1*trig[iGrid*4+3];
	ylm[ind+1] = temp1*trig[iGrid*4+2];
      }
      //if(l==1)printf("m=0 ylm %.16lg\n",ylm[1]);
      //printf("m=1 ylm Re %.16lg im %.16lg\n",ylm[numGrid+2],ylm[numGrid+3]);
      break;
    case 2:
      for(iGrid=0;iGrid<numGrid;iGrid++){
	ylm[iGrid] = pre20*(3.0*trig[iGrid*4+1]*trig[iGrid*4+1]-1.0);
      }
      for(iGrid=0;iGrid<numGrid;iGrid++){
	temp1 = pre21*trig[iGrid*4]*trig[iGrid*4+1];
	ind = numGrid+2*iGrid;
	ylm[ind] = temp1*trig[iGrid*4+3];
	ylm[ind+1] = temp1*trig[iGrid*4+2];
      }
      for(iGrid=0;iGrid<numGrid;iGrid++){
	temp1 = pre22*trig[iGrid*4]*trig[iGrid*4];
	ind = numGrid*3+iGrid*2;
	ylm[ind] = temp1*(2.0*trig[iGrid*4+3]*trig[iGrid*4+3]-1.0);
	ylm[ind+1] = temp1*(2.0*trig[iGrid*4+3]*trig[iGrid*4+2]);
      }
      break;
    case 3:
      for(iGrid=0;iGrid<numGrid;iGrid++){
	ylm[iGrid] = pre30*(5.0*trig[iGrid*4+1]*trig[iGrid*4+1]*trig[iGrid*4+1]-3.0*trig[iGrid*4+1]);
      }
      for(iGrid=0;iGrid<numGrid;iGrid++){
	temp1 = pre31*trig[iGrid*4]*(5.0*trig[iGrid*4+1]*trig[iGrid*4+1]-1.0);
	ind = numGrid+iGrid*2;
	ylm[ind] = temp1*trig[iGrid*4+3];
	ylm[ind+1] = temp1*trig[iGrid*4+2];
      }
      for(iGrid=0;iGrid<numGrid;iGrid++){
        temp1 = pre32*trig[iGrid*4]*trig[iGrid*4]*trig[iGrid*4+1];
	ind = numGrid*3+iGrid*2;
	ylm[ind] = temp1*(2.0*trig[iGrid*4+3]*trig[iGrid*4+3]-1.0);
	ylm[ind+1] = temp1*(2.0*trig[iGrid*4+3]*trig[iGrid*4+2]);
      }
      for(iGrid=0;iGrid<numGrid;iGrid++){
	temp1 = pre33*trig[iGrid*4]*trig[iGrid*4]*trig[iGrid*4];
	ind = numGrid*5+iGrid*2;
	ylm[ind] = temp1*(4.0*trig[iGrid*4+3]*trig[iGrid*4+3]*trig[iGrid*4+3]-3.0*trig[iGrid*4+3]);
	ylm[ind+1] = temp1*(3.0*trig[iGrid*4+2]-4.0*trig[iGrid*4+2]*trig[iGrid*4+2]*trig[iGrid*4+2]);
      }
      break;
	
  }

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcRadFun(double *gridAtomNbhd,int radIndex,PSEUDO_REAL *pseudoReal,
	double *radFun,int numGrid)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* KS potential from local pseudo pp can be calculated before SCF        */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int numInterpGrid = pseudoReal->numInterpGrid;
  int interpGridSt = radIndex*numInterpGrid;
  int gridInd;
  int interpInd;
  int iGrid;

  double rMin = 0.0;
  double dr = pseudoReal->dr;
  double h;
  double x,y,z,r,r0;

  double *vps0 = pseudoReal->vpsReal0;
  double *vps1 = pseudoReal->vpsReal1;
  double *vps2 = pseudoReal->vpsReal2;
  double *vps3 = pseudoReal->vpsReal3;


  for(iGrid=0;iGrid<numGrid;iGrid++){
    x = gridAtomNbhd[iGrid*3];
    y = gridAtomNbhd[iGrid*3+1];
    z = gridAtomNbhd[iGrid*3+2];
    r = sqrt(x*x+y*y+z*z);
    gridInd = (int)((r-rMin)/dr)+1;
    r0 = (gridInd-1)*dr+rMin;
    h = r-r0;
    interpInd = interpGridSt+gridInd;
    radFun[iGrid] = ((vps3[interpInd]*h+vps2[interpInd])*h+vps1[interpInd])*h+vps0[interpInd];
    //printf("rrrrrrrrrrrrrr %i %.16lg %.16lg\n",iGrid,r,radFun[iGrid]);
  }
  
/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcDotNlpp(double *wfNbhd,double *radFun,double *ylm,
	 double *dotRe,double *dotIm,int numGrid)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* KS potential from local pseudo pp can be calculated before SCF        */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  int iGrid;
  double dotReLocal = 0.0;
  double dotImLocal = 0.0;
  
  for(iGrid=0;iGrid<numGrid;iGrid++){
    dotReLocal += wfNbhd[iGrid]*radFun[iGrid]*ylm[2*iGrid];
    dotImLocal += -wfNbhd[iGrid]*radFun[iGrid]*ylm[2*iGrid+1];
  }
  *dotRe = dotReLocal;
  *dotIm = dotImLocal;

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/


