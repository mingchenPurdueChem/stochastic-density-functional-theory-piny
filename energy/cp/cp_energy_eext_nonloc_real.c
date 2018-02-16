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
			   double *zfft_tmp,double *zfft,int flag)
/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
  PSEUDO *pseudo = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);

  PARA_FFT_PKG3D *cpParaFftPkg3dLgBigBox = &(cp->cp_para_fft_pkg3d_lg);
  int nfft = cpParaFftPkg3dLgBigBox->nfft;
  int numGrid = nfft/2;
  int iGrid;
  int forceCalcFlag = pseudoReal->forceCalcFlag;
  double *wfReal,*wfForceReal;

  wfReal = (double*)cmalloc(numGrid*sizeof(double));
  wfForceReal = (double*)cmalloc(numGrid*sizeof(double));
  printf("forceCalcFlag %i\n",forceCalcFlag);

  for(iGrid=0;iGrid<numGrid;iGrid++){
    wfReal[iGrid] = zfft_tmp[iGrid*2+1];
    wfForceReal[iGrid] = 0.0;
  }
  nlppKBRealEnergy(cp,class,generalData,wfReal,wfForceReal);
  if(forceCalcFlag==1)nlppKBRealEnergyForce(cp,class,generalData,wfReal);
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
    nlppKBRealEnergy(cp,class,generalData,wfReal,wfForceReal);
    if(forceCalcFlag==1)nlppKBRealEnergyForce(cp,class,generalData,wfReal);
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

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void nlppKBRealEnergy(CP *cp,CLASS *class,GENERAL_DATA *generalData,
		      double *wfReal,double *forceRealNlpp)
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
  PARA_FFT_PKG3D *cpParaFftPkg3dLgBigBox = &(cp->cp_para_fft_pkg3d_lg);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CELL *cell = &(generalData->cell);
  CLATOMS_POS *clatoms_pos = &(class->clatoms_pos[1]);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  ATOMMAPS *atommaps = &(class->atommaps);
  STAT_AVG *stat_avg = &(generalData->stat_avg);

  int numAtomType = atommaps->natm_typ;
  int numAtom = clatoms_info->natm_tot;
  int numGrid;
  int numGridMax;
  int countRad,countNlppRe,countNlppIm;
  int iPart,iRad,l,m,iType,iAtom,iGrid;
  int radIndex,gridIndex;
  int aIndex,bIndex,cIndex;
  int atomType;
  int nkc = cpParaFftPkg3dLgBigBox->nkf3;
  int nkb = cpParaFftPkg3dLgBigBox->nkf2;
  int nka = cpParaFftPkg3dLgBigBox->nkf1;
  int numGridTot = nka*nkb*nkc;
  int div,res;
  int gridShiftNowRe,gridShiftNowIm;
  int energyCalcFlag = pseudoReal->energyCalcFlag;

  int *gridMap;
  int *iAtomAtomType = atommaps->iatm_atm_typ;
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
  double *wfNbhd;
  double *radFun;
  double *forceTemp;
  double *pseudoFunTemp;
  double *vnlPhiAtomGridRe = pseudoReal->vnlPhiAtomGridRe;
  double *vnlPhiAtomGridIm = pseudoReal->vnlPhiAtomGridIm;

  double **dotReAll = pseudoReal->dotReAll;
  double **dotImAll = pseudoReal->dotImAll;

/*======================================================================*/
/* I) Allocate local memory                                             */

  numGridMax = numGridNlppMap[0];
  for(iPart=0;iPart<numAtom;iPart++){
    if(numGridNlppMap[iPart]>numGridMax)numGridMax = numGridNlppMap[iPart];
  }
  forceTemp = (double*)cmalloc(numGridMax*sizeof(double));
  wfNbhd = (double*)cmalloc(numGridMax*sizeof(double));
  vol = getdeth(hmat);
  volInv = 1.0/vol;
  volElem = vol/numGridTot;

/*======================================================================*/
/* II) Loop over iPart/irad/ipart/m                                 */

  gridShiftNowRe = 0;
  gridShiftNowIm = 0;
  for(iAtom=0;iAtom<numAtom;iAtom++){
    atomType = iAtomAtomType[iAtom+1]-1;
    numGrid = numGridNlppMap[iAtom];
    countRad = 0;
    countNlppRe = 0;
    countNlppIm = 0;
    /* cpy the wave function */
    if(numGrid>0){ //if numGrid=0, only local pp will be calculated
      for(iGrid=0;iGrid<numGrid;iGrid++){
	gridIndex = gridNlppMap[iAtom][iGrid];
	wfNbhd[iGrid] = wfReal[gridIndex];
	forceTemp[iGrid] = 0.0;
      }
      for(l=0;l<atomLMax[atomType];l++){
	for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
	  radIndex = atomRadMap[atomType][countRad+iRad];
	  energyl = 0.0;
	  for(m=0;m<=l;m++){
	    //calcDotNlpp(wfNbhd,radFun,&ylm[ylmShift],&dotRe,&dotIm);
	    if(m!=0){
	      dotRe = ddotBlasWrapper(numGrid,wfNbhd,1,
		      &vnlPhiAtomGridRe[gridShiftNowRe],1)*volElem;
	      dotIm = ddotBlasWrapper(numGrid,wfNbhd,1,
                      &vnlPhiAtomGridIm[gridShiftNowIm],1)*volElem;
	      dotReAll[iAtom][countNlppRe+m] = dotRe;
	      dotImAll[iAtom][countNlppIm+m-1] = dotIm;
	      energyl += 2.0*(dotRe*dotRe+dotIm*dotIm)*vpsNormList[radIndex]*volInv;
	      //printf("energy %lg\n",energyl);
	      //printf("m %i dotRe %lg dotIm %lg vpsNormList[radIndex] %lg\n",m,dotRe,dotIm,vpsNormList[radIndex]);
	      dotRe *= 2.0*vpsNormList[radIndex];
	      dotIm *= 2.0*vpsNormList[radIndex];
	      daxpyBlasWrapper(numGrid,dotRe,&vnlPhiAtomGridRe[gridShiftNowRe],1,
			       forceTemp,1);
	      daxpyBlasWrapper(numGrid,dotIm,&vnlPhiAtomGridIm[gridShiftNowIm],1,
                               forceTemp,1);
	      /*
	      for(iGrid=0;iGrid<numGrid;iGrid++){
		forceTemp[iGrid] += 2.0*radFun[iGrid]*(ylm[ylmShift+iGrid]*dotRe-
				    ylm[ylmShift+numGrid+iGrid]*dotIm)*vpsNormList[radIndex];
	      }//endfor iGrid
	      */
	      gridShiftNowRe += numGrid;
	      gridShiftNowIm += numGrid;
	    }
	    else{
              dotRe = ddotBlasWrapper(numGrid,wfNbhd,1,
                      &vnlPhiAtomGridRe[gridShiftNowRe],1)*volElem;
	      dotReAll[iAtom][countNlppRe] = dotRe;
	      //printf("numGrid %i\n",numGrid);
	      /*
	      for(iGrid=0;iGrid<numGrid;iGrid++){
		printf("griddddddddddd %lg %lg\n",wfNbhd[iGrid],vnlPhiAtomGridRe[gridShiftNowRe+iGrid]);
	      }
	      */
	      
	      energyl += dotRe*dotRe*vpsNormList[radIndex]*volInv;
	      //printf("energy %lg\n",energyl);
	      /*
	      for(iGrid=0;iGrid<numGrid;iGrid++){
		forceTemp[iGrid] += radFun[iGrid]*ylm[ylmShift+iGrid]*dotRe*vpsNormList[radIndex];
	      }//endfor iGrid
	      */
              //printf("m %i dotRe %lg volElem %lg vpsNormList[radIndex] %lg\n",m,dotRe*volElem,volElem,vpsNormList[radIndex]);

	      dotRe *= vpsNormList[radIndex];
              //printf("m %i dotRe %lg volElem %lg vpsNormList[radIndex] %lg\n",m,dotRe,volElem,vpsNormList[radIndex]);

	      daxpyBlasWrapper(numGrid,dotRe,&vnlPhiAtomGridRe[gridShiftNowRe],1,
			       forceTemp,1);
	      gridShiftNowRe += numGrid;
	      //printf("forcetemppppppp %lg\n",forceTemp[0]);
	    }//endif m
	  }//endfor m
	  //printf("energyl %lg\n",energyl);
	  energy += energyl;
	}//endfor iRad
        countRad += atomLRadNum[atomType][l];
	countNlppRe += m+1;
	countNlppIm += m;
      }//endfor l
      for(iGrid=0;iGrid<numGrid;iGrid++){
	gridIndex = gridNlppMap[iAtom][iGrid];
	forceRealNlpp[gridIndex] += forceTemp[iGrid];
      }//endfor iGrid
    }//endif numGrid
  }//endfor iAtom

  if(energyCalcFlag==1){
    stat_avg->cp_enl += energy;
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
  PARA_FFT_PKG3D *cpParaFftPkg3dLgBigBox = &(cp->cp_para_fft_pkg3d_lg);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CELL *cell = &(generalData->cell);
  CLATOMS_POS *clatoms_pos = &(class->clatoms_pos[1]);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  ATOMMAPS *atommaps = &(class->atommaps);

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
  int nkc = cpParaFftPkg3dLgBigBox->nkf3;
  int nkb = cpParaFftPkg3dLgBigBox->nkf2;
  int nka = cpParaFftPkg3dLgBigBox->nkf1;
  int div,res;
  int ylmShift;

  int *gridMap;
  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *atomLMax = pseudoReal->numLMax; //max L for each atom
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
  double *ylm;
  double a[3],b[3],c[3];
  double aGrid[3],bGrid[3],cGrid[3];
  double reduceCoord[3];
  double nucleiCoord[3];
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
  int nfft = cpParaFftPkg3dLgBigBox->nfft;
  int numGridTot = nfft/2;
  
  double *testwfReal = (double*)calloc(numGridTot,sizeof(double));

  numGridMax = numGridNlppMap[0];
  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    if(numGridNlppMap[iAtom]>numGridMax)numGridMax = numGridNlppMap[iAtom];
  }
  printf("numGridMax %i\n",numGridMax);
  gridAtomNbhd = (double*)cmalloc(numGridMax*3*sizeof(double));
  radFun = (double*)cmalloc(numGridMax*sizeof(double));
  trig = (double*)cmalloc(4*numGridMax*sizeof(double));
  a[0] = hmat[1];a[1] = hmat[2];a[2] = hmat[3];
  b[0] = hmat[4];b[1] = hmat[5];b[2] = hmat[6];
  c[0] = hmat[7];c[1] = hmat[8];c[2] = hmat[9];
  aGrid[0] = a[0]/nka;aGrid[1] = a[1]/nka;aGrid[2] = a[2]/nka;
  bGrid[0] = b[0]/nkb;bGrid[1] = b[1]/nkb;bGrid[2] = b[2]/nkb;
  cGrid[0] = c[0]/nkc;cGrid[1] = c[1]/nkc;cGrid[2] = c[2]/nkc;

  countWfTot = 0;
  gridShiftRe = 0;
  gridShiftIm = 0;
  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    atomType = iAtomAtomType[iAtom+1]-1;
    numGrid = numGridNlppMap[iAtom];
    nucleiCoord[0] = x[iAtom+1];
    nucleiCoord[1] = y[iAtom+1];
    nucleiCoord[2] = z[iAtom+1];
    // calculate difference between nuclei and grid
    for(iGrid=0;iGrid<numGrid;iGrid++){
      gridIndex = gridNlppMap[iAtom][iGrid];
      cIndex = gridIndex/(nka*nkb);
      res = gridIndex%(nka*nkb);
      bIndex = res/nkb;
      aIndex = res%nkb;
      xTemp = aGrid[0]*aIndex+bGrid[0]*bIndex+cGrid[0]*cIndex-nucleiCoord[0];
      yTemp = aGrid[1]*aIndex+bGrid[1]*bIndex+cGrid[1]*cIndex-nucleiCoord[1];
      zTemp = aGrid[2]*aIndex+bGrid[2]*bIndex+cGrid[2]*cIndex-nucleiCoord[2];

      //fold the differences into the box
      xTemp2 = xTemp*hmati[1]+yTemp*hmati[4]+zTemp*hmati[7];
      yTemp2 = xTemp*hmati[2]+yTemp*hmati[5]+zTemp*hmati[8];
      zTemp2 = xTemp*hmati[3]+yTemp*hmati[6]+zTemp*hmati[9];
      xTemp = xTemp2-NINT(xTemp2);
      yTemp = yTemp2-NINT(yTemp2);
      zTemp = zTemp2-NINT(zTemp2);
      gridAtomNbhd[iGrid*3]   = xTemp*hmat[1]+yTemp*hmat[4]+zTemp*hmat[7];
      gridAtomNbhd[iGrid*3+1] = xTemp*hmat[2]+yTemp*hmat[5]+zTemp*hmat[8];
      gridAtomNbhd[iGrid*3+2] = xTemp*hmat[3]+yTemp*hmat[6]+zTemp*hmat[9];
    }
    calcTrig(gridAtomNbhd,numGrid,trig);
    countRad = 0;
    for(l=0;l<atomLMax[atomType];l++){
      /* calculate sherical harmonic */
      ylm = (double*)cmalloc((2*l+1)*numGrid*sizeof(double));
      calcSpHarm(ylm,l,gridAtomNbhd,numGrid,trig);
      for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
	radIndex = atomRadMap[atomType][countRad+iRad];
	/* calculate radial fun u*phi */
	calcRadFun(gridAtomNbhd,radIndex,pseudoReal,radFun,numGrid);
	for(m=0;m<=l;m++){
	  if(m!=0){
	    ylmShift = (m*2-1)*numGrid;
	    for(iGrid=0;iGrid<numGrid;iGrid++){
	      vnlPhiAtomGridRe[gridShiftRe+iGrid] = radFun[iGrid]*ylm[ylmShift+iGrid*2];
	      vnlPhiAtomGridIm[gridShiftIm+iGrid] = radFun[iGrid]*ylm[ylmShift+iGrid*2+1];
	    }//endfor iGrid
	    gridShiftRe += numGrid;
	    gridShiftIm += numGrid;
	  }else{
	    for(iGrid=0;iGrid<numGrid;iGrid++){
	      vnlPhiAtomGridRe[gridShiftRe+iGrid] = radFun[iGrid]*ylm[iGrid];
	      gridIndex = gridNlppMap[iAtom][iGrid];
	      testwfReal[gridIndex] = vnlPhiAtomGridRe[gridShiftRe+iGrid];
	    }
	    gridShiftRe += numGrid;
	  }//endif m
	}//endfor m
      }//endfor iRad
      countRad += atomLRadNum[atomType][l];
      free(ylm);
    }//endfor l
  }//endfor iAtom
  /*
  for(iGrid=0;iGrid<numGridTot;iGrid++){
    printf("55555 %.16lg\n",testwfReal[iGrid]);
  }
  fflush(stdout);
  exit(0);
  */
  free(radFun);
  free(trig);
  free(gridAtomNbhd);

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
    /*
    debugFlag = 0;
    for(j=0;j<4;j++){
      if(isinf(trig[iGrid*4+j])==1||isnan(trig[iGrid*4+j])==1){
	printf("iGrid %i j %i x %lg y %lg z %lg trig %lg\n",iGrid,j,x,y,z,trig[iGrid*4+j]);
      }
    }
    */
  }
  
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
    //if(r<1.0e-10)printf("rrrrrrrrrrrrrr %i %.16lg %.16lg\n",iGrid,r,radFun[iGrid]);
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
  double dotReLocal;
  double dotImLocal;
  
  for(iGrid=0;iGrid<numGrid;iGrid++){
    dotReLocal += wfNbhd[iGrid]*radFun[iGrid]*ylm[2*iGrid];
    dotImLocal += -wfNbhd[iGrid]*radFun[iGrid]*ylm[2*iGrid+1];
  }
  *dotRe = dotReLocal;
  *dotIm = dotImLocal;

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/


