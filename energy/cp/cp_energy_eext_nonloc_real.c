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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void nlppKBRealFilter(CP *cp,CLASS *class,GENERAL_DATA *general_data,double *wfReal,
	      double *forceRealNlpp)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Real space nlpp, only used in filtering	             */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  PSEUDO *pseudo = &(cp->pseudo);
  PARA_FFT_PKG3D *cpParaFftPkg3dLgBigBox = &(cp->cp_para_fft_pkg3d_lg);
  PSEUDO *pseudoReal = &(pseudo->pseudoReal);
  CELL *cell = &(generalData->cell);
  CLATOMS_POS *clatoms_pos = &(class->clatoms_pos[1]);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  ATOMMAPS *atommaps = &(class->atommaps);

  int numAtomType = atommaps->natm_typ;
  int numGrid;
  int numGridMax;
  int countRad;
  int iPart,iRad,l,m,iType;
  int radIndex,gridIndex;
  int aIndex,bIndex,cIndex;
  int atomType;
  int nkc = cpParaFftPkg3dLgBigBox->nkf3;
  int nkb = cpParaFftPkg3dLgBigBox->nkf2;
  int nka = cpParaFftPkg3dLgBigBox->nkf1;
  int numGridTot = nka*nkb*nkc;
  int div,res;

  int *gridMap;
  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *atomLMax = pseudoReal->numLMax; //max L for each atom
  int *atomGridTot = pseudoReal->atomGridTot;
  int **gridShift = pseudoReal->gridShift;
  int **atomLRadNum; //num of radical functions for each atom and each L
  int **atomRadMap;  //map of radical functions for each atom, starting from l=0 to l=lmax
  int **numRadFun = pseudoReal->numRadFun;
  int *numGridNlppMap = pseudoReal->numgridNlppMap;
  int **gridNlppMap = pseudoReal->gridNlppMap;

  double vpsNorm;
  double vol;
  double volElem;
  double dotRe,dotIm;

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

  double **vnlPhiAtomGridRe = pseudoReal->vnlPhiAtomGridRe;
  double **vnlPhiAtomGridIm = pseudoReal->vnlPhiAtomGridIm;

/*======================================================================*/
/* I) Allocate local memory                                             */

  numGridMax = numGridNlppMap[0];
  for(iPart=0;iPart<numAtom;iPart++){
    if(numGridNlppMap[iPart]>numGridMax)numGridMax = numGridNlppMap[iPart];
  }
  wfNbhd = (double*)cmalloc(numGridMax*sizeof(double));
  forceTemp = (double*)cmalloc(numGridMax*sizeof(double));
  volElem = getdeth(hmat)/numGridTot;

/*======================================================================*/
/* II) Loop over iPart/irad/ipart/m                                 */

  gridShiftNowRe = 0;
  gridShiftNowIm = 0;
  for(iAtom=0;iAtom<numAtom;iAtom++){
    atomType = iAtomAtomType[iAtom+1]-1;
    numGrid = numGridNlppMap[atomType];
    countRad = 0;
    /* cpy the wave function */
    if(numGrid>0){ //if numGrid=0, only local pp will be calculated
      for(iGrid=0;iGrid<numGrid;iGrid++){
	gridIndex = gridNlppMap[iAtom][iGrid];
	wfNbhd[iGrid] = wfReal[gridIndex];
	forceTemp[iGrid] = 0.0;
      }
      for(l=0;l<numLMax;l++){
	for(iRad=0;iRad<atomLRadNum[atomType][l]){
	  radIndex = atomRadMap[atomType][countRad+iRad];
	  for(m=0;m<=l;m++){
	    //calcDotNlpp(wfNbhd,radFun,&ylm[ylmShift],&dotRe,&dotIm);
	    if(m!=0){
	      dotRe = ddotBlasWrapper(numGrid,WfNhbd,1,
		      &vnlPhiAtomGridRe[gridShiftNowRe],1);
	      dotIm = ddotBlasWrapper(numGrid,WfNhbd,1,
                      &vnlPhiAtomGridIm[gridShiftNowIm],1);
	      dotRe *= 2.0*vpsNormList[radIndex];
	      dotIm *= -2.0*vpsNormList[radIndex];
	      daxpyBlasWrapper(numGrid,dotRe,&vnlPhiAtomGridRe[gridShiftNowRe],
			       forceTemp);
	      daxpyBlasWrapper(numGrid,dotIm,&vnlPhiAtomGridIm[gridShiftNowIm],
                               forceTemp);
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
              dotRe = ddotBlasWrapper(numGrid,WfNhbd,1,
                      &vnlPhiAtomGridRe[gridShiftNowRe],1);
	      /*
	      for(iGrid=0;iGrid<numGrid;iGrid++){
		forceTemp[iGrid] += radFun[iGrid]*ylm[ylmShift+iGrid]*dotRe*vpsNormList[radIndex];
	      }//endfor iGrid
	      */
	      dotRe *= vpsNormList[radIndex];
	      daxpyBlasWrapper(numGrid,dotRe,&vnlPhiAtomGridRe[gridShiftNowRe],
			       forceTemp);
	      gridShiftNowRe += numGrid;
	    }//endif m
	  }//endfor m
	}//endfor iRad
        countRad += atomLRadNum[iAtom][l];
      }//endfor l
      for(iGrid=0;iGrid<numGrid;iGrid++){
	gridIndex = gridNlppMap[iAtom][iGrid];
	forceRealNlpp[gridIndex] += forceTemp[iGrid];
      }//endfor iGrid
    }//endif numGrid
  }//endfor iAtom

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
void calcPseudoWf(CP *cp,CLASS *class,GENERAL_DATA *general_data,
		  double *nlWf)
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
  PSEUDO *pseudoReal = &(pseudo->pseudoReal);
  CELL *cell = &(generalData->cell);
  CLATOMS_POS *clatoms_pos = &(class->clatoms_pos[1]);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  ATOMMAPS *atommaps = &(class->atommaps);

  int numAtomType = atommaps->natm_typ;
  int iAtom;
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

  int *gridMap;
  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *atomLMax = pseudoReal->numLMax; //max L for each atom
  int *atomGridTot = pseudoReal->atomGridTot;
  int **atomLRadNum = pseudoReal->atomLRadNum; //num of radical functions for each atom and each L
  int **atomRadMap = pseudoReal->atomRadMap;  //map of radical functions for each atom, starting from l=0 to l=lmax
  int **numRadFun = pseudoReal->numRadFun;
  int *numGridNlppMap = pseudoReal->numgridNlppMap;
  int **gridNlppMap = pseudoReal->gridNlppMap;

  double vpsNorm;
  double dotRe,dotIm;

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
  double *radFun = pseudoReal->radFun;
  double *forceTemp;
  double *pseudoFunTemp;

  double **vnlPhiAtomGridRe = pseudoReal->vnlPhiAtomGridRe;
  double **vnlPhiAtomGridIm = pseudoReal->vnlPhiAtomGridIm;

  numGridMax = numGridNlppMap[0];
  for(iPart=0;iPart<numAtom;iPart++){
    if(numGridNlppMap[iPart]>numGridMax)numGridMax = numGridNlppMap[iPart];
  }
  gridAtomNbhd = (double*)cmalloc(numGridMax*3*sizeof(double));
  wfNbhd = (double*)cmalloc(numGridMax*sizeof(double));
  radFun = (double*)cmalloc(numGridMax*sizeof(double));
  forceTemp = (double*)cmalloc(numGridMax*sizeof(double));
  trig = (double*)cmalloc(4*numGridMax*sizeof(double));
  a[0] = hmat[1];a[1] = hmat[2];a[2] = hmat[3];
  b[0] = hmat[4];b[1] = hmat[5];b[2] = hmat[6];
  c[0] = hmat[7];c[1] = hmat[8];c[2] = hmat[9];
  aGrid[0] = a[0]/nka;aGrid[1] = a[1]/nka;aGrid[2] = a[2]/nka;
  bGrid[0] = b[0]/nkb;bGrid[1] = b[1]/nkb;bGrid[2] = b[2]/nkb;
  cGrid[0] = c[0]/nkc;cGrid[1] = c[1]/nkc;cGrid[2] = c[2]/nkc;


  for(iGrid=0;iGrid<numGrid;iGrid++){
    gridIndex = gridNlppMap[iAtom][iGrid];
    cIndex = gridIndex/(nka*nkb);
    res = gridIndex%(nka*nkb);
    bIndex = res/nkb;
    cIndex = res%nkb;
    gridAtomNbhd[iGrid*3] = aGrid[0]*aIndex+bGrid[0]*bIndex+cGrid[0]*cIndex;
    gridAtomNbhd[iGrid*3+1] = aGrid[1]*aIndex+bGrid[1]*bIndex+cGrid[1]*cIndex;
    gridAtomNbhd[iGrid*3+2] = aGrid[2]*aIndex+bGrid[2]*bIndex+cGrid[2]*cIndex;
    wfNbhd[iGrid] = wfReal[gridIndex];
    forceTemp[iGrid] = 0.0;
  }
  calcTrig(gridAtomNbhd,numGrid,trig);
  countWfTot = 0;
  gridShiftRe = 0;
  gridShiftIm = 0;
  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    atomType = iAtomAtomType[iAtom+1]-1;
    for(l=0;l<numLMax;l++){
      /* calculate sherical harmonic */
      ylm = (double*)cmalloc((2*l+1)*numGrid*sizeof(double));
      calcSpHarm(ylm,l,gridAtomNbhd,numGrid,&nucleiCoord[0]);
      for(iRad=0;iRad<atomLRadNum[atomType][l]){
	radIndex = atomRadMap[atomType][countRad+iRad];
	/* calculate radial fun u*phi */
	calcRadFun(gridAtomNbhd,radIndex,pseudoReal,radFun);
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
	    vnlPhiAtomGridRe[gridShiftRe+iGrid] = radFun[iGrid]*ylm[iGrid];
	    gridShiftRe += numGrid;
	  }//endif m
	}//endfor m
      }//endfor iRad
      countRad += atomLRadNum[iAtom][l];
      free(ylm);
    }//endfor l
  }//endfor iAtom

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
  double x,y,z,rInv,rProj2,rProj,rProjInv;
  
  for(iGrid=0;iGrid<numGrid;iGrid++){
    x = gridAtomNbhd[iGrid*3];
    y = gridAtomNbhd[iGrid*3+1];
    z = gridAtomNbhd[iGrid*3+2];
    rProj2 = x*x+y*y;
    rInv = 1.0/sqrt(rProj2+z*z);
    rProj = sqrt(rProj2);
    rProjInv = 1.0/rProj;
    trig[iGrid*4] = rProj*rInv; //sin(theta)
    trig[iGrid*4+1] = z*rInv; //cos(theta)
    trig[iGrid*4+2] = y*rProjInv;
    trig[iGrid*4+3] = x*rProjInv;
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
/* KS potential from local pseudo pp can be calculated before SCF        */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int iGrid;
  int ind;
  double temp1;
  double piInv = 1.0/M_PI;

  switch(l){
    case 0:
      double pre00 = 0.5*sqrt(piInv);
      for(iGrid=0;iGrid<numGrid;iGrid++)ylm[iGrid] = pre00;
      break;
    case 1:
      double pre10 = 0.5*sqrt(3.0*piInv);
      double pre11 = -0.5*sqrt(1.5*piInv);
      for(iGrid=0;iGrid<numGrid;iGrid++)ylm[iGrid] = pre10*trig[iGrid*4+1];
      for(iGrid=0;iGrid<numGrid;iGrid++){
	temp1 = pre11*trig[iGrid*4];
	ind = numGrid+iGrid*2;
	ylm[ind] = temp1*trig[iGrid*4+3];
	ylm[ind+1] = temp1*trig[iGrid*4+2];
      }
      break;
    case 2:
      double pre20 = 0.25*sqrt(5.0*piInv);
      double pre21 = -0.5*sqrt(7.5*piInv);
      double pre22 = 0.25*sqrt(7.5*piInv);
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
	double *radFun,int numGrid,double *nucleiCoord)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* KS potential from local pseudo pp can be calculated before SCF        */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int numInterpGrid = pseudoReal->numInterpGrid[radIndex];
  int interpGridSt = pseudoReal->interpGridSt[radIndex];
  int gridInd;
  int interpInd;

  double rMin = pseudoReal->rMin;
  double dr = pseudoReal->dr;
  double h;
  double x,y,z,r,r0;

  double *vps0 = pseudoReal->vps0;
  double *vps1 = pseudoReal->vps1;
  double *vps2 = pseudoReal->vps2;
  double *vps3 = pseudoReal->vps3;


  for(iGrid=0;iGrid<numGrid;iGrid++numGrid;iGrid++){
    x = gridAtomNbhd[iGrid*3]-nucleiCoord[0];
    y = gridAtomNbhd[iGrid*3+1]-nucleiCoord[1];
    z = gridAtomNbhd[iGrid*3+2]-nucleiCoord[2];
    r = sqrt(x*x+y*y+z*z);
    gridInd = int((r-rMin)/dr)+1;
    r0 = (gridInd-1)*dr+rMin;
    h = r-r0;
    interpInd = interpGridSt+gridInd;
    radFun[iGrid] = ((vps3[interpInd]*h+vps2[interpInd])*h+vps1[interpInd])*h+vps0[interpInd];
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


