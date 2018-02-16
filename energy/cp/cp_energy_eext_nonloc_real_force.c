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
void nlppKBRealEnergyForce(CP *cp,CLASS *class,GENERAL_DATA *generalData,
                           double *wfReal)
/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*************************************************************************/
/* Calculate nuclei force						 */
/*************************************************************************/
/*==========================================================================*/
/*               Local variable declarations                                */

  PARA_FFT_PKG3D *cpParaFftPkg3dLgBigBox = &(cp->cp_para_fft_pkg3d_lg);
  PSEUDO *pseudo = &(cp->pseudo);
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
  int countRad;
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
  int countNlppRe,countNlppIm;

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
  double dotRe,dotIm,dotDevRe,dotDevIm;
  double energy = 0.0;
  double energyl;
  double forceNlX,forceNlY,forceNlZ;

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
  double *fx = clatoms_pos->fx;
  double *fy = clatoms_pos->fy;
  double *fz = clatoms_pos->fz;
  double *wfNbhd;
  double *radFun;
  double *forceTemp;
  double *pseudoFunTemp;

  //double *forceNlX = pseudoReal->forceNlX;
  //double *forceNlY = pseudoReal->forceNlY;
  //double *forceNlZ = pseudoReal->forceNlZ;
  double **dotReAll = pseudoReal->dotReAll;
  double **dotImAll = pseudoReal->dotImAll;

  double *vnlPhiAtomGridRe = pseudoReal->vnlPhiAtomGridRe;
  double *vnlPhiAtomGridIm = pseudoReal->vnlPhiAtomGridIm;
  double *vnlPhiDxAtomGridRe = pseudoReal->vnlPhiDxAtomGridRe;
  double *vnlPhiDxAtomGridIm = pseudoReal->vnlPhiDxAtomGridIm;
  double *vnlPhiDyAtomGridRe = pseudoReal->vnlPhiDyAtomGridRe;
  double *vnlPhiDyAtomGridIm = pseudoReal->vnlPhiDyAtomGridIm;
  double *vnlPhiDzAtomGridRe = pseudoReal->vnlPhiDzAtomGridRe;
  double *vnlPhiDzAtomGridIm = pseudoReal->vnlPhiDzAtomGridIm;


  gridShiftNowRe = 0;
  gridShiftNowIm = 0;

  for(iAtom=0;iAtom<numAtom;iAtom++){ //
    atomType = iAtomAtomType[iAtom+1]-1;
    numGrid = numGridNlppMap[iAtom];
    countRad = 0;
    forceNlX = 0.0;
    forceNlY = 0.0;
    forceNlZ = 0.0;
    countNlppRe = 0;
    countNlppIm = 0;
    if(numGrid>0){ //if numGrid=0, only local pp will be calculated
      for(iGrid=0;iGrid<numGrid;iGrid++){
        gridIndex = gridNlppMap[iAtom][iGrid];
        wfNbhd[iGrid] = wfReal[gridIndex];
        forceTemp[iGrid] = 0.0;
      }//endfor iGrid
      for(l=0;l<atomLMax[atomType];l++){
        for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
          radIndex = atomRadMap[atomType][countRad+iRad];
          energyl = 0.0;
          for(m=0;m<=l;m++){
            //calcDotNlpp(wfNbhd,radFun,&ylm[ylmShift],&dotRe,&dotIm);
            if(m!=0){
              dotRe = dotReAll[iAtom][countNlppRe+m];
              dotIm = dotReAll[iAtom][countNlppIm+m-1];
	      // x      
              dotDevRe = ddotBlasWrapper(numGrid,wfNbhd,1,
					 &vnlPhiDxAtomGridRe[gridShiftNowRe],1)*volElem;
              dotDevIm = ddotBlasWrapper(numGrid,wfNbhd,1,
				         &vnlPhiDxAtomGridIm[gridShiftNowIm],1)*volElem;
	      forceNlX += (dotRe*dotDevRe-dotIm*dotDevIm)*4.0*vpsNormList[radIndex]*volInv;
	      // y
	      dotDevRe = ddotBlasWrapper(numGrid,wfNbhd,1,
                                         &vnlPhiDyAtomGridRe[gridShiftNowRe],1)*volElem;
              dotDevIm = ddotBlasWrapper(numGrid,wfNbhd,1,
                                         &vnlPhiDyAtomGridIm[gridShiftNowIm],1)*volElem;
	      forceNlY += (dotRe*dotDevRe-dotIm*dotDevIm)*4.0*vpsNormList[radIndex]*volInv;
	      // z
              dotDevRe = ddotBlasWrapper(numGrid,wfNbhd,1,
                                         &vnlPhiDzAtomGridRe[gridShiftNowRe],1)*volElem;
              dotDevIm = ddotBlasWrapper(numGrid,wfNbhd,1,
                                         &vnlPhiDzAtomGridIm[gridShiftNowIm],1)*volElem;
              forceNlZ += (dotRe*dotDevRe-dotIm*dotDevIm)*4.0*vpsNormList[radIndex]*volInv;
              gridShiftNowRe += numGrid;
              gridShiftNowIm += numGrid;
            }
            else{
	      dotRe = dotReAll[iAtom][countNlppRe];
	      // x
              dotDevRe = ddotBlasWrapper(numGrid,wfNbhd,1,
					 &vnlPhiDxAtomGridRe[gridShiftNowRe],1)*volElem;
	      forceNlX += dotRe*dotDevRe*2.0*vpsNormList[radIndex]*volInv;
	      // y
              dotDevRe = ddotBlasWrapper(numGrid,wfNbhd,1,
                                         &vnlPhiDyAtomGridRe[gridShiftNowRe],1)*volElem;
              forceNlY += dotRe*dotDevRe*2.0*vpsNormList[radIndex]*volInv;
	      // z
              dotDevRe = ddotBlasWrapper(numGrid,wfNbhd,1,
                                         &vnlPhiDzAtomGridRe[gridShiftNowRe],1)*volElem;
	      forceNlZ += dotRe*dotDevRe*2.0*vpsNormList[radIndex]*volInv;
              gridShiftNowRe += numGrid;
            }//endif m
          }//endfor m
        }//endfor iRad
        countRad += atomLRadNum[atomType][l];
        countNlppRe += m+1;
        countNlppIm += m;
      }//endfor l
    }//endif numGrid
    fx[iAtom+1] += forceNlX;
    fy[iAtom+1] += forceNlY;
    fz[iAtom+1] += forceNlZ;
  }//endfor iAtom
/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcPseudoWfDev(CP *cp,CLASS *class,GENERAL_DATA *generalData)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate pseudo wave function derivative with respect to nuclei	 */
/* coordinate. Only update this before you decide to calculate force     */
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
  double rProj;

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
  double *radFun,*radDevFun;
  double *pseudoFunTemp;
  double *ylmTheta,*ylmPhi;
  double *dTheta,*dPhi,*ylmDx,*ylmDy,*ylmDz;

  double *vnlPhiAtomGridRe = pseudoReal->vnlPhiAtomGridRe;
  double *vnlPhiAtomGridIm = pseudoReal->vnlPhiAtomGridIm;
  double *vnlPhiDxAtomGridRe = pseudoReal->vnlPhiDxAtomGridRe;
  double *vnlPhiDxAtomGridIm = pseudoReal->vnlPhiDxAtomGridIm;
  double *vnlPhiDyAtomGridRe = pseudoReal->vnlPhiDyAtomGridRe;
  double *vnlPhiDyAtomGridIm = pseudoReal->vnlPhiDyAtomGridIm;
  double *vnlPhiDzAtomGridRe = pseudoReal->vnlPhiDzAtomGridRe;
  double *vnlPhiDzAtomGridIm = pseudoReal->vnlPhiDzAtomGridIm;
    
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
      // remove singularity
      rProj = gridAtomNbhd[iGrid*3]*gridAtomNbhd[iGrid*3]+
	      gridAtomNbhd[iGrid*3+1]*gridAtomNbhd[iGrid*3+1];
      if(rProj<1.0e-14){
	gridAtomNbhd[iGrid*3] = 1.0e-14;
      }
    }
    calcTrig(gridAtomNbhd,numGrid,trig);
    countRad = 0;
    dTheta = (double*)cmalloc(3*numGrid*sizeof(double));
    dPhi = (double*)cmalloc(3*numGrid*sizeof(double));
    calcAngleDeriv(trig,gridAtomNbhd,dTheta,dPhi,numGrid);
    for(l=0;l<atomLMax[atomType];l++){
      /* calculate sherical harmonic */
      ylm = (double*)cmalloc((2*l+1)*numGrid*sizeof(double));
      ylmTheta = (double*)cmalloc((2*l+1)*numGrid*sizeof(double));
      ylmPhi = (double*)cmalloc((2*l+1)*numGrid*sizeof(double));
      ylmDx = (double*)cmalloc((2*l+1)*numGrid*sizeof(double));
      ylmDy = (double*)cmalloc((2*l+1)*numGrid*sizeof(double));
      ylmDz = (double*)cmalloc((2*l+1)*numGrid*sizeof(double));
      calcSpHarm(ylm,l,gridAtomNbhd,numGrid,trig);
      calcSpHarmDeriv(ylm,ylmTheta,ylmPhi,l,gridAtomNbhd,numGrid,trig);
      for(m=0;m<=l;m++){
	if(m!=0){
	  ylmShift = (m*2-1)*numGrid;
	  for(iGrid=0;iGrid<numGrid;iGrid++){
	    ylmDx[ylmShift+iGrid*2] = ylmTheta[ylmShift+iGrid*2]*dTheta[iGrid*3]+
				      ylmPhi[ylmShift+iGrid*2]*dPhi[iGrid*3];
            ylmDx[ylmShift+iGrid*2+1] = ylmTheta[ylmShift+iGrid*2+1]*dTheta[iGrid*3]+
                                      ylmPhi[ylmShift+iGrid*2+1]*dPhi[iGrid*3];
            ylmDy[ylmShift+iGrid*2] = ylmTheta[ylmShift+iGrid*2]*dTheta[iGrid*3+1]+
                                      ylmPhi[ylmShift+iGrid*2]*dPhi[iGrid*3+1];
            ylmDy[ylmShift+iGrid*2+1] = ylmTheta[ylmShift+iGrid*2+1]*dTheta[iGrid*3+1]+
                                      ylmPhi[ylmShift+iGrid*2+1]*dPhi[iGrid*3+1];
            ylmDz[ylmShift+iGrid*2] = ylmTheta[ylmShift+iGrid*2]*dTheta[iGrid*3+2]+
                                      ylmPhi[ylmShift+iGrid*2]*dPhi[iGrid*3+2];
            ylmDz[ylmShift+iGrid*2+1] = ylmTheta[ylmShift+iGrid*2+1]*dTheta[iGrid*3+2]+
                                      ylmPhi[ylmShift+iGrid*2+1]*dPhi[iGrid*3+2];
	  }//endfor iGrid
	}
	else{
          ylmDx[iGrid] = ylmTheta[iGrid]*dTheta[iGrid*3]+ylmPhi[iGrid]*dPhi[iGrid*3];  
          ylmDy[iGrid] = ylmTheta[iGrid]*dTheta[iGrid*3+1]+ylmPhi[iGrid]*dPhi[iGrid*3+1];
          ylmDz[iGrid] = ylmTheta[iGrid]*dTheta[iGrid*3+2]+ylmPhi[iGrid]*dPhi[iGrid*3+2];
	}//endif m
      }//endfor m
      for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
        radIndex = atomRadMap[atomType][countRad+iRad];
        /* calculate radial fun u*phi */
        calcRadFun(gridAtomNbhd,radIndex,pseudoReal,radFun,numGrid);
        calcRadFunDev(gridAtomNbhd,radIndex,pseudoReal,radDevFun,numGrid);
        for(m=0;m<=l;m++){
          if(m!=0){
            ylmShift = (m*2-1)*numGrid;
            for(iGrid=0;iGrid<numGrid;iGrid++){
	      vnlPhiDxAtomGridRe[gridShiftRe+iGrid] = ylmDx[ylmShift+iGrid*2]*radFun[iGrid]+ylm[ylmShift+iGrid*2]*radDevFun[3*iGrid];
              vnlPhiDxAtomGridIm[gridShiftIm+iGrid] = ylmDx[ylmShift+iGrid*2+1]*radFun[iGrid]+ylm[ylmShift+iGrid*2+1]*radDevFun[3*iGrid];
              vnlPhiDyAtomGridRe[gridShiftRe+iGrid] = ylmDy[ylmShift+iGrid*2]*radFun[iGrid]+ylm[ylmShift+iGrid*2]*radDevFun[3*iGrid+1];
              vnlPhiDyAtomGridIm[gridShiftIm+iGrid] = ylmDy[ylmShift+iGrid*2+1]*radFun[iGrid]+ylm[ylmShift+iGrid*2+1]*radDevFun[3*iGrid+1];
              vnlPhiDzAtomGridRe[gridShiftRe+iGrid] = ylmDz[ylmShift+iGrid*2]*radFun[iGrid]+ylm[ylmShift+iGrid*2]*radDevFun[3*iGrid+2];
              vnlPhiDzAtomGridIm[gridShiftIm+iGrid] = ylmDz[ylmShift+iGrid*2+1]*radFun[iGrid]+ylm[ylmShift+iGrid*2+1]*radDevFun[3&iGrid+2];
            }//endfor iGrid
            gridShiftRe += numGrid;
            gridShiftIm += numGrid;
          }else{
            for(iGrid=0;iGrid<numGrid;iGrid++){
	      vnlPhiDxAtomGridRe[gridShiftRe+iGrid] = ylmDx[ylmShift+iGrid*2]*radFun[iGrid]+ylm[ylmShift+iGrid*2]*radDevFun[3*iGrid];
              vnlPhiDyAtomGridRe[gridShiftRe+iGrid] = ylmDy[ylmShift+iGrid*2]*radFun[iGrid]+ylm[ylmShift+iGrid*2]*radDevFun[3*iGrid+1];
              vnlPhiDzAtomGridRe[gridShiftRe+iGrid] = ylmDz[ylmShift+iGrid*2]*radFun[iGrid]+ylm[ylmShift+iGrid*2]*radDevFun[3*iGrid+2];
            }
            gridShiftRe += numGrid;
          }//endif m
        }//endfor m
      }//endfor iRad
      countRad += atomLRadNum[atomType][l];
      free(ylm);
      free(ylmTheta);
      free(ylmPhi);
      free(ylmDx);
      free(ylmDy);
      free(ylmDz);
    }//endfor l
    free(dTheta);
    free(dPhi);
  }//endfor iAtom
  /*
  for(iGrid=0;iGrid<numGridTot;iGrid++){
    printf("55555 %.16lg\n",testwfReal[iGrid]);
  }
  fflush(stdout);
  exit(0);
  */
  free(radFun);
  free(radDevFun);
  free(trig);
  free(gridAtomNbhd);

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcSpHarmDeriv(double *ylm,double *ylmTheta, double *ylmPhi,int l,
		     double *gridAtomNbhd,int numGrid,double *trig)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Spheircal harmonic values at each grid point.                         */
/* For each l, array ylmTheta(Phi) stores yl0,yl1(Re),yl1(Im)...         */
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
  // trig 0:sin(theta) 1:cos(theta) 2:sin(phi) 3:cos(phi)
  

  switch(l){
    case 0: //Y00
      for(iGrid=0;iGrid<numGrid;iGrid++){
	ylmTheta[iGrid] = 0.0;
	ylmPhi[iGrid] = 0.0;
      }
      break;
    case 1:
      for(iGrid=0;iGrid<numGrid;iGrid++){ //Y10
	ylmTheta[iGrid] = -pre10*trig[iGrid*4];
	ylmPhi[iGrid] = 0.0;
      }
      for(iGrid=0;iGrid<numGrid;iGrid++){ //Y11
        ind = numGrid+iGrid*2;
	ylmTheta[ind] = pre11*trig[iGrid*4+1]*trig[iGrid*4+3];
	ylmTheta[ind+1] = pre11*trig[iGrid*4+1]*trig[iGrid*4+2];
	ylmPhi[ind] = -pre11*trig[iGrid*4]*trig[iGrid*4+2];
	ylmPhi[ind+1] = pre11*trig[iGrid*4]*trig[iGrid*4+3];
      }
      break;
    case 2:
      for(iGrid=0;iGrid<numGrid;iGrid++){ //Y20
        ylmTheta[iGrid] = -pre20*6.0*trig[iGrid*4]*trig[iGrid*4+1];
	ylmPhi[iGrid] = 0.0;
      }
      for(iGrid=0;iGrid<numGrid;iGrid++){ //Y21
        ind = numGrid+2*iGrid;
	ylmTheta[ind] = pre21*(trig[iGrid*4+1]*trig[iGrid*4+1]-
			trig[iGrid*4]*trig[iGrid*4])*trig[iGrid*4+3];
	ylmTheta[ind+1] = pre21*(trig[iGrid*4+1]*trig[iGrid*4+1]-
			  trig[iGrid*4]*trig[iGrid*4])*trig[iGrid*4+2];
	ylmPhi[ind] = -pre21*trig[iGrid*4+1]*trig[iGrid*4]*trig[iGrid*4+2];
	ylmPhi[ind+1] = pre21*trig[iGrid*4+1]*trig[iGrid*4]*trig[iGrid*4+3];
      }
      for(iGrid=0;iGrid<numGrid;iGrid++){
        ind = numGrid*3+iGrid*2;
	ylmTheta[ind] = pre22*2.0*trig[iGrid*4]*trig[iGrid*4+1]*trig[iGrid*4+3];
	ylmTheta[ind] = pre22*2.0*trig[iGrid*4]*trig[iGrid*4+1]*trig[iGrid*4+2];
	ylmPhi[ind] = -pre22*trig[iGrid*4]*trig[iGrid*4]*trig[iGrid*4+2];
	ylmPhi[ind+1] = pre22*trig[iGrid*4]*trig[iGrid*4]*trig[iGrid*4+3];
      }
      break;
  }
/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcAngleDeriv(double *trig,double *gridAtomNbhd,double *dTheta,double *dPhi,
		    int numGrid)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Spheircal harmonic values at each grid point.                         */
/* For each l, array ylmTheta(Phi) stores yl0,yl1(Re),yl1(Im)...         */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int iGrid;
  double r,r2,rProj2,rProj;
  double x,y,z;
  double pre1,pre2;

  for(iGrid=0;iGrid<numGrid;iGrid++){
    x = gridAtomNbhd[3*iGrid];
    y = gridAtomNbhd[3*iGrid+1];
    z = gridAtomNbhd[3*iGrid+2];
    rProj2 = x*x+y*y;
    if(rProj2<1.0e-30){
      x += 1.0e-14;
      rProj2 = x*x+y*y;
    }
    r2 = x*x+y*y+z*z;
    r = sqrt(r2);
    pre1 = -trig[4*iGrid+1]/trig[4*iGrid]/r2;
    dTheta[3*iGrid] = pre1*x;
    dTheta[3*iGrid+1] = pre1*y;
    dTheta[3*iGrid+2] = trig[4*iGrid]/r;
    dPhi[3*iGrid] = y/rProj2;
    dPhi[3*iGrid+1] = -x/rProj2;
    dPhi[3*iGrid+2] = 0.0;
  }
/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcRadFunDev(double *gridAtomNbhd,int radIndex,PSEUDO_REAL *pseudoReal,
		   double *radDevFun,int numGrid)
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

  double *vps0 = pseudoReal->vpsDevReal0;
  double *vps1 = pseudoReal->vpsDevReal1;
  double *vps2 = pseudoReal->vpsDevReal2;
  double *vps3 = pseudoReal->vpsDevReal3;

  for(iGrid=0;iGrid<numGrid;iGrid++){
    x = gridAtomNbhd[iGrid*3];
    y = gridAtomNbhd[iGrid*3+1];
    z = gridAtomNbhd[iGrid*3+2];
    r = sqrt(x*x+y*y+z*z);
    gridInd = (int)((r-rMin)/dr)+1;
    r0 = (gridInd-1)*dr+rMin;
    h = r-r0;
    interpInd = interpGridSt+gridInd;
    radDevFun[iGrid] = ((vps3[interpInd]*h+vps2[interpInd])*h+vps1[interpInd])*h+vps0[interpInd];
    //if(r<1.0e-10)printf("rrrrrrrrrrrrrr %i %.16lg %.16lg\n",iGrid,r,radFun[iGrid]);
  }

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/




