/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: init-coord-hmat-fft.c                          */
/*                                                                          */
/* This routine pass atomic coordinates to fragments and initialize FFT	    */
/* grids and h-matrix.				                            */
/* 1) Pass atomic coordinates to fragments and remove pbc of the coords     */
/* 2) Find the FFT grid point that is closest to center of atomic positions */
/* 3) Use the point as reference point, generate a box (parallel to the big */
/*    h-matrix that contains all fragments with skin, with integer number   */
/*    FFT grid points on each dimension.				    */
/* 4) Generate FFT grids in real space, and use this to initialize all FFT  */
/*    variables.							    */
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
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initCoordHmatFFT(GENERAL_DATA *generalData,CLASS *class,CP *cp,
		      GENERAL_DATA *generalDataMini,CLASS *classMini,CP *cpMini)
/*========================================================================*/
/*             Begin Routine                                              */
/*************************************************************************/
{/*Begin subprogram: */
/*************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
  double geoCnt[3];
  passAtomCoord(generalData,class,cp,generalDataMini,classMini,cpMini,1,geoCnt);
  
  initFFTMap(generalData,class,cp,generalDataMini,classMini,cpMini,1,geoCnt);


/*------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Parse: note there is a noncommuting order of calls                      */
/*==========================================================================*/
void passAtomCoord(GENERAL_DATA *generalData,CLASS *class,CP *cp,
		   GENERAL_DATA *generalDataMini,CLASS *classMini,CP *cpMini,
		   int ip_now,double *geoCnt)
/*========================================================================*/
/*             Begin Routine                                              */
/*************************************************************************/
{/*Begin subprogram: */
/*************************************************************************/
/* This routine takes the coords of a fragment, remove the pbc and       */
/* center the fragment at the geometric center.				 */
/*************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
  STODFTINFO *stodftInfo	= cp->stodftInfo;
  FRAGINFO *fragInfo		= stodftInfo->fragInfo;
  CLATOMS_POS *clatomsPos	= &(cp->clatoms_pos[ip_now]);
  CLATOMS_POS *clatomsPosMini	= &(cpMini->clatoms_pos[1]);
  CELL *cell			= &(generalData->cell);
  

  int numFragProc	= fragInfo->numFragProc;
  int iFrag		= fragInfo->iFrag;
  int iAtom,iProj,iDim;
  
  int numAtomFrag	= fragInfo->numAtomFragProc[iFrag];
  int *atomFragMap	= fragInfo->atomFragMapProc[iFrag];

  double xCnt,yCnt,zCnt;
  double *x = clatomsPos->x;
  double *y = clatomsPos->y;
  double *z = clatomsPos->z;
  double *xMini = clatomsPosMini->x;
  double *yMini	= clatomsPosMini->y;
  double *zMini = clatomsPosMini->z;
  double *xDiff,*yDiff,*zDiff;
  double *hmat	= cell->hmat;
  double *hmati = cell->hmati;
  double xRef,yRef,zRef;

/*======================================================================*/
/* I) Pass the coordinates                                              */
  
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xMini[iAtom] = x[atomFragMap[iAtom]];
    yMini[iAtom] = y[atomFragMap[iAtom]];
    zMini[iAtom] = z[atomFragMap[iAtom]];
  }
  
/*======================================================================*/
/* II) Remove PBC and shift the geometric center to 0                   */

  // We assume our fragment is small
  xRef = xMini[1];
  yRef = yMini[1];
  zRef = zMini[1];
  
  double *xDiff = (double*)cmalloc((numAtomFrag+1)*sizeof(double));
  double *yDiff = (double*)cmalloc((numAtomFrag+1)*sizeof(double));
  double *zDiff = (double*)cmalloc((numAtomFrag+1)*sizeof(double));

  // Shift the first atom to the center of the box
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xDiff[iAtom] = xMini[iAtom]-xRef[iAtom];
    yDiff[iAtom] = yMini[iAtom]-yRef[iAtom];
    zDiff[iAtom] = zMini[iAtom]-zRef[iAtom];
  }
  // Scale to cubic box
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xMini[iAtom] = xDiff[iAtom]*hmati[1]+yDiff[iAtom]*hmati[4]+zDiff[iAtom]*hmati[7];
    yMini[iAtom] = xDiff[iAtom]*hmati[2]+yDiff[iAtom]*hmati[5]+zDiff[iAtom]*hmati[8];
    zMini[iAtom] = xDiff[iAtom]*hmati[3]+yDiff[iAtom]*hmati[6]+zDiff[iAtom]*hmati[9];
  }
  // Remove the PBC. I assume the fragment is no bigger then half of the box.
  // If it is, then a better way is to make it into even smaller pieces and remove pbc
  // for those small pieces. 
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    if(xMini[iAtom]>0.5)xMini[iAtom] -= 1.0;
    if(yMini[iAtom]>0.5)yMini[iAtom] -= 1.0;
    if(zMini[iAtom]>0.5)zMini[iAtom] -= 1.0;
    if(xMini[iAtom]<-0.5)xMini[iAtom] += 1.0;
    if(yMini[iAtom]<-0.5)yMini[iAtom] += 1.0;
    if(zMini[iAtom]<-0.5)zMini[iAtom] += 1.0;
  }

  // Rescale to Big box
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xDiff[iAtom] = xMini[iAtom]*hmat[1]+yMini[iAtom]*hmat[4]+zMini[iAtom]*hmat[7];
    yDiff[iAtom] = xMini[iAtom]*hmat[2]+yMini[iAtom]*hmat[5]+zMini[iAtom]*hmat[8];
    zDiff[iAtom] = xMini[iAtom]*hmat[3]+yMini[iAtom]*hmat[6]+zMini[iAtom]*hmat[9];
  }
  // Copy x/y/zDiff to x/y/zMini
  memcpy(&xMini[1],&xDiff[1],numAtomFrag*sizeof(double));
  memcpy(&yMini[1],&yDiff[1],numAtomFrag*sizeof(double));
  memcpy(&zMini[1],&zDiff[1],numAtomFrag*sizeof(double));

  // Get the geometric center
  xCnt = 0.0;
  yCnt = 0.0;
  zCnt = 0.0;
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xCnt += xMini[iAtom];
    yCnt += yMini[iAtom];
    zCnt += zMini[iAtom];
  }
  xCnt /= numAtomFrag;
  yCnt /= numAtomFrag;
  zCnt /= numAtomFrag;
  geoCnt[0] = xCnt;
  geoCnt[1] = yCnt;
  geoCnt[2] = zCnt;
  // Shift the center to 0
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xMini[iAtom] -= xCnt;
    yMini[iAtom] -= yCnt;
    zMini[iAtom] -= zCnt;
  }
  
/*======================================================================*/
/* II) Free local array				                        */
  free(xDiff);
  free(yDiff);
  free(zDiff); 
/*------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initFFTMap(GENERAL_DATA *generalData,CLASS *class,CP *cp,
                   GENERAL_DATA *generalDataMini,CLASS *classMini,CP *cpMini,
                   int ip_now,double *geoCnt)
/*========================================================================*/
/*             Begin Routine                                              */
{/*Begin subprogram: */
/*************************************************************************/
/* This routine find the FFT grid point in big box that is closest to    */
/* the geometric center of the fragment and set that point as center     */
/*************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
/*------------------------------------------------------------------------*/
  PARA_FFT_PKG3D *cpParaFftPkg3dLgBigBox = &(cp->cp_para_fft_pkg3d_lg);
  CELL *cell = &(generalData->cell);
  CELL *cellMini = &(generalDataMini->cell);
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo            = stodftInfo->fragInfo; 
  CLATOMS_POS *clatomsPosMini   = &(cpMini->clatoms_pos[1]);

  int iAtom,iProj,iDim;
  int iGrid,jGrid,kGrid;
  int iFrag = fragInfo->iFrag;
  int numAtomFrag       = fragInfo->numAtomFragProc[iFrag];
  int numGridBigBoxC = cpParaFftPkg3dLgBigBox->nkf3;
  int numGridBigBoxB = cpParaFftPkg3dLgBigBox->nkf2;
  int numGridBigBoxA = cpParaFftPkg3dLgBigBox->nkf1;
  int numGridBigBox[3],numGridMiniBox[3];
  int indexGrid[3],zeroGrid[3];
  int negativeGridNum;
  int numGridTotMiniBox;
  int index,indexa,indexb,indexc;

  double geoCntBox[3],gridSize[3],zeroShift[3];
  double zeroShift[3] = {0};
  double aBig[3],bBig[3],cBig[3];
  double aGrid[3],bGrid[3],cGrid[3];
  double aNorm[3],bNorm[3],cNorm[3];
  double aGridLen,bGridLen,cGridLen;
  double crossProd[3];
  double norm,dotProd;
  double projMin,projMax;
  double distVert,distProjAxis;
  double negativeLength;

  double *hmat  = cell->hmat;
  double *hmati = cell->hmati;
  double *skinFrag = fragInfo->skinFragBox[iFrag];
  double *xMini = clatomsPosMini->x;
  double *yMini = clatomsPosMini->y;
  double *zMini = clatomsPosMini->z;
  double *hmatMini = cellMini->hmat;
  double *projList;

/*======================================================================*/
/* I) Put the genmetric center into the box if necessary                */

  geoCntBox[0] = geoCnt[0]*hmati[1]+geoCnt[1]*hmati[4]+geoCnt[2]*hmati[7];
  geoCntBox[1] = geoCnt[0]*hmati[2]+geoCnt[1]*hmati[5]+geoCnt[2]*hmati[8];
  geoCntBox[2] = geoCnt[0]*hmati[3]+geoCnt[1]*hmati[6]+geoCnt[2]*hmati[9];

  geoCntBox[0] -= floor(geoCntBox[0]);
  geoCntBox[1] -= floor(geoCntBox[1]);
  geoCntBox[2] -= floor(geoCntBox[2]);

/*======================================================================*/
/* II) Find the closest FFT grid point in the big box.                  */
 
  // Get the Bin index
  numGridBigBox[0] = cpParaFftPkg3dLgBigBox->nkf1; //a
  numGridBigBox[1] = cpParaFftPkg3dLgBigBox->nkf2; //b
  numGridBigBox[2] = cpParaFftPkg3dLgBigBox->nkf3; //c

  //gridSize[0] = 1.0/numGridBigBox[0];
  //gridSize[1] = 1.0/numGridBigBox[1];
  //gridSize[2] = 1.0/numGridBigBox[2];

  indexGrid[0] = NINT(geoCntBox[0]*numGridBigBox[0]);
  indexGrid[1] = NINT(geoCntBox[1]*numGridBigBox[1]);
  indexGrid[2] = NINT(geoCntBox[2]*numGridBigBox[2]);


/*======================================================================*/
/* II) Find the number of grid points on each dimension.                */
/*     e.g. We want to calculate # grid points along c direction. We	*/
/*     first shift our molecule to the center FFT grid we find in the   */
/*     last step. Then we calculate the signed distances between all	*/
/*     atoms and the <ab> surface. This is done by project the atom	*/
/*     positions along the aXb direction. Then for each projection, we  */
/*     add/substract skin value. Now we have 2*(atom number) values and */
/*     we pick the max(positive) and min(negative) value. The max-min   */
/*     is the vertical distance between box top and box botom. Finally, */
/*     we calculate the c length by using (max-min)/cos(<c,aXb>).       */
/*     The # of grid along c is a integer multiplication of grid	*/
/*     that just larger then c length calculated before.		*/
  aBig[0] = hmat[1];aBig[1] = hmat[2];aBig[2] = hmat[3];
  bBig[0] = hmat[4];bBig[1] = hmat[5];bBig[2] = hmat[6];
  cBig[0] = hmat[7];cBig[1] = hmat[8];cBig[2] = hmat[9];
  aGrid[0] = aBig[0]/numGridBigBox[0];
  aGrid[1] = aBig[1]/numGridBigBox[0];
  aGrid[2] = aBig[2]/numGridBigBox[0];
  bGrid[0] = bBig[0]/numGridBigBox[1];
  bGrid[1] = bBig[1]/numGridBigBox[1];
  bGrid[2] = bBig[2]/numGridBigBox[1];
  cGrid[0] = cBig[0]/numGridBigBox[2];
  cGrid[1] = cBig[1]/numGridBigBox[2];
  cGrid[2] = cBig[2]/numGridBigBox[2];
  
  for(iDim=0;iDim<3;iDim++){
    aNorm[iDim] = aGrid[iDim];
    bNorm[iDim] = bGrid[iDim];
    cNorm[iDim] = cGrid[iDim];
  } 
  aGridLen = normalized3d(aNorm);
  bGridLen = normalized3d(bNorm);
  cGridLen = normalized3d(cNorm);

  //update the center to the grid point
  geoCntBox[0] = aGrid[0]*indexGrid[0]+bGrid[0]*indexGrid[1]+cGrid[0]*indexGrid[2];
  geoCntBox[1] = aGrid[1]*indexGrid[0]+bGrid[1]*indexGrid[1]+cGrid[1]*indexGrid[2];
  geoCntBox[2] = aGrid[2]*indexGrid[0]+bGrid[2]*indexGrid[1]+cGrid[2]*indexGrid[2];

  // Reshift mini coords
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xMini[iAtom] -= geoCntBox[0];
    yMini[iAtom] -= geoCntBox[1];
    zMini[iAtom] -= geoCntBox[2];
  }
/*--------------------------------------------------------------------------*/
/*  Along c direction					                    */
  negativeLength = 0.0;
  distProjAxis = calcMiniBoxLength(numAtomFrag,&xMini[1],&yMini[1],&zMini[1],
				    aNorm,bNorm,cNorm,skinFrag,&negativeLength);
  
  numGridMiniBox[2] = (int)(distProjAxis/cGridLen);
  if(numGridMiniBox[2]%2!=0)numGridMinBox[2] += 1;
  numGridMiniBox[2] += 2; 
  negativeGridNum = numGridMiniBox[2]/2;
  zeroGrid[2] = indexGrid[2]-negativeGridNum;
  for(iDim=0;iDim<3;iDim++)zeroShift[iDim] -= negativeGridNum*cGrid[iDim];
  if(zeroGrid[2]<0)zeroGrid[2] += numGridBigBoxC;
  // round it if necessary

/*--------------------------------------------------------------------------*/
/*  Along b direction                                                       */

  negativeLength = 0.0;
  distProjAxis = calcMiniBoxLength(numAtomFrag,&xMini[1],&yMini[1],&zMini[1],
                                    cNorm,aNorm,bNorm,skinFrag,&negativeLength);
  numGridMiniBox[1] = (int)distProjAxis/bGridLen;
  if(numGridMiniBox[1]%2!=0)numGridMiniBox[1] += 1;
  numGridMiniBox[1] += 2;
  negativeGridNum = numGridMiniBox[1]/2;
  zeroGrid[1] = indexGrid[1]-negativeGridNum;
  for(iDim=0;iDim<3;iDim++)zeroShift[iDim] -= negativeGridNum*bGrid[iDim];
  if(zeroGrid[1]<0)zeroGrid[1] += numGridBigBoxB;
  // round it if necessary
/*--------------------------------------------------------------------------*/
/*  Along a direction                                                       */

  negativeLength = 0.0;
  distProjAxis = calcMiniBoxLength(numAtomFrag,&xMini[1],&yMini[1],&zMini[1],
                                    bNorm,cNorm,aNorm,skinFrag,&negativeLength);
  numGridMiniBox[0] = (int)distProjAxis/aGridLen;
  if(numGridMiniBox[0]%2!=0)numGridMiniBox[0] += 1;
  numGridMiniBox[0] += 2;
  negativeGridNum = numGridMiniBox[0]/2;
  zeroGrid[0] = indexGrid[0]-negativeGridNum;
  for(iDim=0;iDim<3;iDim++)zeroShift[iDim] -= negativeGridNum*aGrid[iDim];
  if(zeroGrid[0]<0)zeroGrid[0] += numGridBigBoxA;
  // round it if necessary

/*--------------------------------------------------------------------------*/
/*  shift the zero point from central Grid to corner Grid                   */

  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xMini[iAtom] -= zeroShift[0];
    yMini[iAtom] -= zeroShift[1];
    zMini[iAtom] -= zeroShift[2];
  }

/*======================================================================*/
/* II) Map the grid.					                */
  
  numGridTotMiniBox = numGridMiniBox[2]*numGridMiniBox[1]*numGridMiniBox[0];

  fragInfo->numGridFragProc[iFrag] = numGridTotMiniBox;
  fragInfo->gridMapProc[iFrag] = (int*)cmalloc(numGridTotMiniBox*sizeof(int));
 
  for(iGrid=0;iGrid<numGridMiniBox[2];iGrid++){//c
    for(jGrid=0;jGrid<numGridMiniBox[1];jGrid++){//b
      for(kGrid=0;kGrid<numGridMiniBox[0];kGrid++){//a
	index = iGrid*numGridMiniBox[1]*numGridMiniBox[0]
		+jGrid*numGridMiniBox[0]+kGrid;    
	indexc = zeroGrid[2]+iGrid;
	indexb = zeroGrid[1]+jGrid;
	indexa = zeroGrid[0]+kGrid;	
      }//endfor kGrid
    }//endfor jGrid
  }//endfor iGrid

/*======================================================================*/
/* III) Get the mini cell matrix                                        */

  hmatMini[1] = numGridMinBox[0]*aGrid[0];
  hmatMini[2] = numGridMinBox[0]*aGrid[1];
  hmatMini[3] = numGridMinBox[0]*aGrid[2];
  hmatMini[4] = numGridMinBox[1]*bGrid[0];
  hmatMini[5] = numGridMinBox[1]*bGrid[1];
  hmatMini[6] = numGridMinBox[1]*bGrid[2];
  hmatMini[7] = numGridMinBox[2]*cGrid[0];
  hmatMini[8] = numGridMinBox[2]*cGrid[1];
  hmatMini[9] = numGridMinBox[2]*cGrid[2];
  
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double calcMiniBoxLength(int numAtomFrag,double *x,double *y,double *z,
			 double *aNorm,double *bNorm,double *cNorm,
			 double *skinFrag,double *negativelength) 
/*========================================================================*/
/*             Begin Routine                                              */
{/*Begin subprogram: */
/*========================================================================*/
/*             Local variable declarations                                */
/*------------------------------------------------------------------------*/
  int iAtom,iProj;

  double projMin,projMax;
  double dotProd;
  double distVert,distProjAxis;
  double crossProd[3];
  double *projList = (double*)cmalloc(2*numAtomFrag*sizeof(double));;
  
  cross_product(aNorm,bNorm,crossProd);
  normalize3d(crossProd);
  // generate atom-plan distance list
  for(iAtom=0;iAtom<numAtomFrag;iAtom++){
    dotProd = x[iAtom]*crossProd[0]+y[iAtom]*crossProd[1]
              +z[iAtom]*crossProd[2];
    projList[iAtom*2] = dotProd+skinFrag[iAtom];
    projList[iAtom*2+1] = dotProd-skinFrag[iAtom];
  }

  // find the max/min element
  projMin = 1.0e10;
  projMax = -1.0e10;

  for(iProj=0;iProj<2*numAtomFrag;iProj++){
    if(projList[iProj]>projMax)projMax = projList[iProj];
    if(projList[iProj]<projMin)projMin = projList[iProj];
  }

  // Get the vertical distance between upper/lower surface
  distVert = projMax-projMin;
  // Get the c direction distance   
  dotProd = dot(cNorm,crossProd);
  distProjAxis = disVert/fabs(dotProd);
  *negativelength = fabs(projMin/dotProd);
  cfree(projList);

  return distProjAxis
  
}/*end routine*/
/*==========================================================================*/

