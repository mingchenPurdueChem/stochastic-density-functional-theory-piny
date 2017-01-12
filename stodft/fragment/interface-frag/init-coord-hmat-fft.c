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
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_interface_frag_local.h"


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
  CLATOMS_POS *clatomsPos	= &(class->clatoms_pos[ip_now]);
  CLATOMS_POS *clatomsPosMini	= &(classMini->clatoms_pos[1]);
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
  double *xTemp,*yTemp,*zTemp;
  double *hmat	= cell->hmat;
  double *hmati = cell->hmati;
  double xRef,yRef,zRef;

/*======================================================================*/
/* I) Pass the coordinates                                              */
  
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xMini[iAtom] = x[atomFragMap[iAtom-1]];
    yMini[iAtom] = y[atomFragMap[iAtom-1]];
    zMini[iAtom] = z[atomFragMap[iAtom-1]];
  }
  
/*======================================================================*/
/* II) Remove PBC and shift the geometric center to 0                   */

  // We assume our fragment is small
  xRef = xMini[1];
  yRef = yMini[1];
  zRef = zMini[1];
  
  xDiff = (double*)cmalloc((numAtomFrag+1)*sizeof(double));
  yDiff = (double*)cmalloc((numAtomFrag+1)*sizeof(double));
  zDiff = (double*)cmalloc((numAtomFrag+1)*sizeof(double));
  xTemp = (double*)cmalloc((numAtomFrag+1)*sizeof(double));
  yTemp = (double*)cmalloc((numAtomFrag+1)*sizeof(double));
  zTemp = (double*)cmalloc((numAtomFrag+1)*sizeof(double));


  // Shift the first atom to the center of the box
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xDiff[iAtom] = xMini[iAtom]-xRef;
    yDiff[iAtom] = yMini[iAtom]-yRef;
    zDiff[iAtom] = zMini[iAtom]-zRef;
  }
  // Scale to cubic box
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xTemp[iAtom] = xDiff[iAtom]*hmati[1]+yDiff[iAtom]*hmati[4]+zDiff[iAtom]*hmati[7];
    yTemp[iAtom] = xDiff[iAtom]*hmati[2]+yDiff[iAtom]*hmati[5]+zDiff[iAtom]*hmati[8];
    zTemp[iAtom] = xDiff[iAtom]*hmati[3]+yDiff[iAtom]*hmati[6]+zDiff[iAtom]*hmati[9];
  }
  // Remove the PBC. I assume the fragment is no bigger then half of the box.
  // If it is, then a better way is to make it into even smaller pieces and remove pbc
  // for those small pieces. 
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    if(xTemp[iAtom]>0.5)xTemp[iAtom] -= 1.0;
    if(yTemp[iAtom]>0.5)yTemp[iAtom] -= 1.0;
    if(zTemp[iAtom]>0.5)zTemp[iAtom] -= 1.0;
    if(xTemp[iAtom]<-0.5)xTemp[iAtom] += 1.0;
    if(yTemp[iAtom]<-0.5)yTemp[iAtom] += 1.0;
    if(zTemp[iAtom]<-0.5)zTemp[iAtom] += 1.0;
  }

  // Rescale to Big box
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xDiff[iAtom] = xTemp[iAtom]*hmat[1]+yTemp[iAtom]*hmat[4]+zTemp[iAtom]*hmat[7];
    yDiff[iAtom] = xTemp[iAtom]*hmat[2]+yTemp[iAtom]*hmat[5]+zTemp[iAtom]*hmat[8];
    zDiff[iAtom] = xTemp[iAtom]*hmat[3]+yTemp[iAtom]*hmat[6]+zTemp[iAtom]*hmat[9];
    //printf("%lg %lg %lg\n",xDiff[iAtom],yDiff[iAtom],zDiff[iAtom]);
  }
  //debug
  /*
  FILE *testxyz = fopen("test.xyz","w");
  fprintf(testxyz,"3\n\n");
  fprintf(testxyz,"O %lg %lg %lg\n",xDiff[1]*BOHR,yDiff[1]*BOHR,zDiff[1]*BOHR);
  fprintf(testxyz,"H %lg %lg %lg\n",xDiff[2]*BOHR,yDiff[2]*BOHR,zDiff[2]*BOHR);
  fprintf(testxyz,"H %lg %lg %lg\n",xDiff[3]*BOHR,yDiff[3]*BOHR,zDiff[3]*BOHR);
  fclose(testxyz);
  */

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
  geoCnt[0] = xCnt+xRef;
  geoCnt[1] = yCnt+yRef;
  geoCnt[2] = zCnt+zRef;
  // Shift the center to 0
  
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xMini[iAtom] -= xCnt;
    yMini[iAtom] -= yCnt;
    zMini[iAtom] -= zCnt;
  }
  
  /*
  FILE *testxyz = fopen("test.xyz","w");
  fprintf(testxyz,"3\n\n");
  fprintf(testxyz,"O %lg %lg %lg\n",xMini[1]*BOHR,yMini[1]*BOHR,zMini[1]*BOHR);
  fprintf(testxyz,"H %lg %lg %lg\n",xMini[2]*BOHR,yMini[2]*BOHR,zMini[2]*BOHR);
  fprintf(testxyz,"H %lg %lg %lg\n",xMini[3]*BOHR,yMini[3]*BOHR,zMini[3]*BOHR);
  fclose(testxyz);
  */

  
/*======================================================================*/
/* II) Free local array				                        */
  free(xDiff);
  free(yDiff);
  free(zDiff); 
  free(xTemp);
  free(yTemp);
  free(zTemp);
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
  CLATOMS_POS *clatomsPosMini   = &(classMini->clatoms_pos[1]);

  int iAtom,iProj,iDim;
  int iGrid,jGrid,kGrid;
  int iFrag = fragInfo->iFrag;
  int numAtomFrag       = fragInfo->numAtomFragProc[iFrag];
  int numGridBigBoxC = cpParaFftPkg3dLgBigBox->nkf3;
  int numGridBigBoxB = cpParaFftPkg3dLgBigBox->nkf2;
  int numGridBigBoxA = cpParaFftPkg3dLgBigBox->nkf1;
  int numGridBigBox[3];
  int indexGrid[3],zeroGrid[3];
  int negativeGridNum;
  int numGridTotMiniBox;
  int index,indexa,indexb,indexc;
  int *numGridMiniBox = fragInfo->numGridFragDim[iFrag];

  double geoCntBox[3],geoCntDiff[3],gridSize[3];
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
/* I) Rescale the geometric center into unit box	                */
  //geoCnt[0] = 24.4033431368*-0.01;

  geoCntBox[0] = geoCnt[0]*hmati[1]+geoCnt[1]*hmati[4]+geoCnt[2]*hmati[7];
  geoCntBox[1] = geoCnt[0]*hmati[2]+geoCnt[1]*hmati[5]+geoCnt[2]*hmati[8];
  geoCntBox[2] = geoCnt[0]*hmati[3]+geoCnt[1]*hmati[6]+geoCnt[2]*hmati[9];


  /*
  geoCntBox[0] -= floor(geoCntBox[0]);
  geoCntBox[1] -= floor(geoCntBox[1]);
  geoCntBox[2] -= floor(geoCntBox[2]);
  */

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
  
  //printf("geoCntBox %lg %lg %lg\n",geoCntBox[0],geoCntBox[1],geoCntBox[2]);


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

  //printf("aGridLen %lg bGridLen %lg cGridLen %lg\n",aGridLen,bGridLen,cGridLen);
  
  //double boxlen = sqrt(hmat[1]*hmat[1]+hmat[2]*hmat[2]+hmat[3]*hmat[3]);
  //printf("genCnt %lg %lg %lg\n",geoCntBox[0]*boxlen,geoCntBox[1]*boxlen,geoCntBox[2]*boxlen);

  //update the center to the grid point
  geoCntBox[0] = aGrid[0]*indexGrid[0]+bGrid[0]*indexGrid[1]+cGrid[0]*indexGrid[2];
  geoCntBox[1] = aGrid[1]*indexGrid[0]+bGrid[1]*indexGrid[1]+cGrid[1]*indexGrid[2];
  geoCntBox[2] = aGrid[2]*indexGrid[0]+bGrid[2]*indexGrid[1]+cGrid[2]*indexGrid[2];

  // Reshift mini coords
  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xMini[iAtom] += geoCnt[0]-geoCntBox[0];
    yMini[iAtom] += geoCnt[1]-geoCntBox[1];
    zMini[iAtom] += geoCnt[2]-geoCntBox[2];
  }
/*--------------------------------------------------------------------------*/
/*  Along c direction					                    */
  negativeLength = 0.0;
  distProjAxis = calcMiniBoxLength(numAtomFrag,&xMini[1],&yMini[1],&zMini[1],
				    aNorm,bNorm,cNorm,skinFrag,&negativeLength);
  
  numGridMiniBox[2] = (int)(distProjAxis/cGridLen);
  if(numGridMiniBox[2]%2!=0)numGridMiniBox[2] += 1;
  numGridMiniBox[2] += 2; 
  numGridMiniBox[2] = 72;
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


  numGridMiniBox[1] = (int)(distProjAxis/bGridLen);
  if(numGridMiniBox[1]%2!=0)numGridMiniBox[1] += 1;
  numGridMiniBox[1] += 2;
  numGridMiniBox[1] = 72;
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

  numGridMiniBox[0] = (int)(distProjAxis/aGridLen);
  if(numGridMiniBox[0]%2!=0)numGridMiniBox[0] += 1;
  numGridMiniBox[0] += 2;
  numGridMiniBox[0] = 72;
  negativeGridNum = numGridMiniBox[0]/2;
  zeroGrid[0] = indexGrid[0]-negativeGridNum;
  for(iDim=0;iDim<3;iDim++)zeroShift[iDim] -= negativeGridNum*aGrid[iDim];
  if(zeroGrid[0]<0)zeroGrid[0] += numGridBigBoxA;
  // round it if necessary

/*--------------------------------------------------------------------------*/
/*  shift the zero point from central Grid to corner Grid                   */

  /*
  double zeroPoint[3];
  zeroPoint[0] = zeroGrid[0]*aGrid[0];
  zeroPoint[1] = zeroGrid[1]*bGrid[1];
  zeroPoint[2] = zeroGrid[2]*cGrid[2];
  printf("zero Point %lg %lg %lg\n",zeroPoint[0],zeroPoint[1],zeroPoint[2]);
  */

  for(iAtom=1;iAtom<=numAtomFrag;iAtom++){
    xMini[iAtom] -= zeroShift[0];
    yMini[iAtom] -= zeroShift[1];
    zMini[iAtom] -= zeroShift[2];
    //printf("Final mini x %lg y %lg z %lg\n",xMini[iAtom],yMini[iAtom],zMini[iAtom]);
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
	if(indexc<0)indexc += numGridBigBox[2];
	if(indexc>=numGridBigBox[2])indexc -= numGridBigBox[2];
	if(indexb<0)indexb += numGridBigBox[1];
	if(indexb>=numGridBigBox[1])indexb -= numGridBigBox[1];
        if(indexa<0)indexa += numGridBigBox[0];
        if(indexa>=numGridBigBox[0])indexa -= numGridBigBox[0];
	fragInfo->gridMapProc[iFrag][index] = indexc*numGridBigBox[1]*numGridBigBox[0]+
		    indexb*numGridBigBox[0]+indexa;
      }//endfor kGrid
    }//endfor jGrid
  }//endfor iGrid

  //debug
  for(iGrid=0;iGrid<numGridTotMiniBox;iGrid++){
    fragInfo->gridMapProc[iFrag][iGrid] = iGrid;
  }

/*======================================================================*/
/* III) Get the mini cell matrix                                        */

  hmatMini[1] = numGridMiniBox[0]*aGrid[0];
  hmatMini[2] = numGridMiniBox[0]*aGrid[1];
  hmatMini[3] = numGridMiniBox[0]*aGrid[2];
  hmatMini[4] = numGridMiniBox[1]*bGrid[0];
  hmatMini[5] = numGridMiniBox[1]*bGrid[1];
  hmatMini[6] = numGridMiniBox[1]*bGrid[2];
  hmatMini[7] = numGridMiniBox[2]*cGrid[0];
  hmatMini[8] = numGridMiniBox[2]*cGrid[1];
  hmatMini[9] = numGridMiniBox[2]*cGrid[2];

  //printf("umGridMiniBox %i %i %i\n",numGridMiniBox[0],numGridMiniBox[1],numGridMiniBox[2]);
  //printf("agrid %lg bgrid %lg cgrid %lg\n",aGrid[0],bGrid[1],cGrid[2]);

  printf("hmatMini %lg %lg %lg\n",hmatMini[1],hmatMini[2],hmatMini[3]);
  printf("hmatMini %lg %lg %lg\n",hmatMini[4],hmatMini[5],hmatMini[6]);
  printf("hmatMini %lg %lg %lg\n",hmatMini[7],hmatMini[8],hmatMini[9]);


  //Let's give a final test before this is over
  
  FILE *testxyz = fopen("test.xyz","w");
  fprintf(testxyz,"3\n\n");
  fprintf(testxyz,"O %lg %lg %lg\n",xMini[1]*BOHR,yMini[1]*BOHR,zMini[1]*BOHR);
  fprintf(testxyz,"H %lg %lg %lg\n",xMini[2]*BOHR,yMini[2]*BOHR,zMini[2]*BOHR);
  fprintf(testxyz,"H %lg %lg %lg\n",xMini[3]*BOHR,yMini[3]*BOHR,zMini[3]*BOHR);
  fclose(testxyz);
  printf("%lg %lg %lg\n",hmatMini[1]*BOHR,hmatMini[2]*BOHR,hmatMini[3]*BOHR);
  printf("%lg %lg %lg\n",hmatMini[4]*BOHR,hmatMini[5]*BOHR,hmatMini[6]*BOHR);
  printf("%lg %lg %lg\n",hmatMini[7]*BOHR,hmatMini[8]*BOHR,hmatMini[9]*BOHR);
  
  //printf("zeroGrid %i %i %i\n",zeroGrid[0],zeroGrid[1],zeroGrid[2]);
  //printf("indexGrid %i %i %i\n",indexGrid[0],indexGrid[1],indexGrid[2]);
  //exit(0);
  

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
  distProjAxis = distVert/fabs(dotProd);
  *negativelength = fabs(projMin/dotProd);
  cfree(projList);

  return distProjAxis;

/*==========================================================================*/  
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double normalized3d(double *x)
/*========================================================================*/
/*             Begin Routine                                              */
{/*Begin subprogram: */
/*========================================================================*/
/*             Local variable declarations                                */
/*------------------------------------------------------------------------*/
  double norm = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  x[0] /= norm;
  x[1] /= norm;
  x[2] /= norm;
  
  return norm;

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/


