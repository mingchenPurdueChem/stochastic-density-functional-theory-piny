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
  int iAtom;
  
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
/*  Parse: note there is a noncommuting order of calls                      */
/*==========================================================================*/
void findCnt(GENERAL_DATA *generalData,CLASS *class,CP *cp,
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
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo            = stodftInfo->fragInfo; 
  CLATOMS_POS *clatomsPosMini   = &(cpMini->clatoms_pos[1]);


  int iAtom;
  int iFrag = fragInfo->iFrag;
  int numAtomFrag       = fragInfo->numAtomFragProc[iFrag];
  int numGridBigBoxC = cpParaFftPkg3dLgBigBox->nkf3;
  int numGridBigBoxB = cpParaFftPkg3dLgBigBox->nkf2;
  int numGridBigBoxA = cpParaFftPkg3dLgBigBox->nkf1;
  int numGridBigBox[3];
  int indexGrid[3];

  double geoCntBox[3],gridSize[3];
  double 

  double *hmat  = cell->hmat;
  double *hmati = cell->hmati;
  double *skinFrag = fragInfo->skinFragBox[iFrag];
  double *xMini = clatomsPosMini->x;
  double *yMini = clatomsPosMini->y;
  double *zMini = clatomsPosMini->z;

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
 
  
  


}/*end routine*/
/*==========================================================================*/

