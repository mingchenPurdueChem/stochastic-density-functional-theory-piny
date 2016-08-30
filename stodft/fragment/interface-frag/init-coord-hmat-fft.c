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
		   int ip_now)
/*========================================================================*/
/*             Begin Routine                                              */
{/*Begin subprogram: */
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

  double *x = clatomsPos->x;
  double *y = clatomsPos->y;
  double *z = clatomsPos->z;
  double *xMini = clatomsPosMini->x;
  double *yMini	= clatomsPosMini->y;
  double *zMini = clatomsPosMini->z;
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


/*------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

