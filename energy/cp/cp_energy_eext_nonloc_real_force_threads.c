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
void nlppKBRealEnergyForceThreads(CP *cp,CLASS *class,GENERAL_DATA *generalData,
                           double *wfReal,PARA_FFT_PKG3D *cpParaFftPkg3d)
/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*************************************************************************/
/* Calculate nuclei force			 */
/*************************************************************************/
/*==========================================================================*/
/*               Local variable declarations                                */

  PSEUDO *pseudo = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CPOPTS *cpopts = &(cp->cpopts);
  CELL *cell = &(generalData->cell);
  CLATOMS_POS *clatoms_pos = &(class->clatoms_pos[1]);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  ATOMMAPS *atommaps = &(class->atommaps);
  STAT_AVG *stat_avg = &(generalData->stat_avg);
  CPEWALD *cpewald = &(cp->cpewald);
  
  int numAtomType = atommaps->natm_typ;
  int numAtom = clatoms_info->natm_tot;
  int numGrid;
  int numGridMax;
  int countRad;
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
  double a[3],b[3],c[3];
  double aGrid[3],bGrid[3],cGrid[3];
  double nucleiCoord[3];
  double *hmat = cell->hmat;
  double *fx = clatoms_pos->fx;
  double *fy = clatoms_pos->fy;
  double *fz = clatoms_pos->fz;
  double *wfNbhd;

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

/*======================================================================*/
/* I) Allocate local memory                                             */

  numGridMax = numGridNlppMap[0];
  for(iPart=0;iPart<numAtom;iPart++){
    if(numGridNlppMap[iPart]>numGridMax)numGridMax = numGridNlppMap[iPart];
  }
  wfNbhd = (double*)cmalloc(numGridMax*sizeof(double));
  vol = getdeth(hmat);
  volInv = 1.0/vol;
  volElem = vol/numGridTot;
/*======================================================================*/
/* II) Loop over iPart/irad/ipart/m                                 */

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
    //printf("numGrid %i\n",numGrid);
    if(numGrid>0){ //if numGrid=0, only local pp will be calculated
      for(iGrid=0;iGrid<numGrid;iGrid++){
        gridIndex = gridNlppMap[iAtom][iGrid];
        wfNbhd[iGrid] = wfReal[gridIndex];
      }//endfor iGrid
      for(l=0;l<atomLMax[atomType];l++){
        for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
          radIndex = atomRadMap[atomType][countRad+iRad];
          for(m=0;m<=l;m++){
            //calcDotNlpp(wfNbhd,radFun,&ylm[ylmShift],&dotRe,&dotIm);
            if(m!=0){
              dotRe = dotReAll[iAtom][countNlppRe+m];
              dotIm = dotImAll[iAtom][countNlppIm+m-1];
              //printf("dottttttttt %i %i %lg %i %lg\n",
              //       iAtom,countNlppRe+m,dotRe,countNlppIm+m-1,dotIm);
              // x    
              //printf("dddddx %lg %lg dy %lg %lg dz %lg %lg\n",vnlPhiDxAtomGridRe[gridShiftNowRe],vnlPhiDxAtomGridIm[gridShiftNowIm],vnlPhiDyAtomGridRe[gridShiftNowRe],vnlPhiDyAtomGridIm[gridShiftNowIm],vnlPhiDzAtomGridRe[gridShiftNowRe],vnlPhiDzAtomGridIm[gridShiftNowIm]);
              dotDevRe = ddotBlasWrapper(numGrid,wfNbhd,1,
		     &vnlPhiDxAtomGridRe[gridShiftNowRe],1)*volElem;
              dotDevIm = ddotBlasWrapper(numGrid,wfNbhd,1,
		         &vnlPhiDxAtomGridIm[gridShiftNowIm],1)*volElem;
              //printf("xxx m %i dotRe %lg dotIm %lg dotDevRe %lg dotDevIm %lg\n",m,dotRe,dotIm,dotDevRe,dotDevIm);
	      forceNlX += (dotRe*dotDevRe+dotIm*dotDevIm)*4.0*vpsNormList[radIndex]*volInv;
              // y
	      dotDevRe = ddotBlasWrapper(numGrid,wfNbhd,1,
                                         &vnlPhiDyAtomGridRe[gridShiftNowRe],1)*volElem;
              dotDevIm = ddotBlasWrapper(numGrid,wfNbhd,1,
                                         &vnlPhiDyAtomGridIm[gridShiftNowIm],1)*volElem;
              //printf("yyy m %i dotRe %lg dotIm %lg dotDevRe %lg dotDevIm %lg\n",m,dotRe,dotIm,dotDevRe,dotDevIm);

              forceNlY += (dotRe*dotDevRe+dotIm*dotDevIm)*4.0*vpsNormList[radIndex]*volInv;
              // z
              dotDevRe = ddotBlasWrapper(numGrid,wfNbhd,1,
                                         &vnlPhiDzAtomGridRe[gridShiftNowRe],1)*volElem;
              dotDevIm = ddotBlasWrapper(numGrid,wfNbhd,1,
                                         &vnlPhiDzAtomGridIm[gridShiftNowIm],1)*volElem;
              //printf("zzz m %i dotRe %lg dotIm %lg dotDevRe %lg dotDevIm %lg\n",m,dotRe,dotIm,dotDevRe,dotDevIm);

              forceNlZ += (dotRe*dotDevRe+dotIm*dotDevIm)*4.0*vpsNormList[radIndex]*volInv;
              gridShiftNowRe += numGrid;
              gridShiftNowIm += numGrid;
            }
            else{
              dotRe = dotReAll[iAtom][countNlppRe];
              //printf("dottttttttt %i %i %lg\n",
              //      iAtom,countNlppRe,dotRe);
              // x
              dotDevRe = ddotBlasWrapper(numGrid,wfNbhd,1,
		     &vnlPhiDxAtomGridRe[gridShiftNowRe],1)*volElem;
              //printf("xxx m %i dotRe %lg dotDevRe %lg\n",m,dotRe,dotDevRe);

              forceNlX += dotRe*dotDevRe*2.0*vpsNormList[radIndex]*volInv;
              // y
              dotDevRe = ddotBlasWrapper(numGrid,wfNbhd,1,
                                         &vnlPhiDyAtomGridRe[gridShiftNowRe],1)*volElem;
              //printf("yyy m %i dotRe %lg dotDevRe %lg\n",m,dotRe,dotDevRe);

              forceNlY += dotRe*dotDevRe*2.0*vpsNormList[radIndex]*volInv;
              // z
              dotDevRe = ddotBlasWrapper(numGrid,wfNbhd,1,
                                         &vnlPhiDzAtomGridRe[gridShiftNowRe],1)*volElem;
              //printf("zzz m %i dotRe %lg dotDevRe %lg\n",m,dotRe,dotDevRe);

              forceNlZ += dotRe*dotDevRe*2.0*vpsNormList[radIndex]*volInv;
              gridShiftNowRe += numGrid;
            }//endif m
          }//endfor m
          countNlppRe += l+1;
          countNlppIm += l;
        }//endfor iRad
        countRad += atomLRadNum[atomType][l];
      }//endfor l
    }//endif numGrid
    //if(iAtom==0)printf("forceNl %.16lg %.16lg %.16lg\n",forceNlX,forceNlY,forceNlZ);
    fx[iAtom+1] -= forceNlX;
    fy[iAtom+1] -= forceNlY;
    fz[iAtom+1] -= forceNlZ;
  }//endfor iAtom

/*======================================================================*/
/* III) free local memory                                               */

  free(wfNbhd);

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

