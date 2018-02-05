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

  double *forceNlX = pseudoReal->forceNlX;
  double *forceNlY = pseudoReal->forceNlY;
  double *forceNlZ = pseudoReal->forceNlZ;

  double *vnlPhiAtomGridRe = pseudoReal->vnlPhiAtomGridRe;
  double *vnlPhiAtomGridIm = pseudoReal->vnlPhiAtomGridIm;

  
  
  gridShiftNowRe = 0;
  gridShiftNowIm = 0;
  for(iAtom=0;iAtom<numAtom;iAtom++){
    forceNlX[iAtom] = 0.0;
    forceNlY[iAtom] = 0.0;
    forceNlZ[iAtom] = 0.0;
  }

  for(iAtom=0;iAtom<numAtom;iAtom++){ //
    atomType = iAtomAtomType[iAtom+1]-1;
    numGrid = numGridNlppMap[atomType];
    countRad = 0;
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
              dotRe = ddotBlasWrapper(numGrid,wfNbhd,1,
                      &vnlPhiAtomGridRe[gridShiftNowRe],1);
              dotIm = ddotBlasWrapper(numGrid,wfNbhd,1,
                      &vnlPhiAtomGridIm[gridShiftNowIm],1);
              energyl += 2.0*(dotRe*dotRe+dotIm*dotIm)*vpsNormList[radIndex]*volElem*volElem*volInv;
              dotRe *= 2.0*vpsNormList[radIndex]*volElem;
              dotIm *= 2.0*vpsNormList[radIndex]*volElem;
              daxpyBlasWrapper(numGrid,dotRe,&vnlPhiAtomGridRe[gridShiftNowRe],1,
                               forceTemp,1);
              daxpyBlasWrapper(numGrid,dotIm,&vnlPhiAtomGridIm[gridShiftNowIm],1,
                               forceTemp,1);
              gridShiftNowRe += numGrid;
              gridShiftNowIm += numGrid;
            }
            else{
              dotRe = ddotBlasWrapper(numGrid,wfNbhd,1,
                      &vnlPhiAtomGridRe[gridShiftNowRe],1);

              energyl += dotRe*dotRe*vpsNormList[radIndex]*volElem*volElem*volInv;
              dotRe *= vpsNormList[radIndex]*volElem;
              daxpyBlasWrapper(numGrid,dotRe,&vnlPhiAtomGridRe[gridShiftNowRe],1,
                               forceTemp,1);
              gridShiftNowRe += numGrid;
            }//endif m
          }//endfor m
          energy += energyl;
        }//endfor iRad
        countRad += atomLRadNum[atomType][l];
      }//endfor l
      for(iGrid=0;iGrid<numGrid;iGrid++){
        gridIndex = gridNlppMap[iAtom][iGrid];
        forceRealNlpp[gridIndex] += forceTemp[iGrid];
      }//endfor iGrid
    }//endif numGrid
  }//endfor iAtom

    }

    

  }//endfor iAtom
  


/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/



