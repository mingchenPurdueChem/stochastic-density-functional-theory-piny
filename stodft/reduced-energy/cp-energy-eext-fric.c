/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: cp-energy-eext-stodft.c                        */
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
#define TIME_CP_OFF
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcNlppRealFriction(CLASS *class,GENERAL_DATA *general_data,CP *cp,
                          double *hDevMat)
/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*************************************************************************/
/* Calculate nuclei force                                                */
/*************************************************************************/
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  CPOPTS *cpopts                = &(cp->cpopts);
  CPSCR *cpscr                  = &(cp->cpscr);
  PSEUDO *pseudo                = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal       = &(pseudo->pseudoReal);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  COMMUNICATE *communicate      = &(cp->communicate);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  CPEWALD *cpewald              = &(cp->cpewald);
  CELL *cell                    = &(general_data->cell);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  EWALD *ewald                  = &(general_data->ewald);
  PTENS *ptens                  = &(general_data->ptens);
  FRAGINFO *fragInfo            = stodftInfo->fragInfo;
  METALLIC *metallic            = stodftInfo->metallic;
  ATOMMAPS *atommaps            = &(class->atommaps);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm = &(cp->cp_sclr_fft_pkg3d_sm);

  int pseudoRealFlag = pseudoReal->pseudoRealFlag;
  int cpLsda         = cpopts->cp_lsda;
  int numStateStoUp  = stodftInfo->numStateStoUp;
  int chemPotOpt     = stodftInfo->chemPotOpt;
  int numStateUpProc = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc = cpcoeffs_info->nstate_dn_proc;
  int numCoeff       = cpcoeffs_info->ncoef;
  int numCoeffUpTotal = numStateUpProc*numCoeff;
  int numCoeffDnTotal = numStateDnProc*numCoeff;
  int numCoeffLarge       = cpcoeffs_info->ncoef_l;
  int numCoeffLargeProc;
  //int numCoeffLargeProc   = cp->cp_para_fft_pkg3d_lg.ncoef_proc;
  int nfft = cp_sclr_fft_pkg3d_sm->nfft;
  int numGrid = nfft/2;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int rhoRealGridNum    = stodftInfo->rhoRealGridNum;
  int numInterpPmeDual = pseudo->n_interp_pme_dual;
  int iperd           = cell->iperd;
  int numChemPot = stodftInfo->numChemPot;
  int occNumber = stodftInfo->occNumber;
  int myidState         = communicate->myid_state;
  int numProcStates = communicate->np_states;
  int numAtomTot = clatoms_info->natm_tot;
  int iState,jState,iCoeff,iChem,iAtom,iGrid;
  int smearOpt = stodftInfo->smearOpt;
  int numStateFric = metallic->numStateFric;
  int numAtomFricProc = metallic->numAtomFricProc;
  int numAtomType   = atommaps->natm_typ;

  int atomType,atomIndex;
  int numNlppTemp;

  int is,i,iupper;
  int ioff,ncoef1,ioff2;
  int iii,iis,nis;
  int iAng,iDim,countNlppRe,countNlppIm;
  int iRad,radIndex,countRad;
  int l,m;
  int myid_state = communicate->myid_state;
  int np_states  = communicate->np_states;

  int *atomFricIndProc = metallic->atomFricIndProc;
  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *numLMax = pseudoReal->numLMax;
  int *numNlppAtom = pseudoReal->numNlppAtom;
  int **atomLRadNum = pseudoReal->atomLRadNum;
  int **atomRadMap = pseudoReal->atomRadMap;

  double tpi = 2.0*M_PI;
  double eke,ekeDn;
  double chemPotTrue = stodftInfo->chemPotTrue;
  double energyKe   = stat_avg->cp_eke;
  double energyPnl  = stat_avg->cp_enl;
  double energyHart = stat_avg->cp_ehart;
  double energyEext = stat_avg->cp_eext;
  double energyExc  = stat_avg->cp_exc;
  double vol        = cell->vol;
  double volInv     = 1.0/vol;
  double energyTotElec,energyTot;
  double energyExtTemp,energyExcTemp,energyHartTemp;
  double vInter;
  double vrecip;
  double vself,vbgr;
  double vrecipLocal;
  double entropy = stodftInfo->entropy;
  double smearTemperature = stodftInfo->smearTemperature;

  double *ksStateChemPotRe = metallic->ksStateChemPotRe;
  double *ksStateChemPotIm = metallic->ksStateChemPotIm;
  double *zfft             = cpscr->cpscr_wave.zfft;
  double *wfReal;
  double *vpsNormList = pseudoReal->vpsNormList;
  double **dotRe,**dotIm,**dotDevRe,**dotDevIm;

  /*
  for(iAtomType=0;iAtomType<numAtomType;iType++){
    numNlppAtom[iAtomType] = 0.0;
    for(iAng=0;iAng<=numLMax[iAtomType];iAng++){
      numNlppAtom[iAtomType] += atomLRadNum[iAtomType][iAng]*(iAng+1);
    }//endfor iAng
  }
  */
  dotRe = (double**)cmalloc(numAtomFricProc*sizeof(double*));
  dotIm = (double**)cmalloc(numAtomFricProc*sizeof(double*));
  dotDevRe = (double**)cmalloc(numAtomFricProc*sizeof(double*));
  dotDevIm = (double**)cmalloc(numAtomFricProc*sizeof(double*));

  for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
    atomIndex = atomFricIndProc[iAtom];
    atomType = iAtomAtomType[atomIndex+1]-1;
    numNlppTemp = numNlppAtom[atomType];
    if(numNlppTemp>=1){   
      dotRe[iAtom] = (double*)cmalloc(numStateFric*numNlppTemp*sizeof(double));
      dotIm[iAtom] = (double*)cmalloc(numStateFric*(numNlppTemp-1)*sizeof(double));
      dotDevRe[iAtom] = (double*)cmalloc(numStateFric*numNlppTemp*3*sizeof(double));
      dotDevIm[iAtom] = (double*)cmalloc(numStateFric*(numNlppTemp-1)*3*sizeof(double));
    }
  }

  wfReal = (double*)cmalloc(numGrid*sizeof(double));

/*======================================================================*/
/* I) Loop all states involved in friction calculation                  */

  iupper = numStateFric;
  if(numStateFric%2!= 0){
    iupper = numStateFric-1;
  }// endif

  for(is=1;is<=iupper;is=is+2){
    ioff   = (is-1)*numCoeff;
    ioff2 = (is)*numCoeff;

/*======================================================================*/
/* II) FFT states back to real space                                    */

/*--------------------------------------------------------------------------*/
/*I) pack the complex zfft array with two wavefunctions (real) 
    the wavefunctions are reperesented in spherically cuttof 
    half g space                                                            */

    dble_pack_coef_fftw3d(&ksStateChemPotRe[ioff],&ksStateChemPotIm[ioff],
                          &ksStateChemPotRe[ioff2],&ksStateChemPotIm[ioff2],
                          zfft,cp_sclr_fft_pkg3d_sm);
     
/*--------------------------------------------------------------------------*/
/*II) fourier transform the two wavefunctions to real space
     convention exp(-igr)                                                   */

    para_fft_gen3d_fwd_to_r_fftw3d_threads(zfft,cp_sclr_fft_pkg3d_sm);

/*======================================================================*/
/* III) Calculate dotRe, dotIm, dotDevRe, dotDevIm                      */

    for(iGrid=0;iGrid<numGrid;iGrid++){
      wfReal[iGrid] = zfft[2*iGrid+1];
    }

    for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
      atomIndex = atomFricIndProc[iAtom];
      atomType = iAtomAtomType[atomIndex+1]-1;
      numNlppTemp = numNlppAtom[atomType];
      calcNlppDot(class,general_data,cp,atomIndex,iAtom,wfReal,&dotRe[iAtom][(is-1)*numNlppTemp],
                  &dotIm[iAtom][(is-1)*(numNlppTemp-1)],&dotDevRe[iAtom][(is-1)*numNlppTemp*3],
                  &dotDevIm[iAtom][(is-1)*(numNlppTemp-1)*3]);
    }

    for(iGrid=0;iGrid<numGrid;iGrid++){
      wfReal[iGrid] = zfft[2*iGrid+2];
    }
    for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
      atomIndex = atomFricIndProc[iAtom];
      atomType = iAtomAtomType[atomIndex+1]-1;
      numNlppTemp = numNlppAtom[atomType];
      calcNlppDot(class,general_data,cp,atomIndex,iAtom,wfReal,&dotRe[iAtom][(is-1)*numNlppTemp],
                  &dotIm[iAtom][(is-1)*(numNlppTemp-1)],&dotDevRe[iAtom][(is-1)*numNlppTemp*3],
                  &dotDevIm[iAtom][(is-1)*(numNlppTemp-1)*3]);
    }
  }//endfor is

  if(numStateFric%2!= 0){
    is = numStateFric;
    ioff = (is -1)*numCoeff;

/*--------------------------------------------------------------------------*/
/*   I) sngl pack                                                           */

    sngl_pack_coef_fftw3d(&ksStateChemPotRe[ioff],&ksStateChemPotIm[ioff],zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* II) fourier transform the wavefunctions to real space                    */
/*      convention exp(-igr)                                                */

    para_fft_gen3d_fwd_to_r_fftw3d_threads(zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* III) Calculate dotRe, dotIm, dotDevRe, dotDevIm                          */

    for(iGrid=0;iGrid<numGrid;iGrid++){
      wfReal[iGrid] = zfft[2*iGrid+1];
    }
    for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
      atomIndex = atomFricIndProc[iAtom];
      atomType = iAtomAtomType[atomIndex+1]-1;
      numNlppTemp = numNlppAtom[atomType];
      calcNlppDot(class,general_data,cp,atomIndex,iAtom,wfReal,&dotRe[iAtom][(is-1)*numNlppTemp],
                  &dotIm[iAtom][(is-1)*(numNlppTemp-1)],&dotDevRe[iAtom][(is-1)*numNlppTemp*3],
                  &dotDevIm[iAtom][(is-1)*(numNlppTemp-1)*3]);
    }
  }


/*======================================================================*/
/* V) Calculate vnlDevMat                                               */

  for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
    atomIndex = atomFricIndProc[iAtom];
    if(numNlppAtom[atomType]>=1){
      numNlppTemp = numNlppAtom[atomType];
      for(iDim=0;iDim<3;iDim++){
        for(iState=0;iState<numAtomFricProc;iState++){
          for(jState=0;jState<numAtomFricProc;jState++){
            countNlppRe = 0;
            countNlppIm = 0;
            for(l=0;l<=numLMax[atomType];l++){
              for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
                radIndex = atomRadMap[atomType][countRad+iRad];
                for(m=0;m<=l;m++){
                  //calcDotNlpp(wfNbhd,radFun,&ylm[ylmShift],&dotRe,&dotIm);
                  if(m!=0){
                    hDevMat[(iAtom*3+iDim)*numStateFric*numStateFric+iState*numStateFric+jState] 
                    += (dotRe[iAtom][iState*numNlppTemp+countNlppRe+m]*
                        dotDevRe[iAtom][iState*numNlppTemp*3+(countNlppRe+m)*3+iDim]+
                        dotIm[iAtom][iState*numNlppTemp+countNlppIm+m-1]*
                        dotDevIm[iAtom][iState*numNlppTemp*3+(countNlppIm+m-1)*3+iDim])*
                        4.0*vpsNormList[radIndex]*volInv;
                    
                  }
                  else{
                    hDevMat[(iAtom*3+iDim)*numStateFric*numStateFric+iState*numStateFric+jState]
                    += (dotRe[iAtom][iState*numNlppTemp+countNlppRe+m]*
                        dotDevRe[iAtom][iState*numNlppTemp*3+(countNlppRe+m)*3+iDim])*
                        2.0*vpsNormList[radIndex]*volInv;
                  }//endif m
                }//endfor m
                countNlppRe += l+1;
                countNlppIm += l;
              }//endfor iRad
            }//endfor l
          }//endfor jState
        }//endfor iState
      }//endfor iDim
    }//endifnumNlppAtom
  }//endfor iAtom

/*======================================================================*/
/* V) Free template array                                               */

  free(wfReal);
  for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
    if(numNlppTemp>=1){
      free(dotRe[iAtom]);
      free(dotIm[iAtom]);
      free(dotDevRe[iAtom]);
      free(dotDevIm[iAtom]);
    }
  }
  free(dotRe);
  free(dotIm);
  free(dotDevRe);
  free(dotDevIm);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcNlppDot(CLASS *class,GENERAL_DATA *general_data,CP *cp,
                      int atomIndex,int iAtom, double *wfReal, 
                      double *dotRe,double *dotIm,double *dotDevRe,double *dotDevIm)
/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*************************************************************************/
/* Calculate nuclei force                                                */
/*************************************************************************/
/*==========================================================================*/
/*               Local variable declarations                                */
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  PSEUDO *pseudo = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CELL *cell = &(general_data->cell);
  CLATOMS_POS *clatoms_pos = &(class->clatoms_pos[1]);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  ATOMMAPS *atommaps = &(class->atommaps);
  STAT_AVG *stat_avg = &(general_data->stat_avg);
  CPEWALD *cpewald = &(cp->cpewald);
  METALLIC *metallic            = stodftInfo->metallic;
  PARA_FFT_PKG3D *cpParaFftPkg3d = &(cp->cp_para_fft_pkg3d_lg);
  COMMUNICATE *communicate      = &(cp->communicate);

  int realSparseOpt = cpewald->realSparseOpt;
  int numAtomType = atommaps->natm_typ;
  int numAtom = clatoms_info->natm_tot;
  int numGrid;
  int numGridMax;
  int countRad,countNlppRe,countNlppIm;
  int iPart,iRad,l,m,iType,iGrid;
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
  int *numLMax = pseudoReal->numLMax; //max L for each atom
  int **atomLRadNum = pseudoReal->atomLRadNum; //num of radical functions for each atom and each L
  int **atomRadMap = pseudoReal->atomRadMap;  //map of radical functions for each atom, starting from l=0 to l=lmax
  int *numGridNlppMap = pseudoReal->numGridNlppMap;
  int *locOpt = pseudo->loc_opt;
  int **gridNlppMap = pseudoReal->gridNlppMap;

  double vpsNorm;
  double vol,volInv;
  double volElem;
  double energy = 0.0;
  double energyl;

  //double *dotRe,*dotIm;
  //double *dotDevRe,*dotDevIm;
  double *vnlDevMat = metallic->vnlDevMat;
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
  double *vnlPhiDxAtomGridRe = pseudoReal->vnlPhiDxAtomGridRe;
  double *vnlPhiDxAtomGridIm = pseudoReal->vnlPhiDxAtomGridIm;
  double *vnlPhiDyAtomGridRe = pseudoReal->vnlPhiDyAtomGridRe;
  double *vnlPhiDyAtomGridIm = pseudoReal->vnlPhiDyAtomGridIm;
  double *vnlPhiDzAtomGridRe = pseudoReal->vnlPhiDzAtomGridRe;
  double *vnlPhiDzAtomGridIm = pseudoReal->vnlPhiDzAtomGridIm;
  
  atomType = iAtomAtomType[atomIndex+1]-1;

  numGrid = numGridNlppMap[atomIndex];
  countRad = 0;
  countNlppRe = 0;
  countNlppIm = 0;
  gridShiftNowRe = gridStIndRe[atomIndex];
  gridShiftNowIm = gridStIndIm[atomIndex];
  /* cpy the wave function */
  //printf("numGrid %i locOpt[atomType+1] %i\n",numGrid,locOpt[atomType+1]);
  if(numGrid>0){ //if numGrid=0, only local pp will be calculated
    for(iGrid=0;iGrid<numGrid;iGrid++){
      gridIndex = gridNlppMap[atomIndex][iGrid];
      wfNbhd[iGrid] = wfReal[gridIndex];
    }
    for(l=0;l<=numLMax[atomType];l++){
      if(locOpt[atomType+1]!=l){
        for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
          radIndex = atomRadMap[atomType][countRad+iRad];
          energyl = 0.0;
          for(m=0;m<=l;m++){
            if(m!=0){
              dotRe[countNlppRe+m] = ddotBlasWrapper(numGrid,&wfNbhd[0],1,
                                        &vnlPhiAtomGridRe[gridShiftNowRe],1)*volElem;
              dotIm[countNlppIm+m-1] = ddotBlasWrapper(numGrid,&wfNbhd[0],1,
                                         &vnlPhiAtomGridIm[gridShiftNowIm],1)*volElem;
              dotDevRe[(countNlppRe+m)*3] = ddotBlasWrapper(numGrid,wfNbhd,1,
                                                &vnlPhiDxAtomGridRe[gridShiftNowRe],1)*volElem;
              dotDevIm[(countNlppIm+m-1)*3] = ddotBlasWrapper(numGrid,wfNbhd,1,
                                                &vnlPhiDxAtomGridIm[gridShiftNowIm],1)*volElem;
              dotDevRe[(countNlppRe+m)*3+1] = ddotBlasWrapper(numGrid,wfNbhd,1,
                                                &vnlPhiDyAtomGridRe[gridShiftNowRe],1)*volElem;
              dotDevIm[(countNlppIm+m-1)*3+1] = ddotBlasWrapper(numGrid,wfNbhd,1,
                                                &vnlPhiDyAtomGridIm[gridShiftNowIm],1)*volElem;
              dotDevRe[(countNlppRe+m)*3+2] = ddotBlasWrapper(numGrid,wfNbhd,1,
                                                &vnlPhiDzAtomGridRe[gridShiftNowRe],1)*volElem;
              dotDevIm[(countNlppIm+m-1)*3+2] = ddotBlasWrapper(numGrid,wfNbhd,1,
                                                &vnlPhiDzAtomGridIm[gridShiftNowIm],1)*volElem;
              gridShiftNowRe += numGrid;
              gridShiftNowIm += numGrid;
            }
            else{
              //dotRe = ddot1(numGrid,&wfNbhd[0],1,
              //          &vnlPhiAtomGridRe[gridShiftNowRe],1)*volElem;
              dotRe[countNlppRe] = ddotBlasWrapper(numGrid,&wfNbhd[0],1,
                                      &vnlPhiAtomGridRe[gridShiftNowRe],1)*volElem;
              dotDevRe[countNlppRe*3] = ddotBlasWrapper(numGrid,wfNbhd,1,
                                           &vnlPhiDxAtomGridRe[gridShiftNowRe],1)*volElem;
              dotDevRe[countNlppRe*3+1] = ddotBlasWrapper(numGrid,wfNbhd,1,
                                           &vnlPhiDyAtomGridRe[gridShiftNowRe],1)*volElem;
              dotDevRe[countNlppRe*3+2] = ddotBlasWrapper(numGrid,wfNbhd,1,
                                           &vnlPhiDzAtomGridRe[gridShiftNowRe],1)*volElem;
              //printf("volElem %lg\n",volElem);
              //printf("11111111 dotRe %.8lg\n",dotRe);
              gridShiftNowRe += numGrid;
            }//endif m
          }//endfor m
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

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcLocalPotFriction(CLASS *class,GENERAL_DATA *general_data,CP *cp,
                          double *vlocDevMat)
/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*************************************************************************/
/* Calculate nuclei force                                                */
/*************************************************************************/
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  EWALD        *ewald        = &(general_data->ewald);
  CELL         *cell         = &(general_data->cell);
  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR        *cpscr        = &(cp->cpscr);
  CPEWALD      *cpewald      = &(cp->cpewald);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  COMMUNICATE   *communicate    = &(cp->communicate);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPCOEFFS_POS  *cpcoeffs_pos   = &(cp->cpcoeffs_pos[1]);
  PSEUDO        *pseudo         = &(cp->pseudo);
  FRAGINFO *fragInfo            = stodftInfo->fragInfo;
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  METALLIC *metallic            = stodftInfo->metallic;
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  CLATOMS_POS *clatoms_pos      = &(class->clatoms_pos[1]);

  int smearOpt        = stodftInfo->smearOpt;
  int numStateFric    = metallic->numStateFric;
  int numAtomFricProc = metallic->numAtomFricProc;
  int numCoeff        = cpcoeffs_info->ncoef;  
  int atomType,atomIndex;
  int rhoRealGridNum    = stodftInfo->rhoRealGridNum;
  int rhoRealGridTot    = stodftInfo->rhoRealGridTot;
  int cpLsda = cpopts->cp_lsda;
  int cpGGA  = cpopts->cp_gga;
  int cpDualGridOptOn = cpopts->cp_dual_grid_opt;
  int myidState         = communicate->myid_state;
  int numProcStates     = communicate->np_states;
  int numInterpPmeDual = pseudo->n_interp_pme_dual;
  int numAtomTot = clatoms_info->natm_tot;
  
  int iState,jState,iGrid,iAtom;

  MPI_Comm comm_states   =    communicate->comm_states;


  int *coefFormUp   = &(cpcoeffs_pos->icoef_form_up);
  int *coefFormDn   = &(cpcoeffs_pos->icoef_form_dn);
  int *coefOrthUp   = &(cpcoeffs_pos->icoef_orth_up);
  int *coefOrthDn   = &(cpcoeffs_pos->icoef_orth_dn);
  int *densityMap   = stodftInfo->densityMap;
  int *rhoRealSendCounts = stodftInfo->rhoRealSendCounts;
  int *rhoRealDispls = stodftInfo->rhoRealDispls;
  int *atomFricIndProc = metallic->atomFricIndProc;

  double *ksStateChemPotRe = metallic->ksStateChemPotRe;
  double *ksStateChemPotIm = metallic->ksStateChemPotIm;
  double *wfReal = (double*)cmalloc(numStateFric*rhoRealGridTot*sizeof(double));
  double *rhoTemp  = (double*)cmalloc(rhoRealGridTot*sizeof(double));
  double *rhoCoeffReUp   = cpscr->cpscr_rho.rhocr_up;
  double *rhoCoeffImUp   = cpscr->cpscr_rho.rhoci_up;
  double *rhoUp          = cpscr->cpscr_rho.rho_up;
  double *rhoCoeffReUpDensCpBox = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
  double *rhoCoeffImUpDensCpBox = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
  double *rhoUpCorrect    = stodftCoefPos->rhoUpCorrect;
  double *divRhoxUp       = cpscr->cpscr_grho.d_rhox_up;
  double *divRhoyUp       = cpscr->cpscr_grho.d_rhoy_up;
  double *divRhozUp       = cpscr->cpscr_grho.d_rhoz_up;
  double *d2RhoUp        = cpscr->cpscr_grho.d2_rho_up;
  double *rhoCoeffReDn   = cpscr->cpscr_rho.rhocr_dn;
  double *rhoCoeffImDn   = cpscr->cpscr_rho.rhoci_dn;
  double *rhoDn          = cpscr->cpscr_rho.rho_dn;
  double *rhoCoeffReDnDensCpBox = cpscr->cpscr_rho.rhocr_dn_dens_cp_box;
  double *rhoCoeffImDnDensCpBox = cpscr->cpscr_rho.rhoci_dn_dens_cp_box;
  double *rhoDnCorrect    = stodftCoefPos->rhoDnCorrect;
  double *divRhoxDn       = cpscr->cpscr_grho.d_rhox_dn;
  double *divRhoyDn       = cpscr->cpscr_grho.d_rhoy_dn;
  double *divRhozDn       = cpscr->cpscr_grho.d_rhoz_dn;
  double *d2RhoDn        = cpscr->cpscr_grho.d2_rho_dn;

  double *fx = clatoms_pos->fx;
  double *fy = clatoms_pos->fy;
  double *fz = clatoms_pos->fz;

  for(iState=0;iState<numStateFric;iState++){
    rhoCalcRealFriction(general_data,cp,class,ksStateChemPotRe,ksStateChemPotIm,
                        wfReal,numStateFric);
  }

	  
  for(iState=0;iState<numStateFric;iState++){
    for(jState=0;jState<numStateFric;jState++){
      /* Calculate Real Space "Density" phi_i^*(r)phi_j(r) */
      if(myidState==0){
        for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
          rhoTemp[iGrid] = wfReal[iState*rhoRealGridTot+iGrid]*wfReal[jState*rhoRealGridTot+iGrid];
        }
      }

      /* Scatter rhoTemp */
      if(numProcStates>1){
        Scatterv(rhoTemp,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
                 &rhoUp[1],rhoRealGridNum,MPI_DOUBLE,0,comm_states);
      }
      else{
        memcpy(&rhoUp[1],rhoTemp,rhoRealGridNum*sizeof(double));
      }
      /* Calculate k-space "Density" */
      calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
                         rhoCoeffReUp,rhoCoeffImUp,rhoUp,rhoCoeffReUpDensCpBox,rhoCoeffImUpDensCpBox,
                         divRhoxUp,divRhoyUp,divRhozUp,d2RhoUp,cpGGA,cpDualGridOptOn,numInterpPmeDual,
                         communicate,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));

      /* Calculate F_local(ij) */
      for(iAtom=1;iAtom<=numAtomTot;iAtom++){
        fx[iAtom] = 0.0;
        fy[iAtom] = 0.0;
        fz[iAtom] = 0.0;
      }
         
      calcLocExtPostScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
      /* Broadcast Force */
      Bcast(&fx[1],numAtomTot,MPI_DOUBLE,0,comm_states);
      Bcast(&fy[1],numAtomTot,MPI_DOUBLE,0,comm_states);
      Bcast(&fz[1],numAtomTot,MPI_DOUBLE,0,comm_states);
      for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
        atomIndex = atomFricIndProc[iAtom];
        vlocDevMat[(iAtom*3)*numStateFric*numStateFric+iState*numStateFric+jState] = fx[atomIndex];
        vlocDevMat[(iAtom*3+1)*numStateFric*numStateFric+iState*numStateFric+jState] = fy[atomIndex];
        vlocDevMat[(iAtom*3+2)*numStateFric*numStateFric+iState*numStateFric+jState] = fz[atomIndex];
      }
    }//endfor jState
  }//endfor iState

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void rhoCalcRealFriction(GENERAL_DATA *general_data,CP *cp,CLASS *class,
                        double *ccreal,double *ccimag,
                        double *wfReal,int nstate)
/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"


  CPEWALD *cpewald = &(cp->cpewald);
  CPSCR *cpscr = &(cp->cpscr);
  CPOPTS *cpopts = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cp->cpcoeffs_info);
  COMMUNICATE *communicate = &(cp->communicate);
  EWALD *ewald = &(general_data->ewald);
  CELL *cell = &(general_data->cell);

  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm = &(cp->cp_sclr_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

  int iii,ioff,ioff2;
  int is,i,j,k,iupper;
  int gridoff,gridoff2;
  int igrid;
  int ncoef = cpcoeffs_info->ncoef;
  int   nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
  int   nfft2_proc       =    nfft_proc/2;
  int nfft = cp_para_fft_pkg3d_lg->nfft;
  int nfft2 = nfft/2;

  double *zfft           =    cpscr->cpscr_wave.zfft;
  double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
  double *hmatCP         =    cell->hmat_cp;
  double volCP           = getdeth(hmatCP);
  double invVolCP        = 1.0/invVolCP;

/*==========================================================================*/
/*2)sum the density in real space two states at a time                      */
/*  This is done at state level and uses the scalar packages!               */

  iupper = nstate;
  if(nstate%2!=0){
     iupper = nstate-1;
  }/* endif */

  for(is=1;is<=iupper;is=is+2){
    ioff   = (is-1)*ncoef;
    ioff2 = (is)*ncoef;
    gridoff = (is-1)*nfft2;
    gridoff2 = is*nfft2;

/*--------------------------------------------------------------------------*/
/*I) pack the complex zfft array with two wavefunctions (real) 
    the wavefunctions are reperesented in spherically cuttof 
    half g space                                                            */

    //printf("creal1 %lg creal2 %lg\n",ccreal[ioff+1],ccreal[ioff2+1]);
    dble_pack_coef(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
                      zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*II) fourier transform the two wavefunctions to real space
     convention exp(-igr)                                                   */

    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
    //printf("zfft %lg\n",zfft[1]);

/*--------------------------------------------------------------------------*/
/* III) Copy the real sapce wave function and add the square of the two     
        wave functions to the density(real space)                           */

    for(igrid=0;igrid<nfft2;igrid++){
      wfReal[gridoff+igrid] = zfft[igrid*2+1];
      wfReal[gridoff2+igrid] = zfft[igrid*2+2];
    }
    //printf("wfReal %lg %lg\n",wfReal[(is-1)*nfft2_proc],wfReal[is*nfft2_proc]);
  }/*endfor is*/

/*--------------------------------------------------------------------------*/
/*IV) if there is an odd number of states, use single pack
      for the last state (i caen't be chaengin' the laws of physics
      captn, i've got to be usin the single pack!!!!!!)                     */

  if(nstate%2!=0){
    ioff = (nstate-1)*ncoef;
    gridoff = (nstate-1)*nfft2;
    sngl_pack_coef(&ccreal[ioff],&ccimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* V) fourier transform the last wavefunction to real space
       convention exp(-igr)                                                 */

    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);


/*--------------------------------------------------------------------------*/
/*VI) Copy the real sapce wave function and add the square of the last wave 
      function to the density(real space)   */

    for(igrid=0;igrid<nfft2;igrid++){
      wfReal[gridoff+igrid] = zfft[igrid*2+1];
    }

  }//endif nstat%2
  /*
  printf("nfft2 %i\n",nfft2);
  for(is=0;is<nstate;is++){
    printf("wwwfReal %lg\n",wfReal[is*nfft2_proc]);
  }
  */

  /*
  for(is=0;is<nstate;is++){
    char fname[100]; 
    sprintf(fname,"wf-stodft-%i",is);
    FILE *fout = fopen(fname,"w");
    for(igrid=0;igrid<nfft2;igrid++)fprintf(fout,"%i %.16lg\n",igrid,wfReal[is*nfft2+igrid]);
    fclose(fout);
  }
  */

/*==============================================================*/
}/*end routine*/
/*==============================================================*/

