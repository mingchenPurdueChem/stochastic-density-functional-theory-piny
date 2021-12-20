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
  int *locOpt = pseudo->loc_opt;
  int **atomLRadNum = pseudoReal->atomLRadNum;
  int **atomRadMap = pseudoReal->atomRadMap;

  double tpi = 2.0*M_PI;
  double vol        = cell->vol;
  double volInv     = 1.0/vol;
  double temp;

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
      calcNlppDot(class,general_data,cp,atomIndex,iAtom,wfReal,&dotRe[iAtom][(is)*numNlppTemp],
                  &dotIm[iAtom][(is)*(numNlppTemp-1)],&dotDevRe[iAtom][(is)*numNlppTemp*3],
                  &dotDevIm[iAtom][(is)*(numNlppTemp-1)*3]);
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

  for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
    for(is=0;is<numStateFric;is++){
      printf("dotRe %lg dotIm %lg dotDevRe %lg dotDevIm %lg\n",dotRe[iAtom][(is)*numNlppTemp],
             dotIm[iAtom][(is)*(numNlppTemp-1)],dotDevRe[iAtom][(is)*numNlppTemp*3],
             dotDevIm[iAtom][(is)*(numNlppTemp-1)*3]);
    }
  }

/*======================================================================*/
/* V) Calculate vnlDevMat                                               */

  for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
    atomIndex = atomFricIndProc[iAtom];
    atomType = iAtomAtomType[atomIndex+1]-1;
    if(numNlppAtom[atomType]>=1){
      numNlppTemp = numNlppAtom[atomType];
      for(iDim=0;iDim<3;iDim++){
        for(iState=0;iState<numStateFric;iState++){
          for(jState=0;jState<numStateFric;jState++){
            countNlppRe = 0;
            countNlppIm = 0;
            countRad = 0;
            for(l=0;l<=numLMax[atomType];l++){
              if(locOpt[atomType+1]!=l){
                for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
                  radIndex = atomRadMap[atomType][countRad+iRad];
                  for(m=0;m<=l;m++){
                    //calcDotNlpp(wfNbhd,radFun,&ylm[ylmShift],&dotRe,&dotIm);
                    if(m!=0){
                      // 2.0 due to the contribution of m and -m
                      temp = 
                          (dotRe[iAtom][iState*numNlppTemp+countNlppRe+m]*
                          dotDevRe[iAtom][jState*numNlppTemp*3+(countNlppRe+m)*3+iDim]+
                          dotIm[iAtom][iState*(numNlppTemp-1)+countNlppIm+m-1]*
                          dotDevIm[iAtom][jState*(numNlppTemp-1)*3+(countNlppIm+m-1)*3+iDim]+
                          dotRe[iAtom][jState*numNlppTemp+countNlppRe+m]*
                          dotDevRe[iAtom][iState*numNlppTemp*3+(countNlppRe+m)*3+iDim]+
                          dotIm[iAtom][jState*(numNlppTemp-1)+countNlppIm+m-1]*
                          dotDevIm[iAtom][iState*(numNlppTemp-1)*3+(countNlppIm+m-1)*3+iDim])*
                          2.0*vpsNormList[radIndex]*volInv;
                      hDevMat[(iAtom*3+iDim)*numStateFric*numStateFric+iState*numStateFric+jState] += temp;
                      printf("iDim %i iState %i jState %i l %i m %i dotRe %lg dotIm %lg dotDevRe %lg dotDevIm %lg vpsNormList %lg volInv %lg temp %lg\n",
                             iDim,iState,jState,l,m,dotRe[iAtom][iState*numNlppTemp+countNlppRe+m],
                             dotIm[iAtom][iState*(numNlppTemp-1)+countNlppIm+m-1],
                             dotDevRe[iAtom][jState*numNlppTemp*3+(countNlppRe+m)*3+iDim],
                             dotDevIm[iAtom][jState*(numNlppTemp-1)*3+(countNlppIm+m-1)*3+iDim],
                             vpsNormList[radIndex],volInv,temp);
                    }
                    else{
                      temp = 
                          (dotRe[iAtom][iState*numNlppTemp+countNlppRe+m]*
                          dotDevRe[iAtom][jState*numNlppTemp*3+(countNlppRe+m)*3+iDim]+
                          dotRe[iAtom][jState*numNlppTemp+countNlppRe+m]*
                          dotDevRe[iAtom][iState*numNlppTemp*3+(countNlppRe+m)*3+iDim])*
                          vpsNormList[radIndex]*volInv;
                      hDevMat[(iAtom*3+iDim)*numStateFric*numStateFric+iState*numStateFric+jState] += temp;
                      printf("iDim %i iState %i jState %i l %i m %i dotRe %lg dotIm 0.0 dotDevRe %lg dotDevIm 0.0 vpsNormList %lg volInv %lg temp %lg\n",
                             iDim,iState,jState,l,m,dotRe[iAtom][iState*numNlppTemp+countNlppRe+m],
                             dotDevRe[iAtom][jState*numNlppTemp*3+(countNlppRe+m)*3+iDim],
                             vpsNormList[radIndex],volInv,temp);
                    }//endif m
                  }//endfor m
                  countNlppRe += l+1;
                  countNlppIm += l;
                }//endfor iRad
                countRad += atomLRadNum[atomType][l];
              }
              else{ // local channel
                for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
                  radIndex = atomRadMap[atomType][countRad+iRad];
                  countNlppRe += l+1;
                  countNlppIm += l;
                }//endfor iRad
                countRad += atomLRadNum[atomType][l];
              }//endif locOpt
            }//endfor l
            //printf("hDevMat %lg\n",hDevMat[(iAtom*3+iDim)*numStateFric*numStateFric+iState*numStateFric+jState]);
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
  double vol = cell->vol;
  double volInv = 1.0/vol;
  double volElem = vol/numGridTot;
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
  wfNbhd = (double*)cmalloc(numGrid*sizeof(double));
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
              //printf("wfNbhd %lg vnlPhiAtomGridRe %lg vnlPhiDxAtomGridRe %lg volElem %lg\n",wfNbhd[0],vnlPhiAtomGridRe[gridShiftNowRe],vnlPhiDxAtomGridRe[gridShiftNowRe],volElem);
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
              printf("l %i m %i dotRe %lg dotIm %lg dotDevRe %lg %lg %lg dotDevIm %lg %lg %lg\n",l,m,dotRe[countNlppRe+m],dotIm[countNlppIm+m-1],dotDevRe[(countNlppRe+m)*3],dotDevRe[(countNlppRe+m)*3+1],dotDevRe[(countNlppRe+m)*3+2],dotDevIm[(countNlppIm+m-1)*3],dotDevIm[(countNlppIm+m-1)*3+1],dotDevIm[(countNlppIm+m-1)*3+2]);

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
              printf("l %i m %i dotRe %lg dotDevRe %lg %lg %lg\n",l,m,dotRe[countNlppRe],dotDevRe[(countNlppRe)*3],dotDevRe[(countNlppRe)*3+1],dotDevRe[(countNlppRe)*3+2]);
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

  free(wfNbhd);

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

  cpewald->ewaldLocalOpt = 1;	  
  for(iState=0;iState<numStateFric;iState++){
    for(jState=0;jState<numStateFric;jState++){
      /* Calculate Real Space "Density" phi_i^*(r)phi_j(r) */
      if(myidState==0){
        for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
          rhoTemp[iGrid] = wfReal[iState*rhoRealGridTot+iGrid]*wfReal[jState*rhoRealGridTot+iGrid];
        }
        printf("rhoTemp %lg\n",rhoTemp[0]);
      }

      /* Scatter rhoTemp */
      if(numProcStates>1){
        Scatterv(rhoTemp,rhoRealSendCounts,rhoRealDispls,MPI_DOUBLE,
                 &rhoUp[1],rhoRealGridNum,MPI_DOUBLE,0,comm_states);
        printf("rhoUp %lg\n",rhoUp[1]);
      }
      else{
        memcpy(&rhoUp[1],rhoTemp,rhoRealGridNum*sizeof(double));
      }
      /* Calculate k-space "Density" */
      calcRhoStoRecipFullg(cpewald,cpscr,cpcoeffs_info,ewald,cell,
                         rhoCoeffReUp,rhoCoeffImUp,rhoUp,rhoCoeffReUpDensCpBox,rhoCoeffImUpDensCpBox,
                         divRhoxUp,divRhoyUp,divRhozUp,d2RhoUp,cpGGA,cpDualGridOptOn,numInterpPmeDual,
                         communicate,&(cp->cp_para_fft_pkg3d_lg),&(cp->cp_para_fft_pkg3d_dens_cp_box));

      printf("rhoCoeffReUp %lg rhoCoeffImUp %lg\n",rhoCoeffReUp[1],rhoCoeffImUp[1]);
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
        printf("Fx %lg atomIndex %i\n",fx[1],atomIndex);
        vlocDevMat[(iAtom*3)*numStateFric*numStateFric+iState*numStateFric+jState] = fx[atomIndex+1];
        vlocDevMat[(iAtom*3+1)*numStateFric*numStateFric+iState*numStateFric+jState] = fy[atomIndex+1];
        vlocDevMat[(iAtom*3+2)*numStateFric*numStateFric+iState*numStateFric+jState] = fz[atomIndex+1];
      }
    }//endfor jState
  }//endfor iState
  cpewald->ewaldLocalOpt = 0;

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


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void genDHPhi(CP *cp,CLASS *class, GENERAL_DATA *general_data,
           double *coeffReUp, double *coeffImUp, 
           double *coeffReDn, double *coeffImDn,
           double *daHChiReUp, double *daHChiImUp, 
           double *daHChiReDn, double *daHChiImDn)
/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"

  CPOPTS       *cpopts       = &(cp->cpopts);
  CPSCR        *cpscr        = &(cp->cpscr);
  CPCOEFFS_POS  *cpcoeffs_pos   = &(cp->cpcoeffs_pos[1]);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  COMMUNICATE   *commCP         = &(cp->communicate);
  STODFTINFO   *stodftInfo   = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;
  
/*==========================================================================*/
/* 1) Non-local pp */

  calcDNLPhi(cp,class,general_data,coeffReUp,coeffImUp,
             coeffReDn,coeffImDn,daHChiReUp,daHChiImUp,
             daHChiReDn,daHChiImDn);

/*==========================================================================*/
/* 1) Local pp */

  calcDLocPhi(cp,class,general_data,coeffReUp,coeffImUp,
             coeffReDn,coeffImDn,daHChiReUp,daHChiImUp,
             daHChiReDn,daHChiImDn);

/*==============================================================*/
}/*end routine*/
/*==============================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcDNLPhi(CP *cp,CLASS *class, GENERAL_DATA *general_data,
           double *coeffReUp, double *coeffImUp,
           double *coeffReDn, double *coeffImDn,
           double *daHChiReUp, double *daHChiImUp,
           double *daHChiReDn, double *daHChiImDn)
/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"

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
  int numGridAll = nfft/2;
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
  int *locOpt = pseudo->loc_opt;
  int **atomLRadNum = pseudoReal->atomLRadNum;
  int **atomRadMap = pseudoReal->atomRadMap;

  double tpi = 2.0*M_PI;
  double vol        = cell->vol;
  double volInv     = 1.0/vol;
  double temp;

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
  dotRe = (double**)cmalloc(numAtomFric*sizeof(double*));
  dotIm = (double**)cmalloc(numAtomFric*sizeof(double*));
  dotDevRe = (double**)cmalloc(numAtomFric*sizeof(double*));
  dotDevIm = (double**)cmalloc(numAtomFric*sizeof(double*));

  for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
    atomIndex = atomFricIndProc[iAtom];
    atomType = iAtomAtomType[atomIndex+1]-1;
    numNlppTemp = numNlppAtom[atomType];
    if(numNlppTemp>=1){   
      dotRe[iAtom] = (double*)cmalloc(numStateUpProc*numNlppTemp*sizeof(double));
      dotIm[iAtom] = (double*)cmalloc(numStateUpProc*(numNlppTemp-1)*sizeof(double));
      dotDevRe[iAtom] = (double*)cmalloc(numStateUpProc*numNlppTemp*3*sizeof(double));
      dotDevIm[iAtom] = (double*)cmalloc(numStateUpProc*(numNlppTemp-1)*3*sizeof(double));
    }
  }

  wfReal = (double*)cmalloc(numGrid*sizeof(double));

/*======================================================================*/
/* I) Loop all states involved in friction calculation                  */

  iupper = numStateUpProc;
  if(numStateUpProc%2!= 0){
    iupper = numStateUpProc-1;
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

    dble_pack_coef_fftw3d(&coeffReUp[ioff],&coeffImUp[ioff],
                          &coeffReUp[ioff2],&coeffImUp[ioff2],
                          zfft,cp_sclr_fft_pkg3d_sm);
     
/*--------------------------------------------------------------------------*/
/*II) fourier transform the two wavefunctions to real space
     convention exp(-igr)                                                   */

    para_fft_gen3d_fwd_to_r_fftw3d_threads(zfft,cp_sclr_fft_pkg3d_sm);

/*======================================================================*/
/* III) Calculate dotRe, dotIm, dotDevRe, dotDevIm                      */

    for(iGrid=0;iGrid<numGridAll;iGrid++){
      wfReal[iGrid] = zfft[2*iGrid+1];
    }

    for(iAtom=0;iAtom<numAtomFric;iAtom++){
      atomIndex = atomFricIndProc[iAtom];
      atomType = iAtomAtomType[atomIndex+1]-1;
      numNlppTemp = numNlppAtom[atomType];
      calcNlppDot(class,general_data,cp,atomIndex,iAtom,wfReal,&dotRe[iAtom][(is-1)*numNlppTemp],
                  &dotIm[iAtom][(is-1)*(numNlppTemp-1)],&dotDevRe[iAtom][(is-1)*numNlppTemp*3],
                  &dotDevIm[iAtom][(is-1)*(numNlppTemp-1)*3]);
    }

    for(iGrid=0;iGrid<numGridAll;iGrid++){
      wfReal[iGrid] = zfft[2*iGrid+2];
    }
    for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
      atomIndex = atomFricIndProc[iAtom];
      atomType = iAtomAtomType[atomIndex+1]-1;
      numNlppTemp = numNlppAtom[atomType];
      calcNlppDot(class,general_data,cp,atomIndex,iAtom,wfReal,&dotRe[iAtom][(is)*numNlppTemp],
                  &dotIm[iAtom][(is)*(numNlppTemp-1)],&dotDevRe[iAtom][(is)*numNlppTemp*3],
                  &dotDevIm[iAtom][(is)*(numNlppTemp-1)*3]);
    }
  }//endfor is

  if(numStateFric%2!= 0){
    is = numStateFric;
    ioff = (is -1)*numCoeff;

/*--------------------------------------------------------------------------*/
/*   I) sngl pack                                                           */

    sngl_pack_coef_fftw3d(&coeffReUp[ioff],&coeffImUp[ioff],zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* II) fourier transform the wavefunctions to real space                    */
/*      convention exp(-igr)                                                */

    para_fft_gen3d_fwd_to_r_fftw3d_threads(zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* III) Calculate dotRe, dotIm, dotDevRe, dotDevIm                          */

    for(iGrid=0;iGrid<numGridAll;iGrid++){
      wfReal[iGrid] = zfft[2*iGrid+1];
    }
    for(iAtom=0;iAtom<numAtomFric;iAtom++){
      atomIndex = atomFricIndProc[iAtom];
      atomType = iAtomAtomType[atomIndex+1]-1;
      numNlppTemp = numNlppAtom[atomType];
      calcNlppDot(class,general_data,cp,atomIndex,iAtom,wfReal,&dotRe[iAtom][(is-1)*numNlppTemp],
                  &dotIm[iAtom][(is-1)*(numNlppTemp-1)],&dotDevRe[iAtom][(is-1)*numNlppTemp*3],
                  &dotDevIm[iAtom][(is-1)*(numNlppTemp-1)*3]);
    }
  }

  for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
    for(is=0;is<numStateFric;is++){
      printf("dotRe %lg dotIm %lg dotDevRe %lg dotDevIm %lg\n",dotRe[iAtom][(is)*numNlppTemp],
             dotIm[iAtom][(is)*(numNlppTemp-1)],dotDevRe[iAtom][(is)*numNlppTemp*3],
             dotDevIm[iAtom][(is)*(numNlppTemp-1)*3]);
    }
  }

/*======================================================================*/
/* V) Calculate D_nlpp|phi>                                             */

  double forceTemp = (double*)calloc(numStateUpProc*numAtomFric*3*numGridMax*sizeof(double));


  for(iState=0;iState<numStateUpProc;iState++){
    for(iAtom=0;iAtom<numAtomFric;iAtom++){
      atomIndex = atomFricIndProc[iAtom];
      atomType = iAtomAtomType[atomIndex+1]-1;
      numGrid = numGridNlppMap[atomIndex];
      countRad = 0;
      countNlppRe = 0;
      countNlppIm = 0;
      gridShiftNowRe = gridStIndRe[atomIndex];
      gridShiftNowIm = gridStIndIm[atomIndex];
      if(numGrid>0){
        numNlppTemp = numNlppAtom[atomType];
        for(l=0;l<=numLMax[atomType];l++){
          if(locOpt[atomType+1]!=l){
            for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
              radIndex = atomRadMap[atomType][countRad+iRad];
              // m=0
              dotReNow = dotRe[iAtom][iState*numNlppTemp+countNlppRe]*vpsNormList[radIndex];
              for(iDim=0;iDim<3;iDim++){
                dotDevReNow = dotDevRe[iAtom][iState*numNlppTemp*3+countNlppRe*3+iDim]*vpsNormList[radIndex];
                gridIndex = iState*numAtomFric*numGridMax*3+
                            iAtom*numGridMax*3+
                            iDim*numGridMax;
                daxpyBlasWrapper(numGrid,dotReNow,
                                 &vnlPhiDxAtomGridRe[gridShiftNowRe],1,
                                 &forceTemp[gridIndex],1);
                daxpyBlasWrapper(numGrid,dotDevReNow,
                                 &vnlPhiAtomGridRe[gridShiftNowRe],1,
                                 &forceTemp[gridIndex],1);
              }//endfor iDim
              gridShiftNowRe += numGrid;
              for(m=1;m<=l;m++){
                dotReNow = dotRe[iAtom][iState*numNlppTemp+countNlppRe+m]*2.0*vpsNormList[radIndex];
                dotImNow = dotIm[iAtom][iState*numNlppTemp+countNlppIm+m]*2.0*vpsNormList[radIndex];
                for(iDim=0;iDim<3;iDim++){
                  dotDevReNow = dotDevRe[iAtom][iState*numNlppTemp*3+(countNlppRe+m)*3+iDim]*vpsNormList[radIndex]*2.0*vpsNormList[radIndex];
                  dotDevImNow = dotDevIm[iAtom][iState*numNlppTemp*3+(countNlppIm+m-1)*3+iDim]*vpsNormList[radIndex]*2.0*vpsNormList[radIndex];

                  gridIndex = iState*numAtomFric*numGridMax*3+
                              iAtom*numGridMax*3+
                              iDim*numGridMax;
                  daxpyBlasWrapper(numGrid,dotReNow,
                                &vnlPhiDxAtomGridRe[gridShiftNowRe],1,
                                &forceTemp[gridIndex],1);
                  daxpyBlasWrapper(numGrid,dotImNow,
                                &vnlPhiDxAtomGridIm[gridShiftNowIm],1,
                                &forceTemp[gridIndex],1);

                  daxpyBlasWrapper(numGrid,dotDevReNow,
                                &vnlPhiAtomGridRe[gridShiftNowRe],1,
                                &forceTemp[gridIndex],1);
                  daxpyBlasWrapper(numGrid,dotDevImNow,
                                &vnlPhiAtomGridIm[gridShiftNowIm],1,
                                &forceTemp[gridIndex],1);
                }
                gridShiftNowRe += numGrid;
                gridShiftNowIm += numGrid;
              }//endfor m
              countNlppRe += l+1;
              countNlppIm += l;
            }//endfor iRad
            countRad += atomLRadNum[atomType][l];  
          }
          else{
            for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
              for(m=0;m<=l;m++){
                if(m!=0){
                  gridShiftNowRe += numGrid;
                  gridShiftNowIm += numGrid;
                }
                else{
                  gridShiftNowRe += numGrid;
                }//endif m
              }//endfor m
            }//endfor iRad
          }//endif locOpt  
        }//endfor l
      }//endif numGrid
    }//endfor iAtom
  }//endfor iState

/*======================================================================*/
/* V) D_nlpp|phi> Real to k-space                                       */

  iupper = numStateUpProc*numAtomFric;
  if((numStateUpProc*numAtomFric)%2!= 0){
    iupper = numStateUpProc*numAtomFric-1;
  }// endif

  for(is=1;is<=iupper;is=is+2){
    iState1 = (is-1)/numAtomFric;
    iAtom1 = (is-1)%numAtomFric;
    iState2 = (is-1)/numAtomFric;
    iAtom2 = (is-1)%numAtomFric;
    atomIndex1 = atomFricIndProc[iAtom1];
    atomIndex2 = atomFricIndProc[iAtom2];

/*======================================================================*/
/* II) FFT states back to real space                                    */
    
    for(iDim=0;iDim<3;iDim++){
      for(iGrid=0;iGrid<numGridAll*2;iGrid++)zfft[iGrid+1] = 0.0;

      numGrid = numGridNlppMap[atomIndex1];
      for(iGrid=0;iGrid<numGrid;iGrid++){
        gridIndex = gridNlppMap[atomIndex1][iGrid];
        zfft[2*gridIndex+1] = forceTemp[iState1*numAtomFric*numGridMax*3+
                                         iAtom1*numGridMax*3+
                                         iDim*numGridMax+iGrid];
      }
      numGrid = numGridNlppMap[atomIndex2];
      for(iGrid=0;iGrid<numGrid;iGrid++){
        gridIndex = gridNlppMap[atomIndex2][iGrid];
        zfft[2*gridIndex+2] = forceTemp[iState2*numAtomFric*numGridMax*3+
                                        iAtom2*numGridMax*3+
                                        iDim*numGridMax+iGrid];
      }

      if(fftw3dFlag==0){
        para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
      }
      else{
        para_fft_gen3d_bck_to_g_fftw3d_threads(zfft,cp_sclr_fft_pkg3d_sm);
      }

      ioff = iState1*numAtomFric*3*numCoeffUpTotal+
             iAtom1*3*numCoeffUpTotal+
             iDim*numCoeffUpTotal;
      ioff2 = iState2*numAtomFric*3*numCoeffUpTotal+
              iAtom2*3*numCoeffUpTotal+
              iDim*numCoeffUpTotal;

      if(fftw3dFlag==0){
        dble_upack_coef_sum(&daHChiReUp[ioff],&daHChiImUp[ioff],
                            &daHChiReUp[ioff2],&daHChiImUp[ioff2],
                            zfft,cp_sclr_fft_pkg3d_sm);
      }
      else{
        dble_upack_coef_sum_fftw3d(&daHChiReUp[ioff],&daHChiImUp[ioff],
                            &daHChiReUp[ioff2],&daHChiImUp[ioff2],
                            zfft,cp_sclr_fft_pkg3d_sm);

      }//endif fftw3dFlag
    }//endfor iDim
  }//endfor is

  if(numStateUpProc*numAtomFric%2!=0){
    is = numStateUpProc*numAtomFric;
    iState1 = (is-1)/numAtomFric;
    iAtom1 = (is-1)%numAtomFric;
    
    for(iDim=0;iDim<3;iDim++){
      for(iGrid=0;iGrid<numGridAll*2;iGrid++)zfft[iGrid+1] = 0.0;

      numGrid = numGridNlppMap[atomIndex1];
      for(iGrid=0;iGrid<numGrid;iGrid++){
        gridIndex = gridNlppMap[atomIndex1][iGrid];
        zfft[2*gridIndex+1] = forceTemp[iState1*numAtomFric*numGridMax*3+
                                         iAtom1*numGridMax*3+
                                         iDim*numGridMax+iGrid];
      }

      if(fftw3dFlag==0){
        para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
      }
      else{
        para_fft_gen3d_bck_to_g_fftw3d_threads(zfft,cp_sclr_fft_pkg3d_sm);
      }

      ioff = iState1*numAtomFric*3*numCoeffUpTotal+
             iAtom1*3*numCoeffUpTotal+
             iDim*numCoeffUpTotal;

      if(fftw3dFlag==0){
        sngl_upack_coef_sum(&daHChiReUp[ioff],&daHChiImUp[ioff],
                            zfft,cp_sclr_fft_pkg3d_sm);
      }
      else{
        sngl_upack_coef_sum_fftw3d(&daHChiReUp[ioff],&daHChiImUp[ioff],
                            zfft,cp_sclr_fft_pkg3d_sm);

      }//endif fftw3dFlag
    }//endfor iDim

  }//endif numStateUpProc*numAtomFric%2
 

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

  free(forceTemp);

    

/*==============================================================*/
}/*end routine*/
/*==============================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcDLocPhi(CP *cp,CLASS *class, GENERAL_DATA *general_data,
           double *coeffReUp, double *coeffImUp,
           double *coeffReDn, double *coeffImDn,
           double *daHChiReUp, double *daHChiImUp,
           double *daHChiReDn, double *daHChiImDn)
/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"


  for(icount=1;icount<=ngo;icount++){
/*======================================================================*/
/* I) Get the k vectors                                                 */
    aka = (double)(kastore[(icount+koff)]);
    akb = (double)(kbstore[(icount+koff)]);
    akc = (double)(kcstore[(icount+koff)]);
    xk = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
    yk = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
    zk = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
    g2 = xk*xk+yk*yk+zk*zk;
    g4 = g2*g2;
    g  = sqrt(g2);
    ak2[icount] = g2;

/*======================================================================*/
/* II) If break point number one or you are just starting out calculate */
/*     the helpful vectors                                              */
 

    if(ibreak1[(icount+koff)]==1||icount==1){
      for(ipart=1;ipart<=natm_use;ipart++){
        atemp = ewd_scr_x[ipart];
        btemp = ewd_scr_y[ipart];
        ctemp = ewd_scr_z[ipart];
        arg = (aka*atemp + akb*btemp + akc*ctemp)*tpi;
        //printf("aka %lg akb %lg akc %lg\n",aka,akb,akc);
        //printf("arg %lg\n",arg);
        helr[ipart] = cos(arg);
        heli[ipart] = sin(arg);
      }/*endfor*/
    }/*endif*/
 

/*======================================================================*/
/* III) Get the external potential                                      */
/*               (interaction of electron with particles)               */
  
    if(ipseud_opt==1){
      for(itype=1;itype<=natm_typ;itype++){
        index_atm[itype] =  (itype-1)*nsplin_g*(n_ang_max+1)*n_rad_max
                         +  loc_opt[itype]*nsplin_g*n_rad_max;
      }/*endfor*/

      get_vpsnow(index_atm,nsplin_g,gmin_spl,dg_spl,g,
                  vps0,vps1,vps2,vps3,vtemp,iatm_typ,natm_typ,natm_use,1);


/*-------------------------*/
/* charge correction */
      for(i=1; i<= natm_use; i++){
        vtemp[i] += -fpi*(q_tmp[i] - q_pseud[iatm_typ[i]])/g2;
      }
    }else{  /*q_temp is q */
      get_vpslong(natm_use,vtemp,g2,q_tmp,alpha_conv_dual,pivol);
    }/*endif*/

 /*----------------------------------------------------------------------*/
 /* Cluster boundary condition correction                                */

    if( (iperd != 3) && (idens_opt==0)){
      for(ipart=1;ipart<=natm_use;ipart++){
        vtemp[ipart] -= q[ipart]*clus_corr_r[icount];
      }/* endfor */
    }/* endif cluster boundary conditions */

/*----------------------------------------------------------------------*/

    vextr[icount]  =  ddot1(natm_use,helr,1,vtemp,1)*rvol;
    vexti[icount]  = -ddot1(natm_use,heli,1,vtemp,1)*rvol;

/*======================================================================*/
/* IV) Get the real and imag parts of the structure factor              */

    if(idens_opt==0){/*create charge weighted structure factor*/
      sumr = ddot1(natm_use,helr,1,q,1);
      sumi = ddot1(natm_use,heli,1,q,1);
      smag = sumr*sumr+sumi*sumi;
    }/*endif*/

/*======================================================================*/
/* VI) Get the force on the particles, checking for CBC                 */


    if((iperd == 0)||(idens_opt==1)){
      for(ipart=1;ipart<=natm_use;ipart++){
        srx = xk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
        sry = yk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
        srz = zk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
        six = xk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
        siy = yk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
        siz = zk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
        fx_tmp[ipart] += (srx*heli[ipart] - six*helr[ipart]);
        fy_tmp[ipart] += (sry*heli[ipart] - siy*helr[ipart]);
        fz_tmp[ipart] += (srz*heli[ipart] - siz*helr[ipart]);
      }/* endfor */
    }else{
      sumr_h = sumr;
      sumi_h = sumi;
      sumr = sumr*preg*2.0;
      sumi = sumi*preg*2.0;
      for(ipart=1;ipart<=natm_use;ipart++){
        srx = 0.0;
        sry = 0.0;
        srz = 0.0;
        if(ewaldLocalOpt==0||ewaldLocalOpt==2){
          pre = sumr*q[ipart];
          srx = xk*pre;
          sry = yk*pre;
          srz = zk*pre;
          pre = sumi*q[ipart];
          six = xk*pre;
          siy = yk*pre;
          siz = zk*pre;
        }

        fxCl[ipart] += (srx*heli[ipart]  - six*helr[ipart]);
        fyCl[ipart] += (sry*heli[ipart]  - siy*helr[ipart]);
        fzCl[ipart] += (srz*heli[ipart]  - siz*helr[ipart]);

        if(ewaldLocalOpt==0||ewaldLocalOpt==1){
          pre = 2.0*rhocr[icount]*vtemp[ipart]*rvol;
          srx += xk*pre;
          sry += yk*pre;
          srz += zk*pre;
          pre = -2.0*rhoci[icount]*vtemp[ipart]*rvol;
          six += xk*pre;
          siy += yk*pre;
          siz += zk*pre;
        }

        /*
        srx = xk*(sumr*q[ipart]+2.0*rhocr[icount]*vtemp[ipart]*rvol);
        sry = yk*(sumr*q[ipart]+2.0*rhocr[icount]*vtemp[ipart]*rvol);
        srz = zk*(sumr*q[ipart]+2.0*rhocr[icount]*vtemp[ipart]*rvol);
        six = xk*(sumi*q[ipart]-2.0*rhoci[icount]*vtemp[ipart]*rvol);
        siy = yk*(sumi*q[ipart]-2.0*rhoci[icount]*vtemp[ipart]*rvol);
        siz = zk*(sumi*q[ipart]-2.0*rhoci[icount]*vtemp[ipart]*rvol);
        */

        //debug k space classical force only
        /*
        srx = xk*(sumr*q[ipart]);
        sry = yk*(sumr*q[ipart]);
        srz = zk*(sumr*q[ipart]);
        six = xk*(sumi*q[ipart]);
        siy = yk*(sumi*q[ipart]);
        siz = zk*(sumi*q[ipart]);
        */
        fx_tmp[ipart] += (srx*heli[ipart]  - six*helr[ipart]);
        fy_tmp[ipart] += (sry*heli[ipart]  - siy*helr[ipart]);
        fz_tmp[ipart] += (srz*heli[ipart]  - siz*helr[ipart]);
        // debug print x,y,z
        //printf("icount %i ipart %i rhocr %lg rhoci %lg srx %lg six %lg heli %lg helr %lg fxtmp %lg\n",icount,ipart,rhocr[icount],rhoci[icount],srx,six,heli[ipart],helr[ipart],fx_tmp[ipart]);

      }/*endfor*/
    } /* endif cluster BC */

/*======================================================================*/
/* VII) If break point two, increment the helpful vectors                 */

    if(ibreak2[(icount+koff)]==1){
      for(ipart=1;ipart<=natm_use;ipart++){
        temp = helr[ipart];
        helr[ipart] = helr[ipart]*cossc[ipart] - heli[ipart]*sinsc[ipart];
        heli[ipart] = heli[ipart]*cossc[ipart] + temp*sinsc[ipart];
      }/*endfor*/
    }/*endif*/
  }

/*======================================================================*/
/* VIII) g=0 term (local pseudopotential) including term for CBCs       */

  if((myid_state+1)==np_states){
    if(ipseud_opt==1){
      ak2[(ngo+1)] = 0.0;
      for(ipart=1;ipart<=natm_use;ipart++){
        vtemp[ipart] = gzvps[iatm_typ[ipart]];
      }/*endfor*/
      vextr[(ngo+1)] =  dsum1(natm_use,vtemp,1)*rvol;
      vexti[(ngo+1)] = 0.0;
    }else{ /*large sparse grid */
      vextr[(ngo+1)] = 0.0;
      bgr  = dsum1(natm_use,q_tmp,1);
      bgr  = bgr*M_PI/(alpha_conv_dual*alpha_conv_dual*vol);
      vextr[(ngo+1)] = bgr;
    }/*endif*/

    if( (iperd!=3) && (idens_opt==0) ) {
      vextr[(ngo+1)] -= dsum1(natm_use,q_tmp,1)*clus_corr_r[(ngo+1)]*rvol;
    }/*endif*/

    *pseud_hess_loc = vextr[(ngo+1)];

    if( (iperd>0) && (iperd!=3) && (idens_opt==0)) {
      q_sum1 = dsum1(natm_use,q,1);
      if(iperd == 2){
         vrecip += 0.5*q_sum1*q_sum1*rvol*(clus_corr_r[(ngo+1)]
                 + M_PI/(alp_clus*alp_clus));
       } else {
         vrecip += 0.5*q_sum1*q_sum1*clus_corr_r[(ngo+1)]*rvol;
       }/* endif iperd */
    }/*endif*/
  }/*endif*/


/*==============================================================*/
}/*end routine*/
/*==============================================================*/

