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
  double *zfft_tmp         = cpscr->cpscr_wave.zfft_tmp;

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

  /*
  for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
    for(is=0;is<numStateFric;is++){
      printf("dotRe %lg dotIm %lg dotDevRe %lg dotDevIm %lg\n",dotRe[iAtom][(is)*numNlppTemp],
             dotIm[iAtom][(is)*(numNlppTemp-1)],dotDevRe[iAtom][(is)*numNlppTemp*3],
             dotDevIm[iAtom][(is)*(numNlppTemp-1)*3]);
    }
  }
  */

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
                      //printf("iDim %i iState %i jState %i l %i m %i dotRe %lg dotIm %lg dotDevRe %lg dotDevIm %lg vpsNormList %lg volInv %lg temp %lg\n",
                      //       iDim,iState,jState,l,m,dotRe[iAtom][iState*numNlppTemp+countNlppRe+m],
                      //       dotIm[iAtom][iState*(numNlppTemp-1)+countNlppIm+m-1],
                      //       dotDevRe[iAtom][jState*numNlppTemp*3+(countNlppRe+m)*3+iDim],
                      //       dotDevIm[iAtom][jState*(numNlppTemp-1)*3+(countNlppIm+m-1)*3+iDim],
                      //       vpsNormList[radIndex],volInv,temp);
                    }
                    else{
                      temp = 
                          (dotRe[iAtom][iState*numNlppTemp+countNlppRe+m]*
                          dotDevRe[iAtom][jState*numNlppTemp*3+(countNlppRe+m)*3+iDim]+
                          dotRe[iAtom][jState*numNlppTemp+countNlppRe+m]*
                          dotDevRe[iAtom][iState*numNlppTemp*3+(countNlppRe+m)*3+iDim])*
                          vpsNormList[radIndex]*volInv;
                      hDevMat[(iAtom*3+iDim)*numStateFric*numStateFric+iState*numStateFric+jState] += temp;
                      //printf("iDim %i iState %i jState %i l %i m %i dotRe %lg dotIm 0.0 dotDevRe %lg dotDevIm 0.0 vpsNormList %lg volInv %lg temp %lg\n",
                      //       iDim,iState,jState,l,m,dotRe[iAtom][iState*numNlppTemp+countNlppRe+m],
                      //       dotDevRe[iAtom][jState*numNlppTemp*3+(countNlppRe+m)*3+iDim],
                      //       vpsNormList[radIndex],volInv,temp);
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
              //printf("l %i m %i dotRe %lg dotIm %lg dotDevRe %lg %lg %lg dotDevIm %lg %lg %lg\n",l,m,dotRe[countNlppRe+m],dotIm[countNlppIm+m-1],dotDevRe[(countNlppRe+m)*3],dotDevRe[(countNlppRe+m)*3+1],dotDevRe[(countNlppRe+m)*3+2],dotDevIm[(countNlppIm+m-1)*3],dotDevIm[(countNlppIm+m-1)*3+1],dotDevIm[(countNlppIm+m-1)*3+2]);

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
              //printf("l %i m %i dotRe %lg dotDevRe %lg %lg %lg\n",l,m,dotRe[countNlppRe],dotDevRe[(countNlppRe)*3],dotDevRe[(countNlppRe)*3+1],dotDevRe[(countNlppRe)*3+2]);
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

  double sum;
  
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
  double *fxLoc,*fyLoc,*fzLoc;

  fxLoc = (double*)cmalloc(numAtomTot*sizeof(double));
  fyLoc = (double*)cmalloc(numAtomTot*sizeof(double));
  fzLoc = (double*)cmalloc(numAtomTot*sizeof(double));


  //for(iState=0;iState<numStateFric;iState++){
  rhoCalcRealFriction(general_data,cp,class,ksStateChemPotRe,ksStateChemPotIm,
                      wfReal,numStateFric);
  //}

  cpewald->ewaldLocalOpt = 1;	  
  for(iState=0;iState<numStateFric;iState++){
    for(jState=0;jState<numStateFric;jState++){
      /* Calculate Real Space "Density" phi_i^*(r)phi_j(r) */
      if(myidState==0){
        //sum = 0.0;
        for(iGrid=0;iGrid<rhoRealGridTot;iGrid++){
          rhoTemp[iGrid] = wfReal[iState*rhoRealGridTot+iGrid]*wfReal[jState*rhoRealGridTot+iGrid];
          //sum += rhoTemp[iGrid];
        }
        //sum /= rhoRealGridTot;
        //printf("rho dot %lg\n",sum);
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

      printf("rhoCoeffReUp %i rhoCoeffImUp %lg\n",rhoCoeffReUp[1],rhoCoeffImUp[1]);
      /* Calculate F_local(ij) */
      for(iAtom=1;iAtom<=numAtomTot;iAtom++){
        fx[iAtom] = 0.0;
        fy[iAtom] = 0.0;
        fz[iAtom] = 0.0;
      }
      calcLocExtPostScf(class,general_data,cp,cpcoeffs_pos,clatoms_pos);
      printf("myid %i force CO %lg %lg %lg\n",myidState,fx[9],fy[9],fz[9]);
      /* Reduce Force    */
      if(numProcStates==1){
        memcpy(&fxLoc[0],&fx[1],numAtomTot*sizeof(double));
        memcpy(&fyLoc[0],&fy[1],numAtomTot*sizeof(double));
        memcpy(&fzLoc[0],&fz[1],numAtomTot*sizeof(double));
      }
      else{
        Reduce(&fx[1],&fxLoc[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,comm_states);
        Reduce(&fy[1],&fyLoc[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,comm_states);
        Reduce(&fz[1],&fzLoc[0],numAtomTot,MPI_DOUBLE,MPI_SUM,0,comm_states);
      }

      /* Broadcast Force */
      Bcast(&fxLoc[0],numAtomTot,MPI_DOUBLE,0,comm_states);
      Bcast(&fyLoc[0],numAtomTot,MPI_DOUBLE,0,comm_states);
      Bcast(&fzLoc[0],numAtomTot,MPI_DOUBLE,0,comm_states);
      for(iAtom=0;iAtom<numAtomFricProc;iAtom++){
        atomIndex = atomFricIndProc[iAtom];
        printf("Fx %lg %lg %lg atomIndex %i\n",fxLoc[atomIndex],fyLoc[atomIndex],fzLoc[atomIndex],atomIndex);
        vlocDevMat[(iAtom*3)*numStateFric*numStateFric+iState*numStateFric+jState] = -fxLoc[atomIndex];
        vlocDevMat[(iAtom*3+1)*numStateFric*numStateFric+iState*numStateFric+jState] = -fyLoc[atomIndex];
        vlocDevMat[(iAtom*3+2)*numStateFric*numStateFric+iState*numStateFric+jState] = -fzLoc[atomIndex];
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
  CPOPTS *cpopts                = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  CPSCR *cpscr                  = &(cp->cpscr);  
  PSEUDO *pseudo                = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal       = &(pseudo->pseudoReal);
  CELL *cell                    = &(general_data->cell);
  CLATOMS_POS *clatoms_pos      = &(class->clatoms_pos[1]);
  CLATOMS_INFO *clatoms_info    = &(class->clatoms_info);
  ATOMMAPS *atommaps            = &(class->atommaps);
  STAT_AVG *stat_avg            = &(general_data->stat_avg);
  CPEWALD *cpewald              = &(cp->cpewald);
  METALLIC *metallic            = stodftInfo->metallic;
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm = &(cp->cp_sclr_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

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
  int nfft              = cp_sclr_fft_pkg3d_sm->nfft;
  int numGridAll        = nfft/2;
  int cpDualGridOptOn   = cpopts->cp_dual_grid_opt;
  int rhoRealGridNum    = stodftInfo->rhoRealGridNum;
  int numInterpPmeDual  = pseudo->n_interp_pme_dual;
  int iperd             = cell->iperd;
  int numChemPot        = stodftInfo->numChemPot;
  int occNumber         = stodftInfo->occNumber;
  int myidState         = communicate->myid_state;
  int numProcStates     = communicate->np_states;
  int numAtomTot        = clatoms_info->natm_tot;
  int iState,jState,iCoeff,iChem,iAtom,iGrid;
  int smearOpt          = stodftInfo->smearOpt;
  int numAtomFricProc   = metallic->numAtomFricProc;
  int numAtomFric       = metallic->numAtomFric;
  int numAtomType       = atommaps->natm_typ;
  int fftw3dFlag        = cpopts->fftw3dFlag;
  int numGridMax,numGrid;

  int atomType,atomIndex;
  int numNlppTemp;
  int gridShiftNowRe,gridShiftNowIm;

  int is,i,iupper;
  int ioff,ncoef1,ioff2;
  int iii,iis,nis;
  int iPart;
  int iAng,iDim,countNlppRe,countNlppIm;
  int iRad,radIndex,countRad;
  int l,m;
  int gridIndex;
  int iState1,iState2,iAtom1,iAtom2;
  int atomIndex1,atomIndex2;
  int myid_state = communicate->myid_state;
  int np_states  = communicate->np_states;

  int *atomFricIndProc = metallic->atomFricIndProc;
  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *numLMax = pseudoReal->numLMax;
  int *numNlppAtom = pseudoReal->numNlppAtom;
  int *locOpt = pseudo->loc_opt;
  int *numGridNlppMap = pseudoReal->numGridNlppMap;
  int *gridStIndRe = pseudoReal->gridStIndRe;
  int *gridStIndIm = pseudoReal->gridStIndIm;
  int **atomLRadNum = pseudoReal->atomLRadNum;
  int **atomRadMap = pseudoReal->atomRadMap;
  int **gridNlppMap = pseudoReal->gridNlppMap;

  double tpi = 2.0*M_PI;
  double vol        = cell->vol;
  double volInv     = 1.0/vol;
  double temp;
  double dotReNow,dotDevReNow;
  double dotImNow,dotDevImNow;

  double *ksStateChemPotRe = metallic->ksStateChemPotRe;
  double *ksStateChemPotIm = metallic->ksStateChemPotIm;
  double *zfft             = cpscr->cpscr_wave.zfft;
  double *zfft_tmp         = cpscr->cpscr_wave.zfft_tmp;
  double *wfReal;
  double *forceTemp;
  double *vpsNormList = pseudoReal->vpsNormList;
  double *vnlPhiAtomGridRe = pseudoReal->vnlPhiAtomGridRe;
  double *vnlPhiAtomGridIm = pseudoReal->vnlPhiAtomGridIm;
  double *vnlPhiDxAtomGridRe = pseudoReal->vnlPhiDxAtomGridRe;
  double *vnlPhiDxAtomGridIm = pseudoReal->vnlPhiDxAtomGridIm;
  double *vnlPhiDyAtomGridRe = pseudoReal->vnlPhiDyAtomGridRe;
  double *vnlPhiDyAtomGridIm = pseudoReal->vnlPhiDyAtomGridIm;
  double *vnlPhiDzAtomGridRe = pseudoReal->vnlPhiDzAtomGridRe;
  double *vnlPhiDzAtomGridIm = pseudoReal->vnlPhiDzAtomGridIm;

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

  wfReal = (double*)cmalloc(numGridAll*sizeof(double));

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
    for(iAtom=0;iAtom<numAtomFric;iAtom++){
      atomIndex = atomFricIndProc[iAtom];
      atomType = iAtomAtomType[atomIndex+1]-1;
      numNlppTemp = numNlppAtom[atomType];
      calcNlppDot(class,general_data,cp,atomIndex,iAtom,wfReal,&dotRe[iAtom][(is)*numNlppTemp],
                  &dotIm[iAtom][(is)*(numNlppTemp-1)],&dotDevRe[iAtom][(is)*numNlppTemp*3],
                  &dotDevIm[iAtom][(is)*(numNlppTemp-1)*3]);
    }
  }//endfor is

  if(numStateUpProc%2!= 0){
    is = numStateUpProc;
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
    for(is=0;is<numStateUpProc;is++){
      printf("dotRe %lg dotIm %lg dotDevRe %lg dotDevIm %lg\n",dotRe[iAtom][(is)*numNlppTemp],
             dotIm[iAtom][(is)*(numNlppTemp-1)],dotDevRe[iAtom][(is)*numNlppTemp*3],
             dotDevIm[iAtom][(is)*(numNlppTemp-1)*3]);
    }
  }

/*======================================================================*/
/* V) Calculate D_nlpp|phi>                                             */

  numGridMax = numGridNlppMap[0];
  for(iPart=0;iPart<numAtomTot;iPart++){
    if(numGridNlppMap[iPart]>numGridMax)numGridMax = numGridNlppMap[iPart];
  }

  forceTemp = (double*)calloc(numStateUpProc*numAtomFric*3*numGridMax,sizeof(double));


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
/* II) FFT states back to k space                                       */
    
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

  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS *clatoms_pos = &(class->clatoms_pos[1]);
  CELL *cell = &(general_data->cell);
  EWALD *ewald = &(general_data->ewald);
  CPEWALD *cpewald = &(cp->cpewald);
  CPSCR *cpscr = &(cp->cpscr);
  PSEUDO *pseudo = &(cp->pseudo);
  EWD_SCR *ewd_scr = &(class->ewd_scr);
  CPOPTS *cpopts = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info  = &(cp->cpcoeffs_info);
  ATOMMAPS *atommaps = &(class->atommaps);
  COMMUNICATE *communicate = &(cp->communicate);
  FOR_SCR *for_scr = &(class->for_scr);
  STODFTINFO *stodftInfo        = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos  = cp->stodftCoefPos;  
  METALLIC *metallic            = stodftInfo->metallic;
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm = &(cp->cp_sclr_fft_pkg3d_sm);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  

  int idens_opt,ipseud_opt;
  int istart,ngo,irem,idiv;
  int ipart,jpart,iii,itype,i,j,k;
  int icount,koff,natm_use;
  int realSparseOpt = cpewald->realSparseOpt;
  int ewaldLocalOpt = cpewald->ewaldLocalOpt;
  int cp_dual_grid_opt = cpopts->cp_dual_grid_opt;
  int idual_switch = 0;
  int cp_lsda = cpopts->cp_lsda;

  double falp2,falp_clus2,vol,rvol,pivol,fpi,arg,q_sum1,bgr;
  double aka,akb,akc,xk,yk,zk,atemp,btemp,ctemp;
  double xtemp,ytemp,ztemp;
  double sumr,sumi,g2,g4,preg,prep,tpi,pi,g;
  double sumr_h,sumi_h;
  double srx,sry,srz,six,siy,siz,temp,smag;
  double dx,dy,dz;
  double sx,sy,sz;
  double asx,asy,asz;
  double phase;
  double argp,fargp,argm,fargm,area;
  double pre,pre2;

/*--------------------------------------------*/
/*         Local Pointer declarations         */

  /*------------------*/
  /* Atom information */
  int natm_tot      = clatoms_info->natm_tot;
  int hess_calc     = clatoms_info->hess_calc;
  double *x         = clatoms_pos->x;
  double *y         = clatoms_pos->y;
  double *z         = clatoms_pos->z;
  double *q         = clatoms_info->q;
  int natm_typ      = atommaps->natm_typ;
  int *index_atm    = for_scr->index_atm;
  int *iatm_typ;            /*Assigned below based on flags */
  int *iatm_typ_full;

  /*--------------------------------*/
  /* Cell and pressure information */
  int iperd                 = cell->iperd;
  int cp_ptens              = cpopts->cp_ptens_calc;
  double *cp_box_center     = cell->cp_box_center;
  double *cp_box_center_rel = cell->cp_box_center_rel;
  double *hmat;             /* Assigned below based on flags */
  double *hmati;
  double *hmat_big          = cell->hmat;
  double *hmati_big         = cell->hmati;

  /*----------------------*/
  /* G-vector information */
  int *kastore,*kbstore,*kcstore; // Assigned below based on flags
  int *ibreak1,*ibreak2;
  double *vextr,*vexti;
  double *dvextr,*dvexti;
  double *rhocr,*rhoci;
  double *ak2;
  int nktot;

  /*----------------------------------------------*/
  /* Pseudo-potential and Reduced Periodicity info*/
  int nsplin_g         = pseudo->nsplin_g;
  int n_rad_max        = pseudo->n_rad_max;
  double *clus_corr_r  = ewald->clus_corr_r;
  double *dclus_corr_r = ewald->dclus_corr_r;
  double alpha_conv_dual = pseudo->alpha_conv_dual;
  double dg_spl        = pseudo->dg_spl;
  double gmin_spl      = pseudo->gmin_spl;
  double *vps0         = pseudo->vps0;
  double *vps1         = pseudo->vps1;
  double *vps2         = pseudo->vps2;
  double *vps3         = pseudo->vps3;
  double *dvps0        = pseudo->dvps0;
  double *dvps1        = pseudo->dvps1;
  double *dvps2        = pseudo->dvps2;
  double *dvps3        = pseudo->dvps3;
  double *gzvps        = pseudo->gzvps;
  double *q_pseud      = pseudo->q_pseud;
  int n_ang_max        = pseudo->n_ang_max;
  int *loc_opt         = pseudo->loc_opt;
  int np_loc_cp_box    = pseudo->np_loc_cp_box;
  int *ip_loc_cp_box   = pseudo->ip_loc_cp_box;

  /*---------------------------------*/
  /* Ewald and ewald scr information */
  double alp_ewald  = ewald->alp_ewd;
  double alp_clus   = ewald->alp_clus;
  double *cossc     = ewd_scr->cossc;
  double *sinsc     = ewd_scr->sinsc;
  double *helr      = ewd_scr->helr;
  double *heli      = ewd_scr->heli;
  double *vtemp     = ewd_scr->temp;
  double *dvtemp    = ewd_scr->vtemp_now;
  double *ewd_scr_x = ewd_scr->x;
  double *ewd_scr_y = ewd_scr->y;
  double *ewd_scr_z = ewd_scr->z;
  double *q_tmp     = ewd_scr->q;

  /*---------------------------*/
  /* Communication information */
  int myid_state    = communicate->myid_state;
  int np_states     = communicate->np_states;
  MPI_Comm comm     = communicate->comm_states;
  //FILE *fvrecip = fopen("vrecip","w");

  /*---------------------------*/
  /* Friction Calculation      */
  int numAtomFric      = metallic->numAtomFric;
  int numCoeff         = cp_sclr_fft_pkg3d_sm->ncoef;
  int *atomFricIndProc = metallic->atomFricIndProc;
  int nfft             = cp_sclr_fft_pkg3d_sm->nfft;
  int nfft2            = nfft/2;
  int atomIndex;
  int iupper,is,ioff,ioff2,iCoeff;
  int numStateUpProc   = cpcoeffs_info->nstate_up_proc;
  int numStateDnProc   = cpcoeffs_info->nstate_dn_proc;
  int numCoeffUpTotal  = numStateUpProc*numCoeff;
  int numCoeffDnTotal  = numStateDnProc*numCoeff;
  int fftw3dFlag       = cpopts->fftw3dFlag;

  double *davextr = (double*)cmalloc(3*numAtomFric*(ngo+1)*sizeof(double))-1;
  double *davexti = (double*)cmalloc(3*numAtomFric*(ngo+1)*sizeof(double))-1;
  double *coeffReTemp = (double*)cmalloc(numCoeff*sizeof(double))-1;
  double *coeffImTemp = (double*)cmalloc(numCoeff*sizeof(double))-1;
  double *coeffReTemp2 = (double*)cmalloc(numCoeff*sizeof(double))-1;
  double *coeffImTemp2 = (double*)cmalloc(numCoeff*sizeof(double))-1;
  double *zfft = cpscr->cpscr_wave.zfft;
  double *zfft_tmp = cpscr->cpscr_wave.zfft_tmp;
  double *v_loc_real_up = (double*)cmalloc(nfft2*sizeof(double))-1;
  double *v_loc_real_dn = (double*)cmalloc(nfft2*sizeof(double))-1;

/*======================================================================*/
/* 0) Assign local pointers                                             */

  //Let's test the following:
  //generating vext on dense k grid, then transform back to r space
  //then pooling it by select the sparse grid
  //realSparseOpt = 0;

  if(cp_dual_grid_opt < 2 || idual_switch == 0){
    /* large sparse grid when cp_dual_grid_opt == 2*/
    idens_opt = 0;
    ipseud_opt= (cp_dual_grid_opt==2 ? 0 : 1);

    //if(realSparseOpt==0){
    kastore   = ewald->kastr;
    kbstore   = ewald->kbstr;
    kcstore   = ewald->kcstr;
    ak2       = cpewald->ak2;
    nktot     = ewald->nktot;
    ibreak1   = ewald->ibrk1;
    ibreak2   = ewald->ibrk2;
    vextr     = cpscr->cpscr_loc.vextr;
    vexti     = cpscr->cpscr_loc.vexti;
    dvextr    = cpscr->cpscr_loc.dvextr;
    dvexti    = cpscr->cpscr_loc.dvexti;
    rhocr     = cpscr->cpscr_rho.rhocr_up;
    rhoci     = cpscr->cpscr_rho.rhoci_up;
    //ak2       = cpewald->ak2;
    //nktot     = ewald->nktot;
    hmat      = cell->hmat;
    hmati     = cell->hmati;
    natm_use  = natm_tot;
    iatm_typ  = atommaps->iatm_atm_typ;
  }else{
    // small dense grid 
    idens_opt = 1;
    ipseud_opt= 1;
    kastore   = cpewald->kastr_dens_cp_box;
    kbstore   = cpewald->kbstr_dens_cp_box;
    kcstore   = cpewald->kcstr_dens_cp_box;
    ibreak1   = cpewald->ibrk1_dens_cp_box; /*DY edit*/
    ibreak2   = cpewald->ibrk2_dens_cp_box; /*DY edit*/
    vextr     = cpscr->cpscr_loc.vextr_dens_cp_box;
    vexti     = cpscr->cpscr_loc.vexti_dens_cp_box;
    rhocr     = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
    rhoci     = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
    ak2       = cpewald->ak2_dens_cp_box;
    nktot     = cpewald->nktot_dens_cp_box;
    hmat      = cell->hmat_cp;
    hmati     = cell->hmati_cp;
    natm_use  = np_loc_cp_box;
    iatm_typ  = for_scr->iexcl;
    iatm_typ_full = atommaps->iatm_atm_typ;
 }/*endif*/


/*======================================================================*/
/* I) Get more useful constants                                         */

  pi    = M_PI;
  tpi   = 2.0*pi;
  fpi   = 4.0*pi;

  vol   = getdeth(hmat);
  rvol  = 1.0/vol;
  pivol = vol/4.0/pi;
  falp2 = 4.0*alp_ewald*alp_ewald;
  falp_clus2 = 4.0*alp_clus*alp_clus;

/*======================================================================*/
/* II) Find cos and sin of sc components of the particles               */
/*    ( hmnati rvec = svec   r=(x,y,z) s=(a,b,c) )                       */

  if(idens_opt==0){
    for(ipart=1;ipart<=natm_use;ipart++){
      xtemp = x[ipart];
      ytemp = y[ipart];
      ztemp = z[ipart];

      q_tmp[ipart]     = q[ipart];
      ewd_scr_x[ipart] = xtemp*hmati[1]
                       + ytemp*hmati[4]
                       + ztemp*hmati[7];
      ewd_scr_y[ipart] = xtemp*hmati[2]
                       + ytemp*hmati[5]
                       + ztemp*hmati[8];
      ewd_scr_z[ipart] = xtemp*hmati[3]
                       + ytemp*hmati[6]
                       + ztemp*hmati[9];
      ctemp            = ewd_scr_z[ipart]*tpi;
      cossc[ipart]     = cos(ctemp);
      sinsc[ipart]     = sin(ctemp);
      //printf("ipart %i ewd_scr %lg %lg %lg\n",ipart,ewd_scr_x[ipart],ewd_scr_y[ipart],ewd_scr_z[ipart]);

    }/*endfor*/
  }else{
    for(ipart=1;ipart<=natm_use;ipart++){
      iatm_typ[ipart]  = iatm_typ_full[ip_loc_cp_box[ipart]];
      q_tmp[ipart]     = q[ip_loc_cp_box[ipart]];
      dx               = x[ip_loc_cp_box[ipart]] - cp_box_center[1];
      dy               = y[ip_loc_cp_box[ipart]] - cp_box_center[2];
      dz               = z[ip_loc_cp_box[ipart]] - cp_box_center[3];

      asx              = dx*hmati_big[1]+dy*hmati_big[4]+dz*hmati_big[7];
      asy              = dx*hmati_big[2]+dy*hmati_big[5]+dz*hmati_big[8];
      asz              = dx*hmati_big[3]+dy*hmati_big[6]+dz*hmati_big[9];
      sx               = asx - NINT(asx);
      sy               = asy - NINT(asy);
      sz               = asz - NINT(asz);
      dx               = sx*hmat_big[1]+sy*hmat_big[4]+sz*hmat_big[7];
      dy               = sx*hmat_big[2]+sy*hmat_big[5]+sz*hmat_big[8];
      dz               = sx*hmat_big[3]+sy*hmat_big[6]+sz*hmat_big[9];

      xtemp            = dx + cp_box_center_rel[1];
      ytemp            = dy + cp_box_center_rel[2];
      ztemp            = dz + cp_box_center_rel[3];
      ewd_scr_x[ipart] = xtemp*hmati[1]
                       + ytemp*hmati[4]
                       + ztemp*hmati[7];
      ewd_scr_y[ipart] = xtemp*hmati[2]
                       + ytemp*hmati[5]
                       + ztemp*hmati[8];
      ewd_scr_z[ipart] = xtemp*hmati[3]
                       + ytemp*hmati[6]
                       + ztemp*hmati[9];
      ctemp            = ewd_scr_z[ipart]*tpi;
      cossc[ipart]     = cos(ctemp);
      sinsc[ipart]     = sin(ctemp);
    }/*endfor*/
  }/*endif*/
  
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
      for(ipart=0;ipart<=natm_use;ipart++){
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

    for(i=0;i<numAtomFric;i++){
      atomIndex = atomFricIndProc[i]+1;
      davextr[i*3*ngo+icount] = heli[atomIndex]*vtemp[atomIndex]*xk;
      davexti[i*3*ngo+icount] = -helr[atomIndex]*vtemp[atomIndex]*xk;
      davextr[i*3*ngo+ngo+icount] = heli[atomIndex]*vtemp[atomIndex]*yk;
      davexti[i*3*ngo+ngo+icount] = -helr[atomIndex]*vtemp[atomIndex]*yk;
      davextr[i*3*ngo+2*ngo+icount] = heli[atomIndex]*vtemp[atomIndex]*zk;
      davexti[i*3*ngo+2*ngo+icount] = -helr[atomIndex]*vtemp[atomIndex]*zk;
    }

    //vextr[icount]  =  ddot1(natm_use,helr,1,vtemp,1)*rvol;
    //vexti[icount]  = -ddot1(natm_use,heli,1,vtemp,1)*rvol;

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


  for(i=0;i<numAtomFric;i++){
    for(j=0;j<3;j++){
/*====================================================================*/
/*  II) single pack the derivative of local pp potential 
        for fourier transform routine    */
      
      sngl_pack_coef(&davextr[i*3*ngo+j*ngo],&davexti[i*3*ngo+j*ngo],zfft,cp_para_fft_pkg3d_lg);

/*====================================================================*/
/* III) fourier transform ks potential to real space exp(-igr)        */

      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_lg);

/*====================================================================*/
/* IV) Contract and unpack rho for dual option, otherwise just upack  */

      sngl_upack_rho(zfft,v_loc_real_up,cp_para_fft_pkg3d_lg);


/*====================================================================*/
/* V) Assign up to down for LSDA because vext is vext        */

      if(cp_lsda==1){
        for(i=1;i<=nfft2;i++){v_loc_real_dn[i]=v_loc_real_up[i];}
      }//endif

/*====================================================================*/
/* VI) Apply v_ks to orbitals                                         */

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
/* III) Calculate Da(V_loc)|Psi>                      */

        cp_vpsi(zfft,v_loc_real_up,nfft);

        if(fftw3dFlag==0){
          para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
        }
        else{
          para_fft_gen3d_bck_to_g_fftw3d_threads(zfft,cp_sclr_fft_pkg3d_sm);
        }

        if(fftw3dFlag==0){
          dble_upack_coef_sum(coeffReTemp,coeffImTemp,
                              coeffReTemp2,coeffImTemp2,
                              zfft,cp_sclr_fft_pkg3d_sm);
        }
        else{
          dble_upack_coef_sum_fftw3d(coeffReTemp,coeffImTemp,
                              coeffReTemp2,coeffImTemp2,
                              zfft,cp_sclr_fft_pkg3d_sm);
        }//endif fftw3dFlag


        ioff = (is-1)*numAtomFric*3*numCoeffUpTotal+
               i*3*numCoeffUpTotal+
               j*numCoeffUpTotal;
        ioff2 = is*numAtomFric*3*numCoeffUpTotal+
                i*3*numCoeffUpTotal+
                j*numCoeffUpTotal;
        for(iCoeff=0;iCoeff<numCoeff;iCoeff++){
          daHChiReUp[ioff+iCoeff] += coeffReTemp[iCoeff+1];
          daHChiImUp[ioff+iCoeff] += coeffImTemp[iCoeff+1];
        }
        for(iCoeff=0;iCoeff<numCoeff;iCoeff++){
          daHChiReUp[ioff2+iCoeff] += coeffReTemp2[iCoeff+1];
          daHChiImUp[ioff2+iCoeff] += coeffImTemp2[iCoeff+1];
        }

      }//endfor is

      if(numStateUpProc%2!= 0){
        is = numStateUpProc;
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

        cp_vpsi(zfft,v_loc_real_up,nfft);

        if(fftw3dFlag==0){
          para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
        }
        else{
          para_fft_gen3d_bck_to_g_fftw3d_threads(zfft,cp_sclr_fft_pkg3d_sm);
        }

        if(fftw3dFlag==0){
          sngl_upack_coef_sum(coeffReTemp,coeffImTemp,zfft,
                              cp_sclr_fft_pkg3d_sm);
        }
        else{
          sngl_upack_coef_sum_fftw3d(coeffReTemp,coeffImTemp,zfft,
                              cp_sclr_fft_pkg3d_sm);
        }

        ioff = (is-1)*numAtomFric*3*numCoeffUpTotal+
               i*3*numCoeffUpTotal+
               j*numCoeffUpTotal;

        for(iCoeff=0;iCoeff<numCoeff;iCoeff++){
          daHChiReUp[ioff+iCoeff] += coeffReTemp[iCoeff+1];
          daHChiImUp[ioff+iCoeff] += coeffImTemp[iCoeff+1];
        }//endfor iCoeff
      }//endif numStateProc
    }//endfor j
  }//endfor i

/*==============================================================*/
}/*end routine*/
/*==============================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void daVlocApply(GENERAL_DATA *general_data,CP *cp,CLASS *class,
                        double *ccreal,double *ccimag,
                        double *wfReal,int nstate)
/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"


/*==============================================================*/
}/*end routine*/
/*==============================================================*/




