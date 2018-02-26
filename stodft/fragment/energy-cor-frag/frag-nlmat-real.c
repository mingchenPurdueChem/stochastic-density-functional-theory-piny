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
void calcNonLocalMatrix(CP *cp, CP *cpMini, CLASS *classMini,
                        GENERAL_DATA *generalDataMini)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
/**************************************************************************/
/* This function calculate non-local pseudopotential matrix w.r.t. frag-  */
/* -ment MO, as well as the force component.                              */
/* The mojority part of the code is copied from control_cp_eext_recip     */
/* Right now only KB form is included.                                    */
/**************************************************************************/
/*======================================================================*/
/* Local Variable declaration                                           */
#include "../typ_defs/typ_mask.h"
  CLATOMS_INFO *clatoms_info = &(classMini->clatoms_info);
  CLATOMS_POS *clatoms_pos = &(classMini->clatoms_pos[1]);
  CPCOEFFS_INFO *cpcoeffs_info = &(cpMini->cpcoeffs_info);
  CPCOEFFS_POS *cpcoeffs_pos = &(cpMini->cpcoeffs_pos[1]);
  CPEWALD *cpewald = &(cpMini->cpewald);
  CPSCR *cpscr = &(cpMini->cpscr);
  CPOPTS *cpopts = &(cpMini->cpopts);
  PSEUDO *pseudo = &(cpMini->pseudo);
  EWD_SCR *ewd_scr = &(classMini->ewd_scr);
  ATOMMAPS *atommaps = &(classMini->atommaps);
  CELL *cell = &(generalDataMini->cell);
  EWALD *ewald = &(generalDataMini->ewald);
  PTENS *ptens = &(generalDataMini->ptens);
  COMMUNICATE *communicate = &(cp->communicate);
  FOR_SCR *for_scr = &(classMini->for_scr);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cpMini->cp_para_fft_pkg3d_lg);
  STAT_AVG *stat_avg = &(generalDataMini->stat_avg);
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  
  int i,j,k;
  int iState,jState;
  int iFrag = fragInfo->iFrag;
  int nstate_up = cpcoeffs_info->nstate_up_proc;
  int nstate_dn = cpcoeffs_info->nstate_dn_proc;
  int cp_lsda = cpopts->cp_lsda;
  int numAtomCalc = fragInfo->numAtomFragVnlCalc[iFrag];

  double *coefReUp = cpcoeffs_pos->cre_up;
  double *coefImUp = cpcoeffs_pos->cim_up;
  double *coefReDn = cpcpeffs_pos->cre_dn;
  double *coefImDn = cpcpeffs_pos->cim_dn;


/*======================================================================*/
/* I) Initialize matrix and force                                       */

  // initialize matrix and force
  for(iState=0;iState<nstate_up;iState++){
    for(jState=0;jState<nstate_up;jState++){
      fragInfo->vnlMatrixUp[iFrag][iState*nstate_up+jState] = 0.0;
    }
  }
  if(cp_lsda==1 && nstate_dn != 0){
    for(iState=0;iState<nstate_dn;iState++){
      for(jState=0;jState<nstate_up;jState++){
        fragInfo->vnlMatrixDn[iFrag][iState*nstate_dn+jState] = 0.0;
      }
    }
  }

  for(i=0;i<numAtomCalc;i++){
    for(iState=0;iState<nstate_up;iState++){
      for(jState=0;jState<nstate_up;jState++){
        fragInfo->vnlFxMatrixUp[iFrag][i*nstate_up*nstate_up+iState*nstate_up+jState] = 0.0;
        fragInfo->vnlFyMatrixUp[iFrag][i*nstate_up*nstate_up+iState*nstate_up+jState] = 0.0;
        fragInfo->vnlFzMatrixUp[iFrag][i*nstate_up*nstate_up+iState*nstate_up+jState] = 0.0;
      }
    }
    if(cp_lsda==1 && nstate_dn != 0){
      for(i=0;i<numAtomCalc;i++){
        for(iState=0;iState<nstate_dn;iState++){
          for(jState=0;jState<nstate_dn;jState++){
            fragInfo->vnlFxMatrixDn[iFrag][i*nstate_dn*nstate_dn+iState*nstate_dn+jState] = 0.0;
            fragInfo->vnlFyMatrixDn[iFrag][i*nstate_dn*nstate_dn+iState*nstate_dn+jState] = 0.0;
            fragInfo->vnlFzMatrixDn[iFrag][i*nstate_dn*nstate_dn+iState*nstate_dn+jState] = 0.0;
          }
        }
      }
    }
  }

  calcVnlMatrixSpinUp(cpMini,classMini,generalDataMini,coefReUp,coefImUp);
  if(cp_lsda==1 && nstate_dn != 0){
    calcVnlMatrixSpinDn(cpMini,classMini,generalDataMini,coefReDn,coefImDn);
  }

/*======================================================================*/
    }/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcVnlMatrixSpinUp(CP *cpMini, CLASS *classMini,GENERAL_DATA *generalDataMini,
			double *ccreal,double *ccimag)
/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

   int is,i,iupper;
   double tpi;
   double aka,akb,akc,xk,yk,zk,cfact;
   double eke;
   int ioff,ncoef1,ioff2;
   int iii,iis,nis;

   int nfft       = cp_sclr_fft_pkg3d_sm->nfft;
   int nfft2      = nfft/2;
   int ncoef      = cp_sclr_fft_pkg3d_sm->ncoef;
   int myid_state = communicate->myid_state;
   int np_states  = communicate->np_states;
   int pseudoRealFlag = cp->pseudo.pseudoReal.pseudoRealFlag;

/*            Local pointers                                       */

   int  *kastore_sm    =  cpewald->kastr_sm;
   int  *kbstore_sm    =  cpewald->kbstr_sm;
   int  *kcstore_sm    =  cpewald->kcstr_sm;

#define DEBUG_OFF
#ifdef DEBUG
      int icount;
      double c_g,g2,anorm,sum,vol,cre_now,cim_now;
      double dx,x_pos,y_pos,z_pos,phase_r,phase_i,arg;
      FILE *fp;
#endif

      double sum_check,sum_check_tmp;
      MPI_Comm comm_states = communicate->comm_states;

  int fftw3dFlag = cpewald->fftw3dFlag;
  int onebodyMatrixFlag = cpewald->onebodyMatrixFlag;
  double *keMatrix = cpewald->keMatrix;

/*=================================================================*/
/*  Find the upper state limit                                     */

  ncoef1 = ncoef - 1;
  iupper = nstate;
  if(nstate % 2 == 1){
     iupper = nstate - 1;
  }

  //debug
#ifdef REAL_PP_DEBUG  
  for(i=1;i<=ncoef;i++){
    fccreal[i] = 0.0;
    fccimag[i] = 0.0;
  }
#endif

/*=================================================================*/
/*  get the forces on the coefs of each state                      */

  for(is=1 ; is<= iupper; is+=2 ){

    ioff = (is-1)*ncoef;
    ioff2 = (is)*ncoef;

/*==========================================================================*/
/* 1) get the wave functions in real space two at a time                    */
/*   I) double pack the complex zfft array with two real wavefunctions      */

    if(fftw3dFlag==0){
      dble_pack_coef(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
                     zfft,cp_sclr_fft_pkg3d_sm);
    }
    else{
      dble_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
                     zfft,cp_sclr_fft_pkg3d_sm);
    }

/*--------------------------------------------------------------------------*/
/* II) fourier transform the wavefunctions to real space                    */
/*      convention exp(-igr)                                                */

    if(fftw3dFlag==0){
      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
    }
    else{
      para_fft_gen3d_fwd_to_r_fftw3d(zfft,cp_sclr_fft_pkg3d_sm);
    }

/*==========================================================================*/
/* 2) get v|psi> in g space and store it in zfft                            */
/*   I) get  v|psi> in real space                                           */

      memcpy(&zfft_tmp[1],&zfft[1],nfft*sizeof(double));
      cp_vpsi(zfft,v_ks,nfft);
      cp->pseudo.pseudoReal.energyCalcFlag = 1;
      controlEnergyNlppReal(cp,class,general_data,zfft_tmp,zfft,1);

    }/*endfor is */

/*==========================================================================*/
/*==========================================================================*/
/* 4) if there is an odd number of states, go through                       */
/*      the same procedure using sng_packs                                  */

  if(nstate % 2 != 0){
     is = nstate;
     ioff = (is -1)*ncoef;

/*--------------------------------------------------------------------------*/
/*   I) sngl pack                                                           */

     if(fftw3dFlag==0){
       sngl_pack_coef(&ccreal[ioff],&ccimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);
     }
     else{
       sngl_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);
     }
/*--------------------------------------------------------------------------*/
/* II) fourier transform the wavefunctions to real space                    */
/*      convention exp(-igr)                                                */

      if(fftw3dFlag==0){
        para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);
      }
      else{
        para_fft_gen3d_fwd_to_r_fftw3d(zfft,cp_sclr_fft_pkg3d_sm);
      }

/*==========================================================================*/
/* 5) get v|psi> in g space and store it in zfft                            */
/*   I) get  v|psi> in real space                                           */

      memcpy(&zfft_tmp[1],&zfft[1],nfft*sizeof(double));
      cp_vpsi(zfft,v_ks,nfft);
      cp->pseudo.pseudoReal.energyCalcFlag = 1;
      controlEnergyNlppReal(cp,class,general_data,zfft_tmp,zfft,0);

    }//endif nstate

/*======================================================================*/
    }/*end routine*/
/*======================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void controlEnergyNlppReal(CP *cp,CLASS *class,GENERAL_DATA *generalData,
	       double *zfft_tmp,double *zfft,int flag)
/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
  PSEUDO *pseudo = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);

  PARA_FFT_PKG3D *cpParaFftPkg3dLgBigBox = &(cp->cp_para_fft_pkg3d_lg);
  int nfft = cpParaFftPkg3dLgBigBox->nfft;
  int numGrid = nfft/2;
  int iGrid;
  int forceCalcFlag = pseudoReal->forceCalcFlag;
  double *wfReal,*wfForceReal;

  wfReal = (double*)cmalloc(numGrid*sizeof(double));
  wfForceReal = (double*)cmalloc(numGrid*sizeof(double));
  //printf("forceCalcFlag %i\n",forceCalcFlag);

  for(iGrid=0;iGrid<numGrid;iGrid++){
    wfReal[iGrid] = zfft_tmp[iGrid*2+1];
    wfForceReal[iGrid] = 0.0;
  }
  nlppKBRealEnergy(cp,class,generalData,wfReal,wfForceReal);
  if(forceCalcFlag==1)nlppKBRealEnergyForce(cp,class,generalData,wfReal);
  for(iGrid=0;iGrid<numGrid;iGrid++){
#ifdef REAL_PP_DEBUG  
    zfft[iGrid*2+1] = wfForceReal[iGrid];
    //printf("forceeeeeeeee grid 1 %.16lg\n",wfForceReal[iGrid]);
#else
    zfft[iGrid*2+1] += wfForceReal[iGrid];
#endif
  }
  if(flag==1){//double
    for(iGrid=0;iGrid<numGrid;iGrid++){
      wfReal[iGrid] = zfft_tmp[iGrid*2+2];
      wfForceReal[iGrid] = 0.0;
    }
    nlppKBRealEnergy(cp,class,generalData,wfReal,wfForceReal);
    if(forceCalcFlag==1)nlppKBRealEnergyForce(cp,class,generalData,wfReal);
    for(iGrid=0;iGrid<numGrid;iGrid++){
#ifdef REAL_PP_DEBUG  
      zfft[iGrid*2+2] = wfForceReal[iGrid];
      //printf("forceeeeeeeee grid 2 %.16lg\n",wfForceReal[iGrid]);
#else
      zfft[iGrid*2+2] += wfForceReal[iGrid];
#endif
    }
  }//endif
  free(wfReal);
  free(wfForceReal);

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void nlppKBRealEnergy(CP *cp,CLASS *class,GENERAL_DATA *generalData,
	      double *wfReal,double *forceRealNlpp)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Real space nlpp, only used in filtering		 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  PSEUDO *pseudo = &(cp->pseudo);
  PARA_FFT_PKG3D *cpParaFftPkg3dLgBigBox = &(cp->cp_para_fft_pkg3d_lg);
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
  int countRad,countNlppRe,countNlppIm;
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
  double *vnlPhiAtomGridRe = pseudoReal->vnlPhiAtomGridRe;
  double *vnlPhiAtomGridIm = pseudoReal->vnlPhiAtomGridIm;

  double **dotReAll = pseudoReal->dotReAll;
  double **dotImAll = pseudoReal->dotImAll;

/*======================================================================*/
/* I) Allocate local memory                                             */

  numGridMax = numGridNlppMap[0];
  for(iPart=0;iPart<numAtom;iPart++){
    if(numGridNlppMap[iPart]>numGridMax)numGridMax = numGridNlppMap[iPart];
  }
  forceTemp = (double*)cmalloc(numGridMax*sizeof(double));
  wfNbhd = (double*)cmalloc(numGridMax*sizeof(double));
  vol = getdeth(hmat);
  volInv = 1.0/vol;
  volElem = vol/numGridTot;

/*======================================================================*/
/* I) Allocate local memory                                             */

  // initialize matrix and force
  for(iState=0;iState<nstate_up;iState++){
    for(jState=0;jState<nstate_up;jState++){
      fragInfo->vnlMatrixUp[iFrag][iState*nstate_up+jState] = 0.0;
    }
  }
  if(cp_lsda==1 && nstate_dn != 0){
    for(iState=0;iState<nstate_dn;iState++){
      for(jState=0;jState<nstate_up;jState++){
        fragInfo->vnlMatrixDn[iFrag][iState*nstate_dn+jState] = 0.0;
      }
    }
  }

  for(i=0;i<numAtomCalc;i++){
    for(iState=0;iState<nstate_up;iState++){
      for(jState=0;jState<nstate_up;jState++){
        fragInfo->vnlFxMatrixUp[iFrag][i*nstate_up*nstate_up+iState*nstate_up+jState] = 0.0;
        fragInfo->vnlFyMatrixUp[iFrag][i*nstate_up*nstate_up+iState*nstate_up+jState] = 0.0;
        fragInfo->vnlFzMatrixUp[iFrag][i*nstate_up*nstate_up+iState*nstate_up+jState] = 0.0;
      }
    }
    if(cp_lsda==1 && nstate_dn != 0){
      for(i=0;i<numAtomCalc;i++){
        for(iState=0;iState<nstate_dn;iState++){
          for(jState=0;jState<nstate_dn;jState++){
            fragInfo->vnlFxMatrixDn[iFrag][i*nstate_dn*nstate_dn+iState*nstate_dn+jState] = 0.0;
            fragInfo->vnlFyMatrixDn[iFrag][i*nstate_dn*nstate_dn+iState*nstate_dn+jState] = 0.0;
            fragInfo->vnlFzMatrixDn[iFrag][i*nstate_dn*nstate_dn+iState*nstate_dn+jState] = 0.0;
          }
        }
      }
    }
  }


/*======================================================================*/
/* II) Loop over iPart/irad/ipart/m                                 */

  gridShiftNowRe = 0;
  gridShiftNowIm = 0;
  for(iAtom=0;iAtom<numAtom;iAtom++){
    atomType = iAtomAtomType[iAtom+1]-1;
    numGrid = numGridNlppMap[iAtom];
    countRad = 0;
    countNlppRe = 0;
    countNlppIm = 0;
    /* cpy the wave function */
    if(numGrid>0){ //if numGrid=0, only local pp will be calculated
      for(iGrid=0;iGrid<numGrid;iGrid++){
    gridIndex = gridNlppMap[iAtom][iGrid];
    wfNbhd[iGrid] = wfReal[gridIndex];
    forceTemp[iGrid] = 0.0;
      }
      for(l=0;l<atomLMax[atomType];l++){
    for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
      radIndex = atomRadMap[atomType][countRad+iRad];
      energyl = 0.0;
      for(m=0;m<=l;m++){
        //calcDotNlpp(wfNbhd,radFun,&ylm[ylmShift],&dotRe,&dotIm);
        if(m!=0){
          dotRe = ddotBlasWrapper(numGrid,wfNbhd,1,
	      &vnlPhiAtomGridRe[gridShiftNowRe],1)*volElem;
          dotIm = ddotBlasWrapper(numGrid,wfNbhd,1,
                      &vnlPhiAtomGridIm[gridShiftNowIm],1)*volElem;
          dotReAll[iAtom][countNlppRe+m] = dotRe;
          dotImAll[iAtom][countNlppIm+m-1] = dotIm;
          //printf("dottttttttt %i %i %lg %i %lg\n",
          //       iAtom,countNlppRe+m,dotRe,countNlppIm+m-1,dotIm);
          energyl += 2.0*(dotRe*dotRe+dotIm*dotIm)*vpsNormList[radIndex]*volInv;
          //printf("energy %lg\n",energyl);
          //printf("m %i dotRe %lg dotIm %lg vpsNormList[radIndex] %lg\n",m,dotRe,dotIm,vpsNormList[radIndex]);
          dotRe *= 2.0*vpsNormList[radIndex];
          dotIm *= 2.0*vpsNormList[radIndex];
          daxpyBlasWrapper(numGrid,dotRe,&vnlPhiAtomGridRe[gridShiftNowRe],1,
	           forceTemp,1);
          daxpyBlasWrapper(numGrid,dotIm,&vnlPhiAtomGridIm[gridShiftNowIm],1,
                               forceTemp,1);
          /*
          for(iGrid=0;iGrid<numGrid;iGrid++){
	forceTemp[iGrid] += 2.0*radFun[iGrid]*(ylm[ylmShift+iGrid]*dotRe-
		    ylm[ylmShift+numGrid+iGrid]*dotIm)*vpsNormList[radIndex];
          }//endfor iGrid
          */
          gridShiftNowRe += numGrid;
          gridShiftNowIm += numGrid;
        }
        else{
              dotRe = ddotBlasWrapper(numGrid,wfNbhd,1,
                      &vnlPhiAtomGridRe[gridShiftNowRe],1)*volElem;
          dotReAll[iAtom][countNlppRe] = dotRe;
              //printf("dottttttttt %i %i %lg\n",
              //       iAtom,countNlppRe,dotRe);

          //printf("numGrid %i\n",numGrid);
          /*
          for(iGrid=0;iGrid<numGrid;iGrid++){
	printf("griddddddddddd %lg %lg\n",wfNbhd[iGrid],vnlPhiAtomGridRe[gridShiftNowRe+iGrid]);
          }
          */
          
          energyl += dotRe*dotRe*vpsNormList[radIndex]*volInv;
          //printf("energy %lg\n",energyl);
          /*
          for(iGrid=0;iGrid<numGrid;iGrid++){
	forceTemp[iGrid] += radFun[iGrid]*ylm[ylmShift+iGrid]*dotRe*vpsNormList[radIndex];
          }//endfor iGrid
          */
              //printf("m %i dotRe %lg volElem %lg vpsNormList[radIndex] %lg\n",m,dotRe*volElem,volElem,vpsNormList[radIndex]);

          dotRe *= vpsNormList[radIndex];
              //printf("m %i dotRe %lg volElem %lg vpsNormList[radIndex] %lg\n",m,dotRe,volElem,vpsNormList[radIndex]);

          daxpyBlasWrapper(numGrid,dotRe,&vnlPhiAtomGridRe[gridShiftNowRe],1,
	           forceTemp,1);
          gridShiftNowRe += numGrid;
          //printf("forcetemppppppp %lg\n",forceTemp[0]);
        }//endif m
      }//endfor m
      //printf("energyl %lg\n",energyl);
      energy += energyl;
          countNlppRe += l+1;
          countNlppIm += l;
    }//endfor iRad
        countRad += atomLRadNum[atomType][l];
      }//endfor l
      for(iGrid=0;iGrid<numGrid;iGrid++){
    gridIndex = gridNlppMap[iAtom][iGrid];
    forceRealNlpp[gridIndex] += forceTemp[iGrid];
      }//endfor iGrid
    }//endif numGrid
  }//endfor iAtom

  if(energyCalcFlag==1){
    stat_avg->cp_enl += energy;
  }

/*======================================================================*/
/* III) free local memory                                               */

  free(wfNbhd);
  free(forceTemp);

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/
