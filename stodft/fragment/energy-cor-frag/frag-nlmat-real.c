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
void calcRealNonLocalMatrix(CP *cp, CP *cpMini, CLASS *classMini,
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
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
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
  
  int i,j,k,idot;
  int iState,jState;
  int iFrag = fragInfo->iFrag;
  int nstate_up = cpcoeffs_info->nstate_up_proc;
  int nstate_dn = cpcoeffs_info->nstate_dn_proc;
  int cp_lsda = cpopts->cp_lsda;
  int numAtomCalc = fragInfo->numAtomFragVnlCalc[iFrag];
  int numNlppAll = 0;
  int numNlppAtom;
  int numAtomTot = clatoms_info->natm_tot;

  int *nlppAtomStartIndex;

  double *coefReUp = cpcoeffs_pos->cre_up;
  double *coefImUp = cpcoeffs_pos->cim_up;
  double *coefReDn = cpcpeffs_pos->cre_dn;
  double *coefImDn = cpcpeffs_pos->cim_dn;

  double *dotReAllStatesUp;
  double *dotImAllStatesUp;
  double *dotReAllDxStatesUp;
  double *dotImAllDxStatesUp;
  double *dotReAllDyStatesUp;
  double *dotImAllDyStatesUp;
  double *dotReAllDzStatesUp;
  double *dotImAllDzStatesUp;
  double *dotReAllStatesDn;
  double *dotImAllStatesDn;
  double *dotReAllDxStatesDn;
  double *dotImAllDxStatesDn;
  double *dotReAllDyStatesDn;
  double *dotImAllDyStatesDn;
  double *dotReAllDzStatesDn;
  double *dotImAllDzStatesDn;

  pseudoReal->nlppAtomStartIndex = (int*)cmalloc(numAtomTot*sizeof(int));
  nlppAtomStartIndex = pseudoReal->nlppAtomStartIndex;
  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    atomType = iAtomAtomType[iAtom+1]-1;   
    numNlppAtom = 0;
    nlppAtomStartIndex[iAtom] = numNlppAll;
    for(l=0;l<atomLMax[atomType];l++){
      numNlppAtom += atomLRadNum[atomType][l]*(l+1);
    }
    numNlppAll += numNlppAtom;
  }

  dotReAllStatesUp = (double*)cmalloc(numNlppAll*nstate_up*sizeof(double));
  dotImAllStatesUp = (double*)cmalloc(numNlppAll*nstate_up*sizeof(double));
  dotReAllDxStatesUp = (double*)cmalloc(numNlppAll*nstate_up*sizeof(double));
  dotImAllDxStatesUp = (double*)cmalloc(numNlppAll*nstate_up*sizeof(double));
  dotReAllDyStatesUp = (double*)cmalloc(numNlppAll*nstate_up*sizeof(double));
  dotImAllDyStatesUp = (double*)cmalloc(numNlppAll*nstate_up*sizeof(double));
  dotReAllDzStatesUp = (double*)cmalloc(numNlppAll*nstate_up*sizeof(double));
  dotImAllDzStatesUp = (double*)cmalloc(numNlppAll*nstate_up*sizeof(double));

  for(idot=0;idot<numNlppAll*nstate_up;idot++){
    dotReAllStatesUp[idot] = 0.0;
    dotImAllStatesUp[idot] = 0.0;
    dotReAllDxStatesUp[idot] = 0.0;
    dotImAllDxStatesUp[idot] = 0.0;
    dotReAllDyStatesUp[idot] = 0.0;
    dotImAllDyStatesUp[idot] = 0.0;
    dotReAllDzStatesUp[idot] = 0.0;
    dotImAllDzStatesUp[idot] = 0.0;
  }
  if(cp_lsda==1&&nstate_dn>0){
    dotReAllStatesDn = (double*)cmalloc(numNlppAll*nstate_dn*sizeof(double));
    dotImAllStatesDn = (double*)cmalloc(numNlppAll*nstate_dn*sizeof(double));
    dotReAllDxStatesDn = (double*)cmalloc(numNlppAll*nstate_dn*sizeof(double));
    dotImAllDxStatesDn = (double*)cmalloc(numNlppAll*nstate_dn*sizeof(double));
    dotReAllDyStatesDn = (double*)cmalloc(numNlppAll*nstate_dn*sizeof(double)); 
    dotImAllDyStatesDn = (double*)cmalloc(numNlppAll*nstate_dn*sizeof(double));
    dotReAllDzStatesDn = (double*)cmalloc(numNlppAll*nstate_dn*sizeof(double));
    dotImAllDzStatesDn = (double*)cmalloc(numNlppAll*nstate_dn*sizeof(double));
    for(idot=0;idot<numNlppAll*nstate_up;idot++){
      dotReAllStatesDn[idot] = 0.0;
      dotImAllStatesDn[idot] = 0.0;
      dotReAllDxStatesDn[idot] = 0.0;
      dotImAllDxStatesDn[idot] = 0.0;
      dotReAllDyStatesDn[idot] = 0.0;
      dotImAllDyStatesDn[idot] = 0.0;
      dotReAllDzStatesDn[idot] = 0.0;
      dotImAllDzStatesDn[idot] = 0.0;
    }
  }

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

/*======================================================================*/
/* II) Calculate Tensor of dot product between MO and			*/
/*     pseudo wf(derivative) */


  calcVnlRealDot(cpMini,classMini,generalDataMini,coefReUp,coefImUp,
		 dotReAllStatesUp,dotImAllStatesUp,dotReAllDxStatesUp,
		 dotImAllDxStatesUp,dotReAllDyStatesUp,dotImAllDyStatesUp,
		 dotReAllDzStatesUp,dotImAllDzStatesUp);
  if(cp_lsda==1&&nstate_dn!=0){
    calcVnlRealDot(cpMini,classMini,generalDataMini,coefReDn,coefImDn,
		   dotReAllStatesDn,dotImAllStatesDn,dotReAllDxStatesDn,
		   dotImAllDxStatesDn,dotReAllDyStatesDn,dotImAllDyStatesDn,
		   dotReAllDzStatesDn,dotImAllDzStatesDn);

  }

/*======================================================================*/
/* III) Calculate vnl matrix						*/

  calcMatrixFromDot(cp,cpMini,classMini,coefReUp,coefImUp,
		    dotReAllStatesUp,dotImAllStatesUp,dotReAllDxStatesUp,
                    dotImAllDxStatesUp,dotReAllDyStatesUp,dotImAllDyStatesUp,
                    dotReAllDzStatesUp,dotImAllDzStatesUp,vnlMatrixUp,
		    vnlFxMatrixUp,vnlFyMatrixUp,vnlFzMatrixUp);
  if(cp_lsda==1&&nstate_dn!=0){
    calcMatrixFromDot(cp,cpMini,classMini,coefReDn,coefImDn,
		      dotReAllStatesDn,dotImAllStatesDn,dotReAllDxStatesDn,
		      dotImAllDxStatesDn,dotReAllDyStatesDn,dotImAllDyStatesDn,
		      dotReAllDzStatesDn,dotImAllDzStatesDn,vnlMatrixDn,
		      vnlFxMatrixDn,vnlFyMatrixDn,vnlFzMatrixDn);
  }


/*======================================================================*/
    }/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcVnlRealDot(CP *cpMini, CLASS *classMini,GENERAL_DATA *generalDataMini,
		    double *ccreal,double *ccimag,double *dotReAllStates,
		    double *dotImAllStates,double *dotReAllDxStates,
		    double *dotImAllDxStates,double *dotReAllDyStates,
		    double *dotImAllDyStates,double *dotReAllDzStates,
		    double *dotImAllDzStates)
/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  CPEWALD *cpewald = &(cpMini->cpewald);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm = &(cpMini->cp_sclr_fft_pkg3d_sm);
  PSEUDO *pseudo = &(cpMini->pseudo);
  PSEUFO_REAL *pseudoReal = &(cpMini->pseudoReal);
  
  int is,i,iupper;
  int ioff,ncoef1,ioff2;
  int iii,iis,nis;
  int nfft       = cp_sclr_fft_pkg3d_sm->nfft;
  int nfft2      = nfft/2;
  int ncoef      = cp_sclr_fft_pkg3d_sm->ncoef;
  int myid_state = communicate->myid_state;
  int np_states  = communicate->np_states;
  int pseudoRealFlag = cp->pseudo.pseudoReal.pseudoRealFlag;
  int fftw3dFlag = cpewald->fftw3dFlag;
  int numNlppAll = pseudoReal->numNlppAll;

  MPI_Comm comm_states = communicate->comm_states;

  int *kastore_sm = cpewald->kastr_sm;
  int *kbstore_sm = cpewald->kbstr_sm;
  int *kcstore_sm = cpewald->kcstr_sm;

  double tpi;
  double aka,akb,akc,xk,yk,zk,cfact;
  double eke;
  double sum_check,sum_check_tmp;

  double *wfReal = (double*)cmalloc(numGrid*sizeof(double));

/*=================================================================*/
/*  Find the upper state limit                                     */

  ncoef1 = ncoef-1;
  iupper = nstate;
  if(nstate%2==1)iupper = nstate-1;

/*=================================================================*/
/*  get the forces on the coefs of each state                      */

  for(is=1;is<=iupper;is+=2){
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

    for(iGrid=0;iGrid<nfft2;iGrid++){
      wfReal[iGrid] = zfft[iGrid*2+1];
    }
    calcVnlRealDotState(cpMini,classMini,generalDataMini,wfReal,
			&dotReAllStates[(is-1)*numNlppAll],
			&dotImAllStates[(is-1)*numNlppAll],
			&dotReAllDxStates[(is-1)*numNlppAll],
                        &dotImAllDxStates[(is-1)*numNlppAll],
                        &dotReAllDyStates[(is-1)*numNlppAll],
                        &dotImAllDyStates[(is-1)*numNlppAll],
                        &dotReAllDzStates[(is-1)*numNlppAll],
                        &dotImAllDzStates[(is-1)*numNlppAll]);
    for(iGrid=0;iGrid<nfft2;iGrid++){
      wfReal[iGrid] = zfft[iGrid*2+2];
    }
    calcVnlRealDotState(cpMini,classMini,generalDataMini,wfReal,
                        &dotReAllStates[is*numNlppAll],
                        &dotImAllStates[is*numNlppAll],
                        &dotReAllDxStates[is*numNlppAll],
                        &dotImAllDxStates[is*numNlppAll],
                        &dotReAllDyStates[is*numNlppAll],
                        &dotImAllDyStates[is*numNlppAll],
                        &dotReAllDzStates[is*numNlppAll],
                        &dotImAllDzStates[is*numNlppAll]);
  }//endfor is

/*==========================================================================*/
/*==========================================================================*/
/* 4) if there is an odd number of states, go through                       */
/*      the same procedure using sng_packs                                  */

  if(nstate%2!=0){
    is = nstate;
    ioff = (is-1)*ncoef;
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

     for(iGrid=0;iGrid<nfft2;iGrid++){
       wfReal[iGrid] = zfft[iGrid*2+1];
     }
     calcVnlRealDotState(cpMini,classMini,generalDataMini,wfReal,
                         &dotReAllStates[(is-1)*numNlppAll],
                         &dotImAllStates[(is-1)*numNlppAll],
                         &dotReAllDxStates[(is-1)*numNlppAll],
                         &dotImAllDxStates[(is-1)*numNlppAll],
                         &dotReAllDyStates[(is-1)*numNlppAll],
                         &dotImAllDyStates[(is-1)*numNlppAll],
                         &dotReAllDzStates[(is-1)*numNlppAll],
                         &dotImAllDzStates[(is-1)*numNlppAll]);
    }//endif nstate

/*==========================================================================*/
/* 6) Free local pointers						    */

  free(&wfReal[0]);

/*======================================================================*/
    }/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcVnlRealDotState(CP *cp,CLASS *class,GENERAL_DATA *generalData,
	      double *wfReal,double *dotRe,double *dotIm,double*dotReDx,
	      double *dotImDx,double *dotReDy,double *dotImDy,double *dotReDz,
	      double *dotImDz)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Calculate dot product for one state					 */
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
  int atomIndSt;

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
  wfNbhd = (double*)cmalloc(numGridMax*sizeof(double));
  vol = getdeth(hmat);
  volInv = 1.0/vol;
  volElem = vol/numGridTot;

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
    atomIndSt = nlppAtomStartIndex[iAtom];
    /* cpy the wave function */
    if(numGrid>0){ //if numGrid=0, only local pp will be calculated
      for(iGrid=0;iGrid<numGrid;iGrid++){
	gridIndex = gridNlppMap[iAtom][iGrid];
	wfNbhd[iGrid] = wfReal[gridIndex];
      }
      for(l=0;l<atomLMax[atomType];l++){
	for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
	  radIndex = atomRadMap[atomType][countRad+iRad];
	  for(m=0;m<=l;m++){
	    if(m!=0){
	      dotRe[atomIndSt+countNlppRe+m] 
				= ddotBlasWrapper(numGrid,wfNbhd,1,
				  &vnlPhiAtomGridRe[gridShiftNowRe],1)*volElem;
	      dotIm[atomIndSt+countNlppIm+m-1] 
				= ddotBlasWrapper(numGrid,wfNbhd,1,
				  &vnlPhiAtomGridIm[gridShiftNowIm],1)*volElem;

              dotReDx[atomIndSt+countNlppRe+m] 
				= ddotBlasWrapper(numGrid,wfNbhd,1,
                                  &vnlPhiDxAtomGridRe[gridShiftNowRe],1)*volElem;
              dotImDx[atomIndSt+countNlppIm+m-1] 
				= ddotBlasWrapper(numGrid,wfNbhd,1,
                                  &vnlPhiDxAtomGridIm[gridShiftNowIm],1)*volElem;
              dotReDy[atomIndSt+countNlppRe+m] 
				= ddotBlasWrapper(numGrid,wfNbhd,1,
                                  &vnlPhiDyAtomGridRe[gridShiftNowRe],1)*volElem;
              dotImDy[atomIndSt+countNlppIm+m-1] 
				= ddotBlasWrapper(numGrid,wfNbhd,1,
                                  &vnlPhiDyAtomGridIm[gridShiftNowIm],1)*volElem;
              dotReDz[atomIndSt+countNlppRe+m] 
				= ddotBlasWrapper(numGrid,wfNbhd,1,
                                  &vnlPhiDzAtomGridRe[gridShiftNowRe],1)*volElem;
              dotImDz[atomIndSt+countNlppIm+m-1] 
				= ddotBlasWrapper(numGrid,wfNbhd,1,
                                  &vnlPhiDzAtomGridIm[gridShiftNowIm],1)*volElem;
	      gridShiftNowRe += numGrid;
	      gridShiftNowIm += numGrid;
	    }
	    else{
  	      dotRe[atomIndSt+countNlppRe] 
				= ddotBlasWrapper(numGrid,wfNbhd,1,
			          &vnlPhiAtomGridRe[gridShiftNowRe],1)*volElem;
              dotReDx[atomIndSt+countNlppRe] 
				= ddotBlasWrapper(numGrid,wfNbhd,1,
                                  &vnlPhiDxAtomGridRe[gridShiftNowRe],1)*volElem;
              dotReDy[atomIndSt+countNlppRe] 
				= ddotBlasWrapper(numGrid,wfNbhd,1,
                                  &vnlPhiDyAtomGridRe[gridShiftNowRe],1)*volElem;
              dotReDz[atomIndSt+countNlppRe] 
				= ddotBlasWrapper(numGrid,wfNbhd,1,
                                  &vnlPhiDzAtomGridRe[gridShiftNowRe],1)*volElem;
	      gridShiftNowRe += numGrid;
	    }//endif m
	  }//endfor m
	  countNlppRe += l+1;
	  countNlppIm += l;
        }//endfor iRad
        countRad += atomLRadNum[atomType][l];
      }//endfor l
    }//endif numGrid
  }//endfor iAtom

/*======================================================================*/
/* III) free local memory                                               */

  free(wfNbhd);

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void calcMatrixFromDot(CP *cp, CP *cpMini,double *dotRe,double *dotIm,double*dotReDx,
              double *dotImDx,double *dotReDy,double *dotImDy,double *dotReDz,
              double *dotImDz,double *vnlMatrix,double *vnlFxMatrix,
	      double *vnlFyMatrix,double *vnlFzMatrix,int numStates)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
/**************************************************************************/
/* This function calculate real sapce non-local pseudopotential matrix    */
/* w.r.t. fragment MO, as well as the force component, using dot product  */
/* between wave function and pseudo wave function calculated in		  */
/* calcVnlRealDotState.							  */
/**************************************************************************/
/*======================================================================*/
/* Local Variable declaration                                           */
 
  CPOPTS *cpopts = &(cpMini->cpopts);
  PSEUDO *pseudo = &(cpMini->pseudo);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CELL *cell = &(generalDataMini->cell);
  CLATOMS_POS *clatoms_pos = &(classMini->clatoms_pos[1]);
  CLATOMS_INFO *clatoms_info = &(classMini->clatoms_info);
  ATOMMAPS *atommaps = &(classMini->atommaps);
  STAT_AVG *stat_avg = &(generalDataMini->stat_avg);
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  
  int cpLsda = cpopts->cp_lsda;
  int iState,jState;
  int iFrag = fragInfo->iFrag;
  int numAtomFragVnlCalc = fragInfo->numAtomFragVnlCalc[iFrag];
  int numAtomTot = clatoms_info->natm_tot;
  int numGrid,atomIndSt,atomType,countRad,countNlppRe,countNlppIm;
  int numNlppAll = pseudoReal->numNlppAll;
  int countAtom;

  int *atomFragVnlCalcMap = fragInfo->atomFragVnlCalcMap[iFrag];
  int *numGridNlppMap = *pseudoReal->numGridNlppMap;
  int *nlppAtomStartIndex = pseudoReal->nlppAtomStartIndex;
  int *atomFragVnlCalcMapInv = fragInfo->atomFragVnlCalcMapInv[iFrag];

  double iTempRe,iTempIm,jTempRe,jTempIm,iTempDRe,iTempDIm,jTempDRe,jTempDIm;
  double *Fx = clatoms_pos->fx;
  double *Fy = clatoms_pos->fy;
  double *Fz = clatoms_pos->fz;

  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    Fx[iAtom+1] = 0.0;
    Fy[iAtom+1] = 0.0;
    Fz[iAtom+1] = 0.0;
  }

  for(iState=0;iState<numStates;iState++){
    for(jState=0;jState<numStates;jState++){
      vnlMatrix[iState*numStates+jState] = 0.0;
      for(iAtom=0;iAtom<numAtomFragVnlCalc;iAtom++){
	vnlFxMatrix[iAtom*numStates*numStates+iState*numStates+jState] = 0.0;
        vnlFyMatrix[iAtom*numStates*numStates+iState*numStates+jState] = 0.0;
        vnlFzMatrix[iAtom*numStates*numStates+iState*numStates+jState] = 0.0;
      }//endfor iAtom
    }//endfor jState
  }//endfor iState

  for(iState=0;iState<numStates;iState++){
    for(jState=iState;jState<numStates;jState++){
      vnlMatrix[iState*numStates+jState] = 0.0;
      countAtom = 0;
      for(iAtom=0;iAtom<numAtomTot;iAtom++){
	//atomIndex = atomFragVnlCalcMap[iAtom];
	numGrid = numGridNlppMap[iAtom];
	atomIndSt = nlppAtomStartIndex[iAtom];
	atomType = iAtomAtomType[iAtom+1]-1;
	countRad = 0;
	countNlppRe = 0;
	countNlppIm = 0;
	if(numGrid>0){
	  for(l=0;l<atomLMax[atomType];l++){
	    for(iRad=0;iRad<atomLRadNum[atomType][l];iRad++){
	      radIndex = atomRadMap[atomType][countRad+iRad];
	      if(atomFragVnlCalcMapInv[iAtom]!=-1){ //I can only put if here since I dont wanna skip iRad/l loop
		atomIndex = atomFragVnlCalcMapInv[iAtom];
		for(m=0;m<=l;m++){
		  if(m!=0){
		    iTempRe = dotRe[iState*numNlppAll+atomIndSt+countNlppRe+m];
		    iTempIm = dotIm[iState*numNlppAll+atomIndSt+countNlppRe+m-1];
		    jTempRe = dotRe[jState*numNlppAll+atomIndSt+countNlppRe+m];
		    jTempIm = dotIm[jState*numNlppAll+atomIndSt+countNlppRe+m-1];
		    vnlMatrix[iState*numStates+jState] += 2.0*(iTempRe*jTempRe+iTempIm*jTempIm)*vpsNormList[radIndex]*volInv;
		    iTempDRe = dotReDx[iState*numNlppAll+atomIndSt+countNlppRe+m];
		    iTempDIm = dotImDx[iState*numNlppAll+atomIndSt+countNlppRe+m-1];
		    jTempDRe = dotReDx[jState*numNlppAll+atomIndSt+countNlppRe+m];
                    jTempDIm = dotImDx[jState*numNlppAll+atomIndSt+countNlppRe+m-1];
		    vnlFxMatrix[countAtom*numStates*numStates+iState*numStates+jState] 
			     += 2.0*(iTempDRe*jTempRe+iTempDIm*jTempIm+
			        iTempRe*jTempDRe+iTempIm*jTempDIm)*
				vpsNormList[radIndex]*volInv;
                    iTempDRe = dotReDy[iState*numNlppAll+atomIndSt+countNlppRe+m];
                    iTempDIm = dotImDy[iState*numNlppAll+atomIndSt+countNlppRe+m-1];
                    jTempDRe = dotReDy[jState*numNlppAll+atomIndSt+countNlppRe+m];
                    jTempDIm = dotImDy[jState*numNlppAll+atomIndSt+countNlppRe+m-1];
                    vnlFyMatrix[countAtom*numStates*numStates+iState*numStates+jState]
                             += 2.0*(iTempDRe*jTempRe+iTempDIm*jTempIm+
                                iTempRe*jTempDRe+iTempIm*jTempDIm)*
                                vpsNormList[radIndex]*volInv;
                    iTempDRe = dotReDz[iState*numNlppAll+atomIndSt+countNlppRe+m];
                    iTempDIm = dotImDz[iState*numNlppAll+atomIndSt+countNlppRe+m-1];
                    jTempDRe = dotReDz[jState*numNlppAll+atomIndSt+countNlppRe+m];
                    jTempDIm = dotImDz[jState*numNlppAll+atomIndSt+countNlppRe+m-1];
                    vnlFzMatrix[countAtom*numStates*numStates+iState*numStates+jState]
                             += 2.0*(iTempDRe*jTempRe+iTempDIm*jTempIm+
                                iTempRe*jTempDRe+iTempIm*jTempDIm)*
                                vpsNormList[radIndex]*volInv;
		  }
		  else{
		    iTempRe = dotRe[iState*numNlppAll+atomIndSt+countNlppRe];
		    jTempRe = dotRe[jState*numNlppAll+atomIndSt+countNlppRe];
		    vnlMatrix[iState*numStates+jState] += iTempRe*jTempRe*vpsNormList[radIndex]*volInv; 
		    iTempDRe = dotReDx[iState*numNlppAll+atomIndSt+countNlppRe];
		    jTempDRe = dotReDx[iState*numNlppAll+atomIndSt+countNlppRe];
		    vnlFxMatrix[countAtom*numStates*numStates+iState*numStates+jState]
			     += (iTempDRe*jTempRe+iTempRe*jTempDRe)*
				vpsNormList[radIndex]*volInv;
                    iTempDRe = dotReDy[iState*numNlppAll+atomIndSt+countNlppRe];
                    jTempDRe = dotReDy[iState*numNlppAll+atomIndSt+countNlppRe];
                    vnlFyMatrix[countAtom*numStates*numStates+iState*numStates+jState]
                             += (iTempDRe*jTempRe+iTempRe*jTempDRe)*
                                vpsNormList[radIndex]*volInv;
                    iTempDRe = dotReDz[iState*numNlppAll+atomIndSt+countNlppRe];
                    jTempDRe = dotReDz[iState*numNlppAll+atomIndSt+countNlppRe];
                    vnlFzMatrix[countAtom*numStates*numStates+iState*numStates+jState]
                             += (iTempDRe*jTempRe+iTempRe*jTempDRe)*
                                vpsNormList[radIndex]*volInv;
		  }//endif m
		}//endfor m
		countAtom += 1;
	      }//endif atomFragVnlCalcMapInv
	      countNlppRe += l+1;
	      countNlppIm += l;
	    }//endfor iRad
	    countRad += atomLRadNum[atomType][l];
	  }//endfor l
	}//endif numGrid
      }//endfor iAtom
    }//endfor jState
    stat_avg->cp_enl += vnlMatrix[iState*numStates+iState];
    countAtom = 0;
    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      if(atomFragVnlCalcMapInv[iAtom]!=-1){
	Fx[countAtom+1] += vnlFxMatrix[countAtom*numStates*numStates+iState*numStates+iState];
        Fy[countAtom+1] += vnlFyMatrix[countAtom*numStates*numStates+iState*numStates+iState];
        Fz[countAtom+1] += vnlFzMatrix[countAtom*numStates*numStates+iState*numStates+iState];
	countAtom += 1;
      }
    }
  }//endfor iState

  // Generate the other half
  for(iState=0;iState<numStates;iState++){
    for(jState=0;jState<iState;jState++){
      vnlMatrix[iState*numStates+jState] = vnlMatrix[jState*numStates+iState];
      for(iAtom=0;iAtom<numAtomFragVnlCalc;iAtom++){
        vnlFxMatrix[iAtom*numStates*numStates+iState*numStates+jState] 
		= vnlFxMatrix[iAtom*numStates*numStates+jState*numStates+iState];
        vnlFyMatrix[iAtom*numStates*numStates+iState*numStates+jState] 
		= vnlFyMatrix[iAtom*numStates*numStates+jState*numStates+iState];
        vnlFzMatrix[iAtom*numStates*numStates+iState*numStates+jState] 
		= vnlFzMatrix[iAtom*numStates*numStates+jState*numStates+iState];
      }//endfor iAtom
    }//endfor jState
  }//endfor iState


/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/

