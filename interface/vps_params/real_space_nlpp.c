/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: real_space_nlpp.c                            */
/*                                                                          */
/* Redo the vps initialization for real space nlpp			    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_vps_params_entry.h"
#include "../proto_defs/proto_search_entry.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_vps_params_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"

#define DEBUG_DKNY_OFF

#define JUERG_FACTOR_ON
#ifdef  JUERG_FACTOR_ON
#define JUERG_FACTOR 0.72
#else
#define JUERG_FACTOR 1.0
#endif
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void controlNlppReal(CP *cp,CLASS *class,GENERAL_DATA *generalData,
                        FILENAME_PARSE *filename_parse)
/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  
  ATOMMAPS *atommaps = &(class->atommaps);
  PSEUDO *pseudo = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CELL *cell = &(generalData->cell);
  VPS_FILE *vpsFile = pseudo->vps_file;
  CPEWALD *cpewald = &(cp->cpewald);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sparse = &(cp->cp_sclr_fft_pkg3d_sparse);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  COMMUNICATE *commCP = &(cp->communicate);

  int myidState = commCP->myid_state;
  int numProcStates = commCP->np_states;
  int iType,iRad,iAng,rGrid,iAtom;
  int smoothOpt	    = pseudoReal->smoothOpt;
  int numAtomType   = atommaps->natm_typ;
  int numAtomTot    = clatoms_info->natm_tot;
  int realSparseOpt = cpewald->realSparseOpt;
  int nkf1,nkf2,nkf3;
  //int nkf1 = cp_para_fft_pkg3d_lg->nkf1;
  //int nkf2 = cp_para_fft_pkg3d_lg->nkf2;
  //int nkf3 = cp_para_fft_pkg3d_lg->nkf3;
  int countRad;
  int numR,angNow,numRadTot;
  int countR = 0;
  int numGridTot = 0;
  int numGSm,numGLg;
  MPI_Comm commStates   =    commCP->comm_states;  

  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *numLMax,*numRadMax;
  int *nrad_0 = pseudo->nrad_0;
  int *nrad_1 = pseudo->nrad_1;
  int *nrad_2 = pseudo->nrad_2;
  int *nrad_3 = pseudo->nrad_3;
  int *ivpsLabel = pseudo->ivps_label;
  int *lMap;
  int *numGridRadSmooth;
  int *numGridNlppMap;

  int **atomLRadNum,**atomRadMap;

  double radCutRatio = pseudoReal->radCutRatio;
  double rCutoffMax;
  double junk1,junk2;
  double z1,z2,alpha1,alpha2;
  double zPol,gamma;
  double dr,rMax,rMin,r,drBf;
  double dg = pseudoReal->dg;
  double aLength,bLength,cLength;
  double gmaxTrueSm = cpewald->gmaxTrueSm;
  double gmaxTrueLg = cpewald->gmaxTrueLg;
  double gmaxTrueLgLg;
  double a[3],b[3],c[3];
  double aiLength,biLength,ciLength;
  double pre = 2.0*M_PI;

  double *vLoc,*vNl,*phiNl;
  double *vNlTot;
  double *hmat = cell->hmat;
  double *hmati = cell->hmati;
  double *ppRealCut;
  double *vpsNormList;
  double *vNlG;

  FILE *fvps;
/*==========================================================================*/
/* I) Initialize radial function and angular channel                        */

  Barrier(commStates);
  if(myidState==0){
    printf("Start smoothing real space non-local pseudopotential...\n");
  }
  
  if(realSparseOpt==0){
    nkf1 = cp_para_fft_pkg3d_lg->nkf1;
    nkf2 = cp_para_fft_pkg3d_lg->nkf2;
    nkf3 = cp_para_fft_pkg3d_lg->nkf3;
  }
  else{
    nkf1 = cp_sclr_fft_pkg3d_sparse->nkf1;
    nkf2 = cp_sclr_fft_pkg3d_sparse->nkf2;
    nkf3 = cp_sclr_fft_pkg3d_sparse->nkf3;
  }

  pseudoReal->numLMax = (int*)cmalloc(numAtomType*sizeof(int));  
  pseudoReal->numRadMax = (int*)cmalloc(numAtomType*sizeof(int));
  pseudoReal->numGridRadSmooth = (int*)cmalloc(numAtomType*sizeof(int));
  pseudoReal->atomLRadNum = (int**)cmalloc(numAtomType*sizeof(int*));
  pseudoReal->atomRadMap = (int**)cmalloc(numAtomType*sizeof(int*));
  pseudoReal->ppRealCut = (double*)cmalloc(numAtomType*sizeof(double));
  numGridRadSmooth = (int*)cmalloc(numAtomType*sizeof(int));
  numLMax = pseudoReal->numLMax;
  numRadMax = pseudoReal->numRadMax;
  numGridRadSmooth = pseudoReal->numGridRadSmooth;
  atomLRadNum = pseudoReal->atomLRadNum;
  atomRadMap = pseudoReal->atomRadMap;
  ppRealCut = pseudoReal->ppRealCut;

  countRad = 0;
  for(iType=0;iType<numAtomType;iType++){
    numLMax[iType] = pseudo->n_ang[iType+1];
    atomLRadNum[iType] = NULL;
    if(numLMax[iType]>0)atomLRadNum[iType] = (int*)cmalloc(numLMax[iType]*sizeof(int));
    switch(numLMax[iType]){
      case 0:
	numRadMax[iType] = 0;
	break;
      case 1:
	numRadMax[iType] = nrad_0[iType+1];
	atomLRadNum[iType][0] = nrad_0[iType+1];
	break;
      case 2:
	numRadMax[iType] = nrad_0[iType+1]+nrad_1[iType+1];
	atomLRadNum[iType][0] = nrad_0[iType+1];
        atomLRadNum[iType][1] = nrad_1[iType+1];
	break;
      case 3:
	numRadMax[iType] = nrad_0[iType+1]+nrad_1[iType+1]+nrad_2[iType+1];
        atomLRadNum[iType][0] = nrad_0[iType+1];
        atomLRadNum[iType][1] = nrad_1[iType+1];
	atomLRadNum[iType][2] = nrad_2[iType+1];
        break;	
      case 4:
	numRadMax[iType] = nrad_0[iType+1]+nrad_1[iType+1]+nrad_2[iType+1]+nrad_3[iType+1];
        atomLRadNum[iType][0] = nrad_0[iType+1];
        atomLRadNum[iType][1] = nrad_1[iType+1];
        atomLRadNum[iType][2] = nrad_2[iType+1];
        atomLRadNum[iType][3] = nrad_3[iType+1];
        break;
    }
    //printf("numRadMax %i %i %i %i\n",numRadMax[0],numLMax[0],atomLRadNum[0][0],atomLRadNum[0][1]);
    atomRadMap[iType] = NULL;
    if(numLMax[iType]>0){
      atomRadMap[iType] = (int*)cmalloc(numRadMax[iType]*sizeof(int));
      for(iRad=0;iRad<numRadMax[iType];iRad++){
        atomRadMap[iType][iRad] = countRad+iRad;
      }
    }
    //printf("atomRadMap %i %i\n",atomRadMap[0][0],atomRadMap[0][1]);
    countRad += numRadMax[iType];
  }
  numRadTot = countRad;
  pseudoReal->numRadTot = numRadTot;

  pseudoReal->vpsNormList = (double*)cmalloc(numRadTot*sizeof(double));
  vpsNormList = pseudoReal->vpsNormList;
  lMap = (int*)cmalloc(numRadTot*sizeof(int));
  

/*==========================================================================*/
/* II) Read and smooth Radial function			                    */
  
  // 1. Calcualte the minimum surface distance of the box
  a[0] = hmat[1];a[1] = hmat[2];a[2] = hmat[3];
  b[0] = hmat[4];b[1] = hmat[5];b[2] = hmat[6];
  c[0] = hmat[7];c[1] = hmat[8];c[2] = hmat[9];
  aLength = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  bLength = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  cLength = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  aiLength = sqrt(hmati[1]*hmati[1]+hmati[4]*hmati[4]+hmati[7]*hmati[7])*pre*nkf1;
  biLength = sqrt(hmati[2]*hmati[2]+hmati[5]*hmati[5]+hmati[8]*hmati[8])*pre*nkf2;
  ciLength = sqrt(hmati[3]*hmati[3]+hmati[6]*hmati[6]+hmati[9]*hmati[9])*pre*nkf3;
  
  gmaxTrueLgLg = MIN(aiLength,biLength);
  gmaxTrueLgLg = MIN(gmaxTrueLgLg,ciLength)-gmaxTrueSm;
  //gmaxTrueLgLg = MIN(gmaxTrueLgLg,ciLength);

  printf("ggggggg %lg %lg %lg\n",gmaxTrueSm,gmaxTrueLg,gmaxTrueLgLg);
  pseudoReal->gMaxSm = gmaxTrueSm*0.75;
  pseudoReal->gMaxLg = gmaxTrueLg;
  //pseudoReal->gMaxLg = 3.0*gmaxTrueSm;
  //pseudoReal->gMaxLg = 12.56060614640884;
  //pseudoReal->gMaxLg = gmaxTrueLgLg;
  numGSm = (int)(pseudoReal->gMaxSm/dg)+1;
  numGLg = (int)(pseudoReal->gMaxLg/dg)+1;
  pseudoReal->numGSm = numGSm;
  pseudoReal->numGLg = numGLg;
  //printf("gMaxSm %.8lg gMaxLg %.8lg\n",pseudoReal->gMaxSm,pseudoReal->gMaxLg);
  pseudoReal->vNlG = (double*)cmalloc(numRadTot*numGLg*sizeof(double));
  vNlG = pseudoReal->vNlG;
 
  // 2. Read the radial functions and determine the cutoff 
  countRad = 0;
  for(iType=0;iType<numAtomType;iType++){
    //printf("ivpsLabel %i\n",ivpsLabel[iType+1]==1);
    //fvps = fopen(vpsFile[iType+1].name,"r");
    rCutoffMax = -10000000.0;
    vLoc = NULL;
    vNl = NULL;
    phiNl = NULL;
    if(ivpsLabel[iType+1]==1){// KB
      if(myidState==0){
	fvps = fopen(vpsFile[iType+1].name,"r");
	fscanf(fvps,"%i %lg %i\n",&numR,&rMax,&angNow);
	fscanf(fvps,"%lg %lg %lg %lg\n",&z1,&alpha1,&z2,&alpha2);
	fscanf(fvps,"%lg %lg\n",&zPol,&gamma);
	//printf("%lg %lg %lg %lg %lg %lg\n",z1,alpha1,z2,alpha2,zPol,gamma);
      }
      if(numProcStates>1){
        printf("ffffffffffffffffuck %i\n",numR);
        Bcast(&numR,1,MPI_INT,0,commStates);
        Bcast(&rMax,1,MPI_DOUBLE,0,commStates);
        Bcast(&angNow,1,MPI_INT,0,commStates);
	Bcast(&z1,1,MPI_DOUBLE,0,commStates);
        Bcast(&alpha1,1,MPI_DOUBLE,0,commStates);
        Bcast(&z2,1,MPI_DOUBLE,0,commStates);
        Bcast(&alpha2,1,MPI_DOUBLE,0,commStates);
        Bcast(&zPol,1,MPI_DOUBLE,0,commStates);
        Bcast(&gamma,1,MPI_DOUBLE,0,commStates);
      }
      vLoc = (double*)cmalloc((numR)*sizeof(double));
      vNl = (double*)cmalloc((numR*numRadMax[iType])*sizeof(double));
      phiNl = (double*)cmalloc((numR*numRadMax[iType])*sizeof(double));
      drBf = dr;
      dr = rMax/(double)numR;
      if(iType>0&&fabs(dr-drBf)>1.0e-10){//dr doesn't match, need reset the pseudopoential file
	printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	printf("Real space grid spacing doesn't match between pseudopotential files.\n");
	printf("Please reset the number of real space grid point!\n");
	printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	fflush(stdout);
	exit(0);
      }
      if(myidState==0){
	// get the non-local part
	for(iAng=0;iAng<angNow;iAng++){
	  for(rGrid=0;rGrid<numR;rGrid++){
	    fscanf(fvps,"%lg",&vNl[iAng*numR+rGrid]);
	    fscanf(fvps,"%lg",&phiNl[iAng*numR+rGrid]);
	  }//endfor rGrid
	}//endfor iAng
	// get the local potential
	for(rGrid=0;rGrid<numR;rGrid++){
	  fscanf(fvps,"%lg",&vLoc[rGrid]);
	  fscanf(fvps,"%lg",&junk1);
	}//endfor rGrid
	fclose(fvps);
      }//endif myidState
      if(numProcStates>1){
        Bcast(&vNl[0],numR*numRadMax[iType],MPI_DOUBLE,0,commStates);
        Bcast(&phiNl[0],numR*numRadMax[iType],MPI_DOUBLE,0,commStates);
        Bcast(&vLoc[0],numR,MPI_DOUBLE,0,commStates);
      }
      // Substract the nonlocal part from local part
      for(iAng=0;iAng<angNow;iAng++){
	lMap[countR+iAng] = iAng;
	vpsNormList[countR+iAng] = 0.0;
	for(rGrid=0;rGrid<numR;rGrid++){
	  vNl[iAng*numR+rGrid] = (vNl[iAng*numR+rGrid]-vLoc[rGrid])*phiNl[iAng*numR+rGrid];
	  vpsNormList[countR+iAng] += vNl[iAng*numR+rGrid]*phiNl[iAng*numR+rGrid];
	  printf("1111111 rGridddddd %lg %lg\n",rGrid*dr,vNl[iAng*numR+rGrid]);
	}//endfor rGrid
	vpsNormList[countR+iAng] *= dr;
	//printf("countR %i vpsNormList %lg\n",countR+iAng,vpsNormList[countR+iAng]);
	vpsNormList[countR+iAng] = 1.0/vpsNormList[countR+iAng];
      }//endfor iAng
      for(iAng=0;iAng<angNow;iAng++){
	rGrid = numR-1;
	while(fabs(vNl[iAng*numR+rGrid])<1.0e-5&&rGrid>0){//double check 1.0e-10
	  //if(myidState==1)printf("rrrrrrrGrid %i\n",rGrid);
	  rGrid -= 1;
	}
	if(rGrid==numR-1){
	  printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	  printf("Your pseudo-potential does not decay enough! Please\n");
	  printf("choose bigger radius cutoff!\n");
	  printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	  fflush(stdout);
	  exit(0);
	}
	r = (rGrid+1.0)*dr*radCutRatio;
	if(r>=0.5*aLength||r>=0.5*bLength||r>=0.5*cLength){
	  printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	  printf("Your box is too small to use real space non-local\n");
	  printf("pseudopotential. Please use the pseudopotential in\n");
	  printf("k space for this small system!\n");
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	  fflush(stdout);
	  exit(0);
	}
	if(r>rCutoffMax)rCutoffMax = r;
      }//endfor iAng
      ppRealCut[iType] = rCutoffMax;
      numGridRadSmooth[iType] = (int)(rCutoffMax/dr)+1;
    }//endif ivpsLabel
    // 3. Smooth the radius function
    
    numGridRadSmooth[iType] = (int)(rCutoffMax/dr)+1;

    for(iRad=0;iRad<numRadMax[iType];iRad++){
      switch(smoothOpt){
	case 1:
	  nlppSmoothKS(pseudo,&vNl[iRad*numR],&vNlG[(countR+iRad)*numGLg],
		       iType,rMax,numR,lMap[countR+iRad]);
	  break;
	case 2:
	  nlppSmoothRoi(pseudo,&vNl[iRad*numR],&vNlG[(countR+iRad)*numGLg],
                        iType,rMax,numR,lMap[countR+iRad]);
	  break;
      }//endswitch    
    }//endfor iRad
    free(vLoc);
    free(vNl);
    free(phiNl);
    countR += numRadMax[iType];
  }//endfor iType

  //pseudoReal->dr = rGridSpacing[0]; //debug, need change in the future
  //dr = pseudoReal->dr;
  pseudoReal->dr = dr;
/*==========================================================================*/
/* III) Real space spline						    */

  interpReal(pseudo,numAtomType,lMap);

/*======================================================================*/
/* IV) Initialize dot product                                           */

  int *numNlppAtom;
  int atomType;
  double **dotReAll,**dotImAll;
  pseudoReal->numNlppAtom = (int*)cmalloc(numAtomType*sizeof(int));
  pseudoReal->dotReAll = (double**)cmalloc(numAtomTot*sizeof(double*));
  pseudoReal->dotImAll = (double**)cmalloc(numAtomTot*sizeof(double*));
  numNlppAtom = pseudoReal->numNlppAtom;
  dotReAll = pseudoReal->dotReAll;
  dotImAll = pseudoReal->dotImAll;
  for(iType=0;iType<numAtomType;iType++){
    numNlppAtom[iType] = 0;
    for(iAng=0;iAng<numLMax[iType];iAng++){
      numNlppAtom[iType] += atomLRadNum[iType][iAng]*(iAng+1);
    }//endfor iAng
  }//endfor iType
  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    atomType = iAtomAtomType[iAtom+1]-1;
    dotReAll[iAtom] = NULL;
    dotImAll[iAtom] = NULL;
    if(numNlppAtom[atomType]>=1){
      //printf("atomType %i %i\n",atomType,numNlppAtom[atomType]);
      dotReAll[iAtom] = (double*)cmalloc(numNlppAtom[atomType]*sizeof(double));
      dotImAll[iAtom] = (double*)cmalloc((numNlppAtom[atomType]-1)*sizeof(double));
    }
  }//endfor iAtom

/*======================================================================*/
/* V) Initialize wave function                                          */

  // I need to set the pointer to NULL here since I need to use realloc
  pseudoReal->pseudoWfCalcFlag = 1; // enforce pseudo wf calculation at begining
  pseudoReal->forceCalcFlag = 1;
  pseudoReal->vnlPhiAtomGridRe = NULL;
  pseudoReal->vnlPhiAtomGridIm = NULL;
  pseudoReal->numGridNlppMap = (int*)cmalloc(numAtomTot*sizeof(int));
  pseudoReal->gridNlppMap = (int**)cmalloc(numAtomTot*sizeof(int*));
  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    pseudoReal->gridNlppMap[iAtom] = NULL;
  }
  pseudoReal->vnlPhiDxAtomGridRe = NULL;
  pseudoReal->vnlPhiDxAtomGridIm = NULL;
  pseudoReal->vnlPhiDyAtomGridRe = NULL;
  pseudoReal->vnlPhiDyAtomGridIm = NULL;
  pseudoReal->vnlPhiDzAtomGridRe = NULL;
  pseudoReal->vnlPhiDzAtomGridIm = NULL;

/*======================================================================*/
/* VI) Initialize other flags                                           */

  pseudoReal->nlppForceOnly = 0;

  if(myidState==0){
    printf("Finish smoothing real space non-local pseudopotential...\n");
    fflush(stdout);
  }
  Barrier(commStates);
  
/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initRealNlppWf(CP *cp,CLASS *class,GENERAL_DATA *generalData)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Real space nlpp, only used in filtering				 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  PSEUDO *pseudo = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  ATOMMAPS *atommaps = &(class->atommaps);
  COMMUNICATE *commCP = &(cp->communicate);

  int numAtomTot = clatoms_info->natm_tot;
  int iAtom,iAng,m,iRad;
  int numGrid;
  int numPsuedoWfRe,numPsuedoWfIm;
  int numGridTotRe = 0;
  int numGridTotIm = 0;
  int atomType;
  int forceCalcFlag = pseudoReal->forceCalcFlag;
  int gridShiftNowRe,gridShiftNowIm;
  int *numGridNlppMap;
  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *numLMax = pseudoReal->numLMax;
  int **atomLRadNum = pseudoReal->atomLRadNum;
  int *gridStIndRe,*gridStIndIm;

/*==========================================================================*/
/* I) Construct Real Space Grids Neighbourhood List                         */

  mapRealSpaceGrid(cp,class,generalData);
  testOverlap(cp,class,generalData);
  //fflush(stdout);
  //exit(0);

/*==========================================================================*/
/* II) Construct Real Space V_nl*Phi_nl*Y_lm	                            */
  numGridNlppMap = pseudoReal->numGridNlppMap;
  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    numGrid = numGridNlppMap[iAtom];
    atomType = iAtomAtomType[iAtom+1]-1;
    for(iAng=0;iAng<numLMax[atomType];iAng++){
      numPsuedoWfRe = atomLRadNum[atomType][iAng]*(iAng+1);
      numPsuedoWfIm = atomLRadNum[atomType][iAng]*iAng;
      numGridTotRe += numGrid*numPsuedoWfRe;
      numGridTotIm += numGrid*numPsuedoWfIm;
    }
  }

  gridShiftNowRe = 0;
  gridShiftNowIm = 0;
  pseudoReal->gridStIndRe = (int*)cmalloc(numAtomTot*sizeof(int));
  pseudoReal->gridStIndIm = (int*)cmalloc(numAtomTot*sizeof(int));
  gridStIndRe = pseudoReal->gridStIndRe;
  gridStIndIm = pseudoReal->gridStIndIm;

  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    gridStIndRe[iAtom] = gridShiftNowRe;
    gridStIndIm[iAtom] = gridShiftNowIm;
    atomType = iAtomAtomType[iAtom+1]-1;
    numGrid = numGridNlppMap[iAtom];
    if(numGrid>0){
      for(iAng=0;iAng<numLMax[atomType];iAng++){
        for(iRad=0;iRad<atomLRadNum[atomType][iAng];iRad++){
          for(m=0;m<=iAng;m++){
            if(m!=0){
              gridShiftNowRe += numGrid;
              gridShiftNowIm += numGrid;
            }
            else{
              gridShiftNowRe += numGrid;
            }//endif m
          }//endfor m
        }//endfor iRad
      }//endfor l
    }//endif numGrid
  }//endfor iAtom


  pseudoReal->vnlPhiAtomGridRe = (double*)realloc(pseudoReal->vnlPhiAtomGridRe,
						  numGridTotRe*sizeof(double));  
  pseudoReal->vnlPhiAtomGridIm = (double*)realloc(pseudoReal->vnlPhiAtomGridIm,
						  numGridTotIm*sizeof(double));
  calcPseudoWf(cp,class,generalData);

/*==========================================================================*/
/* III) Calculate d|V_nl*Phi_nl*Y_lm>/dR                                    */

  if(forceCalcFlag==1){
    pseudoReal->vnlPhiDxAtomGridRe = (double*)realloc(pseudoReal->vnlPhiDxAtomGridRe,
						      numGridTotRe*sizeof(double));
    pseudoReal->vnlPhiDxAtomGridIm = (double*)realloc(pseudoReal->vnlPhiDxAtomGridIm,
						      numGridTotIm*sizeof(double));
    pseudoReal->vnlPhiDyAtomGridRe = (double*)realloc(pseudoReal->vnlPhiDyAtomGridRe,
						      numGridTotRe*sizeof(double));
    pseudoReal->vnlPhiDyAtomGridIm = (double*)realloc(pseudoReal->vnlPhiDyAtomGridIm,
						      numGridTotIm*sizeof(double));
    pseudoReal->vnlPhiDzAtomGridRe = (double*)realloc(pseudoReal->vnlPhiDzAtomGridRe,
						      numGridTotRe*sizeof(double));
    pseudoReal->vnlPhiDzAtomGridIm = (double*)realloc(pseudoReal->vnlPhiDzAtomGridIm,
						      numGridTotIm*sizeof(double));
    calcPseudoWfDev(cp,class,generalData);    
  }
  //fflush(stdout);
  //exit(0);
  //printf("finish initializing nl pp wf\n");
  //Barrier(commCP->comm_states);
  
/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void nlppSmoothKS(PSEUDO *pseudo,double *vNl,double *vNlG,
		  int iType,double rMax,int numR,int l)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Real space nlpp, only used in filtering                   */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);  

  int numGSm = pseudoReal->numGSm;
  int numGLg = pseudoReal->numGLg;
  int ig,ir;
  int numGridRadSmooth = pseudoReal->numGridRadSmooth[iType];

  double dg = pseudoReal->dg;
  double dr = rMax/((double)numR);
  double *gGrid,*r;

  gGrid = (double*)cmalloc((numGLg)*sizeof(double));

  for(ig=0;ig<numGLg;ig++){
    gGrid[ig] = ig*dg;
    vNlG[ig] = 0.0;
  }
  //for(ir=0;ir<numR;ir++)r[ir] = ir*dr;

  // Bessel transform to g space from g=0 to gMaxSm
  //printf("numR %i %lg\n",numR,vNl[numR-1]);
  //bessTransform(vNl,numR,dr,l,vNlG,numGSm,gGrid);
  bessTransform(vNl,numR,dr,l,vNlG,numGLg,gGrid);
  
  /*
  for(ig=0;ig<numGLg;ig++){
    printf("111111111 igggggg %lg %lg\n",gGrid[ig],vNlG[ig]);
  }
  for(ig=numGSm;ig<numGLg;ig++)vNlG[ig] = 0.0;
  */

  // Optimize g space coeffcient from gMaxSm to gMaxLg

  optGCoeff(pseudoReal,numGLg,numGSm,numR,dr,dg,numGridRadSmooth,l,vNlG);

  //printf("rMax %lg\n",numGridRadSmooth*dr);
  /*
  for(ig=0;ig<numGLg;ig++){
    printf("222222 igggggg %lg %lg\n",gGrid[ig],vNlG[ig]);
  }
  */

  for(ig=0;ig<numGLg;ig++)vNlG[ig] *= gGrid[ig];

  //fflush(stdout);
  //exit(0);
 
/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void nlppSmoothRoi(PSEUDO *pseudo,double *vNl,double *vNlG,
                  int iType,double rMax,int numR,int l)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Real space nlpp, only used in filtering                   */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void bessTransform(double *funIn,int numIn,double dx,int l,double *funOut,
		   int numOut,double *yList)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Bessel transform */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  int iGrid,jGrid;

  double x,y;
  double arg;
  //double fpidx = 4.0*M_PI*dx;
  //check the output argument
  for(iGrid=0;iGrid<numOut;iGrid++){
    y = yList[iGrid];
    if(y<0.0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Output argument are |g| or |r| values. They can not\n"); 
      printf("be negative. The value you provide is %lg.\n",y);
      printf("Please check your code!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }//endif y
  }//endfor iGrid

  switch(l){
    case 0:
      for(iGrid=0;iGrid<numOut;iGrid++){
	funOut[iGrid] = 0.0;
	y = yList[iGrid];
	if(y>1.0e-100){
	  for(jGrid=1;jGrid<numIn;jGrid++){
	    x = jGrid*dx;
	    arg = x*y;
	    funOut[iGrid] += j0(arg)*x*funIn[jGrid]*dx;
	  }//endfor jGrid
	}//endif y
	if(y<=1.0e-100&&y>=0.0){
	  for(jGrid=1;jGrid<numIn;jGrid++){
	    x = jGrid*dx;
	    funOut[iGrid] += x*funIn[jGrid]*dx;
	  }//endfor jGrid
	}//endif y
      }//endfor iGrid
      break;
    case 1:
      for(iGrid=0;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = yList[iGrid];
	if(y>1.0e-100){
	  for(jGrid=1;jGrid<numIn;jGrid++){
	    x = jGrid*dx;
	    arg = x*y;
	    funOut[iGrid] += j1(arg)*x*funIn[jGrid]*dx;
	  }//endfor jGrid
	}//endif y
      }//endfor iGrid
      break;
    case 2:
      for(iGrid=0;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = yList[iGrid];
        if(y>1.0e-100){
	  for(jGrid=1;jGrid<numIn;jGrid++){
	    x = jGrid*dx;
	    arg = x*y;
	    funOut[iGrid] += j2(arg)*x*funIn[jGrid]*dx;
	  }//endfor jGrid
	}//endif y
      }//endfor iGrid
      break;
  }//endswitch

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void optGCoeff(PSEUDO_REAL *pseudoReal,int numGLg,int numGSm,int numR,
		double dr,double dg,int numRCutoff, int l,double *vNlG)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Bessel transform */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int iGridG,jGridG,kGridR;
  int numGSolve = numGLg-numGSm;
  int gridShift;

  double qi,qj;
  double r,arg1,arg2;

  double *A;
  double *B;
  double *subA;
  //Construct the A matrix

  A = (double*)calloc(numGLg*numGLg,sizeof(double));
  B = (double*)calloc(numGSolve,sizeof(double));
  subA = (double*)calloc(numGSolve*numGSolve,sizeof(double));
  
  //iGridG=0
  for(jGridG=0;jGridG<numGLg;jGridG++){
    A[jGridG] = 0.0;
    A[jGridG*numGLg] = 0.0;
  }

  for(iGridG=1;iGridG<numGLg;iGridG++){
    //printf("iGridG %i\n",iGridG);
    qi = iGridG*dg;
    for(jGridG=iGridG;jGridG<numGLg;jGridG++){
      qj = jGridG*dg; 
      switch(l){
        case 0:
	  for(kGridR=1;kGridR<numRCutoff;kGridR++){
	    r = kGridR*dr;
	    arg1 = qi*r;
	    arg2 = qj*r;
	    A[iGridG*numGLg+jGridG] += j0(arg1)*j0(arg2)*r*r;
	  }
	  A[iGridG*numGLg+jGridG] *= dr*qi*qi*qj*qj;
	  break;
	case 1:
	  for(kGridR=1;kGridR<numRCutoff;kGridR++){
            r = kGridR*dr;
            arg1 = qi*r;
            arg2 = qj*r;
            A[iGridG*numGLg+jGridG] += j1(arg1)*j1(arg2)*r*r;
          }
          A[iGridG*numGLg+jGridG] *= dr*qi*qi*qj*qj;
          break;	  
	case 2:
          for(kGridR=1;kGridR<numRCutoff;kGridR++){
            r = kGridR*dr;
            arg1 = qi*r;
            arg2 = qj*r;
            A[iGridG*numGLg+jGridG] += j2(arg1)*j2(arg2)*r*r;
          }
          A[iGridG*numGLg+jGridG] *= dr*qi*qi*qj*qj;
          break;
      }//endswitch l
      A[jGridG*numGLg+iGridG] = A[iGridG*numGLg+jGridG];
    }//endfor jGridG
  }//endfor iGridG
  
  //Construct the linear system
  // B=diag{pi/2*q^2}*kai-subA*kai
  // B vector

  for(iGridG=0;iGridG<numGSolve;iGridG++){
    for(jGridG=0;jGridG<numGSm;jGridG++){
      B[iGridG] += A[(iGridG+numGSm)*numGLg+jGridG]*vNlG[jGridG];
    }//endfor jGridG
    B[iGridG] *= dg;
  }//endfor iGridG  
  for(iGridG=0;iGridG<numGSolve;iGridG++){
    qi = (iGridG+numGSm)*dg;
    for(jGridG=0;jGridG<numGSolve;jGridG++){
      subA[iGridG*numGSolve+jGridG] = -A[(iGridG+numGSm)*numGLg+jGridG+numGSm]*dg;
    }//endfor jGridG
    subA[iGridG*numGSolve+iGridG] += 0.5*M_PI*qi*qi;
  }//endfor iGridG

  // Solving the linear system
  dsysvWrapper(subA,B,numGSolve);

  for(iGridG=0;iGridG<numGSolve;iGridG++){
    vNlG[iGridG+numGSm] = B[iGridG];
  }
  free(A);
  free(B);
  free(subA);

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void bessTransformGrad(double *funIn,int numIn,double dx,int l,double *funOut,
                   int numOut,double *yList)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Bessel transform */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  int iGrid,jGrid;

  double x,y;
  double arg;
  double djl0; //derivative of spherical bessel function at x=0
  //double fpidx = 4.0*M_PI*dx;
  double test1,test2;

  //check the output argument
  for(iGrid=0;iGrid<numOut;iGrid++){
    y = yList[iGrid];
    if(y<0.0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Output argument are |g| or |r| values. They can not\n");
      printf("be negative. The value you provide is %lg.\n",y);
      printf("Please check your code!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }//endif y
  }//endfor iGrid

  /*
  // r=0
  switch(l){
    case 0: djl0 = 0.0; break;
    case 1: djl0 = 1.0/3.0; break;
    case 2: djl0 = 0.0; break;
  }
  funOut[0] = 0.0;
  for(jGrid=1;jGrid<numIn;jGrid++){
    x = jGrid*dx;
    funOut[0] += x*x*funIn[jGrid];
  }
  funOut[0] *= djl0*dx;
  */

  switch(l){
    case 0:
      for(iGrid=0;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = yList[iGrid];
	if(y>1.0e-100){
	  for(jGrid=1;jGrid<numIn;jGrid++){
	    x = jGrid*dx;
	    arg = x*y;
	    test1 = dj0(arg);
	    test2 = funIn[jGrid];
	    funOut[iGrid] += dj0(arg)*x*x*funIn[jGrid]*dx;
	  }//endfor jGrid
	}//endif
      }//endfor iGrid
      break;
    case 1:
      for(iGrid=0;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = yList[iGrid];
	if(y>1.0e-100){
	  for(jGrid=1;jGrid<numIn;jGrid++){
	    x = jGrid*dx;
	    arg = x*y;
	    funOut[iGrid] += dj1(arg)*x*x*funIn[jGrid]*dx;
	  }//endfor jGrid
	}//endif
        if(y<=1.0e-100&&y>=0.0){
	  for(jGrid=1;jGrid<numIn;jGrid++){
	    x = jGrid*dx;
	    funOut[iGrid] += x*x*funIn[jGrid];
	  }//endfor
	  funOut[iGrid] *= dx/3.0;
	}//endif y
      }//endfor iGrid
      break;
    case 2:
      for(iGrid=0;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = yList[iGrid];
	if(y>1.0e-100){
	  for(jGrid=1;jGrid<numIn;jGrid++){
	    x = jGrid*dx;
	    arg = x*y;
	    funOut[iGrid] += dj2(arg)*x*x*funIn[jGrid]*dx;
	  }//endfor jGrid
	}//endif y
      }//endfor iGrid
      break;
  }//endswitch

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void bessTransformGradGrad(double *funIn,int numIn,double dx,int l,double *funOut,
                   int numOut,double *yList)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Bessel transform */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  int iGrid,jGrid;

  double x,y;
  double arg;
  double djl0; //derivative of spherical bessel function at x=0
  //double fpidx = 4.0*M_PI*dx;
  double test1,test2;

  //check the output argument
  for(iGrid=0;iGrid<numOut;iGrid++){
    y = yList[iGrid];
    if(y<0.0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Output argument are |g| or |r| values. They can not\n");
      printf("be negative. The value you provide is %lg.\n",y);
      printf("Please check your code!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }//endif y
  }//endfor iGrid

  /*
  // r=0
  switch(l){
    case 0: djl0 = 0.0; break;
    case 1: djl0 = 1.0/3.0; break;
    case 2: djl0 = 0.0; break;
  }
  funOut[0] = 0.0;
  for(jGrid=1;jGrid<numIn;jGrid++){
    x = jGrid*dx;
    funOut[0] += x*x*funIn[jGrid];
  }
  funOut[0] *= djl0*dx;
  */

  switch(l){
    case 0:
      for(iGrid=0;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = yList[iGrid];
        if(y>1.0e-100){
          for(jGrid=1;jGrid<numIn;jGrid++){
            x = jGrid*dx;
            arg = x*y;
            test1 = dj0(arg);
            test2 = funIn[jGrid];
            funOut[iGrid] += dj0(arg)*x*x*x*funIn[jGrid]*dx;
          }//endfor jGrid
        }//endif
      }//endfor iGrid
      break;
    case 1:
      for(iGrid=0;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = yList[iGrid];
        if(y>1.0e-100){
          for(jGrid=1;jGrid<numIn;jGrid++){
            x = jGrid*dx;
            arg = x*y;
            funOut[iGrid] += dj1(arg)*x*x*x*funIn[jGrid]*dx;
          }//endfor jGrid
        }//endif
        if(y<=1.0e-100&&y>=0.0){
          for(jGrid=1;jGrid<numIn;jGrid++){
            x = jGrid*dx;
            funOut[iGrid] += x*x*x*funIn[jGrid];
          }//endfor
          funOut[iGrid] *= dx/3.0;
        }//endif y
      }//endfor iGrid
      break;
    case 2:
      for(iGrid=0;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = yList[iGrid];
        if(y>1.0e-100){
          for(jGrid=1;jGrid<numIn;jGrid++){
            x = jGrid*dx;
            arg = x*y;
            funOut[iGrid] += dj2(arg)*x*x*x*funIn[jGrid]*dx;
          }//endfor jGrid
        }//endif y
      }//endfor iGrid
      break;
  }//endswitch

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void interpReal(PSEUDO *pseudo,int numAtomType,int *lMap)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Map the real space grids                                              */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
   
  int numInterpGrid = -1;
  int gridShift;
  int interpFlag = 0;
  int gridIndTest;
  int interpInd,interpGridSt;
  int numRadTot = pseudoReal->numRadTot;
  int rGrid,iRad,iType;
  int countBadSpline = 1;
  int numGLg = pseudoReal->numGLg;

  int *numGridRadSmooth = pseudoReal->numGridRadSmooth;

  double rTest,r0,h,pseudoTest,pseudoDevTest;
  double dr = pseudoReal->dr;
  double dg = pseudoReal->dg;
  double rMax;
  double rMin = 0.0;
  double pre = 2.0/M_PI;
  double diff;

  double *rList = NULL;
  double *rListTest = NULL;
  double *de = NULL;
  double *vNlG = pseudoReal->vNlG;
  double *ppRealCut = pseudoReal->ppRealCut;
  double *vNlSmooth = NULL;
  double *dvNlSmooth = NULL;
  double *ddvNlSmooth = NULL;
  double *vNlSmoothTest = NULL;
  double *dvNlSmoothTest = NULL;
  double *vpsReal0 = NULL;
  double *vpsReal1 = NULL;
  double *vpsReal2 = NULL;
  double *vpsReal3 = NULL;
  double *vpsDevReal0 = NULL;
  double *vpsDevReal1 = NULL;
  double *vpsDevReal2 = NULL;
  double *vpsDevReal3 = NULL;

  for(iType=0;iType<numAtomType;iType++){
    rGrid = numGridRadSmooth[iType];
    if(numInterpGrid<rGrid){
      numInterpGrid = rGrid;
    }
  }
  numInterpGrid += 5;
  rMax = dr*(numInterpGrid-1);

  for(iType=0;iType<numAtomType;iType++){
    ppRealCut[iType] = rMax-4.0*dr; // actura cutoff is smaller then interpolation cutoff
  }
  while(countBadSpline>0){
    // transform to real space
    vNlSmooth = (double*)realloc(vNlSmooth,numRadTot*numInterpGrid*sizeof(double));
    dvNlSmooth = (double*)realloc(dvNlSmooth,numRadTot*numInterpGrid*sizeof(double));
    ddvNlSmooth = (double*)realloc(ddvNlSmooth,numRadTot*numInterpGrid*sizeof(double));
    vNlSmoothTest = (double*)realloc(vNlSmoothTest,numRadTot*(numInterpGrid-1)*sizeof(double));
    dvNlSmoothTest = (double*)realloc(dvNlSmoothTest,numRadTot*(numInterpGrid-1)*sizeof(double));
    rList = (double*)realloc(rList,(numInterpGrid+1)*sizeof(double));
    rListTest = (double*)realloc(rListTest,numInterpGrid*sizeof(double));
    for(rGrid=0;rGrid<numInterpGrid;rGrid++){
      rList[rGrid+1] = rGrid*dr;
      rListTest[rGrid] = (rGrid+0.5)*dr;
    }

    for(iRad=0;iRad<numRadTot;iRad++){
      bessTransform(&vNlG[iRad*numGLg],numGLg,dg,lMap[iRad],
		    &vNlSmooth[iRad*numInterpGrid],numInterpGrid,&rList[1]);
      bessTransformGrad(&vNlG[iRad*numGLg],numGLg,dg,lMap[iRad],
			&dvNlSmooth[iRad*numInterpGrid],numInterpGrid,&rList[1]);
      bessTransformGradGrad(&vNlG[iRad*numGLg],numGLg,dg,lMap[iRad],
			    &ddvNlSmooth[iRad*numInterpGrid],numInterpGrid,&rList[1]);
      bessTransform(&vNlG[iRad*numGLg],numGLg,dg,lMap[iRad],
		    &vNlSmoothTest[iRad*(numInterpGrid-1)],numInterpGrid-1,&rListTest[0]);
      bessTransformGrad(&vNlG[iRad*numGLg],numGLg,dg,lMap[iRad],
                    &dvNlSmoothTest[iRad*(numInterpGrid-1)],numInterpGrid-1,&rListTest[0]);
      /*
      for(rGrid=0;rGrid<numInterpGrid;rGrid++){
        printf("222222222 iRad %i rGriddddddd %lg %lg %lg\n",iRad,rGrid*dr,rGrid*dr*vNlSmooth[iRad*numInterpGrid+rGrid],dvNlSmooth[iRad*numInterpGrid+rGrid]);
      }
      */
      
    }
    for(rGrid=0;rGrid<numRadTot*numInterpGrid;rGrid++){
      vNlSmooth[rGrid] *= pre;
      dvNlSmooth[rGrid] *= pre;
      ddvNlSmooth[rGrid] *= pre;
    }
    for(rGrid=0;rGrid<numRadTot*(numInterpGrid-1);rGrid++){
      vNlSmoothTest[rGrid] *= pre;
      dvNlSmoothTest[rGrid] *= pre;
    }
    // spline
    //pseudoReal->numInterpGrid = numInterpGrid;
    //rList = (double*)cmalloc((numInterpGrid+1)*sizeof(double));
    de = (double*)realloc(de,(numInterpGrid*numRadTot+1)*sizeof(double));
    vpsReal0 = (double*)realloc(vpsReal0,(numInterpGrid*numRadTot+1)*sizeof(double));
    vpsReal1 = (double*)realloc(vpsReal1,(numInterpGrid*numRadTot+1)*sizeof(double));
    vpsReal2 = (double*)realloc(vpsReal2,(numInterpGrid*numRadTot+1)*sizeof(double));
    vpsReal3 = (double*)realloc(vpsReal3,(numInterpGrid*numRadTot+1)*sizeof(double));
    vpsDevReal0 = (double*)realloc(vpsDevReal0,(numInterpGrid*numRadTot+1)*sizeof(double));
    vpsDevReal1 = (double*)realloc(vpsDevReal1,(numInterpGrid*numRadTot+1)*sizeof(double));
    vpsDevReal2 = (double*)realloc(vpsDevReal2,(numInterpGrid*numRadTot+1)*sizeof(double));
    vpsDevReal3 = (double*)realloc(vpsDevReal3,(numInterpGrid*numRadTot+1)*sizeof(double));

    for(rGrid=0;rGrid<numInterpGrid*numRadTot+1;rGrid++){
      de[rGrid] = 0.0;
      vpsReal0[rGrid] = 0.0;
      vpsReal1[rGrid] = 0.0;
      vpsReal2[rGrid] = 0.0;
      vpsReal3[rGrid] = 0.0;
      vpsDevReal0[rGrid] = 0.0;
      vpsDevReal1[rGrid] = 0.0;
      vpsDevReal2[rGrid] = 0.0;
      vpsDevReal3[rGrid] = 0.0;
    }

    //for(rGrid=0;rGrid<numInterpGrid;rGrid++)rList[rGrid+1] = rGrid*dr;
    memcpy(&vpsReal0[1],&vNlSmooth[0],numInterpGrid*numRadTot*sizeof(double));
    memcpy(&de[1],&dvNlSmooth[0],numInterpGrid*numRadTot*sizeof(double));

    for(iRad=0;iRad<numRadTot;iRad++){
      gridShift = iRad*numInterpGrid;
      splineFitWithDerivative(&(vpsReal0[gridShift]),&(vpsReal1[gridShift]),
			      &(vpsReal2[gridShift]),&(vpsReal3[gridShift]),rList,
			      &(de[gridShift]),numInterpGrid);


    }//endfor iRad
    memcpy(&vpsDevReal0[1],&dvNlSmooth[0],numInterpGrid*numRadTot*sizeof(double));
    memcpy(&de[1],&ddvNlSmooth[0],numInterpGrid*numRadTot*sizeof(double));
    for(iRad=0;iRad<numRadTot;iRad++){
      gridShift = iRad*numInterpGrid;
      splineFitWithDerivative(&(vpsDevReal0[gridShift]),&(vpsDevReal1[gridShift]),
                              &(vpsDevReal2[gridShift]),&(vpsDevReal3[gridShift]),
			      rList,&(de[gridShift]),numInterpGrid);


    }//endfor iRad
    // test middle point at each interval

    countBadSpline = 0;
    //double diffMax = 0.0;
    for(iRad=0;iRad<numRadTot;iRad++){
      interpGridSt = iRad*numInterpGrid;
      for(rGrid=0;rGrid<numInterpGrid-1;rGrid++){
	rTest = rListTest[rGrid];
	gridIndTest = (int)((rTest-rMin)/dr)+1;
	r0 = (gridIndTest-1)*dr+rMin;
	h = rTest-r0;
	interpInd = interpGridSt+gridIndTest;
	pseudoTest = ((vpsReal3[interpInd]*h+vpsReal2[interpInd])*h+
			vpsReal1[interpInd])*h+vpsReal0[interpInd];
	pseudoDevTest = ((vpsDevReal3[interpInd]*h+vpsDevReal2[interpInd])*h+
			    vpsDevReal1[interpInd])*h+vpsDevReal0[interpInd];
	//printf("111111 %lg %lg %lg\n",rTest,pseudoTest,vNlSmoothTest[iRad*(numInterpGrid-1)+rGrid]);
	diff = pseudoTest-vNlSmoothTest[iRad*(numInterpGrid-1)+rGrid];
	if(diff*diff>1.0e-12){
	  countBadSpline += 1;
	  //printf("111111 %lg %lg %lg\n",rTest,pseudoTest,vNlSmoothTest[iRad*(numInterpGrid-1)+rGrid]);
        }
	diff = pseudoDevTest-dvNlSmoothTest[iRad*(numInterpGrid-1)+rGrid];
	if(diff*diff>1.0e-12){
	  countBadSpline += 1;
	  //printf("111111 %lg %lg %lg\n",rTest,pseudoDevTest,dvNlSmoothTest[iRad*(numInterpGrid-1)+rGrid]);
	}
	//if(fabs(diff)>diffMax)diffMax = fabs(diff);
      }//endfor rGrid
    }//endfor iRad
    /*
    printf("countBadSpline %i numInterpGrid %i diffMax %lg\n",countBadSpline,numInterpGrid,diffMax);
    fflush(stdout);
    exit(0);
    */

    // calculate flag
    if(countBadSpline>0){
      numInterpGrid *= 2;
      //printf("New number of interpolation point for rational %i\n",numInterpGrid);
      dr = rMax/((double)(numInterpGrid-1));
    }
  }//endwhile

  pseudoReal->numInterpGrid = numInterpGrid;
  pseudoReal->dr = dr;
  pseudoReal->vpsReal0 = (double*)cmalloc((numInterpGrid*numRadTot+1)*sizeof(double));
  pseudoReal->vpsReal1 = (double*)cmalloc((numInterpGrid*numRadTot+1)*sizeof(double));
  pseudoReal->vpsReal2 = (double*)cmalloc((numInterpGrid*numRadTot+1)*sizeof(double));
  pseudoReal->vpsReal3 = (double*)cmalloc((numInterpGrid*numRadTot+1)*sizeof(double));
  pseudoReal->vpsDevReal0 = (double*)cmalloc((numInterpGrid*numRadTot+1)*sizeof(double));
  pseudoReal->vpsDevReal1 = (double*)cmalloc((numInterpGrid*numRadTot+1)*sizeof(double));
  pseudoReal->vpsDevReal2 = (double*)cmalloc((numInterpGrid*numRadTot+1)*sizeof(double));
  pseudoReal->vpsDevReal3 = (double*)cmalloc((numInterpGrid*numRadTot+1)*sizeof(double));

  memcpy(&pseudoReal->vpsReal0[1],&vpsReal0[1],numInterpGrid*numRadTot*sizeof(double));
  memcpy(&pseudoReal->vpsReal1[1],&vpsReal1[1],numInterpGrid*numRadTot*sizeof(double));
  memcpy(&pseudoReal->vpsReal2[1],&vpsReal2[1],numInterpGrid*numRadTot*sizeof(double));
  memcpy(&pseudoReal->vpsReal3[1],&vpsReal3[1],numInterpGrid*numRadTot*sizeof(double));
  memcpy(&pseudoReal->vpsDevReal0[1],&vpsDevReal0[1],numInterpGrid*numRadTot*sizeof(double));
  memcpy(&pseudoReal->vpsDevReal1[1],&vpsDevReal1[1],numInterpGrid*numRadTot*sizeof(double));
  memcpy(&pseudoReal->vpsDevReal2[1],&vpsDevReal2[1],numInterpGrid*numRadTot*sizeof(double));
  memcpy(&pseudoReal->vpsDevReal3[1],&vpsDevReal3[1],numInterpGrid*numRadTot*sizeof(double));


  free(vNlSmooth);
  free(dvNlSmooth);
  free(ddvNlSmooth);
  free(vpsReal0);
  free(vpsReal1);
  free(vpsReal2);
  free(vpsReal3);
  free(vpsDevReal0);
  free(vpsDevReal1);
  free(vpsDevReal2);
  free(vpsDevReal3);
  free(de);
  free(rList);
  free(rListTest);

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void mapRealSpaceGrid(CP *cp, CLASS *class, GENERAL_DATA *generalData)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Map the real space grids						 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  ATOMMAPS *atommaps = &(class->atommaps);
  CPOPTS *cpopts = &(cp->cpopts);
  PSEUDO *pseudo = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CELL *cell = &(generalData->cell);
  VPS_FILE *vpsFile = pseudo->vps_file;
  CPEWALD *cpewald = &(cp->cpewald);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sparse = &(cp->cp_sclr_fft_pkg3d_sparse);  
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS *clatoms_pos = &(class->clatoms_pos[1]);

  int numAtomTot = clatoms_info->natm_tot;
  int atomType,l;
  int iAtom,i,j,k,iGrid;
  int numGridCount;
  int realSparseOpt = cpewald->realSparseOpt;
  int nkc,nkb,nka;
  //int nkc = cp_para_fft_pkg3d_lg->nkf3;
  //int nkb = cp_para_fft_pkg3d_lg->nkf2;
  //int nka = cp_para_fft_pkg3d_lg->nkf1;
  int aNum,bNum,cNum;
  int xInd,yInd,zInd,xGridInd,yGridInd,zGridInd;
  int gridInd;
  int numInterpGrid = pseudoReal->numInterpGrid;  
  int fftw3dFlag = cpopts->fftw3dFlag;

  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *numLMax = pseudoReal->numLMax;
  int *numGridNlppMap = pseudoReal->numGridNlppMap;
  int **gridNlppMap = pseudoReal->gridNlppMap;

  double dr = pseudoReal->dr;
  double rCutMax = numInterpGrid*dr;
  double aLength,bLength,cLength,adotb;
  double acrossbLength;
  double x,y,z,xScale,yScale,zScale;
  double xTemp,yTemp,zTemp,xDiff,yDiff,zDiff;
  double xTemp2,yTemp2,zTemp2;
  double xGrid,yGrid,zGrid;
  double distSq,cutOffSq;
  double debug[3];
  
  double *ppRealCut = pseudoReal->ppRealCut;
  double *hmati = cell->hmati;
  double *hmat = cell->hmat;
  double *xList = clatoms_pos->x;
  double *yList = clatoms_pos->y;
  double *zList = clatoms_pos->z;
 
  double aElem[3],bElem[3],cElem[3];
  double aOrth[3],bOrth[3],cOrth[3],acrossb[3];

  if(realSparseOpt==0){
    nkc = cp_para_fft_pkg3d_lg->nkf3;
    nkb = cp_para_fft_pkg3d_lg->nkf2;
    nka = cp_para_fft_pkg3d_lg->nkf1;
  }
  else{
    nkc = cp_sclr_fft_pkg3d_sparse->nkf3;
    nkb = cp_sclr_fft_pkg3d_sparse->nkf2;
    nka = cp_sclr_fft_pkg3d_sparse->nkf1;
  }

  aElem[0] = hmat[1]/nka;aElem[1] = hmat[2]/nka;aElem[2] = hmat[3]/nka;
  bElem[0] = hmat[4]/nkb;bElem[1] = hmat[5]/nkb;bElem[2] = hmat[6]/nkb;
  cElem[0] = hmat[7]/nkc;cElem[1] = hmat[8]/nkc;cElem[2] = hmat[9]/nkc;
  //printf("aElem %lg %lg %lg\n",aElem[0],aElem[1],aElem[2]);
  //printf("bElem %lg %lg %lg\n",bElem[0],bElem[1],bElem[2]);
  //printf("cElem %lg %lg %lg\n",cElem[0],cElem[1],cElem[2]);

  // orthogonalize a,b,c
  aOrth[0] = aElem[0];aOrth[1] = aElem[1];aOrth[2] = aElem[2];
  aLength = sqrt(aElem[0]*aElem[0]+aElem[1]*aElem[1]+aElem[2]*aElem[2]);
  adotb = aOrth[0]*bElem[0]+aOrth[1]*bElem[1]+aOrth[2]*bElem[2];
  bOrth[0] = bElem[0]-adotb*aOrth[0]/(aLength*aLength);
  bOrth[1] = bElem[1]-adotb*aOrth[1]/(aLength*aLength);
  bOrth[2] = bElem[2]-adotb*aOrth[2]/(aLength*aLength);
  bLength = sqrt(bOrth[0]*bOrth[0]+bOrth[1]*bOrth[1]+bOrth[2]*bOrth[2]);
  cross_product(&aOrth[0],&bOrth[0],&acrossb[0]);
  acrossbLength = sqrt(acrossb[0]*acrossb[0]+acrossb[1]*acrossb[1]+acrossb[2]*acrossb[2]);
  acrossb[0] /= acrossbLength;
  acrossb[1] /= acrossbLength;
  acrossb[2] /= acrossbLength;
  cLength = cElem[0]*acrossb[0]+cElem[1]*acrossb[1]+cElem[2]*acrossb[2];
  //printf("aLength %lg bLength %lg cLength %lg\n",aLength,bLength,cLength);

  aNum = (int)(rCutMax/aLength)+2;
  bNum = (int)(rCutMax/bLength)+2;
  cNum = (int)(rCutMax/cLength)+2;
  //printf("rCutMax %lg\n",rCutMax);
  //printf("Num %i %i %i\n",aNum,bNum,cNum);
  
  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    atomType = iAtomAtomType[iAtom+1]-1;
    //numGridNlppMap[iAtom] = 0;
    numGridCount = 0;
    if(numLMax[atomType]>0){
      gridNlppMap[iAtom] = (int *)cmalloc(100*sizeof(int)); 
      cutOffSq = ppRealCut[atomType]*ppRealCut[atomType];
      x = xList[iAtom+1];
      y = yList[iAtom+1];
      z = zList[iAtom+1];
      // determine the nearest grid point around nuclei position
      //printf("x %lg y %lg z %lg\n",x,y,z);
      xScale = (x*hmati[1]+y*hmati[4]+z*hmati[7])*nka;
      yScale = (x*hmati[2]+y*hmati[5]+z*hmati[8])*nkb;
      zScale = (x*hmati[3]+y*hmati[6]+z*hmati[9])*nkc;
      //printf("xScale %lg yScale %lg zScale %lg\n",xScale,yScale,zScale);
      xInd = NINT(xScale);
      yInd = NINT(yScale);
      zInd = NINT(zScale);
      //printf("xInd %i yInd %i zInd %i\n",xInd,yInd,zInd);    
      if(xInd>=nka)xInd -= nka;
      if(yInd>=nkb)yInd -= nkb;
      if(zInd>=nkc)zInd -= nkc;
      // search the grid points within the parallelogram
      for(i=-cNum;i<=cNum;i++){
	zGridInd = zInd+i;
	if(zGridInd<0)zGridInd += nkc;
	if(zGridInd>=nkc)zGridInd -= nkc;
	for(j=-bNum;j<=bNum;j++){
	  yGridInd = yInd+j;
	  if(yGridInd<0)yGridInd += nkb;
	  if(yGridInd>=nkb)yGridInd -= nkb;
	  for(k=-aNum;k<=aNum;k++){
	    xGridInd = xInd+k;
	    if(xGridInd<0)xGridInd += nka;
	    if(xGridInd>=nka)xGridInd -= nka;
	    xGrid = xGridInd*aElem[0]+yGridInd*bElem[0]+zGridInd*cElem[0];
            yGrid = xGridInd*aElem[1]+yGridInd*bElem[1]+zGridInd*cElem[1];
            zGrid = xGridInd*aElem[2]+yGridInd*bElem[2]+zGridInd*cElem[2];
	    xDiff = xGrid-x;
	    yDiff = yGrid-y;
	    zDiff = zGrid-z;
	    xTemp = xDiff*hmati[1]+yDiff*hmati[4]+zDiff*hmati[7];
            yTemp = xDiff*hmati[2]+yDiff*hmati[5]+zDiff*hmati[8];
            zTemp = xDiff*hmati[3]+yDiff*hmati[6]+zDiff*hmati[9];
	    xTemp2 = xTemp-NINT(xTemp);
            yTemp2 = yTemp-NINT(yTemp);
            zTemp2 = zTemp-NINT(zTemp);
            //debug[0] = xTemp2;
            //debug[1] = yTemp2;
            //debug[2] = zTemp2;

	    xDiff = xTemp2*hmat[1]+yTemp2*hmat[4]+zTemp2*hmat[7];
            yDiff = xTemp2*hmat[2]+yTemp2*hmat[5]+zTemp2*hmat[8];
            zDiff = xTemp2*hmat[3]+yTemp2*hmat[6]+zTemp2*hmat[9];
	    distSq = xDiff*xDiff+yDiff*yDiff+zDiff*zDiff;
	    if(distSq<cutOffSq){
	      //printf("distsq %lg cutoffsq %lg\n",distSq,cutOffSq);
	      //printf("grid inddddd %i %i %i\n",xGridInd,yGridInd,zGridInd);
	      //printf("xxxxxx %lg %lg %lg\n",debug[0],debug[1],debug[2]);
	      if(fftw3dFlag==0)gridInd = zGridInd*nkb*nka+yGridInd*nka+xGridInd;
	      else gridInd = xGridInd*nkb*nkc+yGridInd*nkc+zGridInd;
	      gridNlppMap[iAtom][numGridCount] = gridInd;
	      numGridCount += 1;  
	      if(numGridCount%100==0){
		gridNlppMap[iAtom] = (int *)realloc(gridNlppMap[iAtom],(numGridCount+100)*sizeof(int));
	      }//endif numGridCount
	    }//endif distSq
	  }//endfor k
	}//endfor j
      }//endfor i
    }//endif numLMax
    numGridNlppMap[iAtom] = numGridCount;
    /*
    for(iGrid=0;iGrid<numGridNlppMap[iAtom];iGrid++){
      printf("11111111mapp %i %i %i\n",iAtom,iGrid,gridNlppMap[iAtom][iGrid]);
    }
    */
    //printf("iAtom %i numGridNlppMap %i %i\n",iAtom,numGridNlppMap[iAtom],numGridCount);
  }//endfor iAtom

  //fflush(stdout);
  //exit(0);

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void testOverlap(CP *cp, CLASS *class, GENERAL_DATA *generalData)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Test whether nlpp grid overlaps                                       */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  ATOMMAPS *atommaps = &(class->atommaps);
  PSEUDO *pseudo = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CELL *cell = &(generalData->cell);
  CPEWALD *cpewald = &(cp->cpewald);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS *clatoms_pos = &(class->clatoms_pos[1]);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sparse = &(cp->cp_sclr_fft_pkg3d_sparse);
  COMMUNICATE *commCP = &(cp->communicate);
  
  int iAtom,jAtom,iGrid;
  int numAtomTot = clatoms_info->natm_tot;
  int realSparseOpt = cpewald->realSparseOpt;
  int atomType1,atomType2;
  int overlapFlag = 0;
  int countOverlap = 0;
  int countTotal = 0;
  int numGrid;
  int numGridAll;
  //int numGridAll = (cp_para_fft_pkg3d_lg->nfft)/2;
  int numGridNlppAll;
  int myidState = commCP->myid_state;
  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *atomNbhdListNum;
  int *gridNlppInd;
  int *numGridNlppMap = pseudoReal->numGridNlppMap;;
  int **gridNlppMap = pseudoReal->gridNlppMap;
  int **atomNbhdList;
  int **gridNlppMapInvTemp; //temperate list, include all grids
  int **gridNlppMapInv; //final list, include grids in gridNlppInd
  MPI_Comm commStates   =    commCP->comm_states;

  double x1,y1,z1,x2,y2,z2;
  double dx,dy,dz,dx2,dy2,dz2;
  double r2;
  double rcut1,rcut2;
  double rtot;

  double *ppRealCut = pseudoReal->ppRealCut;
  double *hmati = cell->hmati;
  double *hmat = cell->hmat;
  double *xList = clatoms_pos->x;
  double *yList = clatoms_pos->y;
  double *zList = clatoms_pos->z; 
  
  if(realSparseOpt==0){
    numGridAll = (cp_para_fft_pkg3d_lg->nfft)/2;
  }
  else{
    numGridAll = (cp_sclr_fft_pkg3d_sparse->nfft)/2;
  }

  gridNlppMapInvTemp = (int**)cmalloc(numGridAll*sizeof(int*));
  for(iGrid=0;iGrid<numGridAll;iGrid++){
    gridNlppMapInvTemp[iGrid] = (int*)cmalloc(sizeof(int));
    gridNlppMapInvTemp[iGrid][0] = 0;
  }

  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    atomType1 = iAtomAtomType[iAtom+1]-1;
    rcut1 = ppRealCut[atomType1];
    x1 = xList[iAtom+1];
    y1 = yList[iAtom+1];
    z1 = zList[iAtom+1];
    for(jAtom=iAtom+1;jAtom<numAtomTot;jAtom++){
      atomType2 = iAtomAtomType[jAtom+1]-1;
      rcut2 = ppRealCut[atomType2];
      x2 = xList[jAtom+1];
      y2 = yList[jAtom+1];
      z2 = zList[jAtom+1];
      dx = x2-x1;
      dy = y2-y1;
      dz = z2-z1;
      dx2 = dx*hmati[1]+dy*hmati[4]+dz*hmati[7];
      dy2 = dx*hmati[2]+dy*hmati[5]+dz*hmati[8];
      dz2 = dx*hmati[3]+dy*hmati[6]+dz*hmati[9];
      dx = dx2-NINT(dx2);
      dy = dy2-NINT(dy2);
      dz = dz2-NINT(dz2);
      dx2 = dx*hmat[1]+dy*hmat[4]+dz*hmat[7];
      dy2 = dx*hmat[2]+dy*hmat[5]+dz*hmat[8];
      dz2 = dx*hmat[3]+dy*hmat[6]+dz*hmat[9];
      r2 = dx2*dx2+dy2*dy2+dz2*dz2;
      rtot = rcut1+rcut2;
      if(rtot*rtot>=r2)overlapFlag = 1;
    }//endfor jAtom
  }//endfor iAtom

  /*
  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    printf("iAtom %i ",iAtom);
    for(jAtom=0;jAtom<atomNbhdListNum[iAtom];jAtom++){
      printf("%i ",atomNbhdList[iAtom][jAtom]);
    }
    printf("\n");
  }
  */

  if(overlapFlag==1&&myidState==0){
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    printf("Non-local pseudopotential regions are overlapped.\n");
    printf("Be careful when you use multithread.\n");
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }
  pseudoReal->overlapFlag;
  

  //Barrier(commStates);
  int gridIndex;
  int numUpdate;
  int count;
  int i,j,k;
  numGridNlppAll = 0;
  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    numGrid = numGridNlppMap[iAtom];
    for(iGrid=0;iGrid<numGrid;iGrid++){
      gridIndex = gridNlppMap[iAtom][iGrid];
      if(gridNlppMapInvTemp[gridIndex][0]==0){
	numGridNlppAll += 1;
      }
      gridNlppMapInvTemp[gridIndex][0] += 1;
      numUpdate = gridNlppMapInvTemp[gridIndex][0];
      gridNlppMapInvTemp[gridIndex] = (int*)crealloc(gridNlppMapInvTemp[gridIndex],(2*numUpdate+1)*sizeof(int));
      gridNlppMapInvTemp[gridIndex][2*numUpdate-1] = iAtom;
      gridNlppMapInvTemp[gridIndex][2*numUpdate] = iGrid;
    }
  }

  pseudoReal->numGridNlppAll = numGridNlppAll;
  pseudoReal->gridNlppInd = (int*)cmalloc(numGridNlppAll*sizeof(int));
  pseudoReal->gridNlppMapInv = (int**)cmalloc(numGridNlppAll*sizeof(int*));
  for(iGrid=0;iGrid<numGridNlppAll;iGrid++){
    pseudoReal->gridNlppMapInv[iGrid] = NULL;
  }
  gridNlppInd = pseudoReal->gridNlppInd;
  gridNlppMapInv = pseudoReal->gridNlppMapInv;
  count = 0;
  for(iGrid=0;iGrid<numGridAll;iGrid++){
    numUpdate = gridNlppMapInvTemp[iGrid][0];
    if(numUpdate>0){
      gridNlppInd[count] = iGrid;
      gridNlppMapInv[count] = (int*)cmalloc((2*numUpdate+1)*sizeof(int));
      memcpy(gridNlppMapInv[count],gridNlppMapInvTemp[iGrid],(2*numUpdate+1)*sizeof(int));
      /*
      printf("iGrid %i girdInd %i ",count,iGrid);
      for(i=0;i<2*numUpdate+1;i++){
	printf("%i ",gridNlppMapInv[count][i]);
      }
      printf("\n");
      */
      count += 1;
    }
  } 

  for(iGrid=0;iGrid<numGridAll;iGrid++){
    free(gridNlppMapInvTemp[iGrid]);
  }
  free(gridNlppMapInvTemp);

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

