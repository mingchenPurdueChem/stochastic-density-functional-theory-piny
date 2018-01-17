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
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);

  int iType,iRad,iAng,rGrid,iAtom;
  int smoothOpt	    = pseudoReal->smoothOpt;
  int numAtomType   = atommaps->natm_typ;
  int numAtomTot    = clatoms_info->natm_tot;
  int nkf1 = cp_para_fft_pkg3d_lg->nkf1;
  int nkf2 = cp_para_fft_pkg3d_lg->nkf2;
  int nkf3 = cp_para_fft_pkg3d_lg->nkf3;
  int countRad;
  int numR,angNow,numRadTot;
  int countR = 0;
  int numGridTot = 0;

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
  double *vpsReal0,*vpsReal1,*vpsReal2,*vpsReal3;
  double **vNlSmooth,**dvNlSmooth;

  FILE *fvps;
/*==========================================================================*/
/* I) Initialize radial function and angular channel                        */
  
  pseudoReal->numLMax = (int*)cmalloc(numAtomType*sizeof(int));  
  pseudoReal->numRadMax = (int*)cmalloc(numAtomType*sizeof(int));
  pseudoReal->atomLRadNum = (int**)cmalloc(numAtomType*sizeof(int*));
  pseudoReal->atomRadMap = (int**)cmalloc(numAtomType*sizeof(int*));
  pseudoReal->ppRealCut = (double*)cmalloc(numAtomType*sizeof(double));
  vNlSmooth = (double**)cmalloc(numAtomType*sizeof(double*));
  dvNlSmooth = (double**)cmalloc(numAtomType*sizeof(double*));
  numGridRadSmooth = (int*)cmalloc(numAtomType*sizeof(int));
  numLMax = pseudoReal->numLMax;
  numRadMax = pseudoReal->numRadMax;
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

  printf("ggggggg %lg %lg %lg\n",gmaxTrueSm,gmaxTrueLg,gmaxTrueLgLg);
  pseudoReal->gMaxSm = gmaxTrueSm;
  pseudoReal->gMaxLg = gmaxTrueLg;
  //pseudoReal->gMaxLg = 3.0*gmaxTrueSm;
 
  // 2. Read the radial functions and determine the cutoff 
  for(iType=0;iType<numAtomType;iType++){
    fvps = fopen(vpsFile[iType+1].name,"r");
    rCutoffMax = -100000.0;
    if(ivpsLabel[iType+1]==1){// KB
      fscanf(fvps,"%i %lg %i\n",&numR,&rMax,&angNow);
      fscanf(fvps,"%lg %lg %lg %lg\n",&z1,&alpha1,&z2,&alpha2);
      fscanf(fvps,"%lg %lg\n",&zPol,&gamma);
      printf("%lg %lg %lg %lg %lg %lg\n",z1,alpha1,z2,alpha2,zPol,gamma);
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
      // get the non-local part
      for(iAng=0;iAng<angNow;iAng++){
	for(rGrid=0;rGrid<numR;rGrid++){
	  fscanf(fvps,"%lg %lg\n",&vNl[iAng*numR+rGrid],&phiNl[iAng*numR+rGrid]);
	}//endfor rGrid
      }//endfor iAng
      // get the local potential
      for(rGrid=0;rGrid<numR;rGrid++){
	fscanf(fvps,"%lg %lg\n",&vLoc[rGrid],&junk1);
      }//endfor rGrid
      // Substract the nonlocal part from local part
      for(iAng=0;iAng<angNow;iAng++){
	lMap[countR] = iAng;
	vpsNormList[countR] = 0.0;
	for(rGrid=0;rGrid<numR;rGrid++){
	  vNl[iAng*numR+rGrid] = (vNl[iAng*numR+rGrid]-vLoc[rGrid])*phiNl[iAng*numR+rGrid];
	  vpsNormList[countR] += vNl[iAng*numR+rGrid]*phiNl[iAng*numR+rGrid];
	  //printf("rGridddddd %lg %lg\n",rGrid*dr,vNl[iAng*numR+rGrid]);
	}//endfor rGrid
	vpsNormList[countR] *= dr;
	printf("countR %i vpsNormList %lg\n",countR,vpsNormList[countR]);
	vpsNormList[countR] = 1.0/vpsNormList[countR];
	countR += 1;
      }//endfor iAng
      for(iAng=0;iAng<angNow;iAng++){
	rGrid = numR-1;
	while(vNl[iAng*numR+rGrid]<1.0e-10&&rGrid>0){//double check 1.0e-10
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
       
      //free(vLoc);
      //free(vNl);
      //free(phiNl);
      ppRealCut[iType] = rCutoffMax;
    }//endif ivpsLabel
  // 3. Smooth the radius function
    
    numGridRadSmooth[iType] = (int)(rCutoffMax/dr)+1;
    vNlSmooth[iType] = (double*)cmalloc(numGridRadSmooth[iType]*numRadMax[iType]*sizeof(double));
    dvNlSmooth[iType] = (double*)cmalloc(numGridRadSmooth[iType]*numRadMax[iType]*sizeof(double));

    for(iRad=0;iRad<numRadMax[iType];iRad++){
      switch(smoothOpt){
	case 1:
	  nlppSmoothKS(pseudo,&vNl[iRad*numR],rCutoffMax,iRad,
		       &vNlSmooth[iType][iRad*numGridRadSmooth[iType]],
		       &dvNlSmooth[iType][iRad*numGridRadSmooth[iType]],
		       numGridRadSmooth[iType],rMax,numR,lMap[iRad]);
	  break;
	case 2:
	  nlppSmoothRoi(pseudo,&vNl[iRad*numR],rCutoffMax,iRad,
			&vNlSmooth[iType][iRad*numGridRadSmooth[iType]],
			&dvNlSmooth[iType][iRad*numGridRadSmooth[iType]],
			numGridRadSmooth[iType],rMax,numR,lMap[iRad]);
	  break;
      }//endswitch    
    }//endfor iRad
  }//endfor iType

  pseudoReal->dr = dr;
/*==========================================================================*/
/* III) Real space spline						    */

  // 1. Unified the real space cutoff, its easier to allocate
  int numInterpGrid = -1;
  
  int gridShift;
  double *rList;
  for(iType=0;iType<numAtomType;iType++){
    rGrid = (int)(ppRealCut[iType]/dr)+1;
    if(numInterpGrid<rGrid){
      numInterpGrid = rGrid;
    }
  }
  pseudoReal->numInterpGrid = numInterpGrid;
  rList = (double*)cmalloc((numInterpGrid+1)*sizeof(double));
  pseudoReal->vpsReal0 = (double*)calloc((numInterpGrid*numRadTot+1),sizeof(double));
  pseudoReal->vpsReal1 = (double*)calloc((numInterpGrid*numRadTot+1),sizeof(double));
  pseudoReal->vpsReal2 = (double*)calloc((numInterpGrid*numRadTot+1),sizeof(double));
  pseudoReal->vpsReal3 = (double*)calloc((numInterpGrid*numRadTot+1),sizeof(double));
  vpsReal0 = pseudoReal->vpsReal0;
  vpsReal1 = pseudoReal->vpsReal1;
  vpsReal2 = pseudoReal->vpsReal2;
  vpsReal3 = pseudoReal->vpsReal3;

  for(rGrid=0;rGrid<numInterpGrid;rGrid++)rList[rGrid+1] = rGrid*dr;
 
  gridShift = 0;
  for(iType=0;iType<numAtomType;iType++){
    for(iRad=0;iRad<numRadMax[iType];iRad++){
      memcpy(&vpsReal0[1+gridShift],&vNlSmooth[iType][iRad*numGridRadSmooth[iType]],
	     numGridRadSmooth[iType]*sizeof(double));
      gridShift += numInterpGrid;
    }
  }
    
  for(iRad=0;iRad<numRadTot;iRad++){
    gridShift = iRad*numInterpGrid;
    spline_fit(&(vpsReal0[gridShift]),&(vpsReal1[gridShift]),
               &(vpsReal2[gridShift]),&(vpsReal3[gridShift]),rList,numInterpGrid);  
  }
  // test spline
  /*
  double rtest;
  rMin = 0.0;
  double r0,h,pseudoTest;
  int gridIndTest;
  int interpInd,interpGridSt;
  for(iRad=0;iRad<numRadTot;iRad++){
    interpGridSt = iRad*numInterpGrid;
    for(rGrid=0;rGrid<numInterpGrid;rGrid++){
      rtest = (rGrid+0.5)*dr;
      gridIndTest = (int)((rtest-rMin)/dr)+1;
      r0 = (gridIndTest-1)*dr+rMin;
      h = rtest-r0;
      interpInd = interpGridSt+gridIndTest;
      pseudoTest = ((vpsReal3[interpInd]*h+vpsReal2[interpInd])*h+vpsReal1[interpInd])*h+vpsReal0[interpInd];
      printf("iiiiiirad %i %.8lg %.8lg\n",iRad,rtest,pseudoTest*rtest*2.0/M_PI);
    }
  }
  */

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
/* Real space nlpp, only used in filtering                   */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  PSEUDO *pseudo = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  ATOMMAPS *atommaps = &(class->atommaps);


  int numAtomTot = clatoms_info->natm_tot;
  int iAtom,iAng;
  int numGrid;
  int numPsuedoWfRe,numPsuedoWfIm;
  int numGridTotRe = 0;
  int numGridTotIm = 0;
  int atomType;
  int *numGridNlppMap;
  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *numLMax = pseudoReal->numLMax;
  int **atomLRadNum = pseudoReal->atomLRadNum;

/*==========================================================================*/
/* IV) Construct Real Space Grids Neighbourhood List                        */

  mapRealSpaceGrid(cp,class,generalData);
  //fflush(stdout);
  //exit(0);

/*==========================================================================*/
/* V) Construct Real Space V_nl*Phi_nl		                            */
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
  pseudoReal->vnlPhiAtomGridRe = (double*)cmalloc(numGridTotRe*sizeof(double));  
  pseudoReal->vnlPhiAtomGridIm = (double*)cmalloc(numGridTotIm*sizeof(double));

  calcPseudoWf(cp,class,generalData);
  //fflush(stdout);
  //exit(0);
  
/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void nlppSmoothKS(PSEUDO *pseudo,double *vNl,double rCutoffMax,int iRad,
		  double *vNlSmooth,double *dvNlSmooth,int numGridRadSmooth,
		  double rMax,int numR,int l)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Real space nlpp, only used in filtering                   */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);  

  int numGSm;
  int numGLg;
  int ig,ir;

  double gMaxSm = pseudoReal->gMaxSm;
  double gMaxLg = pseudoReal->gMaxLg;
  double dg = pseudoReal->dg;
  double dr = rMax/((double)numR);
  double *gGrid,*r;

  double *vNlG,*dvNlG;
  
  // initialize
  numGSm = (int)(gMaxSm/dg)+1;
  numGLg = (int)(gMaxLg/dg)+1;
  printf("gMaxSm %lg gMaxLg %lg\n",gMaxSm,gMaxLg);

  vNlG = (double*)calloc((numGLg),sizeof(double)); //fv_rphi
  dvNlG = (double*)cmalloc((numGLg)*sizeof(double)); //dfv_rphi
  r = (double*)cmalloc((numR)*sizeof(double));
  gGrid = (double*)cmalloc((numGLg)*sizeof(double));

  for(ig=0;ig<numGLg;ig++)gGrid[ig] = ig*dg;
  for(ir=0;ir<numR;ir++)r[ir] = ir*dr;
 

  printf("2222222 l %i\n",l);
  // Bessel transform to g space from g=0 to gMaxSm
  bessTransform(vNl,numR,dr,l,vNlG,numGSm,dg);
  //bessTransform(vNl,numR,dr,l,vNlG,numGLg,dg);
  
  /*
  for(ig=0;ig<numGLg;ig++){
    printf("111111111 vnlg %lg %lg\n",ig*dg,vNlG[ig]);
  }
  printf("numGridRadSmooth %i\n",numGridRadSmooth);
  */
  
  //fflush(stdout);
  //exit(0);

  // Optimize g space coeffcient from gMaxSm to gMaxLg

  optGCoeff(pseudoReal,numGLg,numGSm,numR,dr,dg,numGridRadSmooth,l,vNlG);
  
  /*
  for(ig=0;ig<numGLg;ig++){
    printf("111111 vNlGTrunc %lg %lg\n",ig*dg,vNlG[ig]);
  }
  */
  

  // Inverse Bassel transform to r space

  for(ig=0;ig<numGLg;ig++){
    vNlG[ig] *= ig*dg;
  }
  bessTransform(vNlG,numGLg,dg,l,vNlSmooth,numGridRadSmooth,dr);  
  bessTransformGrad(vNlG,numGLg,dg,l,dvNlSmooth,numGridRadSmooth,dr);
  for(ir=0;ir<numGridRadSmooth;ir++){
    vNlSmooth[ir] *= 2.0/M_PI;
    dvNlSmooth[ir] *= 2.0/M_PI;
  }
  /*
  for(ir=0;ir<numGridRadSmooth;ir++){
    printf("111111 vnlr %lg %lg\n",ir*dr,vNlSmooth[ir]*ir*dr*2.0/M_PI);
  }
  */
  //fflush(stdout);
  //exit(0);

  // we may need to resacle the potential. Let's check this


/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void nlppSmoothRoi(PSEUDO *pseudo,double *vNl,double rCutoffMax,int iRad,
                  double *vNlSmooth,double *dvNlSmooth,int numGridRadSmooth,
		  double rMax,int numR,int l)
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
		   int numOut,double dy)
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

  funOut[0] = 0.0;
  if(l==0){
    for(jGrid=1;jGrid<numIn;jGrid++){
      x = jGrid*dx;
      funOut[0] += x*funIn[jGrid]*dx;
    }
  }

  switch(l){
    case 0:
      for(iGrid=1;iGrid<numOut;iGrid++){
	funOut[iGrid] = 0.0;
	y = iGrid*dy;
	for(jGrid=1;jGrid<numIn;jGrid++){
	  x = jGrid*dx;
	  arg = x*y;
	  funOut[iGrid] += j0(arg)*x*funIn[jGrid]*dx;
	}//endfor jGrid
      }//endfor iGrid
      break;
    case 1:
      for(iGrid=1;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = iGrid*dy;
        for(jGrid=1;jGrid<numIn;jGrid++){
          x = jGrid*dx;
          arg = x*y;
          funOut[iGrid] += j1(arg)*x*funIn[jGrid]*dx;
        }//endfor jGrid
      }//endfor iGrid
      break;
    case 2:
      for(iGrid=1;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = iGrid*dy;
        for(jGrid=1;jGrid<numIn;jGrid++){
          x = jGrid*dx;
          arg = x*y;
          funOut[iGrid] += j2(arg)*x*funIn[jGrid]*dx;
        }//endfor jGrid
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
  
  //iGrid=0
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
  /*
  for(iGridG=0;iGridG<numGSolve;iGridG++){
    printf("BBBBBBBBBBB %i %lg\n",iGridG,B[iGridG]);
  }
  */
  for(iGridG=0;iGridG<numGSolve;iGridG++){
    qi = (iGridG+numGSm)*dg;
    for(jGridG=0;jGridG<numGSolve;jGridG++){
      subA[iGridG*numGSolve+jGridG] = -A[(iGridG+numGSm)*numGLg+jGridG+numGSm]*dg;
    }//endfor jGridG
    subA[iGridG*numGSolve+iGridG] += 0.5*M_PI*qi*qi;
  }//endfor iGridG
  /*
  for(iGridG=0;iGridG<numGSolve;iGridG++){
    for(jGridG=0;jGridG<numGSolve;jGridG++){
      printf("subAAAAAAA %i %i %lg\n",iGridG,jGridG,subA[iGridG*numGSolve+jGridG]);
    }
    printf("diag subAAAA %i %lg\n",iGridG,subA[iGridG*numGSolve+iGridG]);
  }
  */
  // Solving the linear system
  dsysvWrapper(subA,B,numGSolve);
  /*
  for(iGridG=0;iGridG<numGSolve;iGridG++){
    printf("solutionnnn %i %lg\n",iGridG,B[iGridG]);
  }
  */
  /*
  double *BTest = (double*)calloc(numGSolve,sizeof(double));
  for(iGridG=0;iGridG<numGSolve;iGridG++){
    for(jGridG=0;jGridG<numGSolve;jGridG++){
      BTest[iGridG] += subA[iGridG*numGSolve+jGridG]*B[jGridG];
    }
  }
  for(iGridG=0;iGridG<numGSolve;iGridG++){
    printf("recoverrrrrr %i %lg\n",iGridG,BTest[iGridG]);
  }
  */

  for(iGridG=0;iGridG<numGSolve;iGridG++){
    vNlG[iGridG+numGSm] = B[iGridG];
  }

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void bessTransformGrad(double *funIn,int numIn,double dx,int l,double *funOut,
                   int numOut,double dy)
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

  switch(l){
    case 0:
      for(iGrid=1;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = iGrid*dy;
        for(jGrid=1;jGrid<numIn;jGrid++){
          x = jGrid*dx;
          arg = x*y;
          funOut[iGrid] += dj0(arg)*x*x*funIn[jGrid]*dx;
        }//endfor jGrid
      }//endfor iGrid
      break;
    case 1:
      for(iGrid=1;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = iGrid*dy;
        for(jGrid=1;jGrid<numIn;jGrid++){
          x = jGrid*dx;
          arg = x*y;
          funOut[iGrid] += dj1(arg)*x*x*funIn[jGrid]*dx;
        }//endfor jGrid
      }//endfor iGrid
      break;
    case 2:
      for(iGrid=1;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = iGrid*dy;
        for(jGrid=1;jGrid<numIn;jGrid++){
          x = jGrid*dx;
          arg = x*y;
          funOut[iGrid] += dj2(arg)*x*x*funIn[jGrid]*dx;
        }//endfor jGrid
      }//endfor iGrid
      break;
  }//endswitch

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
  PSEUDO *pseudo = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  CELL *cell = &(generalData->cell);
  VPS_FILE *vpsFile = pseudo->vps_file;
  CPEWALD *cpewald = &(cp->cpewald);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS *clatoms_pos = &(class->clatoms_pos[1]);

  int numAtomTot = clatoms_info->natm_tot;
  int atomType,l;
  int iAtom,i,j,k;
  int numGridCount;
  int nkc = cp_para_fft_pkg3d_lg->nkf3;
  int nkb = cp_para_fft_pkg3d_lg->nkf2;
  int nka = cp_para_fft_pkg3d_lg->nkf1;
  int aNum,bNum,cNum;
  int xInd,yInd,zInd,xGridInd,yGridInd,zGridInd;
  int gridInd;
  int numInterpGrid = pseudoReal->numInterpGrid;  

  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *numLMax = pseudoReal->numLMax;
  int *numGridNlppMap;
  int **gridNlppMap;

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
  
  pseudoReal->numGridNlppMap = (int*)cmalloc(numAtomTot*sizeof(int));
  pseudoReal->gridNlppMap = (int**)cmalloc(numAtomTot*sizeof(int*));
  numGridNlppMap = pseudoReal->numGridNlppMap;
  gridNlppMap = pseudoReal->gridNlppMap;
  
  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    atomType = iAtomAtomType[iAtom+1]-1;
    gridNlppMap[iAtom] = NULL;
    numGridNlppMap[iAtom] = 0;
    if(numLMax[atomType]>0){
      numGridCount = 0;
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
	      gridInd = zGridInd*nkb*nka+yGridInd*nka+xGridInd;
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
    //printf("iAtom %i numGridNlppMap %i %i\n",iAtom,numGridNlppMap[iAtom],gridNlppMap[iAtom][0]);
  }//endfor iAtom

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

