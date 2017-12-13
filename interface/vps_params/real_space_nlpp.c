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
void controlNlppRealSpline(CP *cp,CLASS *class,GENERAL_DATA *generalData,
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
  CELL *cell = &(general_data->cell);
  VPS_FILE *vpsFile = pseudo->vps_file;

  int iType,iRad,iAng,rGrid;
  int smoothOpt	    = pseudoReal->smoothOpt;
  int numAtomType   = atommaps->natm_typ;
  int countRad;
  int numR,angNow,numRadTot;
  int numGridRadSmooth;
  int *numRadGridSpline,numRadGridSplineTot;

  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *numLMax,*numRadMax;
  int *nrad_0 = pseudo->nrad_0;
  int *nrad_1 = pseudo->nrad_1;
  int *nrad_2 = pseudo->nrad_2;
  int *nrad_3 = pseudo->nrad_3;
  int *nrad_4 = pseudo->nrad_4;
  int *ivpsLabel = pseudo->ivps_label;

  int **atomLRadNum,**atomRadMap;

  double radCutRatio = pseudoReal->radCutRatio;
  double rCutoffMax;
  double junk1,junk2;
  double z1,z2,alpha1,alpha2;
  double zPol,gamma;
  double dr,rMax,rMin,r;
  double aLength,bLength,cLength;
  double a[3],b[3],c[3];
  double *vLoc,*vNl,*phiNl,*vNlSmooth;
  double *vNlTot;
  double *hmat = cell->hmat;
  double *ppRealCut;
  double *vpsNormList;
  
  FILE *fvps;
/*==========================================================================*/
/* I) Initialize radial function and angular channel                        */
  
  pseudoReal->numLMax = (int*)cmalloc(numAtomType*sizeof(int));  
  pseudoReal->numRadMax = (int*)cmalloc(numAtomType*sizeof(int));
  pseudoReal->atomLRadNum = (int**)cmalloc(numAtomType*sizeof(int*));
  pseudoReal->atomRadMap = (int**)cmalloc(numAtomType*sizeof(int*));
  pseudoReal->ppRealCut = (double*)cmalloc(numAtomType*sizeof(double));
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
    atomRadMap[iType] = NULL;
    if(numLMax[iType]>0){
      atomRadMap[iType] = (int*)cmalloc(numRadMax[iType]*sizeof(int));
      for(iRad=0;iRad<numRadMax[iType];iRad++){
        atomRadMap[iType][iRad] = countRad+iRad;
      }
    }
    countRad += numRadMax[iType];
  }
  numRadTot = countRad;
  pseudoReal->numRadTot = numRadTot;

  pseudoReal->vpsNormList = (double*)cmalloc(numRadTot*sizeof(double));
  

/*==========================================================================*/
/* I) Read and smooth Radial function			                    */
  
  // 1. Calcualte the minimum surface distance of the box
  a[0] = hmat[1];a[1] = hmat[2];a[2] = hmat[3];
  b[0] = hmat[4];b[1] = hmat[5];b[2] = hmat[6];
  c[0] = hmat[7];c[1] = hmat[8];c[2] = hmat[9];
  aLength = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  bLength = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  cLength = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
 
  // 2. Read the radial functions and determine the cutoff 
  for(iType=0;iType<numAtomType;iType++){
    fvps = fopen(vpsFile[iType+1].name,"r");
    rCutoffMax = -100000.0;
    if(ivpsLabel[iType+1]==1){// KB
      fscanf(fvps,"%i %lg %i\n",&numR,&rMax,&angNow);
      fscanf(fvps,"%lg %lg %lg %lg\n",&z1,&alpha1,&z2,&alpha2);
      fscanf(fvps,"%lg %lg\n",&zPol,&gamma);
      vLoc = (double*)cmalloc((numR+1)*sizeof(double));
      vNl = (double*)cmalloc((numR*numRadMax[iType]+1)*sizeof(double));
      phiNl = (double*)cmalloc((numR*numRadMax[iType]+1)*sizeof(double));
      dr = rMax/(double)numR;
      // get the non-local part
      for(iAng=0;iAng<angNow;iAng++){
	for(rGrid=0;rGrid<numR;rGrid++){
	  fscanf(fvps,"%lg %lg\n",&vNl[iAng*numR+rGrid+1],phiNl[iAng*numR+rGrid+1]);
	}//endfor rGrid
      }//endfor iAng
      // get the local potential
      for(rGrid=0;rGrid<numR;rGrid++){
	fscanf(fvps,"%lg %lg\n",&vLoc[rGrid+1],&junk1);
      }//endfor rGrid
      // Substract the nonlocal part from local part
      for(iAng=0;iAng<angNow;iAng++){
	for(rGrid=0;rGrid<numR;rGrid++){
	  vNl[iAng*numR+rGrid+1] = (vNl[iAng*numR+rGrid+1]-vLoc[rGrid+1])*phiNl[iAng*numR+rGrid+1];
	}//endfor rGrid
      }//endfor iAng
      for(iAng=0;iAng<angNow;iAng++){
	rGrid = 0;
	while(vNl[iAng*numR+rGrid+1]>1.0e-10&&rGrid<numR){//double check 1.0e-10
	  rGrid += 1;
	}
	if(rGrid==numR){
	  printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	  printf("Your pseudo-potential does not decay enough! Please\n");
	  printf("choose bigger radius cutoff!\n");
	  printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	  fflush(stdout);
	  exit(0);
	}
	r = rGrid*dr*radCutRatio;
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
    
    numGridRadSmooth = (int)(rCutoffMax/dr)+1;
    vNlSmooth = (double*)cmalloc(numGridRadSmooth*numRadMax[iType]*sizeof(double));
    for(iRad=0;iRad<numRadMax[iType];iRad++){
      switch(smoothOpt){
	case 1:
	  nlppSmoothKS(pseudo,vNl,rCutoffMax,iRad,vNlSmooth,numGridRadSmooth,
		       rMax,numR);
	  break;
	case 2:
	  nlppSmoothRoi(pseudo,vNl,rCutoffMax,iRad);
	  break;
      }//endswitch    
    }//endfor iRad
  }//endfor iType

/*==========================================================================*/
/* I) Allocate real space spline                                            */

  // 1. Calculate the real space radius grid number
  
  
/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void nlppSmoothKS(PSEUDO *pseudo,double *vNl,double rCutoffMax,int iRad,
		  double *vNlSmooth,int numGridRadSmooth,double rMax,int numR,
		  int l)
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
  int numRCutoff;

  double gMaxSm = pseudoReal->gMaxSm;
  double gMaxLg = pseudoReal->gMaxLg;
  double dg = pseudoRal->dg;
  double dr = rMax/((double)numR);
  double *gGrid;

  double *vNlG,*dvNlG;
  double *
  
  // initialize
  numGSm = (int)(gMaxSm/dg)+1;
  numGLg = (int)(gMaxLg/dg)+1;

  vNlG = (double*)cmalloc((numGLg)*sizeof(double)); //fv_rphi
  dvNlG = (double*)cmalloc((numGlg)*sizeof(double)); //dfv_rphi
  r = (double*)cmalloc((numR)*sizeof(double));
  gGrid = (double*)cmalloc((numGLg)*sizeof(double));

  for(ig=0;ig<numGLg;ig++)gGrid[ig] = ig*dg;
  for(ir=0;ir<numR;ir++)r[ir] = ir*dr;
 

  // Bessel transform to g space from g=0 to gMaxSm
  bessTransform(vNl,numR,dr,l,vNlG,numGSm,dg);

  // Optimize g space coeffcient from gMaxSm to gMaxLg

  optGCoeff(pseudoReal,numGLg,numGSm,numR,dr,dg,rMax,l,vNlG);

  // Inverse Bassel transform to r space

  bessTransform(vNlG,numGLg,dg,l,vNlSmooth,numRCutoff,dr);  
  // we may need to resacle the potential. Let's check this


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
  for(jGrid=1;jGrid<numIn;jGrid++){
    x = jGrid*dx;
    funOut[0] += x*funIn[jGrid]*dx;
  }

  switch(l){
    case 0:
      for(iGrid=1;iGrid<numOut;iGrid++){
	funOut[iGrid] = 0.0;
	y = iGrid*dy;
	for(jGrid=1;jGrid<numIn,jGrid++){
	  x = jGrid*dx;
	  arg = x*y;
	  funOut[iGrid] += j0(arg)*x*funIn*dx;
	}//endfor jGrid
      }//endfor iGrid
      break;
    case 1:
      for(iGrid=1;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = iGrid*dy;
        for(jGrid=1;jGrid<numIn,jGrid++){
          x = jGrid*dx;
          arg = x*y;
          funOut[iGrid] += j1(arg)*x*funIn*dx;
        }//endfor jGrid
      }//endfor iGrid
      break;
    case 2:
      for(iGrid=1;iGrid<numOut;iGrid++){
        funOut[iGrid] = 0.0;
        y = iGrid*dy;
        for(jGrid=1;jGrid<numIn,jGrid++){
          x = jGrid*dx;
          arg = x*y;
          funOut[iGrid] += j2(arg)*x*funIn*dx;
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
		double dr,double dg,double rCutoff, int l,double *vNlG)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Bessel transform */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int iGridG,jGridG,kGridR;
  int numRCutoff = rCutoff/dr;
  int numGSolve = numGLg-numGSm;
  int gridShift;

  double qi,qj;
  double r,arg1,arg2;

  double *A;
  double *B;
  double *subA;
  //Construct the A matrix

  A = (double*)calloc(numGLg*numGLg*sizeof(double));
  B = (double*)calloc(numGSolve*sizeof(double));
  subA = (double*)calloc(numGSolve*sizeof(double));
  
  for(iGridG=0;iGridG<numGLg;iGridG++){
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
      subA[iGridG*numGSolve+jGridG] = -A[(iGridG+numGSm)+jGridG+numGSm]*dg;
    }//endfor jGridG
    subA[iGridG*numGSolve+iGridG] += 0.5*M_PI*qi*qi;
  }//endfor iGridG

  // Solving the linear system
  dsysvWrapper(subA,B,numGSolve);

  for(iGridG=0;iGridG<numGSolve;iGridG++){
    vNlG[iGridG+numGSm] = B[iGridG];
  }

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/




