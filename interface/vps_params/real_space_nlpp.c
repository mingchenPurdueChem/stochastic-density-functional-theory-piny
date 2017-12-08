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
  int numAtomType = atommaps->natm_typ;
  int countRad;
  int numR,angNow

  int *iAtomAtomType = atommaps->iatm_atm_typ;
  int *numLMax,*numRadMax;
  int *nrad_0 = pseudo->nrad_0;
  int *nrad_1 = pseudo->nrad_1;
  int *nrad_2 = pseudo->nrad_2;
  int *nrad_3 = pseudo->nrad_3;
  int *nrad_4 = pseudo->nrad_4;
  int *ivpsLabel = pseudo->ivps_label;

  int **atomLRadNum,**atomRadMap;

  double junk1,junk2;
  double z1,z2,alpha1,alpha2;
  double zPol,gamma;
  double *vLoc,*vNl,*phiNl;
  
  FILE *fvps;
/*==========================================================================*/
/* I) Initialize radial function and angular channel                        */
  
  pseudoReal->numLMax = (int*)cmalloc(numAtomType*sizeof(int));  
  pseudoReal->numRadMax = (int*)cmalloc(numAtomType*sizeof(int));
  pseudoReal->atomLRadNum = (int**)cmalloc(numAtomType*sizeof(int*));
  pseudoReal->atomRadMap = (int**)cmalloc(numAtomType*sizeof(int*));
  numLMax = pseudoReal->numLMax;
  numRadMax = pseudoReal->numRadMax;
  atomLRadNum = pseudoReal->atomLRadNum;
  atomRadMap = pseudoReal->atomRadMap;

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

/*==========================================================================*/
/* I) Read and smooth Radial function			                    */
 

  /* Read the radial functions and determine the cutoff */
  for(iType=0;iType<numAtomType;iType++){
    fvps = fopen(vpsFile[iType+1].name,"r");
    if(ivpsLabel[iType+1]==1){// KB
      fscanf(fvps,"%i %lg %i\n",&numR,&rMax,&angNow);
      fscanf(fvps,"%lg %lg %lg %lg\n",&z1,&alpha1,&z2,&alpha2);
      fscanf(fvps,"%lg %lg\n",&zPol,&gamma);
      vLoc = (double*)cmalloc((numR+1)*sizeof(double));
      vNl = (double*)cmalloc((numR*angNow+1)*sizeof(double));
      phiNl = (double*)cmalloc((numR*angNow+1)*sizeof(double));
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
	  vNl[iAng*numR+rGrid+1] -= vLoc[rGrid+1];
	}
      }
      
      //free(vLoc);
      //free(vNl);
      //free(phiNl);
    }
    
  }

  pseudoReal->vpsReal0 = (double*)cmalloc();
  
  
/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/






