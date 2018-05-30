/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: frag-scf.c                                   */
/*                                                                          */
/* Driver Routine for fragment SCF calculation                              */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
//#include "../proto_defs/proto_stodft_local.h"
#include "../proto_defs/proto_integrate_cpmin_entry.h"
#include "../proto_defs/proto_interface_frag_entry.h"
#include "../proto_defs/proto_frag_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void fragScf(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
            CP *cp,ANALYSIS *analysis,GENERAL_DATA *generalDataMini,
	    CP *cpMini,CLASS *classMini,ANALYSIS *analysisMini,BONDED *bondedMini,
	    int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo	 = stodftInfo->fragInfo;
  COMMUNICATE   *communicate = &(cp->communicate);
  
  int numFragProc = fragInfo->numFragProc;
  int iFrag;
  //debug
  char fileNameFragMO[100];
  FILE *fileFragMO;
  int iState,iGrid;
  int ncoef,icoef;
  int myidState		= communicate->myid_state;
  int numProcStates     = communicate->np_states;
  int fragOpt		= stodftInfo->fragOpt;
  MPI_Comm commStates   = communicate->comm_states;

  int *numGridFragProc = fragInfo->numGridFragProc;
  int *numElecUpFragProc = fragInfo->numElecUpFragProc;

  //debug
  
    
  sprintf(fileNameFragMO,"frag-MO-%i",myidState);
  if(numFragProc>0){
    fileFragMO = fopen(fileNameFragMO,"r");
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      //printf("%p\n",fileFragMO);
      ncoef = cpMini[iFrag].cpcoeffs_info.ncoef;
      for(iState=0;iState<numElecUpFragProc[iFrag];iState++){
        for(icoef=1;icoef<=ncoef;icoef++){
          fscanf(fileFragMO,"%lg",&(cpMini[iFrag].cpcoeffs_pos[1].cre_up[iState*ncoef+icoef]));
          fscanf(fileFragMO,"%lg",&(cpMini[iFrag].cpcoeffs_pos[1].cim_up[iState*ncoef+icoef]));
        }
      }
    }
    fclose(fileFragMO);
  }
  
  
  for(iFrag=0;iFrag<numFragProc;iFrag++){
/*======================================================================*/
/* I) Allocate Mini Structures                                          */

/*======================================================================*/
/* I) Initialize Fragment SCF					        */
    
    /*
    parseFrag(class,bonded,general_data,cp,analysis,&classMini[iFrag],&bondedMini[iFrag],
           &generalDataMini[iFrag],&cpMini[iFrag],&analysisMini[iFrag]);
    */

/*======================================================================*/
/* II) SCF LOOP					                        */
  
    controlCpMinFrag(&classMini[iFrag],&bondedMini[iFrag],&generalDataMini[iFrag],
                     &cpMini[iFrag],&analysisMini[iFrag]);      
    
  }//endfor iFrag
  
  
  if(numProcStates>1)Barrier(commStates);
  
  /* 
  sprintf(fileNameFragMO,"frag-MO-%i",myidState);
  fileFragMO = fopen(fileNameFragMO,"w");
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    ncoef = cpMini[iFrag].cpcoeffs_info.ncoef;
    for(iState=0;iState<numElecUpFragProc[iFrag];iState++){
      for(icoef=1;icoef<=ncoef;icoef++){
	fprintf(fileFragMO,"%.16lg %.16lg\n",cpMini[iFrag].cpcoeffs_pos[1].cre_up[iState*ncoef+icoef],
	       cpMini[iFrag].cpcoeffs_pos[1].cim_up[iState*ncoef+icoef]);
      }
    }
  }
  fclose(fileFragMO);    
  */
  
  if(numProcStates>1)Barrier(commStates);
  //exit(0);

/*======================================================================*/
/* II) Transfer Data and Free Memory                                    */

    
  sprintf(fileNameFragMO,"frag-MO-%i",myidState);
  if(numFragProc>0){
    fileFragMO = fopen(fileNameFragMO,"r");
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      //printf("%p\n",fileFragMO);
      ncoef = cpMini[iFrag].cpcoeffs_info.ncoef;
      for(iState=0;iState<numElecUpFragProc[iFrag];iState++){
	for(icoef=1;icoef<=ncoef;icoef++){
	  fscanf(fileFragMO,"%lg",&(cpMini[iFrag].cpcoeffs_pos[1].cre_up[iState*ncoef+icoef]));
	  fscanf(fileFragMO,"%lg",&(cpMini[iFrag].cpcoeffs_pos[1].cim_up[iState*ncoef+icoef]));
	}
      }
    }
    fclose(fileFragMO);
  }
  

  if(fragOpt==1){
    projRhoMiniMol(cp,general_data,class,cpMini,generalDataMini,classMini,ip_now);
  }
  if(fragOpt==3){
    projRhoMiniUnitCell(cp,general_data,class,cpMini,generalDataMini,classMini,ip_now);
  }
  //fflush(stdout);
  //exit(0);
  energyCorrect(cpMini,generalDataMini,classMini,cp,class,ip_now);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/











