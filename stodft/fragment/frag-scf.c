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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void fragScf(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                    CP *cp,ANALYSIS *analysis,int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo	 = stodftInfo->fragInfo;
  GENERAL_DATA *generalDataMini;
  CP *cpMini;
  CLASS *classMini;
  ANALYSIS *analysisMini;
  BONDED *bondedMini;
  
  int numFragProc = fragInfo->numFragProc;
  int iFrag;
  


  for(iFrag=0;iFrag<numFragProc;iFrag++){
/*======================================================================*/
/* I) Allocate Mini Structures                                          */
    generalDataMini = (GENERAL_DATA*)cmalloc(sizeof(GENERAL_DATA));
    cpMini = (CP*)cmalloc(sizeof(CP));
    classMini = (CLASS*)cmalloc(sizeof(CLASS));
    analysisMini = (ANALYSIS*)cmalloc(sizeof(ANALYSIS));
    bondedMini = (BONDED*)cmalloc(sizeof(BONDED));

/*======================================================================*/
/* I) Initialize Fragment SCF					        */
    parseFrag(class,bonded,general_data,cp,analysis,classMini,bondedMini,
           generalDataMini,cpMini,analysisMini);


/*======================================================================*/
/* II) SCF LOOP					                        */


/*======================================================================*/
/* II) Transfer Data and Free Memory				        */

  }
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/











