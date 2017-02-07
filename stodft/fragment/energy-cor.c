/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: energy-cor.c                                 */
/*                                                                          */
/* This routine calculate the kinectic energy, non-local pseudo-potential   */
/* energy and nuclei force correction.					    */
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
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_frag_local.h"
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void energyCorrect(CP *cp,GENERAL_DATA *general_data,CLASS *class,int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  STODFTINFO *stodftInfo = cp->stodftInfo;
  FRAGINFO *fragInfo = stodftInfo->fragInfo;
  CPOPTS *cpOpts = &(cp->cpopts);
  CPCOEFFS_INFO *cpcoeffs_info = &(cp->cpcoeffs_info);
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);
  COMMUNICATE *commCP = &(cp->communicate);

  double **keMatrixUp = fragInfo->keMatrixUp;
  double **keMatrixDn = fragInfo->keMatrixDn;

/*======================================================================*/
/* I) Kinetic energy	                                                */

/*======================================================================*/
/* II) Non-local pseudo potential energy                                */


/*======================================================================*/
/* III) Non-local pseudo potential force                                */


/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


