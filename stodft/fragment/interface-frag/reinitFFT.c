/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: reinitFFT.c	                            */
/*                                                                          */
/*  This routine re-initialize the FFT after updating nuclei positions.	    */
/*  We need to change all malloc to realloc.				    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_interface_frag_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void reInitFFT(GENERAL_DATA *generalData,CLASS *class,CP *cp,
               GENERAL_DATA *generalDataMini,CLASS *classMini,CP *cpMini,
               int ip_now)
/*========================================================================*/
/*             Begin Routine                                              */
/*************************************************************************/
{/*Begin subprogram: */
/*========================================================================*/
/*             Local variable declarations                                */




/*------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


