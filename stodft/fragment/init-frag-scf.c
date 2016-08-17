/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: init-frag.c                                    */
/*                                                                          */
/* This routine initialize the fragmentation calculation (scf part).        */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_parse_entry.h"
#include "../proto_defs/proto_parse_local.h"
#include "../proto_defs/proto_sim_params_entry.h"
#include "../proto_defs/proto_mol_params_entry.h"
#include "../proto_defs/proto_intra_params_entry.h"
#include "../proto_defs/proto_coords_entry.h"
#include "../proto_defs/proto_coords_local.h"
#include "../proto_defs/proto_cp_ewald_entry.h"
#include "../proto_defs/proto_surf_params_entry.h"
#include "../proto_defs/proto_inter_params_entry.h"
#include "../proto_defs/proto_inter_params_local.h"
#include "../proto_defs/proto_vps_params_entry.h"
#include "../proto_defs/proto_lists_entry.h"
#include "../proto_defs/proto_scratch_entry.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_vel_sampl_class_entry.h"
#include "../proto_defs/proto_vel_sampl_cp_entry.h"
#include "../proto_defs/proto_vel_sampl_cp_local.h"
#include "../proto_defs/proto_coords_cp_entry.h"
#include "../proto_defs/proto_coords_cp_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_entry.h"
#include "../proto_defs/proto_communicate_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_integrate_cpmin_entry.h"
#include "../proto_defs/proto_energy_cp_entry.h"
#include "../proto_defs/proto_stodft_local.h"
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initFragSCF(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp,
                 CLASS *classMini,BONDED *bondedMini,ENERAL_DATA *generalDataMini,
                 CP *cpMini)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Fragmentation initialization SCF part                                 */
/* This should be done after each time updating atomic coordinates       */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS  *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  ATOMMAPS     *atommaps     = &(class->atommaps);
  CELL         *cell         = &(general_data->cell);
  FOR_SCR      *for_scr      = &(class->for_scr);
  EWD_SCR      *ewd_scr      = &(class->ewd_scr);
  PTENS        *ptens        = &(general_data->ptens);

  CPOPTS        *cpopts           = &(cp->cpopts);
  PSEUDO        *pseudo           = &(cp->pseudo);
  CPCOEFFS_INFO *cpcoeffs_info    = &(cp->cpcoeffs_info);
  COMMUNICATE   *communicate      = &(cp->communicate);
  CPCOEFFS_POS  *cpcoeffs_pos     = &(cp->cpcoeffs_pos[ip_now]);
  STODFTINFO    *stodftInfo       = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos    = cp->stodftCoefPos;
  CPSCR         *cpscr            = &(cp->cpscr);
  NEWTONINFO    *newtonInfo;
  FRAGINFO      *fragInfo         = stodftInfo->fragInfo;
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);


     



/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

