/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: all-read-frag.c                                */
/*                                                                          */
/* This file provide all modified malloc modules for fragmentation in parse */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_coords_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_coords_cp_entry.h"
#include "../proto_defs/proto_coords_cp_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_interface_frag_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void readCoefFrag(CP *cp,GENERAL_DATA *general_data,CLASS *class,
               FILENAME_PARSE *filename_parse,
               CP_PARSE *cp_parse,double *tot_memory)

/*==========================================================================*/
/*               Begin subprogram:                                          */
  {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  FILE *fp_dnameci;   

/* Local pointers */
  SIMOPTS *simopts       = &(general_data->simopts);
  int istart             = cp_parse->istart_cp;
  int myid               = cp->communicate.myid;
  MPI_Comm world         = cp->communicate.world;
  int cp_wave_min        = simopts->cp_wave_min;
  int cp_wave_min_pimd   = simopts->cp_wave_min_pimd;
  int cp_min             = simopts->cp_min;
  int cp_min_on;
  int initial_spread_opt = simopts->initial_spread_opt;
  int np_states          = cp->communicate.np_states;
  int num_proc           = cp->communicate.np;
  char *dnameci          = filename_parse->dnameci;

  int ibinary            = cp->cpopts.iread_coef_binary;
  int iii;

  int iopt_cp_pw         = cp->cpcoeffs_info.iopt_cp_pw;
  int iopt_cp_dvr        = cp->cpcoeffs_info.iopt_cp_dvr;

  cp_min_on = cp_wave_min + cp_wave_min_pimd + cp_min;

/*========================================================================*/
/*  V) Allocate and initialize coefficient arrays                         */

  if(iopt_cp_pw){
    read_coef_alloc_init(cp,cp_min_on,tot_memory);
    if(num_proc>1){Barrier(world);}
  }

/*========================================================================*/
/*  VII) Read/Spread the coefficients  */

  if(iopt_cp_pw==1){
    gen_wave_frag(class,general_data,cp,cp_parse,filename_parse->vps_name);
  }

  
  read_coef_init_nhc(cp);

/*-----------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void readHmatFrag(CLASS *classMini,GENERAL_DATA *generalDataMini,
                CP *cpMini,int cp_dual_grid_opt_on, double *dbox_rat,
                int *box_rat)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */

  SIMOPTS *simopts              = &(generalDataMini->simopts);
  CELL *cell                    = &(generalDataMini->cell);
  CLATOMS_INFO *clatoms_info    = &(classMini->clatoms_info);
  ATOMMAPS *atommaps            = &(classMini->atommaps);
  GENERAL_DATA *general_data = generalDataMini;

  int ii,iii,upper;
  int i,j,ip;              /* Num: For loop counter                 */

  char dnamei[100];

  double sx,sy,sz;
  double vol,vol_cp,area;

/*  Local pointers */
  double *hmat           = cell->hmat;
  double *hmati          = cell->hmati;
  double *hmat_ewd       = cell->hmat_ewd;
  double *hmat_ewd_cp    = cell->hmat_ewd_cp;
  double *hmat_cp        = cell->hmat_cp;
  double *hmati_cp       = cell->hmati_cp;
  double *cp_box_center  = cell->cp_box_center;
  double *cp_box_center_rel  = cell->cp_box_center_rel;

  int hmat_int_typ       = cell->hmat_int_typ;
  int hmat_cons_typ      = cell->hmat_cons_typ;
  int iperd              = cell->iperd;

  int ensemble_flag;
  ensemble_flag          = generalDataMini->ensopts.nve
                         + generalDataMini->ensopts.nvt
                         + generalDataMini->ensopts.npt_i;


/*========================================================================*/
/*  I) Assign Volumes */
  for(i=1;i<=9;i++)hmat_cp[i] = hmat[i];

  gethinv(hmat,hmati,&(vol),iperd);
  gethinv(hmat_cp,hmati_cp,&(vol_cp),iperd);

  general_data->cell.vol        = vol;
  general_data->cell.vol_cp     = vol_cp;
  general_data->cell.vol0       = vol;
  general_data->baro.vol        = vol;
  general_data->par_rahman.vol  = vol;
  general_data->stat_avg.vol    = vol;
  general_data->baro.x_lnv      = log(vol)/3.0;

   for(i=1;i<=9;i++){hmat_ewd_cp[i] = hmat_cp[i];}
   for(i=1;i<=9;i++){hmat_ewd[i] = hmat[i];}

  area = hmat[1]*hmat[5] - hmat[2]*hmat[4];
  general_data->cell.area       = area;
  general_data->baro.area       = area;
  general_data->par_rahman.area = area;

  *dbox_rat = hmat[1]/hmat_cp[1];
  *box_rat  = (int)(*dbox_rat);

/*========================================================================*/
/*  II) Determine if Box is Cubic for Fast Imaging in Period.c */

  if( (hmat[4] == 0.0) && (hmat[7] == 0.0) &&
      (hmat[2] == 0.0) && (hmat[8] == 0.0) &&
      (hmat[3] == 0.0) && (hmat[6] == 0.0) &&
      (ensemble_flag != 0)){
    cell->cubic_box_flag = 1;
  }else{
    cell->cubic_box_flag = 0;
  }/*endif*/

/*========================================================================*/
/*  III) Check cell                                                      */

  check_cell(&(general_data->cell),cp_dual_grid_opt_on,*dbox_rat,dnamei);

/*========================================================================*/
/* convert cp_box_center from xtal coordinates to cartesian coordinates   */

  sx = cp_box_center[1];
  sy = cp_box_center[2];
  sz = cp_box_center[3];

  cp_box_center[1] =sx*hmat[1]+sy*hmat[4]+sz*hmat[7];
  cp_box_center[2] =sx*hmat[2]+sy*hmat[5]+sz*hmat[8];
  cp_box_center[3] =sx*hmat[3]+sy*hmat[6]+sz*hmat[9];

  sx = 0.5;
  sy = 0.5;
  sz = 0.5;

  cp_box_center_rel[1] =sx*hmat_cp[1]+sy*hmat_cp[4]+sz*hmat_cp[7];
  cp_box_center_rel[2] =sx*hmat_cp[2]+sy*hmat_cp[5]+sz*hmat_cp[8];
  cp_box_center_rel[3] =sx*hmat_cp[3]+sy*hmat_cp[6]+sz*hmat_cp[9];

/*========================================================================*/
   }/* end routine */
/*==========================================================================*/

