/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: cp_wannier_control                           */
/*                                                                          */
/* This subprogram performs on the fly analysis of CP data                  */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_wannier_cpcon_local.h"
#include "../proto_defs/proto_wannier_cpcon_entry.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0
#define its_max 200


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_wannier_control(CLASS *class,GENERAL_DATA *general_data,CP *cp)

/*=======================================================================*/
{ /*begin routine*/

/*=======================================================================*/
/*            Local variable declarations                                */
/*=======================================================================*/

#include "../typ_defs/typ_mask.h"


  int myid           = class->communicate.myid;
  int num_proc       = class->communicate.np;
  MPI_Comm comm_stat = class->communicate.comm_states;
  MPI_Comm world     = class->communicate.world;

  int itime     = general_data->timeinfo.itime;
  int cp_lsda   = cp->cpopts.cp_lsda;

  int cp_wan_min_on     = cp->cpopts.cp_wan_min_opt;
  int cp_wan_on         = cp->cpopts.cp_wan_opt;
  int cp_wannier_on     = cp->cpscr.cpscr_wannier.cp_wannier_on;
  int cp_wan_init_on    = cp->cpopts.cp_wan_init_opt;
  int cp_wan_calc_frq   = cp->cp_wannier.cp_wan_calc_frq;
  int cp_dip_calc_frq   = cp->cpcoeffs_info.cp_dip_calc_frq;

  double *weight    = cp->electronic_properties.weight;
  double ***Z_real  = cp->electronic_properties.Z_real;
  double ***Z_imag  = cp->electronic_properties.Z_imag;
  double **wan_cent = cp->electronic_properties.wannier_cent;

  int iopt_cp_pw    = cp->cpcoeffs_info.iopt_cp_pw;
  int iopt_cp_dvr   = cp->cpcoeffs_info.iopt_cp_dvr;

  int cp_on,dip_calc_flag,wan_calc_flag;
  int static ifirst=1;

/*--------------------------------------------------------------------------*/
/* 0) Simulation option Check */

  cp_on=(general_data->simopts.cp+general_data->simopts.cp_wave+
         general_data->simopts.debug_cp);

  if(cp_lsda==1){
    if(myid==0){
      printf("@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("LSDA option is not allowed for Wannier dynamics yet.\n");
      printf("@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@\n");
    }
    fflush(stdout);
    exit(1);
  }/*endif*/

  if(cp_on != 1){
    if(myid==0){
      printf("@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("CP option has to be turned on for Wannier dynamics.\n");
      printf("@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@\n");
    }
    fflush(stdout);
    exit(1);
  }/*endif*/

  if((cp_wan_min_on==1 || cp_wan_on==1) && cp_wan_calc_frq < 1){
    if(myid==0){
      printf("@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("wan_calc_freq should be greater than zero for Wannier\n");
      printf("dynamics (cp_wan_opt=1 or cp_wan_min_opt=1).\n");
      printf("@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@\n");
    }
    fflush(stdout);
    exit(1);
  }/*endif*/


/*-------------------------------------------------------------------------*/
/* I) first time set up  */

  if(ifirst == 1 && cp_wan_init_on ==0){
    comp_miller_weight(general_data,weight);
  }

/*-------------------------------------------------------------------------*/
/* II) Time to compute Dipole moment   */

  dip_calc_flag=0;
  if(cp_dip_calc_frq > 0 && itime%cp_dip_calc_frq ==0){
 
    if(num_proc > 1 && iopt_cp_pw==1){
      control_coef_transpose_bck(cp,1);
    }
    calcul_dipole(class, general_data,cp, weight, Z_real, Z_imag);
    dip_calc_flag=1;
  }

/*-------------------------------------------------------------------------*/
/* III) Time to rotate orbitals onto Wannier gauge */

  wan_calc_flag=0;
  if(cp_wan_min_on==1 || cp_wan_on==1){
    if(cp_wan_calc_frq > 0 && itime%cp_wan_calc_frq == 0){
      wan_calc_flag=1;
      
      if(iopt_cp_pw==1){
        if(dip_calc_flag==0 && num_proc > 1){
           control_coef_transpose_bck(cp,1);
        }
        calcul_wannier(general_data,cp,Z_real,Z_imag,weight,wan_cent,
                       dip_calc_flag);
      }
      
      if(iopt_cp_dvr==1){
        calcul_wannier_dvr(general_data,cp,Z_real,Z_imag,weight,wan_cent,
                           dip_calc_flag);
      }
    }
  }

  if(num_proc> 1 && iopt_cp_pw==1 && wan_calc_flag==0 && dip_calc_flag==1){
    control_coef_transpose_fwd(cp,1);
  }

  if(num_proc>1){Barrier(world);}

  ifirst += 1;

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/
