/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: all-mall-frag.c                                */
/*                                                                          */
/* This file provide all modified malloc modules for fragmentation in parse */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_scratch_entry.h"
#include "../proto_defs/proto_scratch_local.h"
#include "../proto_defs/proto_lists_entry.h"
#include "../proto_defs/proto_lists_local.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_interface_frag_local.h"
#include "../proto_defs/proto_intra_params_entry.h"
#include "../proto_defs/proto_intra_params_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void mallMakeListsFrag(CLASS *class,GENERAL_DATA *general_data,
                                    BONDED *bonded,int error_check_on)

/*==========================================================================*/
{/*begin routine */

int isave,iii;
int np         = class->communicate.np;
int rank       = class->communicate.myid;
MPI_Comm world = class->communicate.world;

/*==========================================================================*/
/* 0) Initialize */

/*==========================================================================*/
/* I) Malloc and create the lnk lists: Must be done first.                  */
/*                                     RESPA lnk lst shifts not             */
/*                                     needed by lnk_ver_update option      */

 if( (class->nbr_list.ilnk==1)||
    ((class->nbr_list.verlist.lnk_ver_update+class->nbr_list.iver)==2)){
   isave = general_data->timeinfo.int_res_ter;
   if(class->nbr_list.iver==1){general_data->timeinfo.int_res_ter=0;}
   if(np>1){Barrier(world);}

   lnk_mall_make(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->nbr_list),&(bonded->excl),
                 &(class->atommaps),&(general_data->cell),
                 &(bonded->intra_scr),
                 &(class->for_scr),&(general_data->timeinfo),
                 &(class->interact),&(class->tot_memory),rank,
                   error_check_on,world,np);

   if(class->nbr_list.iver==1&&rank==0){printf("\n");}
   if(class->nbr_list.iver==1){general_data->timeinfo.int_res_ter=isave;}
 /*endif*/}

/*==========================================================================*/
/* II) Malloc and create the ver lists */

 if(class->nbr_list.iver==1){
   (general_data->stat_avg.itime_update) = 0;
  if(np>1){Barrier(world);}
   ver_mall_make(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->nbr_list),&(bonded->excl),
                 &(class->atommaps),&(general_data->cell),
                 &(bonded->intra_scr),
                 &(class->for_scr),&(general_data->timeinfo),
                 &(class->interact),&(class->tot_memory),rank,
                 error_check_on,world,&(class->class_comm_forc_pkg),np);
 /*endif*/}

/*==========================================================================*/
/* III) Done */

/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void controlMallScratchFrag(CLASS *class,BONDED *bonded,CP *cp,
                          GENERAL_DATA *general_data)

/*======================================================================*/
  {/*begin routine*/
/*======================================================================*/
/*         Local variable declarations                                  */
#include "../typ_defs/typ_mask.h"

 int pimd_on,cp_on,extsys_on;

 /* Local Pointers */
 int myid           = class->communicate.myid;
 int natm_tot       = class->clatoms_info.natm_tot;
 int int_res_ter    = general_data->timeinfo.int_res_ter;
 MPI_Comm world     = class->communicate.world;
 double *tot_memory = &(class->tot_memory);
 int nfft_size      = general_data->pme_fft_pkg.nfft_size;
 int ncoef_pme      = general_data->pme_fft_pkg.ncoef_proc;
 int pme_on         = class->part_mesh.pme_on;
 int ilnk_lst       = class->nbr_list.ilnk;
 int iver_lst       = class->nbr_list.iver;
 int lnk_ver_update = class->nbr_list.verlist.lnk_ver_update;
 int np_forc        = class->communicate.np_forc;
 int cp_dual_grid_opt_on     = cp->cpopts.cp_dual_grid_opt;

/*=======================================================================*/
/* I) Output */
  
/*=======================================================================*/
/* II) Useful constants */

  pimd_on     =  general_data->simopts.pimd 
               + general_data->simopts.cp_pimd 
               + general_data->simopts.cp_wave_pimd 
               + general_data->simopts.cp_wave_min_pimd 
               + general_data->simopts.debug_pimd  
               + general_data->simopts.debug_cp_pimd;
  cp_on       =  general_data->simopts.cp_min 
               + general_data->simopts.cp_wave_min
               + general_data->simopts.cp     
               + general_data->simopts.cp_wave
               + general_data->simopts.cp_pimd
               + general_data->simopts.cp_wave_pimd
               + general_data->simopts.debug_cp
               + general_data->simopts.debug_cp_pimd
               + general_data->simopts.cp_wave_min_pimd;
  extsys_on   =  general_data->ensopts.npt_f
               + general_data->ensopts.npt_i
               + general_data->ensopts.nvt;
  class->for_scr.nlen = bonded->intra_scr.nlen;

/*=======================================================================*/
/* III) Integrator scratch */

  mall_integrator_scr_frag(pimd_on,extsys_on,&(class->clatoms_info),
                      &(class->therm_info_class),&(class->therm_info_bead),
                      &(class->int_scr),tot_memory,myid,world);
  
/*=======================================================================*/
/* IV) Intra scratch */

  mall_intra_scr_frag(&(bonded->intra_scr),tot_memory,myid,world);

/*=======================================================================*/
/* V) Intra scratch */

  mall_ewald_scr_frag(cp_on,natm_tot,int_res_ter,nfft_size,ncoef_pme,
                 &(class->part_mesh),&(class->ewd_scr),tot_memory,myid,world);

/*======================================================================*/
/* VI) CP force scr */

  if(cp_on==1){
    cp->cpscr.cpscr_nonloc.natm_nls_max = natm_tot;
    if( cp->cpcoeffs_info.iopt_cp_pw == 1 ){
      mall_cp_scr_frag(&(cp->cptherm_info),&(cp->cpopts),&(cp->cpewald),&(cp->cpscr),
                  &(cp->cpcoeffs_info),&(cp->pseudo),
                  &(cp->cp_para_fft_pkg3d_dens_cp_box),
                  &(cp->cp_para_fft_pkg3d_lg),
                  &(cp->cp_comm_state_pkg_up),&(cp->cp_comm_state_pkg_dn),
                  class->clatoms_info.hess_calc,
                  tot_memory,myid,cp_dual_grid_opt_on,world);
    }
  }/*endif*/

/*======================================================================*/
/* VII) Atom force scratch */

  mall_atm_forc_scr_frag(natm_tot,&(class->for_scr),pme_on,ilnk_lst,
                    iver_lst,lnk_ver_update,tot_memory,np_forc,myid,world);

/*=======================================================================*/
/* VIII) Output */

/*----------------------------------------------------------------------*/
  }/*end routine */
/*======================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void mall_integrator_scr_frag(int pimd_on,int extsys_on,
                        CLATOMS_INFO *clatoms_info,
                        THERM_INFO *therm_info_class,
                        THERM_INFO *therm_info_bead,
                        INT_SCR *int_scr,double *tot_memory,
                        int myid,MPI_Comm world)

/*===========================================================================*/
  {/*begin routine*/
/*===========================================================================*/
/*         Local variable declarations                                       */

  int num_nhc1_max,num_nhc1,natm_mall;
  int num_nhc_tot,ip;
  double now_memory;

  /* Local Pointers */
  int natm_tot      = clatoms_info->natm_tot;
  int pi_beads      = clatoms_info->pi_beads;
  int num_nhc       = therm_info_class->num_nhc;
  int len_nhc       = therm_info_class->len_nhc;
  int num_nhc_bead  = therm_info_bead->num_nhc;
  if(pimd_on==0){num_nhc_bead=0;}

/*===========================================================================*/
/* I) Malloc size calculation */

  natm_mall = natm_tot;
  if((natm_mall % 2)==0){natm_mall += 1;}

  num_nhc_tot  = num_nhc*len_nhc;
  num_nhc1     = num_nhc+1; 
  if((num_nhc1 % 2)==0){num_nhc1+=1;}
  num_nhc1_max = MAX(num_nhc1,num_nhc_bead+1);
  if((num_nhc1_max % 2)==0){num_nhc1_max+=1;}

/*===========================================================================*/
/* II) Malloc the vectors  */

  now_memory = 0.0;

 /*  i) Cartesian                                               */
 /*-------------------------------------------------------------*/
  int_scr->x     = (double *)cmalloc(natm_mall*sizeof(double))-1;
  int_scr->y     = (double *)cmalloc(natm_mall*sizeof(double))-1;
  int_scr->z     = (double *)cmalloc(natm_mall*sizeof(double))-1;

  int_scr->vx    = (double *)cmalloc(natm_mall*sizeof(double))-1;
  int_scr->vy    = (double *)cmalloc(natm_mall*sizeof(double))-1;
  int_scr->vz    = (double *)cmalloc(natm_mall*sizeof(double))-1;
  now_memory    += (natm_mall*6)*(sizeof(double));


 /*  ii) Extended system                                        */
 /*-------------------------------------------------------------*/

  int_scr->x_nhc  = cmall_mat(1,len_nhc,1,num_nhc1);
  int_scr->v_nhc  = cmall_mat(1,len_nhc,1,num_nhc1);
  int_scr->f_nhc  = cmall_mat(1,len_nhc,1,num_nhc1);
  int_scr->atm_kin= (double *)cmalloc(num_nhc1_max*sizeof(double))-1;
  int_scr->sc     = (double *)cmalloc(num_nhc1_max*sizeof(double))-1;
  now_memory     += (3*len_nhc*num_nhc1+2*num_nhc1_max)*(sizeof(double));

  if(pimd_on==1){
    int_scr->therm_scr = (THERM_SCR *)cmalloc(pi_beads*sizeof(THERM_SCR))-1;
    for(ip=1;ip<=pi_beads;ip++){
      int_scr->therm_scr[ip].atm_kin = 
                     (double *)cmalloc(num_nhc1_max*sizeof(double))-1;
      int_scr->therm_scr[ip].sc = 
                     (double *)cmalloc(num_nhc1_max*sizeof(double))-1;
    }/*endfor*/
    now_memory     += (pi_beads*num_nhc1_max)*(sizeof(double));
  }/*endif*/

  int_scr->x_vol_nhc = (double *)cmalloc(len_nhc*sizeof(double))-1;
  int_scr->v_vol_nhc = (double *)cmalloc(len_nhc*sizeof(double))-1;
  int_scr->f_vol_nhc = (double *)cmalloc(len_nhc*sizeof(double))-1;
  int_scr->vgmat     = (double *)cmalloc(9*sizeof(double))-1;
  now_memory        += (9 + 3*len_nhc)*(sizeof(double));

/*===========================================================================*/
/* III)  Output */

  now_memory *= 1.0e-06;
  *tot_memory += now_memory;

/*---------------------------------------------------------------------------*/
  }/*end routine */
/*===========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void mall_intra_scr_frag(INTRA_SCR *intra_scr,double *tot_memory, 
	    int myid,MPI_Comm world)

/*===========================================================================*/
  {/*begin routine*/
/*===========================================================================*/
/*         Local variable declarations                                       */

  double now_memory;
  int nlen_mall    = intra_scr->nlen;

/*===========================================================================*/
/* I) Malloc size calculation */
  
  if((nlen_mall % 2)==0){nlen_mall++;} 

/*===========================================================================*/
/* II) Malloc the vectors  */

  intra_scr->vpot  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->p11   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->p12   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->p13   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->p21   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->p22   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->p23   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->p31   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->p32   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->p33   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  
  /*-----------------------------------------------------------------------*/

  intra_scr->iatm_typ = (int *)cmalloc(nlen_mall*sizeof(int))-1;
  intra_scr->jatm_typ = (int *)cmalloc(nlen_mall*sizeof(int))-1;

  intra_scr->x1    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->y1    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->z1    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->x2    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->y2    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->z2    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->x3    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->y3    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->z3    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->x4    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->y4    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->z4    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->x5    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->y5    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->z5    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->x6    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->y6    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->z6    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->x7    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->y7    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->z7    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->x8    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->y8    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->z8    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->rmin  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dr    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dri   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->r     = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->r3    = (double *)cmalloc(nlen_mall*sizeof(double))-1;

  /*-----------------------------------------------------------------------*/

  intra_scr->fx1   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->fy1   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->fz1   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->fx2   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->fy2   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->fz2   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->fx3   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->fy3   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->fz3   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->fx4   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->fy4   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->fz4   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  
  /*-----------------------------------------------------------------------*/

  intra_scr->c_0   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->c_1   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->c_2   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->c_3   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->c_4   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->c_5   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->c_6   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->s_0   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->s_1   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->s_2   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->s_3   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->s_4   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->s_5   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->s_6   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dc_0  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dc_1  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dc_2  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dc_3  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dc_4  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dc_5  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dc_6  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->ds_0  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->ds_1  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->ds_2  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->ds_3  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->ds_4  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->ds_5  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->ds_6  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->eq    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->q1    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->q2    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->s6    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->feps  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->sc    = (double *)cmalloc(nlen_mall*sizeof(double))-1;

  /*------------------------------------------------------------------------*/

  intra_scr->dx12  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dy12  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dz12  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dx13  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dy13  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dz13  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dx23  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dy23  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dz23  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dx42  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dy42  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dz42  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dx43  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dy43  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dz43  = (double *)cmalloc(nlen_mall*sizeof(double))-1;

  intra_scr->dx56  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dy56  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dz56  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dx57  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dy57  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dz57  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dx67  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dy67  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dz67  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dx86  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dy86  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dz86  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dx87  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dy87  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->dz87  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  
  /*-----------------------------------------------------------------------*/

  intra_scr->del_r    = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->del_rc   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->swit     = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->cut      = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->spl_tmp  = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->cutoff   = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->cutoff_res=(double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->vcut_coul= (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->m1       = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->m2       = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->vx1      = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->vy1      = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->vz1      = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->vx2      = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->vy2      = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->vz2      = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->sc_1     = (double *)cmalloc(nlen_mall*sizeof(double))-1;
  intra_scr->sc_2     = (double *)cmalloc(nlen_mall*sizeof(double))-1;

/*===========================================================================*/
/* III)  Output */

  now_memory   = (nlen_mall*(sizeof(double)*134 + sizeof(int)*2 ))*1.0e-06;;
  *tot_memory += (now_memory);

/*---------------------------------------------------------------------------*/
  }/*end routine */
/*===========================================================================*/


  

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void mall_ewald_scr_frag(int cp_on,int natm_tot,int int_res_ter,
                    int nfft_size,int ncoef_pme,PART_MESH *part_mesh, 
                    EWD_SCR *ewd_scr,double *tot_memory,
                    int myid,MPI_Comm world)

/*==========================================================================*/
  {/*begin routine*/
/*==========================================================================*/
/*         Local variable declarations                                      */

  double now_memory;
  int num_re=0,num_int=0;
  int pme_on       = part_mesh->pme_on;
  int nlen_pme     = part_mesh->nlen_pme;
  int n_interp     = part_mesh->n_interp;
  int ngrid_c      = part_mesh->ngrid_c;
  
/*===========================================================================*/
/* I) Malloc size calculation */

  int natm_mall = natm_tot;
  if((natm_mall % 2)==0){natm_mall += 1;}

/*===========================================================================*/
/* II) Malloc the vectors  */

  ewd_scr->q    = (double *)cmalloc(natm_mall*sizeof(double))-1;
  ewd_scr->x    = (double *)cmalloc(natm_mall*sizeof(double))-1;
  ewd_scr->y    = (double *)cmalloc(natm_mall*sizeof(double))-1;
  ewd_scr->z    = (double *)cmalloc(natm_mall*sizeof(double))-1;
  ewd_scr->fx   = (double *)cmalloc(natm_mall*sizeof(double))-1;
  ewd_scr->fy   = (double *)cmalloc(natm_mall*sizeof(double))-1;
  ewd_scr->fz   = (double *)cmalloc(natm_mall*sizeof(double))-1;
  num_re+=7*natm_mall;

  if(pme_on==0||cp_on==1){
    ewd_scr->heli = (double *)cmalloc(natm_mall*sizeof(double))-1;
    ewd_scr->helr = (double *)cmalloc(natm_mall*sizeof(double))-1;
    ewd_scr->cossc= (double *)cmalloc(natm_mall*sizeof(double))-1;
    ewd_scr->sinsc= (double *)cmalloc(natm_mall*sizeof(double))-1;
    num_re+=4*natm_mall;
  }/*endif*/

  if(cp_on==1 || (pme_on==0 && int_res_ter==1)){
    ewd_scr->fx2  = (double *)cmalloc(natm_mall*sizeof(double))-1;
    ewd_scr->fy2  = (double *)cmalloc(natm_mall*sizeof(double))-1;
    ewd_scr->fz2  = (double *)cmalloc(natm_mall*sizeof(double))-1;
    num_re += 3*natm_mall;
  }/*endif*/

  if(cp_on==1){
    ewd_scr->temp = (double *)cmalloc(natm_mall*sizeof(double))-1;
    ewd_scr->helr_now  = (double *)cmalloc(natm_mall*sizeof(double))-1;
    ewd_scr->heli_now  = (double *)cmalloc(natm_mall*sizeof(double))-1;
    ewd_scr->vtemp_now = (double *)cmalloc(natm_mall*sizeof(double))-1;
    num_re += 4*natm_mall;
  }/*endif*/

 
  if(pme_on == 1 && cp_on==0){

    part_mesh->qgrid          = (double *)cmalloc(nfft_size*sizeof(double))-1;
    part_mesh->qgrid_scr      = (double *)cmalloc(nfft_size*sizeof(double))-1;
    part_mesh->qgrid_tmp_real = (double *) cmalloc(ncoef_pme*sizeof(double))-1;
    part_mesh->qgrid_tmp_imag = (double *) cmalloc(ncoef_pme*sizeof(double))-1;
    num_re                   += (2*nfft_size + 2*ncoef_pme);

    part_mesh->iatemp    = (int *)cmalloc(nlen_pme*sizeof(int))-1;
    part_mesh->ibtemp    = (int *)cmalloc(nlen_pme*sizeof(int))-1;
    part_mesh->ictemp    = (int *)cmalloc(nlen_pme*sizeof(int))-1;
    part_mesh->frac_a    = (double *)cmalloc(nlen_pme*sizeof(double))-1;
    part_mesh->frac_b    = (double *)cmalloc(nlen_pme*sizeof(double))-1;
    part_mesh->frac_c    = (double *)cmalloc(nlen_pme*sizeof(double))-1;
    num_re              += (3*nlen_pme); 
    num_int             += (3*nlen_pme);

    part_mesh->ua        = cmall_mat(1,n_interp,1,nlen_pme);
    part_mesh->ub        = cmall_mat(1,n_interp,1,nlen_pme);
    part_mesh->uc        = cmall_mat(1,n_interp,1,nlen_pme);
    part_mesh->mn_a      = cmall_mat(1,n_interp,1,nlen_pme);
    part_mesh->mn_b      = cmall_mat(1,n_interp,1,nlen_pme);
    part_mesh->mn_c      = cmall_mat(1,n_interp,1,nlen_pme);
    part_mesh->dmn_a     = cmall_mat(1,n_interp,1,nlen_pme);
    part_mesh->dmn_b     = cmall_mat(1,n_interp,1,nlen_pme);
    part_mesh->dmn_c     = cmall_mat(1,n_interp,1,nlen_pme);
    part_mesh->qgrid_now = cmall_mat(1,n_interp,1,nlen_pme);
    part_mesh->igrid_a   = cmall_int_mat(1,n_interp,1,nlen_pme);
    part_mesh->igrid_b   = cmall_int_mat(1,n_interp,1,nlen_pme);
    part_mesh->igrid_c   = cmall_int_mat(1,n_interp,1,nlen_pme);
    part_mesh->igrid_now = cmall_int_mat(1,n_interp,1,nlen_pme);
    num_re              += (14*nlen_pme*n_interp);

    part_mesh->ioff_c   = (int *)cmalloc(ngrid_c*sizeof(int))-1;
    part_mesh->nc       = (int *)cmalloc(ngrid_c*sizeof(int))-1;
    num_int            +=(3*ngrid_c);

  }/*endif*/

/*=========================================================================*/
/* III) Output */
  
  now_memory   = (num_re*sizeof(double)+num_int*sizeof(int))*1.e-06;
  *tot_memory += now_memory;

/*--------------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void mall_cp_scr_frag(CPTHERM_INFO *cptherm_info,CPOPTS *cpopts,CPEWALD *cpewald,
                 CPSCR *cpscr,CPCOEFFS_INFO *cpcoeffs_info,PSEUDO *pseudo,
                 PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box,
                 PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg,
                 CP_COMM_STATE_PKG *cp_comm_state_pkg_up,
                 CP_COMM_STATE_PKG *cp_comm_state_pkg_dn,int atm_hess_calc,
                 double *tot_memory,
                 int myid,int cp_dual_grid_opt_on, MPI_Comm world)

/*==========================================================================*/
  {/*begin routine*/
/*==========================================================================*/
/*         Local pointer declarations                                      */

  int nstate_up         = cpcoeffs_info->nstate_up_proc;
  int nstate_dn         = cpcoeffs_info->nstate_dn_proc;
  int nstate_up_tot     = cpcoeffs_info->nstate_up;
  int nstate_dn_tot     = cpcoeffs_info->nstate_dn;
  int ncoef             = cpcoeffs_info->ncoef;
  int ncoef_l           = cpcoeffs_info->ncoef_l;
  int pi_beads          = cpcoeffs_info->pi_beads;
  int cp_laplacian_on   = cpcoeffs_info->cp_laplacian_on;
  int cp_tau_functional = cpcoeffs_info->cp_tau_functional;
  int num_c_nhc_proc    = cptherm_info->num_c_nhc_proc;
  int massiv_flag       = cptherm_info->massiv_flag;
  int n_ang_max         = pseudo->n_ang_max;
  int n_rad_max         = pseudo->n_rad_max;
  int natm_nls_max      = cpscr->cpscr_nonloc.natm_nls_max;
  int cp_lsda           = cpopts->cp_lsda;
  int cp_ptens_calc     = cpopts->cp_ptens_calc;
  int cp_hess_calc      = cpopts->cp_hess_calc;
  int cp_gga            = cpopts->cp_gga;
  int cp_ke_dens_on     = cpcoeffs_info->cp_ke_dens_on;
  int cp_elf_calc_frq   = cpcoeffs_info->cp_elf_calc_frq;
  int cp_norb           = cpopts->cp_norb;
  int cp_para_opt       = cpopts->cp_para_opt;
  int np_states         = cp_comm_state_pkg_up->num_proc;
  int nfft_up_proc      = cp_para_fft_pkg3d_lg->nfft_proc;
  int nfft_up           = cp_para_fft_pkg3d_lg->nfft;
  int nstate_max_up     = cp_comm_state_pkg_up->nstate_max;
  int nstate_ncoef_proc_max_up = cp_comm_state_pkg_up->nstate_ncoef_proc_max;
  int nstate_max_dn     = cp_comm_state_pkg_up->nstate_max;
  int nstate_ncoef_proc_max_dn = cp_comm_state_pkg_up->nstate_ncoef_proc_max;
  int num_c_nhc1        = num_c_nhc_proc+1;

  int ncoef_l_pme_dual,ncoef_l_pme_dual_proc;
  int ncoef_l_proc_max_mall;
  int ncoef_l_proc_max_mall_ke;

  int ncoef_l_dens_cp_box;
  int ncoef_l_proc_max_dens_cp_box;
  int nfft_up_proc_dens_cp_box,nfft_up_dens_cp_box,nfft2_up_dens_cp_box;
  int nfft_dn_proc_dens_cp_box,nfft_dn_dens_cp_box ;
  int ncoef_l_proc_max_dn;

  int nkf1_cp_box = cp_para_fft_pkg3d_dens_cp_box->nkf1;
  int nkf2_cp_box = cp_para_fft_pkg3d_dens_cp_box->nkf2;
  int nkf3_cp_box = cp_para_fft_pkg3d_dens_cp_box->nkf3;
  int n_interp_pme_dual = pseudo->n_interp_pme_dual;
/*--------------------------------------------------------------------------*/
/*         Local variable declarations                                      */

  double now_memory;
  int i,iii;
  int nlscr_up,nlscr_dn,nlscr_up_pv,nlscr_dn_pv,ncoef_l_pv,ncoef_l_proc_max;
  int ncoef_l_proc_max_mall_cp_box,ncoef_l_proc_max_mall_cp_box_dn;
  int nfft2_mall_up,nfft2_mall_dn,nfft2_mall_up_proc,nfft2_mall_dn_proc;
  int nfft2_up,nfft2_dn,irem;
  int nfft2_up_proc,nfft2_dn_proc,nfft_dn;
  int nfft2_up_gga,nfft2_dn_gga,nlap_up,nlap_dn,nlap_up_ptens;
  int nfft2_up_ke_dens,nfft2_dn_ke_dens;
  int nlap_dn_ptens,ngga_up,ngga_dn,nlap_g_up,nlap_g_dn;
  int nlap_g_up_ptens, nlap_g_dn_ptens;
  int ncoef_up,ncoef_dn;
  int par_size_up,par_size_dn;
  int ncoef2_up_c,ncoef2_dn_c;
  int ncoef2_up,ncoef2_dn,ncoef2_up_spec,ncoef2_dn_spec;
  int ncoef2_up_par,nstate,nstate2,nstate_tot,nstate2_tot;
  int num=0;
  int zero=0;
  int ncoef_l_mall_proc_max_dual,ncoef_l_mall_proc_max_dual_dn;
  int map_count;
  int mtemp;
  int nlen_pme,pme_nkf3,ninterp_pme,nmall;
  int mall_size;

  int cp_wan_opt        = cpopts->cp_wan_opt;
  int cp_wan_min_opt    = cpopts->cp_wan_min_opt;
  int cp_wan_init_opt   = cpopts->cp_wan_init_opt;
  int ndim_wannier;
  int mm=5;

/*=========================================================================*/
/* I) Malloc size calculation */

 /*-------------------------------------------------------------------------*/
 /* i) Dual grid CP : Define the small dense grid sizes */

  if(cp_dual_grid_opt_on >= 1){

    nfft_up_proc_dens_cp_box   = cp_para_fft_pkg3d_dens_cp_box->nfft_proc;
    nfft_up_dens_cp_box        = cp_para_fft_pkg3d_dens_cp_box->nfft;
    nfft2_up_dens_cp_box       = nfft_up_dens_cp_box/2;

    nfft_dn_proc_dens_cp_box   = (cp_lsda == 1 ? nfft_up_proc_dens_cp_box : 0);
    nfft_dn_dens_cp_box        = (cp_lsda == 1 ? nfft_up_dens_cp_box : 0);

    ncoef_l_dens_cp_box          = cpcoeffs_info->ncoef_l_dens_cp_box;
    ncoef_l_proc_max_dens_cp_box = ncoef_l_dens_cp_box/np_states;
    irem                         = (ncoef_l_proc_max_dens_cp_box % np_states);
    if(irem>0){ncoef_l_proc_max_dens_cp_box++;}

  }/*endif cp_dual_grid_opt_on */

 if(cp_dual_grid_opt_on == 2){
   ncoef_l_pme_dual = cp_para_fft_pkg3d_lg->ncoef;
 }/*endif cp_dual_grid_opt_on*/

 /*-------------------------------------------------------------------------*/
 /* ii) Normal CP : Define the grid size              */
 /*     Dual   CP : Define the large sparse grid size */

  nfft2_up      = nfft_up/2;
  nfft2_up_proc = nfft_up_proc/2;

  nfft_dn       = (cp_lsda == 1 ? nfft_up : 0);
  nfft2_dn      = (cp_lsda == 1 ? nfft2_up : 0);
  nfft2_dn_proc = (cp_lsda == 1 ? nfft2_up_proc : 0);

 /*-------------------------------------------------------------------------*/
 /* iii) Choose the correct size for your application                       */
 /*      This is always the small dense grid                                */

  nfft2_mall_up      = (cp_dual_grid_opt_on >= 1 ? 
                        nfft_up_dens_cp_box/2 : nfft2_up);
  nfft2_mall_up_proc = (cp_dual_grid_opt_on >= 1 ? 
                        nfft_up_proc_dens_cp_box/2 : nfft2_up_proc);

  nfft2_mall_dn      = (cp_dual_grid_opt_on >= 1 ?
                        nfft_dn_dens_cp_box/2 : nfft2_dn);
  nfft2_mall_dn_proc = (cp_dual_grid_opt_on >= 1 ? 
                        nfft_dn_proc_dens_cp_box/2 : nfft2_dn_proc);

 /*-------------------------------------------------------------------------*/
 /* iv) Wave function size (spherically cutoff small g-space for dense box) */

  ncoef_up  = ncoef;
  ncoef_dn  = (cp_lsda == 1 ? ncoef : 0);

 /*-------------------------------------------------------------------------*/
 /* v) Normal CP: The sphere cut large g-space for the dense box            */
 /*    Dual CP  : The sphere cut large g-space for the large sparse box     */

  ncoef_l_proc_max = ncoef_l/np_states;
  irem = (ncoef_l % np_states);
  if(irem>0){ncoef_l_proc_max++;}
 
 if(cp_dual_grid_opt_on == 2){
   ncoef_l_pme_dual_proc = ncoef_l_pme_dual/np_states;
   irem = (ncoef_l_pme_dual_proc % np_states);
   if(irem>0){ncoef_l_pme_dual_proc++;}
  }/*endif cp_dual_grid_opt_on */

  ncoef_l_proc_max_mall = (cp_dual_grid_opt_on == 2 ? ncoef_l_pme_dual_proc 
                                                    : ncoef_l_proc_max);
  ncoef_l_proc_max_mall_ke = (cp_ke_dens_on == 1 ? ncoef_l_proc_max_mall:0);
  ncoef_l_proc_max_dn      = (cp_lsda == 1 ? ncoef_l_proc_max_mall : 0);

 /*-------------------------------------------------------------------------*/
 /* vi) Choose the large g-space malloc size based on the dual or normal opt*/
 /*     The malloc size is the always the dense grid.                       */
 /*     Set special dual malloc sizes to zero to avoid mallocing extra      */
 /*     memory during normal CP.                                            */

  if(cp_dual_grid_opt_on==0){

    ncoef_l_proc_max_mall_cp_box = ncoef_l_proc_max;

  }else{

    ncoef_l_proc_max_mall_cp_box = ncoef_l_dens_cp_box/np_states;
    irem = (ncoef_l_dens_cp_box % np_states);
    if(irem>0){ncoef_l_proc_max_mall_cp_box++;}
    ncoef_l_proc_max_mall_cp_box_dn =
                   (cp_lsda == 1 ? ncoef_l_proc_max_mall_cp_box : 0);

  }/*endif cp_dual_grid_opt_on*/

  ncoef_l_mall_proc_max_dual    = (cp_dual_grid_opt_on >= 1 ? 
                                   ncoef_l_proc_max_mall_cp_box : 0);
  ncoef_l_mall_proc_max_dual_dn = (cp_dual_grid_opt_on >= 1 ? 
                                   ncoef_l_proc_max_mall_cp_box_dn : 0);

 /*-------------------------------------------------------------------*/
 /* vii) Wavefunction scratch sizes : always on small dense grid     */

  par_size_up = nstate_max_up*nstate_ncoef_proc_max_up;
  par_size_dn = nstate_max_dn*nstate_ncoef_proc_max_dn;

  if(massiv_flag==0){
    ncoef2_up_c = MAX(ncoef_up*nstate_up,num_c_nhc_proc);
  }else{
    ncoef2_up_c = MAX(ncoef_up*nstate_up,2*ncoef_up);
  }/*endif*/  

  if(np_states>1){ncoef2_up_c = MAX(ncoef2_up_c,par_size_up);}
  ncoef2_up      = ncoef_up*nstate_up;
  ncoef2_dn      = ncoef_dn*nstate_dn;
  ncoef2_up      = MAX(ncoef2_up,par_size_up);
  ncoef2_dn      = MAX(ncoef2_dn,par_size_dn);
  ncoef2_up_spec = MAX(ncoef2_up,ncoef2_dn);
  ncoef2_dn_spec = ncoef2_up_spec;
  if(np_states == 1 &&  cp_norb==0){ncoef2_up_spec = 1;}
  ncoef2_up_par = ncoef2_up_spec;
  if(np_states == 1){ncoef2_up_par = 1;}

  nstate        = MAX(nstate_up, nstate_dn);
  nstate2       = nstate*nstate;
  nstate_tot    = MAX(nstate_up_tot,nstate_dn_tot);
  nstate2_tot   = nstate_tot*nstate_tot;

 /*-------------------------------------------------------------------*/
 /* viii) Nonlocal sizes : always on small dense grid */

   nlscr_up  = nstate_up*(n_ang_max+1)*(n_ang_max+1)
                        *(n_rad_max)*natm_nls_max;
  nlscr_dn  = 0;
  if(cp_lsda==1){
     nlscr_dn = nstate_dn*(n_ang_max+1)*(n_ang_max+1)
                         *(n_rad_max)*(natm_nls_max);
  }/*endif*/
  nlscr_up_pv = 0;
  nlscr_dn_pv = 0;
  ncoef_l_pv  = 0;
  if(cp_ptens_calc == 1){
    nlscr_up_pv = nlscr_up;
    nlscr_dn_pv = nlscr_dn;
    ncoef_l_pv  = ncoef_l_proc_max_mall;
  }/* endif */
  if(atm_hess_calc == 3){
    nlscr_up_pv = nlscr_up;
    nlscr_dn_pv = nlscr_dn;
  }/* endif */

 /*-------------------------------------------------------------------*/
 /* ix) GGA sizes : Always the small dense grid                       */

  nfft2_up_gga  = ((cp_gga == 1 || cp_elf_calc_frq > 0) ? nfft2_mall_up_proc:0);
  nfft2_dn_gga  = (((cp_gga == 1 || cp_elf_calc_frq > 0) && cp_lsda == 1) 
                ? nfft2_mall_dn_proc:0);
  nfft2_up_ke_dens = (cp_ke_dens_on == 1 ? nfft2_mall_up_proc:0);
  nfft2_dn_ke_dens = ((cp_ke_dens_on == 1 && cp_lsda == 1) ? nfft2_mall_dn_proc:0);
  nlap_up       = (cp_laplacian_on == 1 ? nfft2_up_gga : 1);
  nlap_dn       = (cp_laplacian_on == 1 ? nfft2_dn_gga : 1);
  nlap_up_ptens = (cp_laplacian_on == 1&&cp_ptens_calc == 1?nfft2_up_gga : 1);
  nlap_dn_ptens = (cp_laplacian_on == 1&&cp_ptens_calc == 1?nfft2_dn_gga : 1);

  nlap_up       = ( cp_laplacian_on == 1 ? nfft2_up_gga : 1);
  nlap_dn       = ( cp_laplacian_on == 1 ? nfft2_dn_gga : 1); 
  nlap_up_ptens = ((cp_laplacian_on == 1 && cp_ptens_calc == 1) ?
                                           nfft2_up_gga : 1);
  nlap_dn_ptens = ((cp_laplacian_on == 1 && cp_ptens_calc == 1) ?
                                           nfft2_dn_gga : 1);  

  ngga_up         = ((cp_gga == 1 || cp_elf_calc_frq > 0)  
                  ? ncoef_l_proc_max_mall_cp_box:0);
  ngga_dn         = (cp_lsda == 1 ? ngga_up:0);
  nlap_g_up       = (cp_laplacian_on == 1 ? ngga_up : 1);
  nlap_g_dn       = (cp_laplacian_on == 1 ? ngga_dn : 1);
  nlap_g_up_ptens = (cp_laplacian_on == 1 && cp_ptens_calc == 1 ? ngga_up : 1);
  nlap_g_dn_ptens = (cp_laplacian_on == 1 && cp_ptens_calc == 1 ? ngga_dn : 1);

  /*------------------------------------------------------------------------*/
  /* x) Wannier scratch size   */

  if(cp_wan_opt ==1 || cp_wan_min_opt==1 || cp_wan_init_opt==1){
    cpscr->cpscr_wannier.cp_wannier_on=1;
    ndim_wannier=nstate_up_tot*nstate_up_tot;
  }else{
    cpscr->cpscr_wannier.cp_wannier_on=0;
    ndim_wannier=0;
  }



/*===========================================================================*/
/* II) Malloc the vectors  */

/*------------------------------------------------------------------*/
/* Non_local  : Always small dense grid                             */

  cpscr->cpscr_nonloc.vnlre_up 
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr->cpscr_nonloc.vnlim_up     
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;

  cpscr->cpscr_nonloc.vnlre_dn     
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr->cpscr_nonloc.vnlim_dn     
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;

  cpscr->cpscr_nonloc.dvnlre_x_up  
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_y_up  
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_z_up  
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_x_up  
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_y_up  
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_z_up  
                   = (double *)cmalloc(nlscr_up*sizeof(double))-1;

  num += 10*nlscr_up;

  cpscr->cpscr_nonloc.dvnlre_x_dn  
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_y_dn  
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_z_dn  
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_x_dn  
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_y_dn  
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_z_dn  
                   = (double *)cmalloc(nlscr_dn*sizeof(double))-1;

  num += 6*nlscr_dn;

  cpscr->cpscr_nonloc.dvnlre_gxgx_up  
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gygy_up  
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gzgz_up  
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gxgy_up  
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gygz_up  
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gxgz_up  
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gxgx_up  
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gygy_up  
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gzgz_up  
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gxgy_up  
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gygz_up  
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gxgz_up  
                   = (double *)cmalloc(nlscr_up_pv*sizeof(double))-1;

  num += 12*nlscr_up_pv;

  cpscr->cpscr_nonloc.dvnlre_gxgx_dn  
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gygy_dn  
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gzgz_dn  
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gxgy_dn  
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gygz_dn  
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_gxgz_dn  
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;

  cpscr->cpscr_nonloc.dvnlim_gxgx_dn  
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gygy_dn  
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gzgz_dn  
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gxgy_dn  
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gygz_dn  
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_gxgz_dn  
                   = (double *)cmalloc(nlscr_dn_pv*sizeof(double))-1;

  num += 12*nlscr_dn_pv;

/*------------------------------------------------------------------*/
/* rho and tau-dependent quantities if necessary */

 if(cp_para_opt==0){
  cpscr->cpscr_rho.v_ks_up = (double *)cmalloc(nfft2_mall_up*sizeof(double))-1;
  cpscr->cpscr_rho.v_ks_dn = (double *)cmalloc(nfft2_mall_dn*sizeof(double))-1;

  if(cp_tau_functional==1){
    cpscr->cpscr_rho.v_ks_tau_up = (double *)cmalloc(nfft2_mall_up*sizeof(double))-1;
    cpscr->cpscr_rho.v_ks_tau_dn = (double *)cmalloc(nfft2_mall_dn*sizeof(double))-1;
  }/* endif tau functional */

   num += (nfft2_mall_up + nfft2_mall_dn);
   if(cp_tau_functional==1){
     num += (nfft2_mall_up + nfft2_mall_dn);
   }
 }else{
   cpscr->cpscr_rho.v_ks_up = (double *)
                               cmalloc(nfft2_mall_up_proc*sizeof(double))-1;
   cpscr->cpscr_rho.v_ks_dn = (double *)
                               cmalloc(nfft2_mall_dn_proc*sizeof(double))-1; 
   if(cp_tau_functional==1){
     cpscr->cpscr_rho.v_ks_tau_up = (double *)
                                 cmalloc(nfft2_mall_up_proc*sizeof(double))-1;
     cpscr->cpscr_rho.v_ks_tau_dn = (double *)
                                 cmalloc(nfft2_mall_dn_proc*sizeof(double))-1; 
   }/* endif tau functional */
   num += (nfft2_mall_up_proc + nfft2_mall_dn_proc);
   if(cp_tau_functional == 1){
      num += (nfft2_mall_up_proc + nfft2_mall_dn_proc);
   }
  }/*endif*/
  cpscr->cpscr_rho.rho_up = (double *)
                            cmalloc(nfft2_mall_up_proc*sizeof(double))-1;
  cpscr->cpscr_rho.rho_dn = (double *)
                            cmalloc(nfft2_mall_dn_proc*sizeof(double))-1; 

  num += (nfft2_mall_up_proc + nfft2_mall_dn_proc);

  cpscr->cpscr_rho.rhocr_up = (double *)
                        cmalloc(ncoef_l_proc_max_mall*sizeof(double))-1;
  cpscr->cpscr_rho.rhoci_up = (double *)
                        cmalloc(ncoef_l_proc_max_mall*sizeof(double))-1;

  cpscr->cpscr_rho.rhocr_scr = (double *)
                        cmalloc(ncoef_l_proc_max_mall_ke*sizeof(double))-1;
  cpscr->cpscr_rho.rhoci_scr = (double *)
                        cmalloc(ncoef_l_proc_max_mall_ke*sizeof(double))-1;

  cpscr->cpscr_rho.rhocr_dn = (double *)
                              cmalloc(ncoef_l_proc_max_dn*sizeof(double))-1;
  cpscr->cpscr_rho.rhoci_dn = (double *)
                              cmalloc(ncoef_l_proc_max_dn*sizeof(double))-1;

  num += 2*ncoef_l_proc_max_mall;
  num += 2*ncoef_l_proc_max_mall_ke;
  num += 2*ncoef_l_proc_max_dn;

  cpscr->cpscr_rho.rhocr_up_dens_cp_box = (double *)
                        cmalloc(ncoef_l_mall_proc_max_dual*sizeof(double))-1;
  cpscr->cpscr_rho.rhoci_up_dens_cp_box = (double *)
                        cmalloc(ncoef_l_mall_proc_max_dual*sizeof(double))-1;

  cpscr->cpscr_rho.rhocr_dn_dens_cp_box = (double *)
                       cmalloc(ncoef_l_mall_proc_max_dual_dn*sizeof(double))-1;
  cpscr->cpscr_rho.rhoci_dn_dens_cp_box = (double *)
                       cmalloc(ncoef_l_mall_proc_max_dual_dn*sizeof(double))-1;

   num += 2*ncoef_l_mall_proc_max_dual;
   num += 2*ncoef_l_mall_proc_max_dual_dn;

/*------------------------------------------------------------------*/
/* local */

  cpscr->cpscr_loc.vextr = (double *)
                        cmalloc(ncoef_l_proc_max_mall*sizeof(double))-1;
  cpscr->cpscr_loc.vexti = (double *)
                        cmalloc(ncoef_l_proc_max_mall*sizeof(double))-1;
  num += 2*ncoef_l_proc_max;

  if(cp_dual_grid_opt_on == 2){
    cpscr->cpscr_loc.vextr_dens_cp_box = (double *)
                            cmalloc(ncoef_l_proc_max_mall_cp_box*sizeof(double))-1;
    cpscr->cpscr_loc.vexti_dens_cp_box = (double *)
                             cmalloc(ncoef_l_proc_max_mall_cp_box*sizeof(double))-1;
    num += 2*ncoef_l_proc_max_mall_cp_box;
  }

  cpscr->cpscr_loc.dvextr = (double *)cmalloc(ncoef_l_pv*sizeof(double))-1;
  cpscr->cpscr_loc.dvexti = (double *)cmalloc(ncoef_l_pv*sizeof(double))-1;

  num += 2*ncoef_l_pv;

/*------------------------------------------------------------------*/
/* gga */

  cpscr->cpscr_grho.d2_rho_up = (double *)cmalloc(nlap_up*sizeof(double))-1;
  cpscr->cpscr_grho.d2_rho_dn = (double *)cmalloc(nlap_dn*sizeof(double))-1;

  num += (nlap_up + nlap_dn);

  cpscr->cpscr_grho.d2_rho_up_store 
                    = (double *)cmalloc(nlap_up_ptens*sizeof(double))-1;
  cpscr->cpscr_grho.d2_rho_dn_store    
                    = (double *)cmalloc(nlap_dn_ptens*sizeof(double))-1;

  num += (nlap_up_ptens + nlap_dn_ptens);

  cpscr->cpscr_grho.d_rhox_up = (double *)
                                cmalloc(nfft2_up_gga*sizeof(double))-1;
  cpscr->cpscr_grho.d_rhoy_up = (double *)
                                cmalloc(nfft2_up_gga*sizeof(double))-1;
  cpscr->cpscr_grho.d_rhoz_up = (double *)
                                cmalloc(nfft2_up_gga*sizeof(double))-1;
  cpscr->cpscr_grho.d_rhox_dn = (double *)
                                cmalloc(nfft2_dn_gga*sizeof(double))-1;
  cpscr->cpscr_grho.d_rhoy_dn = (double *)
                                cmalloc(nfft2_dn_gga*sizeof(double))-1;
  cpscr->cpscr_grho.d_rhoz_dn = (double *)
                                cmalloc(nfft2_dn_gga*sizeof(double))-1;

  cpscr->cpscr_grho.elec_ke_dens_up = (double *)
                                cmalloc(nfft2_up_ke_dens*sizeof(double))-1;
  cpscr->cpscr_grho.elec_ke_dens_dn = (double *)
                                cmalloc(nfft2_dn_ke_dens*sizeof(double))-1;
  num += 3*(nfft2_up_gga + nfft2_dn_gga) 
       +   (nfft2_up_ke_dens + nfft2_dn_ke_dens);

  cpscr->cpscr_grho.g_rhor_x  = (double *)cmalloc(ngga_up*sizeof(double))-1;
  cpscr->cpscr_grho.g_rhor_y  = (double *)cmalloc(ngga_up*sizeof(double))-1;
  cpscr->cpscr_grho.g_rhor_z  = (double *)cmalloc(ngga_up*sizeof(double))-1;
  cpscr->cpscr_grho.g_rhoi_x  = (double *)cmalloc(ngga_up*sizeof(double))-1;
  cpscr->cpscr_grho.g_rhoi_y  = (double *)cmalloc(ngga_up*sizeof(double))-1;
  cpscr->cpscr_grho.g_rhoi_z  = (double *)cmalloc(ngga_up*sizeof(double))-1;


  num += 6*ngga_up;

  cpscr->cpscr_grho.g2_rhor   = (double *)cmalloc(nlap_g_up*sizeof(double))-1;
  cpscr->cpscr_grho.g2_rhoi   = (double *)cmalloc(nlap_g_up*sizeof(double))-1;
  cpscr->cpscr_grho.rhocr_up_st = (double *)
                                  cmalloc(nlap_g_up_ptens*sizeof(double))-1;
  cpscr->cpscr_grho.rhoci_up_st = (double *)
                                  cmalloc(nlap_g_up_ptens*sizeof(double))-1;
  cpscr->cpscr_grho.rhocr_dn_st = (double *)
                                  cmalloc(nlap_g_dn_ptens*sizeof(double))-1;
  cpscr->cpscr_grho.rhoci_dn_st = (double *)
                                  cmalloc(nlap_g_dn_ptens*sizeof(double))-1;

  num += 2*nlap_g_up + 2*(nlap_g_up_ptens + nlap_g_dn_ptens);

/*------------------------------------------------------------------*/
/* cpewald */

  cpewald->ak2    =  (double *)
              cmalloc((ncoef_l_proc_max_mall)*sizeof(double))-1;
  cpewald->ak2_sm = (double *)cmalloc(ncoef*sizeof(double))-1;
  num += (ncoef_l_proc_max + ncoef);
  if(cp_dual_grid_opt_on == 2){
   cpewald->ak2_dens_cp_box = (double *)cmalloc(ncoef_l_proc_max_mall_cp_box*
                                 sizeof(double))-1;
   num += ncoef_l_proc_max_mall_cp_box;
  }

/*------------------------------------------------------------------*/
/* wave */
  cpscr->cpscr_wave.cre_up       
                    = (double *)cmalloc(ncoef2_up_c*sizeof(double))-1;
  cpscr->cpscr_wave.cim_up       
                    = (double *)cmalloc(ncoef2_up*sizeof(double))-1;
  cpscr->cpscr_wave.cre_dn       
                    = (double *)cmalloc(ncoef2_dn*sizeof(double))-1;
  cpscr->cpscr_wave.cim_dn       
                    = (double *)cmalloc(ncoef2_dn*sizeof(double))-1;

  num += 2*(ncoef2_up + ncoef2_dn);

  if(cp_para_opt==0){
   if(cp_dual_grid_opt_on >= 1){mtemp = MAX(nfft_up,nfft_up_dens_cp_box);}
   if(cp_dual_grid_opt_on == 0){mtemp = nfft_up;}

   cpscr->cpscr_wave.zfft         
                    = (double *)cmalloc(mtemp*sizeof(double))-1;
   cpscr->cpscr_wave.zfft_tmp
                    = (double *)cmalloc(mtemp*sizeof(double))-1;
   num += 2*(mtemp);
  }else{
   if(cp_dual_grid_opt_on >= 1){mtemp = 
                          MAX(nfft_up_proc,nfft_up_proc_dens_cp_box);}
   if(cp_dual_grid_opt_on == 0){mtemp = nfft_up_proc;}

   cpscr->cpscr_wave.zfft         
                    = (double *)cmalloc(mtemp*sizeof(double))-1;
   cpscr->cpscr_wave.zfft_tmp
                    = (double *)cmalloc(mtemp*sizeof(double))-1;

   num += 2*(mtemp);
 }/*endif*/
    cpscr->cpscr_wave.zfft_mall_size = mtemp;

  if(cp_dual_grid_opt_on >= 1 && np_states >= 1){
    map_count = nfft_up_dens_cp_box/nkf1_cp_box;

    cpscr->cpscr_rho.map_dual       = (int *) cmalloc(map_count*sizeof(int))-1;
    cpscr->cpscr_rho.map_dual_upack = (int *) cmalloc(map_count*sizeof(int))-1;

   num += 2*(map_count);
  }else{
   cpscr->cpscr_rho.map_dual       = (int *) cmalloc(zero*sizeof(int))-1;
   cpscr->cpscr_rho.map_dual_upack = (int *) cmalloc(zero*sizeof(int))-1;

  }/*endif*/

/*------------------------------------------------------------------*/
/* ovmat */
  cpscr->cpscr_ovmat.ovlap1       
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.ovlap2       
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.ovlap3       
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.ovlap4       
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.ovlap5       
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.ovlap6       
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.ovlap7       
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.ovlap8       
                    = (double *)cmalloc(nstate2_tot*sizeof(double))-1;

  num += 8*nstate2_tot;

  cpscr->cpscr_ovmat.state_vec1       
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.state_vec2       
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.state_vec3       
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.state_vec4       
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.state_vec5       
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.state_vec6       
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;
  cpscr->cpscr_ovmat.state_vec7       
                    = (double *)cmalloc(nstate_tot*sizeof(double))-1;

  num += 7*nstate_tot;

/*------------------------------------------------------------------*/
/* thermostat */
  if(massiv_flag==0){
    cpscr->cpscr_therm.sc_cp   
                     = (double *)cmalloc(num_c_nhc1*sizeof(double))-1;
    cpscr->cpscr_therm.coef_kin
                     = (double *)cmalloc(num_c_nhc1*sizeof(double))-1;
    num += 2*num_c_nhc1;
  }/*endif*/

/*------------------------------------------------------------------*/
/* dual_pme map  */

  if(cp_dual_grid_opt_on == 2){

   cpscr->cpscr_dual_pme.iatemp  = (int *) 
                                 cmalloc(nkf1_cp_box*sizeof(int))-1;
   cpscr->cpscr_dual_pme.ibtemp  = (int *) 
                                 cmalloc(nkf2_cp_box*sizeof(int))-1;
   cpscr->cpscr_dual_pme.ictemp  = (int *) 
                                 cmalloc(nkf3_cp_box*sizeof(int))-1;

   num += (nkf1_cp_box + nkf2_cp_box + nkf3_cp_box);
/*---*/
   cpscr->cpscr_dual_pme.igrid_a = (int **) 
                                 cmalloc(n_interp_pme_dual*sizeof(int *))-1;
   cpscr->cpscr_dual_pme.igrid_b = (int **) 
                                 cmalloc(n_interp_pme_dual*sizeof(int *))-1;
   cpscr->cpscr_dual_pme.igrid_c = (int **) 
                                 cmalloc(n_interp_pme_dual*sizeof(int *))-1;

   for(i=1; i<= n_interp_pme_dual; i++){
   cpscr->cpscr_dual_pme.igrid_a[i] = (int *) 
                                 cmalloc(nkf1_cp_box*sizeof(int))-1;
   cpscr->cpscr_dual_pme.igrid_b[i] = (int *) 
                                 cmalloc(nkf2_cp_box*sizeof(int))-1;
   cpscr->cpscr_dual_pme.igrid_c[i] = (int *) 
                                 cmalloc(nkf3_cp_box*sizeof(int))-1;
   }

   num += ((nkf1_cp_box + nkf2_cp_box + nkf3_cp_box)*n_interp_pme_dual);
/*---*/

   cpscr->cpscr_dual_pme.igrid_now  = (int **)
                             cmalloc(n_interp_pme_dual*sizeof(int *))-1;

   for(i=1 ; i<= n_interp_pme_dual; i++){
   cpscr->cpscr_dual_pme.igrid_now[i]  = (int *)
                             cmalloc(nkf1_cp_box*sizeof(int ))-1;
   }
 
   num += (n_interp_pme_dual*nkf1_cp_box);
/*---*/ 

   cpscr->cpscr_dual_pme.a_pme = (double *) 
                               cmalloc(nkf1_cp_box*sizeof(double))-1;
   cpscr->cpscr_dual_pme.b_pme = (double *) 
                               cmalloc(nkf2_cp_box*sizeof(double))-1;
   cpscr->cpscr_dual_pme.c_pme = (double *) 
                               cmalloc(nkf3_cp_box*sizeof(double))-1;

   cpscr->cpscr_dual_pme.frac_a = (double *) 
                               cmalloc(nkf1_cp_box*sizeof(double))-1;
   cpscr->cpscr_dual_pme.frac_b = (double *) 
                               cmalloc(nkf2_cp_box*sizeof(double))-1;
   cpscr->cpscr_dual_pme.frac_c = (double *) 
                               cmalloc(nkf3_cp_box*sizeof(double))-1;

  num += ((nkf1_cp_box + nkf2_cp_box + nkf3_cp_box)*2);
/*---*/ 

   cpscr->cpscr_dual_pme.aj  = (double *)
                             cmalloc(n_interp_pme_dual*sizeof(double ))-1;
   cpscr->cpscr_dual_pme.rn  = (double *)
                             cmalloc(n_interp_pme_dual*sizeof(double ))-1;
   cpscr->cpscr_dual_pme.rn1 = (double *)
                             cmalloc(n_interp_pme_dual*sizeof(double ))-1;

  num += (n_interp_pme_dual*3);
/*---*/ 

   cpscr->cpscr_dual_pme.mn_a = (double **)
                             cmalloc(n_interp_pme_dual*sizeof(double *))-1;
   cpscr->cpscr_dual_pme.mn_b = (double **)
                             cmalloc(n_interp_pme_dual*sizeof(double *))-1;
   cpscr->cpscr_dual_pme.mn_c = (double **)
                             cmalloc(n_interp_pme_dual*sizeof(double *))-1;

   for(i=1; i<= n_interp_pme_dual; i++){
   cpscr->cpscr_dual_pme.mn_a[i] = (double *)
                             cmalloc(nkf1_cp_box*sizeof(double ))-1;
   cpscr->cpscr_dual_pme.mn_b[i] = (double *)
                             cmalloc(nkf2_cp_box*sizeof(double ))-1;
   cpscr->cpscr_dual_pme.mn_c[i] = (double *)
                             cmalloc(nkf3_cp_box*sizeof(double ))-1;
   }
   num += ((nkf1_cp_box + nkf2_cp_box + nkf3_cp_box)*n_interp_pme_dual);
/*---*/ 

   cpscr->cpscr_dual_pme.ua = (double **)
                             cmalloc(n_interp_pme_dual*sizeof(double *))-1;
   cpscr->cpscr_dual_pme.ub = (double **)
                             cmalloc(n_interp_pme_dual*sizeof(double *))-1;
   cpscr->cpscr_dual_pme.uc = (double **)
                             cmalloc(n_interp_pme_dual*sizeof(double *))-1;


  for(i=1; i<= n_interp_pme_dual; i++){
   cpscr->cpscr_dual_pme.ua[i] = (double *)
                             cmalloc(nkf1_cp_box*sizeof(double ))-1;
   cpscr->cpscr_dual_pme.ub[i] = (double *)
                             cmalloc(nkf2_cp_box*sizeof(double ))-1;
   cpscr->cpscr_dual_pme.uc[i] = (double *)
                             cmalloc(nkf3_cp_box*sizeof(double ))-1;
  }

   num += ((nkf1_cp_box + nkf2_cp_box + nkf3_cp_box)*n_interp_pme_dual);
/*---*/ 

   cpscr->cpscr_dual_pme.bw_r = (double *) 
                              cmalloc(ncoef_l_proc_max_mall*sizeof(double))-1;
   cpscr->cpscr_dual_pme.bw_i = (double *) 
                              cmalloc(ncoef_l_proc_max_mall*sizeof(double))-1;

   num += (3*ncoef_l_proc_max);
/*---*/     
  }/*endif cp_dual_grid_opt*/ 

/*-----------------------------------------------------------------------*/
/* Wannier scratch                                                       */

  if(cpscr->cpscr_wannier.cp_wannier_on==1){
    cpscr->cpscr_wannier.U_final = (double *) cmalloc(ndim_wannier*sizeof(double))-1;
    cpscr->cpscr_wannier.g       = (double *) cmalloc(ndim_wannier*sizeof(double))-1;
    cpscr->cpscr_wannier.gg      = (double *) cmalloc(ndim_wannier*sizeof(double))-1;
    cpscr->cpscr_wannier.diag    = (double *) cmalloc(ndim_wannier*sizeof(double))-1;
    cpscr->cpscr_wannier.scr     = 
                   (double *) cmalloc((ndim_wannier*(2*mm+1)+2*mm)*sizeof(double))-1;
    cpscr->cpscr_wannier.iprint  = (int *)cmalloc(2*sizeof(int))-1;


    num += 4*ndim_wannier+ndim_wannier*(2*mm+1)+2*mm;

    cpscr->cpscr_wannier.A       = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.R_real  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.R_imag  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.U_real  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.U_imag  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.U_tmp1  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.U_tmp2  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.U_tmp3  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.U_tmp4  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.Bt_real = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.Bt_imag = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    cpscr->cpscr_wannier.M_real  = cmall_mat(1,nstate_up_tot,1,nstate_up_tot);
    
    num += 12*ndim_wannier;

    cpscr->cpscr_wannier.real    = (double *)cmalloc(nstate_up_tot*sizeof(double))-1;
    cpscr->cpscr_wannier.imag    = (double *)cmalloc(nstate_up_tot*sizeof(double))-1;
    cpscr->cpscr_wannier.D       = (double *)cmalloc(nstate_up_tot*sizeof(double))-1;
    cpscr->cpscr_wannier.norm    = (double *)cmalloc(nstate_up_tot*sizeof(double))-1;

    cpscr->cpscr_wannier.phi     = cmall_mat(1,nstate_up_tot,1,3);
    cpscr->cpscr_wannier.HMatrix = cmall_mat(1,3,1,3);
    printf("CP_WANNIER memory allocation \n");

    num += 7*nstate_up_tot;
  }




/*=========================================================================*/
/* II)Malloc PME variables for structure factor                            */

 if(cpscr->cpscr_atom_pme.pme_on == 1 && cp_dual_grid_opt_on == 2){

  nlen_pme = cpscr->cpscr_atom_pme.nlen_pme;

  cpscr->cpscr_atom_pme.nkf1     = cp_para_fft_pkg3d_lg->nkf1;
  cpscr->cpscr_atom_pme.nkf2     = cp_para_fft_pkg3d_lg->nkf2;
  cpscr->cpscr_atom_pme.nkf3     = cp_para_fft_pkg3d_lg->nkf3;

  cpscr->cpscr_atom_pme.ngrid_a  = cp_para_fft_pkg3d_lg->nkf1;
  cpscr->cpscr_atom_pme.ngrid_b  = cp_para_fft_pkg3d_lg->nkf2;
  cpscr->cpscr_atom_pme.ngrid_c  = cp_para_fft_pkg3d_lg->nkf3;

  cpscr->cpscr_atom_pme.iatemp   = (int *) cmalloc(nlen_pme*sizeof(int))-1;
  cpscr->cpscr_atom_pme.ibtemp   = (int *) cmalloc(nlen_pme*sizeof(int))-1;
  cpscr->cpscr_atom_pme.ictemp   = (int *) cmalloc(nlen_pme*sizeof(int))-1;

  num += 3*nlen_pme;
  
  pme_nkf3 = cpscr->cpscr_atom_pme.nkf3;

  cpscr->cpscr_atom_pme.nc      = (int *) cmalloc(pme_nkf3*sizeof(int))-1;
  cpscr->cpscr_atom_pme.ioff_c  = (int *) cmalloc(pme_nkf3*sizeof(int))-1;

  num += 2*pme_nkf3;

  ninterp_pme = cpscr->cpscr_atom_pme.n_interp;

  cpscr->cpscr_atom_pme.igrid_a    = cmall_int_mat(1,ninterp_pme,1,nlen_pme);
  cpscr->cpscr_atom_pme.igrid_b    = cmall_int_mat(1,ninterp_pme,1,nlen_pme);
  cpscr->cpscr_atom_pme.igrid_c    = cmall_int_mat(1,ninterp_pme,1,nlen_pme);
  cpscr->cpscr_atom_pme.igrid_now  = cmall_int_mat(1,ninterp_pme,1,nlen_pme);

  num += 4*ninterp_pme*nlen_pme;

  cpscr->cpscr_atom_pme.frac_a  = (double *) cmalloc(nlen_pme*sizeof(double))-1;
  cpscr->cpscr_atom_pme.frac_b  = (double *) cmalloc(nlen_pme*sizeof(double))-1;
  cpscr->cpscr_atom_pme.frac_c  = (double *) cmalloc(nlen_pme*sizeof(double))-1;

  num += 3*nlen_pme;

  nmall = ninterp_pme*nlen_pme;

  cpscr->cpscr_atom_pme.aj   = (double *) cmalloc(nmall*sizeof(double))-1;
  cpscr->cpscr_atom_pme.rn   = (double *) cmalloc(nmall*sizeof(double))-1;
  cpscr->cpscr_atom_pme.rn1  = (double *) cmalloc(nmall*sizeof(double))-1;

  num += 3*nmall;

  cpscr->cpscr_atom_pme.ua  = cmall_mat(1,ninterp_pme,1,nlen_pme);
  cpscr->cpscr_atom_pme.ub  = cmall_mat(1,ninterp_pme,1,nlen_pme);
  cpscr->cpscr_atom_pme.uc  = cmall_mat(1,ninterp_pme,1,nlen_pme);

  cpscr->cpscr_atom_pme.mn_a  = cmall_mat(1,ninterp_pme,1,nlen_pme);
  cpscr->cpscr_atom_pme.mn_b  = cmall_mat(1,ninterp_pme,1,nlen_pme);
  cpscr->cpscr_atom_pme.mn_c  = cmall_mat(1,ninterp_pme,1,nlen_pme);

  cpscr->cpscr_atom_pme.dmn_a  = cmall_mat(1,ninterp_pme,1,nlen_pme);
  cpscr->cpscr_atom_pme.dmn_b  = cmall_mat(1,ninterp_pme,1,nlen_pme);
  cpscr->cpscr_atom_pme.dmn_c  = cmall_mat(1,ninterp_pme,1,nlen_pme);

  cpscr->cpscr_atom_pme.qgrid_now  = cmall_mat(1,ninterp_pme,1,nlen_pme);

  num += (10*ninterp_pme*nlen_pme);

  mall_size = 2*(cp_para_fft_pkg3d_lg->nkf1*
                 cp_para_fft_pkg3d_lg->nkf2*
                 cp_para_fft_pkg3d_lg->nkf3);
  
  cpscr->cpscr_atom_pme.qgrid  = (double *) cmalloc(mall_size*sizeof(double))-1;
  cpscr->cpscr_atom_pme.qgrid_scr  = (double *)
                                   cmalloc(mall_size*sizeof(double))-1;

  /*Lst: Lth: nktot                     */
  cpscr->cpscr_atom_pme.qgrid_tmp_real = (double *) 
                                     cmalloc((nfft_up+1)*sizeof(double))-1;
  cpscr->cpscr_atom_pme.qgrid_tmp_imag = (double *) 
                                     cmalloc((nfft_up+1)*sizeof(double))-1;

  cpscr->cpscr_atom_pme.bw_r = (double *) 
                              cmalloc((nfft_up+1)*sizeof(double))-1;
  cpscr->cpscr_atom_pme.bw_i = (double *) 
                              cmalloc((nfft_up+1)*sizeof(double))-1;
  cpscr->cpscr_atom_pme.bweight_tot = (double *) 
                              cmalloc((nfft_up+1)*sizeof(double))-1;

  num += 2*(nfft_up+1);
 }/*endif*/


/*===========================================================================*/
/* III) Output */


  now_memory = num*sizeof(double)*(1.0e-6);
  *tot_memory += now_memory;

/*--------------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void mall_atm_forc_scr_frag(int natm_tot,FOR_SCR *for_scr,int pme_on,
                       int ilnk_lst,int iver_lst, int lnk_ver_update,
                       double *tot_memory,int np_forc,int myid,MPI_Comm world)

/*=========================================================================*/
  {/*begin routine*/
/*=========================================================================*/
/*         Local variable declarations                                     */

  double now_memory;
  int nlen_mall    = for_scr->nlen;
  int natm_mall    = natm_tot;

/*=========================================================================*/
/* I) Malloc size calculation */

  if((nlen_mall % 2)==0){nlen_mall++;} 
  if((natm_mall % 2)==0){natm_mall += 1;}

/*===========================================================================*/
/* II) Malloc the vectors  */

  for_scr->wght_lnk  = 1.0;
  for_scr->j_index   = (int *)cmalloc(nlen_mall*sizeof(int))-1;
  for_scr->i_index   = (int *)cmalloc(nlen_mall*sizeof(int))-1;
  for_scr->j_indext  = (int *)cmalloc(nlen_mall*sizeof(int))-1;
  for_scr->i_indext  = (int *)cmalloc(nlen_mall*sizeof(int))-1;
  for_scr->i_indext  = (int *)cmalloc(nlen_mall*sizeof(int))-1;
  for_scr->i_indext2 = (int *)cmalloc(nlen_mall*sizeof(int))-1;
  for_scr->num_brk_i = (int *)cmalloc(nlen_mall*sizeof(int))-1;
  for_scr->iexcl     = (int *)cmalloc(natm_mall*sizeof(int))-1;
  for_scr->index_atm = (int *)cmalloc(natm_mall*sizeof(int))-1;

  now_memory         = ( (9*nlen_mall+natm_mall)*(sizeof(int)) );

  if(pme_on==1&&np_forc>1){
    if((ilnk_lst==0)&&((lnk_ver_update+iver_lst)!=2)){
      for_scr->i_lnk   = (list_int *)cmalloc(natm_mall*sizeof(list_int))-1;
      for_scr->j_lnk   = (list_int *)cmalloc(natm_mall*sizeof(list_int))-1;

      now_memory      += ( (2*natm_mall)*(sizeof(list_int)) );
    }/*endif*/
  }/*endif*/

/*===========================================================================*/
/* III) Output */

  now_memory  *= 1.0e-06;
  *tot_memory += now_memory;

/*-------------------------------------------------------------------------*/
  }/*end routine */
/*=========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void mall_pressure_frag(CLASS *class,GENERAL_DATA *general_data)

/*======================================================================*/
     {/*begin routine*/
/*======================================================================*/
/*          Local variable declarations                                */

  double now_memory;
  int myid = class->communicate.myid;

/*======================================================================*/
/*         Output                                                       */


  now_memory   = (sizeof(int)*0 + sizeof(double)*29*9)*1.e-06;
  class->tot_memory  += now_memory;

  /*=======================================================================*/
  /* 0) Average pressure Memory */

  general_data->stat_avg.apten_out = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->stat_avg.apten     = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->stat_avg.aipten    = (double *)cmalloc((size_t)9*sizeof(double))-1;

  general_data->ptens.tvten        = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten        = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_tot    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_tmp    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_tmp_res = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_t1 = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_t2 = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_t3 = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_t4 = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.count_inc    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_a  = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_old = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_glob = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_whol = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_std  = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_mode     = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_mode_tot = (double *)cmalloc((size_t)9*sizeof(double))-1;

  /*=======================================================================*/
  /* 0) Baro/cell/par_rahman Memory */

  general_data->par_rahman.vgmat    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.vgmat_glob= (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.vgmat_g  = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.fgmat_p  = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.fgmat_v  = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.roll_mtv = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.roll_mtx = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.roll_mtvv = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.vtemps   = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.veigv    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.veig     = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.vexpdt   = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.vsindt   = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.vtempx   = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.vtempv   = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.hmat_t   = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.hmato    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.fv1    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.fv2    = (double *)cmalloc((size_t)9*sizeof(double))-1;

  general_data->baro.hmato       = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->cell.hmat        = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->cell.hmati       = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->cell.hmat_ewd    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->cell.hmat_ewd_cp = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->cell.hmat_cp     = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->cell.hmati_cp    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->cell.cp_box_center = (double *) cmalloc((size_t)3*sizeof(double ))-1;
  general_data->cell.cp_box_center_rel = (double *) cmalloc((size_t)3*sizeof(double ))-1;
  general_data->cell.cp_vbox_center = (double *) cmalloc((size_t)3*sizeof(double ))-1;
  general_data->cell.cp_fbox_center = (double *) cmalloc((size_t)3*sizeof(double ))-1;

  /*=======================================================================*/
  /* 0) More output */

/*------------------------------------------------------------------------*/
} /*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void close_intra_params_frag(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                        GHOST_ATOMS *ghost_atoms,
                        ATOMMAPS *atommaps,
                        BUILD_INTRA *build_intra,
                        BONDED *bonded,
                        NULL_INTER_PARSE *null_inter_parse,
                        double *tot_memory,SIMOPTS *simopts,int np_forc)

/*=======================================================================*/
/*        Begin routine                                                 */
{/*begin routine*/
/*=======================================================================*/
/*           Local Variables                                             */

  int i,ntyp_pow,iii,nsec,inow,know;
  double now_memory;
  int natm_typ_mall;
  int natm_mall;
  int nghost_mall;
  int nfreeze_mall;
  int ncomp_mall;
  int nchrg_mall;
  int ngrp_21_mall;
  int ngrp_typ_21_mall;
  int ngrp_33_mall;
  int ngrp_typ_33_mall;
  int ngrp_watt_33_mall;
  int ngrp_typ_watt_33_mall;
  int ngrp_43_mall;
  int ngrp_typ_43_mall;
  int ngrp_23_mall;
  int ngrp_typ_23_mall;
  int ngrp_46_mall;
  int ngrp_typ_46_mall;
  int nbond_pow_mall;
  int nbond_typ_pow_mall;
  int nbond_con_mall;
  int nbond_typ_con_mall;
  int nbend_pow_mall;
  int nbend_typ_pow_mall;
  int nbend_con_mall;
  int nbend_typ_con_mall;
  int nbend_bnd_mall;
  int nbend_bnd_typ_mall;
  int ntors_pow_mall;
  int ntors_typ_pow_mall;
  int ntors_con_mall;
  int ntors_typ_con_mall;
  int nonfo_mall;
  int nonfo_typ_mall;
  int pimd_on;
  int nghost_old,nfreeze_old,ncomp_old;
  int pi_beads = clatoms_info->pi_beads;
  int pimd_freez_typ = atommaps->pimd_freez_typ;

/*======================================================================*/
/* 0) Make the length of each vector an odd number                      */

  nghost_old          = build_intra->nghost_tot_max;
  nfreeze_old         = build_intra->nfreeze_max;
  ncomp_old           = NCOEF_GHOST_MAX;
  natm_mall           = clatoms_info->natm_tot;
  natm_typ_mall       = atommaps->natm_typ;
  nghost_mall         = ghost_atoms->nghost_tot;
  nfreeze_mall        = atommaps->nfreeze;
  ncomp_mall          = ghost_atoms->natm_comp_max;
  ngrp_21_mall        = bonded->grp_bond_con.num_21;
  ngrp_33_mall        = bonded->grp_bond_con.num_33;
  ngrp_watt_33_mall   = bonded->grp_bond_watts.num_33;
  ngrp_43_mall        = bonded->grp_bond_con.num_43;
  ngrp_23_mall        = bonded->grp_bond_con.num_23;
  ngrp_46_mall        = bonded->grp_bond_con.num_46;
  ngrp_typ_21_mall    = bonded->grp_bond_con.ntyp_21;
  ngrp_typ_33_mall    = bonded->grp_bond_con.ntyp_33;
  ngrp_typ_watt_33_mall  = bonded->grp_bond_watts.ntyp_33;
  ngrp_typ_43_mall    = bonded->grp_bond_con.ntyp_43;
  ngrp_typ_23_mall    = bonded->grp_bond_con.ntyp_23;
  ngrp_typ_46_mall    = bonded->grp_bond_con.ntyp_46;
  nbond_pow_mall      = bonded->bond.npow;
  nbond_typ_pow_mall  = bonded->bond.ntyp_pow;
  nbond_con_mall      = bonded->bond.ncon;
  nbond_typ_con_mall  = bonded->bond.ntyp_con;
  nbend_pow_mall      = bonded->bend.npow;
  nbend_typ_pow_mall  = bonded->bend.ntyp_pow;
  nbend_con_mall      = bonded->bend.ncon;
  nbend_typ_con_mall  = bonded->bend.ntyp_con;
  nbend_bnd_mall      = bonded->bend_bnd.num;
  nbend_bnd_typ_mall  = bonded->bend_bnd.ntyp;
  ntors_pow_mall      = bonded->tors.npow;
  ntors_typ_pow_mall  = bonded->tors.ntyp_pow;
  ntors_con_mall      = bonded->tors.ncon;
  ntors_typ_con_mall  = bonded->tors.ntyp_con;
  nonfo_mall          = bonded->onfo.num;
  nonfo_typ_mall      = bonded->onfo.ntyp;
  if((natm_mall           % 2)==0){natm_mall      += 1;}
  nsec=natm_mall/np_forc;
  if((natm_mall%np_forc)!=0){nsec++;}
  natm_mall = nsec*np_forc;
  if((natm_typ_mall       % 2)==0){natm_typ_mall  += 1;}
  if(((nghost_mall        % 2)==0)
    &&(nghost_mall        !=0)){nghost_mall++;}
  if(((nfreeze_mall        % 2)==0)
    &&(nfreeze_mall        !=0)){nfreeze_mall++;}
  if(((ncomp_mall      % 2)==0)
    &&(ncomp_mall     !=0)){ncomp_mall++;}
  if(((ngrp_21_mall        % 2)==0)
    &&(ngrp_21_mall       !=0)){ngrp_21_mall++;}
  if(((ngrp_33_mall        % 2)==0)
    &&(ngrp_33_mall       !=0)){ngrp_33_mall++;}
  if(((ngrp_watt_33_mall        % 2)==0)
    &&(ngrp_watt_33_mall       !=0)){ngrp_watt_33_mall++;}
  if(((ngrp_43_mall        % 2)==0)
    &&(ngrp_43_mall       !=0)){ngrp_43_mall++;}
  if(((ngrp_23_mall        % 2)==0)
    &&(ngrp_23_mall       !=0)){ngrp_23_mall++;}
  if(((ngrp_46_mall        % 2)==0)
    &&(ngrp_46_mall       !=0)){ngrp_46_mall++;}
  if(((ngrp_typ_21_mall    % 2)==0)
    &&(ngrp_typ_21_mall   !=0)){ngrp_typ_21_mall++;}
  if(((ngrp_typ_33_mall    % 2)==0)
    &&(ngrp_typ_33_mall   !=0)){ngrp_typ_33_mall++;}
  if(((ngrp_typ_watt_33_mall    % 2)==0)
    &&(ngrp_typ_watt_33_mall   !=0)){ngrp_typ_watt_33_mall++;}
  if(((ngrp_typ_43_mall    % 2)==0)
    &&(ngrp_typ_43_mall   !=0)){ngrp_typ_43_mall++;}
  if(((ngrp_typ_23_mall    % 2)==0)
    &&(ngrp_typ_23_mall   !=0)){ngrp_typ_23_mall++;}
  if(((ngrp_typ_46_mall    % 2)==0)
    &&(ngrp_typ_46_mall   !=0)){ngrp_typ_46_mall++;}
  if(((nbond_pow_mall     % 2)==0)
    &&(nbond_pow_mall     !=0)){nbond_pow_mall    += 1;}
  if(((nbond_typ_pow_mall % 2)==0)
    &&(nbond_typ_pow_mall !=0)){nbond_typ_pow_mall+= 1;}
  if(((nbond_con_mall     % 2)==0)
    &&(nbond_con_mall     !=0)){nbond_pow_mall    += 1;}
  if(((nbond_typ_con_mall % 2)==0)
    &&(nbond_typ_con_mall !=0)){nbond_typ_pow_mall+= 1;}
  if(((nbend_pow_mall     % 2)==0)
    &&(nbend_pow_mall     !=0)){nbend_pow_mall    += 1;}
  if(((nbend_typ_pow_mall % 2)==0)
    &&(nbend_typ_pow_mall !=0)){nbend_typ_pow_mall+= 1;}
  if(((nbend_con_mall     % 2)==0)
    &&(nbend_con_mall     !=0)){nbend_pow_mall    += 1;}
  if(((nbend_typ_con_mall % 2)==0)
    &&(nbend_typ_con_mall !=0)){nbend_typ_pow_mall+= 1;}
  if(((nbend_bnd_mall     % 2)==0)
    &&(nbend_bnd_mall     !=0)){nbend_bnd_mall    += 1;}
  if(((nbend_bnd_typ_mall % 2)==0)
    &&(nbend_bnd_typ_mall !=0)){nbend_bnd_typ_mall+= 1;}
  if(((ntors_pow_mall     % 2)==0)
    &&(ntors_pow_mall     !=0)){ntors_pow_mall    += 1;}
  if(((ntors_typ_pow_mall % 2)==0)
    &&(ntors_typ_pow_mall !=0)){ntors_typ_pow_mall+= 1;}
  if(((ntors_con_mall     % 2)==0)
    &&(ntors_con_mall     !=0)){ntors_pow_mall    += 1;}
  if(((ntors_typ_con_mall % 2)==0)
    &&(ntors_typ_con_mall !=0)){ntors_typ_pow_mall+= 1;}
  if(((nonfo_mall         % 2)==0)
    &&(nonfo_mall         !=0)){nonfo_mall        += 1;}
  if(((nonfo_typ_mall     % 2)==0)
    &&(nonfo_typ_mall     !=0)){nonfo_typ_mall    += 1;}


/*======================================================================*/
/* I) Atm stuff: Note constrain flag set here                           */
/*               and charged atoms are counted                          */
  bonded->constrnt.iconstrnt=0;
  if(bonded->bond.ncon > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->grp_bond_con.num_21 > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->grp_bond_con.num_23 > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->grp_bond_con.num_33 > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->grp_bond_watts.num_33 > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->grp_bond_con.num_43 > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->grp_bond_con.num_46 > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->bend.ncon > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->tors.ncon > 0) bonded->constrnt.iconstrnt=1;


  clatoms_info->ichrg      = (int *)cmalloc(natm_mall*sizeof(int))-1;
  clatoms_info->nfree      = 0;
  clatoms_info->nfree_pimd = 0;
  for(i=1;i<=atommaps->nmol_typ;i++){
    inow                 = atommaps->nfree_1mol_jmol_typ[i]
                          *atommaps->nmol_jmol_typ[i];
    know                 = 3*atommaps->natm_1mol_jmol_typ[i]
                            *atommaps->nmol_jmol_typ[i];
    clatoms_info->nfree += inow;
    if(pimd_freez_typ==2){
      clatoms_info->nfree_pimd += inow*pi_beads;
    }else{
      clatoms_info->nfree_pimd += (inow + know*(pi_beads-1));
    }/*endif*/
  }/*endfor*/

  clatoms_info->nchrg = 0;
  for(i=1;i<=clatoms_info->natm_tot;i++){
    clatoms_info->mass[i] *= PROT_MASS;
    if(clatoms_info->q[i] != 0) {
     clatoms_info->nchrg++;
     clatoms_info->ichrg[(clatoms_info->nchrg)] = i;
   }
  }
  nchrg_mall = clatoms_info->nchrg;
  if((nchrg_mall % 2)==0){nchrg_mall += 1;}


  now_memory = ((double)(((natm_mall)*sizeof(double)*13+sizeof(int)*4)+
                         natm_typ_mall*(sizeof(double)*0+sizeof(int)*25))
                *1.e-06);
  *tot_memory += now_memory;

  atommaps->atm_typ = (NAME *)crealloc(&(atommaps->atm_typ)[1],
                                       natm_typ_mall*sizeof(NAME))-1;
  clatoms_info->mass     = (double *)crealloc(&(clatoms_info->mass)[1],
                                         natm_mall*sizeof(double))-1;
  clatoms_info->q        = (double *)crealloc(&(clatoms_info->q)[1],
                                         natm_mall*sizeof(double))-1;
  clatoms_info->cp_vlnc_up  = (int *)crealloc(&(clatoms_info->cp_vlnc_up)[1],
                                         natm_mall*sizeof(int))-1;
  clatoms_info->cp_vlnc_dn  = (int *)crealloc(&(clatoms_info->cp_vlnc_dn)[1],
                                         natm_mall*sizeof(int))-1;
  clatoms_info->cp_atm_flag  = (int *)crealloc(&(clatoms_info->cp_atm_flag)[1],
                                         natm_mall*sizeof(int))-1;
  clatoms_info->ichrg  = (int *)crealloc(&(clatoms_info->ichrg)[1],
                                    nchrg_mall*sizeof(int))-1;
  clatoms_info->alp_pol  = (double *)crealloc(&(clatoms_info->alp_pol)[1],
                                         natm_mall*sizeof(double))-1;
  clatoms_info->b_neut   = (double *)crealloc(&(clatoms_info->b_neut)[1],
                                         natm_mall*sizeof(double))-1;
  clatoms_info->text_atm = (double *)crealloc(&(clatoms_info->text_atm)[1],
                                         natm_mall*sizeof(double))-1;
  atommaps->iatm_mol_typ  = (int *)crealloc(&(atommaps->iatm_mol_typ)[1],
                                            natm_mall*sizeof(int))-1;
  atommaps->iatm_atm_typ  = (int *)crealloc(&(atommaps->iatm_atm_typ)[1],
                                            natm_mall*sizeof(int))-1;
  atommaps->iatm_res_typ  = (int *)crealloc(&(atommaps->iatm_res_typ)[1],
                                            natm_mall*sizeof(int))-1;
  atommaps->iatm_mol_num  = (int *)crealloc(&(atommaps->iatm_mol_num)[1],
                                            natm_mall*sizeof(int))-1;
  atommaps->iatm_res_num  = (int *)crealloc(&(atommaps->iatm_res_num)[1],
                                            natm_mall*sizeof(int))-1;
  pimd_on = simopts->pimd + simopts->cp_pimd
          + simopts->cp_wave_pimd + simopts->cp_wave_min_pimd
          + simopts->debug_pimd + simopts->debug_cp_pimd;

  clatoms_info->xold= (double *)cmalloc(natm_mall*sizeof(double))-1;
  clatoms_info->yold= (double *)cmalloc(natm_mall*sizeof(double))-1;
  clatoms_info->zold= (double *)cmalloc(natm_mall*sizeof(double))-1;
  if(pi_beads>1||pimd_on==1){
   clatoms_info->xmod = (double *)cmalloc(natm_mall*sizeof(double))-1;
   clatoms_info->ymod = (double *)cmalloc(natm_mall*sizeof(double))-1;
   clatoms_info->zmod = (double *)cmalloc(natm_mall*sizeof(double))-1;
  }/*endif*/
  if(clatoms_info->pi_beads>1||pimd_on==1){
   clatoms_info->prekf = (double *)cmalloc(natm_mall*sizeof(double))-1;
  }
  clatoms_info->roll_sc   = (double *)cmalloc(natm_mall*sizeof(double))-1;
/*========================================================================*/
/* I) Ghost_atm stuff */

  ghost_atoms->ighost_map     = (int *) crealloc(&(ghost_atoms->ighost_map)[1],
                                                 nghost_mall*sizeof(int))-1;
  ghost_atoms->natm_comp      = (int *) crealloc(&(ghost_atoms->natm_comp)[1],
                                                 nghost_mall*sizeof(int))-1;
  atommaps->ighost_flag       = (int *)   crealloc(&(atommaps->ighost_flag)[1],
                                                 natm_mall*sizeof(int))-1;

  /* need to set up reallocate matrix */
  ghost_atoms->iatm_comp      = creall_int_mat(ghost_atoms->iatm_comp,
                                               1,ncomp_old ,1,nghost_old,
                                               1,ncomp_mall,1,nghost_mall);
  ghost_atoms->coef           = creall_mat(ghost_atoms->coef,
                                           1,ncomp_old ,1,nghost_old,
                                           1,ncomp_mall,1,nghost_mall);
/*========================================================================*/
/* I) Freeze_atm stuff */

  atommaps->freeze_map     = (int *) crealloc(&(atommaps->freeze_map)[1],
                                                 nfreeze_mall*sizeof(int))-1;
  atommaps->freeze_flag       = (int *)   crealloc(&(atommaps->freeze_flag)[1],
                                                 natm_mall*sizeof(int))-1;
  atommaps->atom_label       = (int *)   crealloc(&(atommaps->atom_label)[1],
                                                 natm_mall*sizeof(int))-1;
/*========================================================================*/
/* II) Grp constrnts */



  bonded->grp_bond_con.j1_21 =
                    (int *)crealloc(&(bonded->grp_bond_con.j1_21)[1],
                                             ngrp_21_mall*sizeof(int))-1;
  bonded->grp_bond_con.j2_21 =
                    (int *)   crealloc(&(bonded->grp_bond_con.j2_21)[1],
                                             ngrp_21_mall*sizeof(int))-1;
  bonded->grp_bond_con.jtyp_21 =
                    (int *) crealloc(&(bonded->grp_bond_con.jtyp_21)[1],
                                             ngrp_21_mall*sizeof(int))-1;
  bonded->grp_bond_con.j1_33 =
                    (int *)crealloc(&(bonded->grp_bond_con.j1_33)[1],
                                             ngrp_33_mall*sizeof(int))-1;
  bonded->grp_bond_con.j2_33 =
                    (int *)   crealloc(&(bonded->grp_bond_con.j2_33)[1],
                                             ngrp_33_mall*sizeof(int))-1;
  bonded->grp_bond_con.j3_33 =
                    (int *)   crealloc(&(bonded->grp_bond_con.j3_33)[1],
                                             ngrp_33_mall*sizeof(int))-1;
  bonded->grp_bond_con.jtyp_33 =
                    (int *) crealloc(&(bonded->grp_bond_con.jtyp_33)[1],
                                             ngrp_33_mall*sizeof(int))-1;
  bonded->grp_bond_watts.j1_33 =
                    (int *)crealloc(&(bonded->grp_bond_watts.j1_33)[1],
                                        ngrp_watt_33_mall*sizeof(int))-1;
  bonded->grp_bond_watts.j2_33 =
                    (int *)   crealloc(&(bonded->grp_bond_watts.j2_33)[1],
                                        ngrp_watt_33_mall*sizeof(int))-1;
  bonded->grp_bond_watts.j3_33 =
                    (int *)   crealloc(&(bonded->grp_bond_watts.j3_33)[1],
                                        ngrp_watt_33_mall*sizeof(int))-1;
  bonded->grp_bond_watts.jtyp_33 =
                    (int *) crealloc(&(bonded->grp_bond_watts.jtyp_33)[1],
                                        ngrp_watt_33_mall*sizeof(int))-1;
  bonded->grp_bond_con.j1_23 =
                    (int *)   crealloc(&(bonded->grp_bond_con.j1_23)[1],
                                             ngrp_23_mall*sizeof(int))-1;
  bonded->grp_bond_con.j2_23 =
                    (int *)   crealloc(&(bonded->grp_bond_con.j2_23)[1],
                                             ngrp_23_mall*sizeof(int))-1;
  bonded->grp_bond_con.j3_23 =
                    (int *)   crealloc(&(bonded->grp_bond_con.j3_23)[1],
                                             ngrp_23_mall*sizeof(int))-1;
  bonded->grp_bond_con.jtyp_23 =
                    (int *) crealloc(&(bonded->grp_bond_con.jtyp_23)[1],
                                             ngrp_23_mall*sizeof(int))-1;
  bonded->grp_bond_con.j1_43 =
                    (int *)crealloc(&(bonded->grp_bond_con.j1_43)[1],
                                             ngrp_43_mall*sizeof(int))-1;
  bonded->grp_bond_con.j2_43 =
                    (int *)   crealloc(&(bonded->grp_bond_con.j2_43)[1],
                                             ngrp_43_mall*sizeof(int))-1;
  bonded->grp_bond_con.j3_43 =
                    (int *)   crealloc(&(bonded->grp_bond_con.j3_43)[1],
                                             ngrp_43_mall*sizeof(int))-1;
  bonded->grp_bond_con.j4_43 =
                    (int *)   crealloc(&(bonded->grp_bond_con.j4_43)[1],
                                             ngrp_43_mall*sizeof(int))-1;
  bonded->grp_bond_con.jtyp_43 =
                    (int *) crealloc(&(bonded->grp_bond_con.jtyp_43)[1],
                                             ngrp_43_mall*sizeof(int))-1;
  bonded->grp_bond_con.j1_46 =
                    (int *)   crealloc(&(bonded->grp_bond_con.j1_46)[1],
                                             ngrp_46_mall*sizeof(int))-1;
  bonded->grp_bond_con.j2_46 =
                    (int *)   crealloc(&(bonded->grp_bond_con.j2_46)[1],
                                              ngrp_46_mall*sizeof(int))-1;
  bonded->grp_bond_con.j3_46 =
                    (int *)   crealloc(&(bonded->grp_bond_con.j3_46)[1],
                                             ngrp_46_mall*sizeof(int))-1;
  bonded->grp_bond_con.j4_46 =
                    (int *)   crealloc(&(bonded->grp_bond_con.j4_46)[1],
                                             ngrp_46_mall*sizeof(int))-1;
  bonded->grp_bond_con.jtyp_46 =
                    (int *) crealloc(&(bonded->grp_bond_con.jtyp_46)[1],
                                             ngrp_46_mall*sizeof(int))-1;
  bonded->grp_bond_con.al_21   = cmall_mat(1,1,1,ngrp_21_mall);
  bonded->grp_bond_con.al_33   = cmall_mat(1,3,1,ngrp_33_mall);
  bonded->grp_bond_con.al_43   = cmall_mat(1,3,1,ngrp_43_mall);
  bonded->grp_bond_con.al_23   = cmall_mat(1,2,1,ngrp_23_mall);
  bonded->grp_bond_con.al_46   = cmall_mat(1,6,1,ngrp_46_mall);
  bonded->grp_bond_con.eq_21   = cmall_mat(1,1,1,ngrp_typ_21_mall);
  bonded->grp_bond_con.eq_33   = cmall_mat(1,3,1,ngrp_typ_33_mall);
  bonded->grp_bond_watts.eq_33 = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_con.eq_43   = cmall_mat(1,3,1,ngrp_typ_43_mall);
  bonded->grp_bond_con.eq_23   = cmall_mat(1,2,1,ngrp_typ_23_mall);
  bonded->grp_bond_con.eq_46   = cmall_mat(1,6,1,ngrp_typ_46_mall);

  bonded->grp_bond_watts.cos_thet0_2   =
                (double *)malloc(ngrp_typ_watt_33_mall*sizeof(int))-1;
  bonded->grp_bond_watts.sin_thet0_2   =
                (double *)malloc(ngrp_typ_watt_33_mall*sizeof(int))-1;
  bonded->grp_bond_watts.c_0_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_1_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_2_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_3_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_4_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_5_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_6_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_0_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_1_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_2_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_3_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_4_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_5_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_6_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);

/*======================================================================*/
/*   II) Bonds  */


  now_memory = ((nbond_pow_mall)*(sizeof(double)*0+sizeof(int)*3)+
                (nbond_con_mall)*(sizeof(double)*1 +  sizeof(int)*3)+
                (nbond_typ_pow_mall)*(sizeof(double)*14+ sizeof(int)*0 )+
                (nbond_typ_con_mall)*(sizeof(double)*1 +  sizeof(int)*0)+
                (nbend_pow_mall)*(sizeof(double)*0 +  sizeof(int)*4)+
                (nbend_typ_pow_mall)*(sizeof(double)*29+  sizeof(int)*0)+
                (nbend_con_mall)*(sizeof(double)*1 +  sizeof(int)*4)+
                (nbend_typ_con_mall)*(sizeof(double)*1 +  sizeof(int)*0)+
                (nbend_bnd_mall)*(sizeof(double)*0 +  sizeof(int)*4)+
                (nbend_bnd_typ_mall)*(sizeof(double)*44 +  sizeof(int)*0)+
                (ntors_pow_mall)*(sizeof(double)*0 +  sizeof(int)*5)+
                (ntors_con_mall)*(sizeof(double)*1 +  sizeof(int)*5)+
                (ntors_typ_pow_mall)*(sizeof(double)*30 +  sizeof(int)*0)+
                (nonfo_mall)*(sizeof(double)*0 +  sizeof(int)*3)+
                (nonfo_typ_mall)*(sizeof(double)*3 +  sizeof(int)*0)
                )*1.e-06;

  *tot_memory += now_memory;
  ntyp_pow = nbond_typ_pow_mall;


  bonded->bond.j1_pow = (int *)    crealloc(&(bonded->bond.j1_pow)[1],
                                            nbond_pow_mall*sizeof(int))-1;
  bonded->bond.j2_pow = (int *)    crealloc(&(bonded->bond.j2_pow)[1],
                                            nbond_pow_mall*sizeof(int))-1;
  bonded->bond.jtyp_pow   = (int *)crealloc(&(bonded->bond.jtyp_pow)[1],
                                            nbond_pow_mall*sizeof(int))-1;
  bonded->bond.j1_con = (int *)    crealloc(&(bonded->bond.j1_con)[1],
                                            nbond_con_mall*sizeof(int))-1;
  bonded->bond.j2_con = (int *)    crealloc(&(bonded->bond.j2_con)[1],
                                            nbond_con_mall*sizeof(int))-1;
  bonded->bond.jtyp_con   = (int *)crealloc(&(bonded->bond.jtyp_con)[1],
                                            nbond_con_mall*sizeof(int))-1;
  bonded->bond.al_con = (double *)cmalloc(nbond_con_mall*sizeof(double))-1;
  bonded->bond.eq_pow = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.eq_pow_res = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.c_0    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.c_1    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.c_2    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.c_3    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.c_4    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.c_5    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.c_6    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.dc_0   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.dc_1   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.dc_2   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.dc_3   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.dc_4   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.dc_5   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.dc_6   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.eq_con = (double *)cmalloc(nbond_typ_con_mall
                                          *sizeof(double))-1;

  null_inter_parse->jbond1_nul= (int *)
    crealloc(&(null_inter_parse->jbond1_nul)[1],
             null_inter_parse->nbond_nul*sizeof(int))-1;
  null_inter_parse->jbond2_nul= (int *)
    crealloc(&(null_inter_parse->jbond2_nul)[1],
             null_inter_parse->nbond_nul*sizeof(int))-1;


/*=======================================================================*/
/*  III) Bends  */

  ntyp_pow = nbend_typ_pow_mall;
  bonded->bend.j1_pow = (int *)   crealloc(&(bonded->bend.j1_pow)[1],
                                           nbend_pow_mall*sizeof(int))-1;
  bonded->bend.j2_pow = (int *)   crealloc(&(bonded->bend.j2_pow)[1],
                                           nbend_pow_mall*sizeof(int))-1;
  bonded->bend.j3_pow = (int *)   crealloc(&(bonded->bend.j3_pow)[1],
                                           nbend_pow_mall*sizeof(int))-1;
  bonded->bend.jtyp_pow = (int *)crealloc(&(bonded->bend.jtyp_pow)[1],
                                          nbend_pow_mall*sizeof(int))-1;
  bonded->bend.j1_con = (int *)   crealloc(&(bonded->bend.j1_con)[1],
                                           nbend_con_mall*sizeof(int))-1;
  bonded->bend.j2_con = (int *)   crealloc(&(bonded->bend.j2_con)[1],
                                           nbend_con_mall*sizeof(int))-1;
  bonded->bend.j3_con = (int *)   crealloc(&(bonded->bend.j3_con)[1],
                                           nbend_con_mall*sizeof(int))-1;
  bonded->bend.jtyp_con   = (int *)crealloc(&(bonded->bend.jtyp_con)[1],
                                            nbend_con_mall*sizeof(int))-1;
  bonded->bend.eq_pow = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_0    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_1    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_2    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_3    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_4    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_5    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_6    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_0    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_1    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_2    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_3    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_4    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_5    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_6    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_0   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_1   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_2   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_3   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_4   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_5   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_6   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_0   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_1   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_2   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_3   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_4   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_5   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_6   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.eq_con = (double *)cmalloc(nbend_typ_con_mall
                                          *sizeof(double))-1;
  bonded->bend.al_con = (double *)cmalloc(nbend_con_mall
                                          *sizeof(double))-1;
  null_inter_parse->jbend1_nul= (int *)
    crealloc(&(null_inter_parse->jbend1_nul)[1],
             null_inter_parse->nbend_nul*sizeof(int))-1;
  null_inter_parse->jbend2_nul= (int *)
    crealloc(&(null_inter_parse->jbend2_nul)[1],
             null_inter_parse->nbend_nul*sizeof(int))-1;
  null_inter_parse->jbend3_nul= (int *)
    crealloc(&(null_inter_parse->jbend3_nul)[1],
             null_inter_parse->nbend_nul*sizeof(int))-1;
/*=======================================================================*/
/*  III) Bend_bnds  */

  bonded->bend_bnd.j1 = (int *)   crealloc(&(bonded->bend_bnd.j1)[1],
                                   nbend_bnd_mall*sizeof(int))-1;
  bonded->bend_bnd.j2 = (int *)   crealloc(&(bonded->bend_bnd.j2)[1],
                                   nbend_bnd_mall*sizeof(int))-1;
  bonded->bend_bnd.j3 = (int *)   crealloc(&(bonded->bend_bnd.j3)[1],
                                   nbend_bnd_mall*sizeof(int))-1;
  bonded->bend_bnd.jtyp = (int *)crealloc(&(bonded->bend_bnd.jtyp)[1],
                                   nbend_bnd_mall*sizeof(int))-1;
  bonded->bend_bnd.eq_bond = (double *)cmalloc(nbend_bnd_typ_mall*
                                               sizeof(double))-1;
  bonded->bend_bnd.eq_bend = (double *)cmalloc(nbend_bnd_typ_mall*
                                               sizeof(double))-1;
  bonded->bend_bnd.cbond_0    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbond_1    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbond_2    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbond_3    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbond_4    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbond_5    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbond_6    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbond_0   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbond_1   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbond_2   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbond_3   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbond_4   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbond_5   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbond_6   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;

  bonded->bend_bnd.cbend_0    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbend_1    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbend_2    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbend_3    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbend_4    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbend_5    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbend_6    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.sbend_0    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.sbend_1    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.sbend_2    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.sbend_3    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.sbend_4    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.sbend_5    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.sbend_6    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbend_0   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbend_1   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbend_2   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbend_3   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbend_4   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbend_5   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbend_6   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dsbend_0   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dsbend_1   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dsbend_2   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dsbend_3   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dsbend_4   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dsbend_5   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dsbend_6   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;

/*======================================================================*/
/*  IV) Tors  */

  ntyp_pow = ntors_typ_pow_mall;
  bonded->tors.j1_pow = (int *)   crealloc(&(bonded->tors.j1_pow)[1],
                                           ntors_pow_mall*sizeof(int))-1;
  bonded->tors.j2_pow = (int *)   crealloc(&(bonded->tors.j2_pow)[1],
                                           ntors_pow_mall*sizeof(int))-1;
  bonded->tors.j3_pow = (int *)   crealloc(&(bonded->tors.j3_pow)[1],
                                           ntors_pow_mall*sizeof(int))-1;
  bonded->tors.j4_pow = (int *)   crealloc(&(bonded->tors.j4_pow)[1],
                                           ntors_pow_mall*sizeof(int))-1;
  bonded->tors.jtyp_pow = (int *)crealloc(&(bonded->tors.jtyp_pow)[1],
                                            ntors_pow_mall*sizeof(int))-1;
  bonded->tors.j1_con = (int *)   crealloc(&(bonded->tors.j1_con)[1],
                                           ntors_con_mall*sizeof(int))-1;
  bonded->tors.j2_con = (int *)   crealloc(&(bonded->tors.j2_con)[1],
                                           ntors_con_mall*sizeof(int))-1;
  bonded->tors.j3_con = (int *)   crealloc(&(bonded->tors.j3_con)[1],
                                           ntors_con_mall*sizeof(int))-1;
  bonded->tors.j4_con = (int *)   crealloc(&(bonded->tors.j4_con)[1],
                                           ntors_con_mall*sizeof(int))-1;
  bonded->tors.jtyp_con = (int *)crealloc(&(bonded->tors.jtyp_con)[1],
                                          ntors_con_mall*sizeof(int))-1;
  bonded->tors.al_con = (double *)cmalloc(ntors_con_mall
                                          *sizeof(double))-1;
  bonded->tors.eq_pow = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_0    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_1    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_2    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_3    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_4    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_5    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_6    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_0    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_1    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_2    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_3    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_4    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_5    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_6    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_0   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_1   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_2   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_3   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_4   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_5   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_6   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_0   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_1   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_2   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_3   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_4   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_5   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_6   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.eq_con = (double *)cmalloc(ntors_typ_pow_mall
                                            *sizeof(double))-1;
  null_inter_parse->jtors1_nul= (int *)
    crealloc(&(null_inter_parse->jtors1_nul)[1],
             null_inter_parse->ntors_nul*sizeof(int))-1;
  null_inter_parse->jtors2_nul= (int *)
    crealloc(&(null_inter_parse->jtors2_nul)[1],
             null_inter_parse->ntors_nul*sizeof(int))-1;
  null_inter_parse->jtors3_nul= (int *)
    crealloc(&(null_inter_parse->jtors3_nul)[1],
             null_inter_parse->ntors_nul*sizeof(int))-1;
  null_inter_parse->jtors4_nul= (int *)
    crealloc(&(null_inter_parse->jtors4_nul)[1],
             null_inter_parse->ntors_nul*sizeof(int))-1;

/*=======================================================================*/
/*   V) Onfo  */

  bonded->onfo.j1     = (int *)   crealloc(&(bonded->onfo.j1)[1],
                                           nonfo_mall*sizeof(int))-1;
  bonded->onfo.j2     = (int *)   crealloc(&(bonded->onfo.j2)[1],
                                           nonfo_mall*sizeof(int))-1;
  bonded->onfo.jtyp   = (int *)   crealloc(&(bonded->onfo.jtyp)[1],
                                           nonfo_mall*sizeof(int))-1;
  bonded->onfo.feps   = (double *)cmalloc(nonfo_typ_mall*sizeof(double))-1;
  bonded->onfo.s6 = (double *)cmalloc(nonfo_typ_mall*sizeof(double))-1;
  bonded->onfo.sc = (double *)cmalloc(nonfo_typ_mall*sizeof(double))-1;
  null_inter_parse->jonfo1_nul= (int *)
    crealloc(&(null_inter_parse->jonfo1_nul)[1],
             null_inter_parse->nonfo_nul*sizeof(int))-1;
  null_inter_parse->jonfo2_nul  = (int *)
    crealloc(&(null_inter_parse->jonfo2_nul)[1],
             null_inter_parse->nonfo_nul*sizeof(int))-1;

/*=======================================================================*/
/*   V) Assign mall variables  */

    clatoms_info->natm_mall = natm_mall;
    clatoms_info->nchrg_mall = nchrg_mall;

    ghost_atoms->nghost_mall = nghost_mall;
    ghost_atoms->ncomp_mall = ncomp_mall;

    atommaps->nfreeze_mall = nfreeze_mall;
    atommaps->natm_typ_mall = natm_typ_mall;

    bonded->bond.nbond_pow_mall = nbond_pow_mall;
    bonded->bond.nbond_typ_pow_mall = nbond_typ_pow_mall;
    bonded->bond.nbond_con_mall = nbond_con_mall;
    bonded->bond.nbond_typ_con_mall = nbond_typ_con_mall;

    bonded->grp_bond_con.ngrp_21_mall = ngrp_21_mall;
    bonded->grp_bond_watts.ngrp_33_mall = ngrp_watt_33_mall;
    bonded->grp_bond_con.ngrp_33_mall = ngrp_33_mall;
    bonded->grp_bond_con.ngrp_43_mall = ngrp_43_mall;
    bonded->grp_bond_con.ngrp_23_mall = ngrp_23_mall;
    bonded->grp_bond_con.ngrp_46_mall = ngrp_46_mall;
    bonded->grp_bond_con.ngrp_typ_21_mall = ngrp_typ_21_mall;
    bonded->grp_bond_con.ngrp_typ_33_mall = ngrp_typ_33_mall;
    bonded->grp_bond_watts.ngrp_typ_33_mall = ngrp_typ_watt_33_mall;
    bonded->grp_bond_con.ngrp_typ_43_mall = ngrp_typ_43_mall;
    bonded->grp_bond_con.ngrp_typ_23_mall = ngrp_typ_23_mall;
    bonded->grp_bond_con.ngrp_typ_46_mall = ngrp_typ_46_mall;

    bonded->bend.nbend_con_mall     = nbend_con_mall;
    bonded->bend.nbend_typ_con_mall = nbend_typ_con_mall;
    bonded->bend.nbend_pow_mall     = nbend_pow_mall;
    bonded->bend.nbend_typ_pow_mall = nbend_typ_pow_mall;

    bonded->bend_bnd.nbend_bnd_mall = nbend_bnd_mall;
    bonded->bend_bnd.nbend_bnd_typ_mall = nbend_bnd_typ_mall;

    bonded->tors.ntors_pow_mall = ntors_pow_mall;
    bonded->tors.ntors_typ_pow_mall = ntors_typ_pow_mall;
    bonded->tors.ntors_con_mall = ntors_con_mall;
    bonded->tors.ntors_typ_con_mall = ntors_typ_con_mall;

    bonded->onfo.nonfo_mall = nonfo_mall;
    bonded->onfo.nonfo_typ_mall = nonfo_typ_mall;

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



