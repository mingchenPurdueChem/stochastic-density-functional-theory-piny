/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                Module: control_communicate_groups (DVR)                  */
/*                                                                          */
/* This subprogram reads in all user inputs                                 */
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
#include "../proto_defs/proto_communicate_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_communicate_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 void build_cp_comm_pkg_dvr(CP *cp,MPI_Comm world)
/*==========================================================================*/
{/* begin routine */
/*==========================================================================*/
/*          Local variable declarations                                     */
#include "../typ_defs/typ_mask.h"
  int irem,idiv,iii;
  int nstate_ncoef_proc_max,nstate_ncoef_proc_min;
  int num_coef,*num_coef_v,ncoef_proc,ncoef_proc_yz;
/*==========================================================================*/
/* I) Up states                                                             */

 /*------------------------------------*/
 /* i) states per processor stuff      */

  idiv =  cp->cpcoeffs_info.nstate_up/cp->communicate.np_states;
  irem = (cp->cpcoeffs_info.nstate_up % cp->communicate.np_states);
  cp->cpcoeffs_info.nstate_up_proc = idiv;
  if(cp->communicate.myid_state < irem) {
     cp->cpcoeffs_info.nstate_up_proc = idiv+1;
  }/*endif*/
  if(cp->communicate.myid_state <= irem) {
    cp->cpcoeffs_info.istate_up_st = cp->communicate.myid_state*(idiv+1)+1;
  } else {
    cp->cpcoeffs_info.istate_up_st = irem*(idiv+1)
                                   + (cp->communicate.myid_state-irem)*idiv+1;
  }/*endif*/
    cp->cpcoeffs_info.istate_up_end = cp->cpcoeffs_info.istate_up_st +
                                    cp->cpcoeffs_info.nstate_up_proc-1;

 /*------------------------------------*/
 /* ii) coefs per processor stuff      */

  cp->cp_comm_state_pkg_dvr_up.num_proc   = cp->communicate.np_states;
  cp->cp_comm_state_pkg_dvr_up.myid       = cp->communicate.myid_state;
  cp->cp_comm_state_pkg_dvr_up.nstate     = cp->cpcoeffs_info.nstate_up;
  cp->cp_comm_state_pkg_dvr_up.ncoef      = cp->cpcoeffs_info.ncoef;
  cp->cp_comm_state_pkg_dvr_up.nstate_proc= cp->cpcoeffs_info.nstate_up_proc;
  cp->cp_comm_state_pkg_dvr_up.world      = world;
  if(cp->communicate.np_states > 1){
    Comm_dup(cp->communicate.comm_states,&(cp->cp_comm_state_pkg_dvr_up.comm));
  } else {
    cp->cp_comm_state_pkg_dvr_up.comm = cp->communicate.comm_states;
  }/* endif */


  irem             = (cp->cp_comm_state_pkg_dvr_up.nstate %
                      cp->cp_comm_state_pkg_dvr_up.num_proc);
  cp->cp_comm_state_pkg_dvr_up.nstate_proc_max  = (irem > 0 ? idiv+1 : idiv);

  cp->cp_comm_state_pkg_dvr_up.nstate_max = (irem > 0 ?
                          ((idiv+1)*cp->communicate.np_states) :
                          (idiv*cp->communicate.np_states)) ;


  /* different from PW code*/

  cp->cp_comm_state_pkg_dvr_up.nstate_proc_min  = idiv; 


  idiv = (cp->cpcoeffs_info.grid_ny)*(cp->cpcoeffs_info.grid_nz)/
         (cp->cp_comm_state_pkg_dvr_up.num_proc);

  irem = (cp->cpcoeffs_info.grid_ny * cp->cpcoeffs_info.grid_nz) %
         cp->cp_comm_state_pkg_dvr_up.num_proc;

  ncoef_proc_yz =  (cp->communicate.myid_state < irem ? idiv+1 : idiv);
  ncoef_proc = ncoef_proc_yz * (cp->cpcoeffs_info.grid_nx);

  cp->cpcoeffs_info.nstate_ncoef_proc_up =  ncoef_proc;

  cp->cp_comm_state_pkg_dvr_up.nstate_ncoef_proc  =
                         cp->cpcoeffs_info.nstate_ncoef_proc_up;

  if(cp->communicate.np_states > 1){
    Allreduce(&(cp->cpcoeffs_info.nstate_ncoef_proc_up),
              &nstate_ncoef_proc_max,
              1,MPI_INT,MPI_MAX,0,world);
    /* Not defined anymore 
    Allreduce(&(cp->cpcoeffs_info.nstate_ncoef_proc_up),
              &nstate_ncoef_proc_min,
              1,MPI_INT,MPI_MIN,0,world); */
  }else{
    nstate_ncoef_proc_max = cp->cpcoeffs_info.nstate_ncoef_proc_up;
    /* nstate_ncoef_proc_min = cp->cpcoeffs_info.nstate_ncoef_proc_up; */
  }

  cp->cpcoeffs_info.nstate_ncoef_proc_max_up          = nstate_ncoef_proc_max;
  cp->cp_comm_state_pkg_dvr_up.nstate_ncoef_proc_max  = nstate_ncoef_proc_max;
  /*cp->cp_comm_state_pkg_dvr_up.nstate_ncoef_proc_min  = nstate_ncoef_proc_min; */


  if(cp->communicate.np_states > 1){
    num_coef   =  cp->cp_comm_state_pkg_dvr_up.nstate_ncoef_proc;
    num_coef_v = (int *) malloc((cp->communicate.np_states)*sizeof(int))-1;
    Allgather(&num_coef,1,MPI_INT,&num_coef_v[1],1,MPI_INT,0,world);

    cp->cpcoeffs_info.icoef_start_up = 1;
    for(iii=1; iii <= cp->communicate.myid; iii++){
     cp->cpcoeffs_info.icoef_start_up += num_coef_v[iii];
    }
    cp->cp_comm_state_pkg_dvr_up.icoef_start =
                                   cp->cpcoeffs_info.icoef_start_up;
  }else{
    cp->cpcoeffs_info.icoef_start_up     = 1;
    cp->cp_comm_state_pkg_dvr_up.icoef_start = 1;
  }

/*==========================================================================*/
/* II) Down states                                                          */

 /*------------------------------------*/
 /* i) states per processor stuff      */
  idiv = cp->cpcoeffs_info.nstate_dn/cp->communicate.np_states;
  irem = (cp->cpcoeffs_info.nstate_dn % cp->communicate.np_states);
  cp->cpcoeffs_info.nstate_dn_proc = idiv;
  if(cp->communicate.myid_state < irem) {
     cp->cpcoeffs_info.nstate_dn_proc = idiv+1;
  }/*endif*/
  if(cp->communicate.myid_state <= irem) {
    cp->cpcoeffs_info.istate_dn_st = cp->communicate.myid_state*(idiv+1)+1;
  } else {
    cp->cpcoeffs_info.istate_dn_st = irem*(idiv+1)
                                   + (cp->communicate.myid_state-irem)*idiv+1;
  }/*endif*/
  cp->cpcoeffs_info.istate_dn_end = cp->cpcoeffs_info.istate_dn_st +
                                    cp->cpcoeffs_info.nstate_dn_proc-1;


 /*------------------------------------*/
 /* ii) coefs per processor stuff      */

  cp->cp_comm_state_pkg_dvr_dn.num_proc   = cp->communicate.np_states;
  cp->cp_comm_state_pkg_dvr_dn.myid       = cp->communicate.myid_state;
  cp->cp_comm_state_pkg_dvr_dn.nstate     = cp->cpcoeffs_info.nstate_dn;
  cp->cp_comm_state_pkg_dvr_dn.ncoef      = cp->cpcoeffs_info.ncoef;
  cp->cp_comm_state_pkg_dvr_dn.nstate_proc= cp->cpcoeffs_info.nstate_dn_proc;
  cp->cp_comm_state_pkg_dvr_dn.world      = world;
  if(cp->communicate.np_states > 1){
    Comm_dup(cp->communicate.comm_states,&(cp->cp_comm_state_pkg_dvr_dn.comm));
  } else {
    cp->cp_comm_state_pkg_dvr_dn.comm = cp->communicate.comm_states;
  }/* endif */

  irem             = (cp->cp_comm_state_pkg_dvr_dn.nstate %
                      cp->cp_comm_state_pkg_dvr_dn.num_proc);
  cp->cp_comm_state_pkg_dvr_dn.nstate_proc_max  = (irem > 0 ? idiv+1 : idiv);
  cp->cp_comm_state_pkg_dvr_dn.nstate_max = (irem > 0 ?
                          ((idiv+1)*cp->communicate.np_states) :
                          (idiv*cp->communicate.np_states)) ;

  cp->cp_comm_state_pkg_dvr_dn.nstate_proc_min  =  idiv; 

  cp->cpcoeffs_info.nstate_ncoef_proc_dn = ncoef_proc;

  cp->cp_comm_state_pkg_dvr_dn.nstate_ncoef_proc  =
              cp->cpcoeffs_info.nstate_ncoef_proc_dn;

  if(cp->communicate.np_states > 1){
    Allreduce(&(cp->cpcoeffs_info.nstate_ncoef_proc_dn),
              &nstate_ncoef_proc_max,
              1,MPI_INT,MPI_MAX,0,world);
    /* Allreduce(&(cp->cpcoeffs_info.nstate_ncoef_proc_dn),
              &nstate_ncoef_proc_min,
              1,MPI_INT,MPI_MIN,0,world); */
  }else{
    nstate_ncoef_proc_max = cp->cpcoeffs_info.nstate_ncoef_proc_dn;
    /* nstate_ncoef_proc_min = cp->cpcoeffs_info.nstate_ncoef_proc_dn; */
  }

  cp->cpcoeffs_info.nstate_ncoef_proc_max_dn = nstate_ncoef_proc_max;
  cp->cp_comm_state_pkg_dvr_dn.nstate_ncoef_proc_max  =
              cp->cpcoeffs_info.nstate_ncoef_proc_max_dn;
  /* cp->cp_comm_state_pkg_dvr_dn.nstate_ncoef_proc_min  = nstate_ncoef_proc_min; */

  if(cp->communicate.np_states > 1){

    for(iii=1; iii <= cp->communicate.np_states; iii++){
      num_coef_v[iii] = 0;
    }
  
    num_coef   =  cp->cp_comm_state_pkg_dvr_dn.nstate_ncoef_proc;
    Allgather(&num_coef,1,MPI_INT,&num_coef_v[1],1,MPI_INT,0,world);

    cp->cpcoeffs_info.icoef_start_dn = 1;
    for(iii=1; iii <= cp->communicate.myid; iii++){
     cp->cpcoeffs_info.icoef_start_dn += num_coef_v[iii];
    }
    cp->cp_comm_state_pkg_dvr_dn.icoef_start =
                                   cp->cpcoeffs_info.icoef_start_dn;
  }else{
    cp->cpcoeffs_info.icoef_start_dn     = 1;
    cp->cp_comm_state_pkg_dvr_dn.icoef_start = 1;
  }

  if(cp->communicate.np_states > 1){
    free(&num_coef_v[1]);
  }

/*==========================================================================*/
   }/* end routine */
/*==========================================================================*/

