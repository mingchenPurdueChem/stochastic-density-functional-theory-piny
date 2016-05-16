/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: mall_scratch.c                               */
/*                                                                          */
/* This subprogram mallocs the scratch memory                               */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_scratch_entry.h"
#include "../proto_defs/proto_scratch_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void mall_cp_scr_dvr(CPTHERM_INFO *cptherm_info,CPOPTS *cpopts,CPEWALD *cpewald,
                 CPSCR *cpscr,CPCOEFFS_INFO *cpcoeffs_info,
                 DVR_MATRIX *dvr_matrix, PSEUDO *pseudo,
                 PARA_FFT_PKG3D *cp_para_fft_pkg3d_sm,
                 CP_COMM_STATE_PKG *cp_comm_state_pkg_up_dvr,
                 CP_COMM_STATE_PKG *cp_comm_state_pkg_dn_dvr,int atm_hess_calc,
                 double *tot_memory,
                 int myid, MPI_Comm world)

/*==========================================================================*/
  {/*begin routine*/
/*==========================================================================*/
/*         Local pointer declarations                                      */

  int nstate_up         = cpcoeffs_info->nstate_up_proc;
  int nstate_dn         = cpcoeffs_info->nstate_dn_proc;
  int nstate_up_tot     = cpcoeffs_info->nstate_up;
  int nstate_dn_tot     = cpcoeffs_info->nstate_dn;

  int ncoef             = cpcoeffs_info->ncoef;   /* total number of rectangular grid */
  int ncoef_l           = cpcoeffs_info->ncoef_l; /* ncoef_l = ncoef in DVR */

  int n_ang_max         = pseudo->n_ang_max;
  int n_rad_max         = pseudo->n_rad_max;
  int natm_nls_max      = cpscr->cpscr_nonloc.natm_nls_max;
  int np_states         = cp_comm_state_pkg_up_dvr->num_proc;
  int num_c_nhc_proc    = cptherm_info->num_c_nhc_proc;
  int nfft_up_proc      = cp_para_fft_pkg3d_sm->nfft_proc;
  int nfft_up           = cp_para_fft_pkg3d_sm->nfft;
  int num_c_nhc1        = num_c_nhc_proc+1;
  int grid_nx           = cpcoeffs_info->grid_nx;
  int grid_ny           = cpcoeffs_info->grid_ny;
  int grid_nz           = cpcoeffs_info->grid_nz;
  int ncoef_l_proc_max_mall;
  int ncoef_l_proc_max_mall_ke;
  int ncoef_l_proc_max_dn;
  int ndim_wannier;
  int mm=5;

  int nstate_ncoef_proc_max_up = cpcoeffs_info->nstate_ncoef_proc_max_up;
  int nstate_proc_max = cp_comm_state_pkg_up_dvr->nstate_proc_max;


  /* Options Turned On/Off */
  int massiv_flag       = cptherm_info->massiv_flag;
  int cp_lsda           = cpopts->cp_lsda;
  int cp_gga            = cpopts->cp_gga;
  int cp_ke_dens_on     = cpcoeffs_info->cp_ke_dens_on;
  int cp_elf_calc_frq   = cpcoeffs_info->cp_elf_calc_frq;
  int cp_para_opt       = cpopts->cp_para_opt;
  int cp_wan_opt        = cpopts->cp_wan_opt;
  int cp_wan_min_opt    = cpopts->cp_wan_min_opt;
  int cp_wan_init_opt   = cpopts->cp_wan_init_opt;

/*--------------------------------------------------------------------------*/
/*         Local variable declarations                                      */

  double now_memory;
  int i,iii,irem;
  int nlscr_up,nlscr_dn,nlscr_mall;
  int ncoef_l_proc_max;
  int nfft2_mall_up,nfft2_mall_dn,nfft2_mall_up_proc,nfft2_mall_dn_proc;
  int nfft2_up,nfft2_dn;
  int nfft2_up_proc,nfft2_dn_proc,nfft_dn;
  int nfft2_up_gga,nfft2_dn_gga,nlap_up,nlap_dn;
  int nfft2_up_ke_dens,nfft2_dn_ke_dens;

  int ncoef_up,ncoef_dn;
  int par_size_up,par_size_dn;
  int ncoef2_up,ncoef2_dn;
  int ncoef2_max_up, ncoef2_max_dn;
  int nstate,nstate2,nstate_tot,nstate2_tot;
  int num=0;

  int mtemp;

/*=========================================================================*/
/* I) Malloc size calculation */

  /* (a) DVR CP : Define the grid size              */

  nfft2_up      = nfft_up/2;
  nfft2_up_proc = nfft_up_proc/2;

  nfft_dn       = (cp_lsda == 1 ? nfft_up : 0);
  nfft2_dn      = (cp_lsda == 1 ? nfft2_up : 0);
  nfft2_dn_proc = (cp_lsda == 1 ? nfft2_up_proc : 0);

  nfft2_mall_up      =  nfft2_up;
  nfft2_mall_up_proc =  nfft2_up_proc;

  nfft2_mall_dn      =  nfft2_dn;
  nfft2_mall_dn_proc =  nfft2_dn_proc;

  ncoef_l_proc_max = ncoef_l/np_states;

  ncoef_up  = nfft2_up;
  ncoef_dn  = (cp_lsda == 1 ? nfft2_up : 0);

  irem = (ncoef_l % np_states);
  if(irem>0){ncoef_l_proc_max++;}

  ncoef_l_proc_max_mall    = ncoef_l_proc_max;
  ncoef_l_proc_max_mall_ke = (cp_ke_dens_on == 1 ? ncoef_l_proc_max_mall:0);
  ncoef_l_proc_max_dn      = (cp_lsda == 1 ? ncoef_l_proc_max_mall : 0);

 /* (b) Wavefunction scratch sizes     */

  ncoef2_up      = MAX(nfft_up, (nfft2_up_proc*nstate_up_tot));
  ncoef2_up      = MAX(ncoef2_up, nfft2_up*nstate_up); 
  ncoef2_max_up  = MAX(ncoef2_up, nstate_ncoef_proc_max_up*nstate_up_tot);

  ncoef2_dn      = MAX(nfft_dn, (nfft2_dn_proc*nstate_dn_tot));
  ncoef2_dn      = MAX(ncoef2_dn, nfft2_dn*nstate_dn);

  nstate        = MAX(nstate_up, nstate_dn);
  nstate2       = nstate*nstate;
  nstate_tot    = MAX(nstate_up_tot,nstate_dn_tot);
  nstate2_tot   = nstate_tot*nstate_tot;

 /* (c) Nonlocal sizes  */

  nlscr_up  = nstate_up_tot*(n_ang_max+1)*(n_ang_max+1)*natm_nls_max;
  nlscr_dn  = 0;

  if(cp_lsda==1){
    nlscr_dn  = nstate_dn_tot*(n_ang_max+1)*(n_ang_max+1)*natm_nls_max;
  }/*endif*/

  nlscr_mall = (nlscr_up > nlscr_dn ? nlscr_up : nlscr_dn);

  /* (d) Wannier scratch size   */

  if(cp_wan_opt ==1 || cp_wan_min_opt==1 || cp_wan_init_opt==1){
    cpscr->cpscr_wannier.cp_wannier_on=1;
    ndim_wannier=nstate_up_tot*nstate_up_tot;
  }else{
    cpscr->cpscr_wannier.cp_wannier_on=0;
    ndim_wannier=0;
  }

/*------------------------------------------------------------------*/
/* II) Local External Potential Arrays                                  */

  cpscr->cpscr_loc.vextr = (double *)
                        cmalloc(nfft2_up_proc*sizeof(double))-1;
  cpscr->cpscr_loc.vexti = (double *)
                        cmalloc(nfft2_up_proc*sizeof(double))-1;

  num += 2*nfft2_up_proc;

/*------------------------------------------------------------------*/
/* III) Non-Local External Potential Arrays                        */

  cpscr->cpscr_nonloc.vnlre_up
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
  cpscr->cpscr_nonloc.vnlim_up
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;

  cpscr->cpscr_nonloc.dvnlre_up
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_up
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;


  cpscr->cpscr_nonloc.dvnlre_x_up
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_y_up
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_z_up
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;

  cpscr->cpscr_nonloc.dvnlim_x_up
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_y_up
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_z_up
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;


  cpscr->cpscr_nonloc.vnlre_tmp
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
  cpscr->cpscr_nonloc.vnlim_tmp
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;

  cpscr->cpscr_nonloc.dvnlre_x_tmp
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_y_tmp
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlre_z_tmp
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;

  cpscr->cpscr_nonloc.dvnlim_x_tmp
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_y_tmp
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
  cpscr->cpscr_nonloc.dvnlim_z_tmp
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;

  if(cp_lsda==1){
    cpscr->cpscr_nonloc.dvnlre_x_dn
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
    cpscr->cpscr_nonloc.dvnlre_y_dn
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
    cpscr->cpscr_nonloc.dvnlre_z_dn
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;

    cpscr->cpscr_nonloc.dvnlim_x_dn
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
    cpscr->cpscr_nonloc.dvnlim_y_dn
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
    cpscr->cpscr_nonloc.dvnlim_z_dn
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;

    cpscr->cpscr_nonloc.vnlre_dn
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
    cpscr->cpscr_nonloc.vnlim_dn
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;

    cpscr->cpscr_nonloc.dvnlre_dn
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
    cpscr->cpscr_nonloc.dvnlim_dn
                   = (double *)cmalloc(nlscr_mall*sizeof(double))-1;
  }
  num += 20*nlscr_mall;

 /*-------------------------------------------------------------------*/
 /* IV) GGA sizes        */

  /* change due to the nstate_proc form of gga routine */
  /* nfft2_up_gga  = ((cp_gga == 1 || cp_elf_calc_frq > 0) ? nfft2_mall_up_proc:0);
     nfft2_dn_gga  = (((cp_gga == 1 || cp_elf_calc_frq > 0) && cp_lsda == 1)
                ? nfft2_mall_dn_proc:0); */

  nfft2_up_gga  = ((cp_gga == 1 || cp_elf_calc_frq > 0) ? nfft2_mall_up:0);
  nfft2_dn_gga  = (((cp_gga == 1 || cp_elf_calc_frq > 0) && cp_lsda == 1)
                ? nfft2_mall_dn:0);

  nfft2_up_ke_dens = (cp_ke_dens_on == 1 ? nfft2_mall_up_proc:0);
  nfft2_dn_ke_dens = ((cp_ke_dens_on == 1 && cp_lsda == 1) ? nfft2_mall_dn_proc:0);

/*------------------------------------------------------------------*/
/* V) DVR DENSITY REAL and RECIPROCAL SPACE                            */

  cpscr->cpscr_rho.rho_up = (double *)
                            cmalloc(nfft2_mall_up_proc*sizeof(double))-1;
  cpscr->cpscr_rho.rho_dn = (double *)
                            cmalloc(nfft2_mall_dn_proc*sizeof(double))-1;

  num += (nfft2_mall_up_proc + nfft2_mall_dn_proc);

  cpscr->cpscr_rho.rhocr_up = (double *)
                        cmalloc(nfft2_up_proc*sizeof(double))-1;

  cpscr->cpscr_rho.rhoci_up = (double *)
                        cmalloc(nfft2_up_proc*sizeof(double))-1;

  cpscr->cpscr_rho.rhocr_scr = (double *)
                        cmalloc(ncoef_l_proc_max_mall_ke*sizeof(double))-1;
  cpscr->cpscr_rho.rhoci_scr = (double *)
                        cmalloc(ncoef_l_proc_max_mall_ke*sizeof(double))-1;

  cpscr->cpscr_rho.rhocr_dn = (double *)
                              cmalloc(ncoef_l_proc_max_dn*sizeof(double))-1;
  cpscr->cpscr_rho.rhoci_dn = (double *)
                              cmalloc(ncoef_l_proc_max_dn*sizeof(double))-1;

  num += 2*nfft2_up_proc;
  num += 2*ncoef_l_proc_max_mall_ke;
  num += 2*ncoef_l_proc_max_dn;

/*------------------------------------------------------------------*/
/* VI) Kinetic Energy Arrays   */

  dvr_matrix->Tsx = (double *) cmalloc((grid_nx*grid_nx)*sizeof(double))-1;
  dvr_matrix->Tsy = (double *) cmalloc((grid_ny*grid_ny)*sizeof(double))-1;
  dvr_matrix->Tsz = (double *) cmalloc((grid_nz*grid_nz)*sizeof(double))-1;

  dvr_matrix->dGx = (double *) cmalloc((grid_nx*grid_nx)*sizeof(double))-1;
  dvr_matrix->dGy = (double *) cmalloc((grid_ny*grid_ny)*sizeof(double))-1;
  dvr_matrix->dGz = (double *) cmalloc((grid_nz*grid_nz)*sizeof(double))-1;

  num += (2*(grid_nx*grid_nx + grid_ny*grid_ny + grid_nz*grid_nz));

/*-------------------------------------------------------------------*/
/* VII) DVR-FBR transformation arrays */

  dvr_matrix->TRx = (double *) cmalloc((grid_nx*grid_nx)*sizeof(double))-1;
  dvr_matrix->TIx = (double *) cmalloc((grid_nx*grid_nx)*sizeof(double))-1;
  dvr_matrix->TRy = (double *) cmalloc((grid_ny*grid_ny)*sizeof(double))-1;
  dvr_matrix->TIy = (double *) cmalloc((grid_ny*grid_ny)*sizeof(double))-1;
  dvr_matrix->TRz = (double *) cmalloc((grid_nz*grid_nz)*sizeof(double))-1;
  dvr_matrix->TIz = (double *) cmalloc((grid_nz*grid_nz)*sizeof(double))-1;

  dvr_matrix->Tfbr_x = (double *) cmalloc((grid_nx)*sizeof(double))-1;
  dvr_matrix->Tfbr_y = (double *) cmalloc((grid_ny)*sizeof(double))-1;
  dvr_matrix->Tfbr_z = (double *) cmalloc((grid_nz)*sizeof(double))-1;

  num += (2*(grid_nx*grid_nx + grid_ny*grid_ny + grid_nz*grid_nz)
         + (grid_nx+grid_ny+grid_nz));

/*------------------------------------------------------------------*/
/* VIII) gga */

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

 /* change due to the nstate_proc form of gga routine */
 /* cpscr->cpscr_grho.d_rhox_up_c = (double *)
                           cmalloc((nfft2_up_gga*nstate_up_tot)*sizeof(double)) -1;
    cpscr->cpscr_grho.d_rhoy_up_c = (double *)
                           cmalloc((nfft2_up_gga*nstate_up_tot)*sizeof(double)) -1;
    cpscr->cpscr_grho.d_rhoz_up_c = (double *)
                           cmalloc((nfft2_up_gga*nstate_up_tot)*sizeof(double)) -1;

    num += (3*nfft2_up_gga*nstate_up_tot); */

  cpscr->cpscr_grho.d_rhox_up_c = (double *)
                           cmalloc((nfft2_up_gga*nstate_up)*sizeof(double)) -1; 
  cpscr->cpscr_grho.d_rhoy_up_c = (double *)
                           cmalloc((nfft2_up_gga*nstate_up)*sizeof(double)) -1; 
  cpscr->cpscr_grho.d_rhoz_up_c = (double *)
                           cmalloc((nfft2_up_gga*nstate_up)*sizeof(double)) -1; 
  num += (3*nfft2_up_gga*nstate_up);

/*------------------------------------------------------------------*/
/* IX) cpewald */

  /* This does not work for DVR !!! */
  /* cpewald->ak2    = (double *)cmalloc((ncoef_l_proc_max_mall)*sizeof(double))-1;*/

  cpewald->ak2    = (double *)cmalloc((nfft2_mall_up_proc)*sizeof(double))-1;

  cpewald->ak2_sm = (double *)cmalloc(ncoef*sizeof(double))-1; /* NEEDED FOR DVR ? */

  num += (nfft2_mall_up_proc + ncoef);

/*------------------------------------------------------------------*/
/* X) wave */

  cpscr->cpscr_wave.cre_up
                    = (double *)cmalloc(ncoef2_max_up*sizeof(double))-1;
  cpscr->cpscr_wave.cre_dn
                    = (double *)cmalloc(ncoef2_dn*sizeof(double))-1;
  cpscr->cpscr_wave.cim_up
                    = (double *)cmalloc(ncoef2_max_up*sizeof(double))-1;
  cpscr->cpscr_wave.cim_dn
                    = (double *)cmalloc(ncoef2_dn*sizeof(double))-1;

  num += 2*(ncoef2_up + ncoef2_dn);

  if(cp_para_opt==0){ /* hybrid*/

     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("No hybrid option for CP-DVR calculation (mall_scratch)\n");
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);

  }else{  /* full g */

    cpscr->cpscr_wave.zfft
                    = (double *)cmalloc(ncoef2_up*sizeof(double))-1;
    cpscr->cpscr_wave.zfft_tmp
                    = (double *)cmalloc(ncoef2_up*sizeof(double))-1;

    num += 2*(ncoef2_up);
  }/*endif*/

  cpscr->cpscr_wave.zfft_mall_size = ncoef2_up;

/*------------------------------------------------------------------*/
/* XI) ovmat */
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
/* XII) thermostat */
  if(massiv_flag==0){
    cpscr->cpscr_therm.sc_cp
                     = (double *)cmalloc(num_c_nhc1*sizeof(double))-1;
    cpscr->cpscr_therm.coef_kin
                     = (double *)cmalloc(num_c_nhc1*sizeof(double))-1;
    num += 2*num_c_nhc1;
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* XIII) Wannier scratch                                                       */

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


/*===========================================================================*/
/* XIV) Output time has arrived */

  now_memory = num*sizeof(double)*(1.0e-6);
  *tot_memory += now_memory;
  if(myid==0){
    printf("CP allocation: %g Mbytes; Total memory: %g Mbytes\n",
            now_memory,*tot_memory);
  }/*endif for myid=0*/

/*--------------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/



