/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: analysis_cp                                  */
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
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_wannier_cpcon_local.h"
#include "../proto_defs/proto_wannier_cpcon_entry.h"


#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0
#define its_max 200


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calcul_wannier_dvr(GENERAL_DATA *general_data,CP *cp, double ***Z_real,
                        double ***Z_imag,double *W, double **wan_cent,  
                        int dip_calc_flag)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/
/*            Local variable declarations                                    */
/*===========================================================================*/
#include "../typ_defs/typ_mask.h"

  int i,j,k;

  int myid         = cp->communicate.myid_state;
  int num_proc     = cp->communicate.np_states;
  MPI_Comm comm    = cp->communicate.comm_states;

  int itime = general_data->timeinfo.itime;

  int icoef_form_up    = cp->cpcoeffs_pos_dvr[1].icoef_form_up;
  int icoef_form_dn    = cp->cpcoeffs_pos_dvr[1].icoef_form_dn;
  int cp_wan_min_on    = cp->cpopts.cp_wan_min_opt;
  int cp_wan_on        = cp->cpopts.cp_wan_opt;
  int *ioff_st         = cp->cpcoeffs_info.ioff_upt;
  int cp_lsda          = cp->cpopts.cp_lsda;
  int ncoef_up         = cp->cp_comm_state_pkg_dvr_up.nstate_ncoef_proc;
  int ncoef_dn         = cp->cp_comm_state_pkg_dvr_dn.nstate_ncoef_proc;
  int nstate_up        = cp->cpcoeffs_info.nstate_up;
  int nstate_dn        = cp->cpcoeffs_info.nstate_dn;
  int nstate=nstate_up;
  int tot_beads        = cp->cpcoeffs_info.pi_beads;
  int ncoef_tot;

  double *dvrc             = cp->cpcoeffs_pos_dvr[1].dvrc_up;
  double *dvrfc            = cp->cpcoeffs_pos_dvr[1].dvrfc_up;
  double *dvrvc            = cp->cpcoeffs_pos_dvr[1].dvrvc_up;
  double *dvrc_tmp         = cp->cpscr.cpscr_wave.cre_up;
  double alpha=1.0, beta=0.0;

  int NDIM=nstate*nstate;
  double *elementA = cp->electronic_properties.elementA;
  double *U_final = cp->cpscr.cpscr_wannier.U_final;

/*=======================================================================*/
/* 0) Parallel checks                                         */

  if(num_proc>1){
    if(icoef_form_up!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up Coef forces are not in transposed form \n");
      printf("on state processor %d in calcul_wannier \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
      if(icoef_form_dn!=1){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Dn Coef forces are not in transposed form \n");
        printf("on state processor %d in calcul_wannier \n",myid);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
    }/*endif*/
  }/*endif*/

  if(tot_beads != 1){
    if(myid==0){
      printf("@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("  No Path-Integral MD yet for Wannier/Diople calculation\n");
      printf("@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    }
    fflush(stdout);
    Finalize();
    exit(1);
  }

/*=============================================================================*/
/* I) Initialization                                                              */

  for (j=1;j<=nstate;j++){
    for (k=1;k<=3;k++){
      wan_cent[j][k]=0.0;
    }
  }
  for (i=1;i<=NDIM;i++){
    elementA[i]=0.0;
  }

/*==============================================================================*/
/* II) Compute z_I,ij tensor elements                                           */

  if(dip_calc_flag==0 ){
    comp_tensor_z_dvr(general_data,cp,Z_real,Z_imag);
  }

/*==============================================================================*/
/* III) Wannier minimization for proc=0                                            */

  if(num_proc > 1){
    Barrier(comm);
  }

  if(myid ==0){

    if(cp_wan_min_on==1 ){
      /*Full minimization: LDFP Minimization*/
      lbfg_min(general_data,cp,elementA,NDIM,Z_real,Z_imag,W);
    }else if(cp_wan_on==1){
      /*on-the_fly wannier: NEWTON Minimization*/
      newton(general_data,cp,elementA,NDIM,Z_real,Z_imag,W);
    }

    comp_wannier_center(general_data,cp,elementA,Z_real,Z_imag,
                        W,wan_cent,U_final);
  }

/*===========================================================================*/
/* IV) Broadcast the optimized A and Wannier centers                         */

   if(num_proc > 1){

     Barrier(comm);

     Bcast(&elementA[1],NDIM,MPI_DOUBLE,0,comm);
     Bcast(&U_final[1],NDIM,MPI_DOUBLE,0,comm);
     Bcast(&wan_cent[1][1],3*nstate,MPI_DOUBLE,0,comm);
   }

/*============================================================================*/
/* V) Rotate orbital, velocities and force                                    */

    ncoef_tot=nstate_up*ncoef_up;

    cp_rotation_prim_dvr(dvrc, icoef_form_up, dvrc_tmp, icoef_form_up, U_final,
                         alpha,beta,ioff_st,&(cp->cp_comm_state_pkg_dvr_up));

    for(i=1;i<=ncoef_tot;i++){
      dvrc[i]=dvrc_tmp[i];
    }

    cp_rotation_prim_dvr(dvrvc, icoef_form_up, dvrc_tmp, icoef_form_up, U_final,
                         alpha,beta,ioff_st,&(cp->cp_comm_state_pkg_dvr_up));

    for(i=1;i<=ncoef_tot;i++){
      dvrvc[i]=dvrc_tmp[i];
    }

    cp_rotation_prim_dvr(dvrfc, icoef_form_up, dvrc_tmp, icoef_form_up, U_final,
                         alpha,beta,ioff_st,&(cp->cp_comm_state_pkg_dvr_up));

    for(i=1;i<=ncoef_tot;i++){
      dvrfc[i]=dvrc_tmp[i];
    }

/*===========================================================================*/
} /*end routine */
/*============================================================================*/

void calcul_initial_wannier_dvr(GENERAL_DATA *general_data,CP *cp)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/
/*            Local variable declarations                                    */
/*===========================================================================*/
#include "../typ_defs/typ_mask.h"

  int i,j,k;

  int myid         = cp->communicate.myid_state;
  int num_proc     = cp->communicate.np_states;
  MPI_Comm comm    = cp->communicate.comm_states;

  int itime = general_data->timeinfo.itime;

  int icoef_form_up    = cp->cpcoeffs_pos_dvr[1].icoef_form_up;
  int icoef_form_dn    = cp->cpcoeffs_pos_dvr[1].icoef_form_dn;
  int cp_wan_min_on = cp->cpopts.cp_wan_min_opt;
  int *ioff_st      = cp->cpcoeffs_info.ioff_upt;
  int cp_lsda        = cp->cpopts.cp_lsda;
  int ncoef_up       = cp->cp_comm_state_pkg_dvr_up.nstate_ncoef_proc;
  int ncoef_dn       = cp->cp_comm_state_pkg_dvr_dn.nstate_ncoef_proc;
  int nstate_up      = cp->cpcoeffs_info.nstate_up;
  int nstate_dn      = cp->cpcoeffs_info.nstate_dn;
  int nstate=nstate_up;
  int ncoef_tot;

  double *dvrc             = cp->cpcoeffs_pos_dvr[1].dvrc_up;
  double *dvrfc            = cp->cpcoeffs_pos_dvr[1].dvrfc_up;
  double *dvrvc            = cp->cpcoeffs_pos_dvr[1].dvrvc_up;
  double *dvrc_tmp         = cp->cpscr.cpscr_wave.cre_up;
  double alpha=1.0, beta=0.0;

  int NDIM=nstate*nstate;
  double *elementA = cp->electronic_properties.elementA;
  double *U_final = cp->cpscr.cpscr_wannier.U_final;

  int I=6;
  double *W         = cp->electronic_properties.weight;
  double ***Z_real  = cp->electronic_properties.Z_real;
  double ***Z_imag  = cp->electronic_properties.Z_imag;
  double **wan_cent = cp->electronic_properties.wannier_cent;;

  int iwrite_init_wcent = cp->cpopts.iwrite_init_wcent;
  int iwrite_init_worb  = cp->cpopts.iwrite_init_worb;

/*=======================================================================*/
/* 0) Parallel checks   */

  if(num_proc>1){
    if(icoef_form_up!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up Coef forces are not in transposed form \n");
      printf("on state processor %d in calcul_wannier_min \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
      if(icoef_form_dn!=1){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Dn Coef forces are not in transposed form \n");
        printf("on state processor %d in calcul_wannier_min \n",myid);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
    }/*endif*/
  }/*endif*/

  if(cp_lsda==1){
    if(myid==0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("LSDA has not been implemented for CP Wannier.\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    }
    fflush(stdout);
    exit(1);
  }
  
/*=============================================================================*/
/* I) Initialization                                                              */

  for (j=1;j<=nstate;j++){
    for (k=1;k<=3;k++){
      wan_cent[j][k]=0.0;
    }
  }
  for (i=1;i<=NDIM;i++){
    elementA[i]=0.0;
  }

/*==============================================================================*/
/* II) Compute Miller weight and z_I,ij tensor elements                         */

   comp_miller_weight(general_data,W);

   comp_tensor_z_dvr(general_data,cp,Z_real,Z_imag);

/*==============================================================================*/
/* III) Wannier minimization for proc=0                                            */

   if(num_proc > 1){
     Barrier(comm);
   }

   if(myid ==0){

    /*Full minimization*/
     lbfg_min(general_data,cp,elementA,NDIM,Z_real,Z_imag,W);

     comp_wannier_center(general_data,cp,elementA,Z_real,Z_imag,
                          W,wan_cent,U_final);

   }

/*===========================================================================*/
/* IV) Broadcast the optimized A and Wannier centers                         */

   if(num_proc > 1){

     Barrier(comm);

     Bcast(&elementA[1],NDIM,MPI_DOUBLE,0,comm);
     Bcast(&wan_cent[1][1],3*nstate,MPI_DOUBLE,0,comm);
     Bcast(&U_final[1],NDIM,MPI_DOUBLE,0,comm);

   }

/*============================================================================*/
/* V) Rotate orbital, velocities and force */

    ncoef_tot=nstate_up*ncoef_up;

    cp_rotation_prim_dvr(dvrc, icoef_form_up, dvrc_tmp, icoef_form_up, U_final,
                         alpha,beta,ioff_st,&(cp->cp_comm_state_pkg_dvr_up));

    for(i=1;i<=ncoef_tot;i++){
      dvrc[i]=dvrc_tmp[i];
    }

    cp_rotation_prim_dvr(dvrvc, icoef_form_up, dvrc_tmp, icoef_form_up, U_final,
                         alpha,beta,ioff_st,&(cp->cp_comm_state_pkg_dvr_up));

    for(i=1;i<=ncoef_tot;i++){
      dvrvc[i]=dvrc_tmp[i];
    }

    cp_rotation_prim_dvr(dvrfc, icoef_form_up, dvrc_tmp, icoef_form_up, U_final,
                         alpha,beta,ioff_st,&(cp->cp_comm_state_pkg_dvr_up));

    for(i=1;i<=ncoef_tot;i++){
      dvrfc[i]=dvrc_tmp[i];
    }

/*============================================================================*/
/* V) print out initial Wannier center and orbital */

  if(iwrite_init_wcent==1){
    if(myid==0){
      printf("\n");
      printf("     WFC X(A)          WFC Y(A)          WFC Z(A)     State\n\n");
      for (i=1;i<=nstate;i++){
        printf("  %- 18.5e%- 18.5e%- 18.5e%-8d\n",
                 wan_cent[i][1], wan_cent[i][2], wan_cent[i][3],i);
      }
    }
  }

  if(iwrite_init_worb==1){
    write_wannier_orb_dvr(general_data,cp); 
  }

/*=======================================================================*/
}/*end routine*/
/*==========================================================================*/

void write_wannier_orb_dvr(GENERAL_DATA *general_data,CP *cp)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/

#include "../typ_defs/typ_mask.h"

  int myid_state          = cp->communicate.myid_state;
  int np_states           = cp->communicate.np_states;
  MPI_Comm comm_states    = cp->communicate.comm_states;

  int iwrite_st       = cp->cpopts.iwrite_init_state;

  double *dvrc_all       =  cp->cpscr.cpscr_wave.zfft;
  double *dvrc_tmp       =  cp->cpscr.cpscr_wave.zfft_tmp;
  double *hmat           =  general_data->cell.hmat;

  double *dvrc       = cp->cpcoeffs_pos_dvr[1].dvrc_up;
  int icoef_orth     = cp->cpcoeffs_pos_dvr[1].icoef_orth_up;
  int icoef_form     = cp->cpcoeffs_pos_dvr[1].icoef_form_up;

  int nfft        = cp->cp_para_fft_pkg3d_sm.nfft;
  int nfft2       = nfft/2;
  int nfft_proc   = cp->cp_para_fft_pkg3d_sm.nfft_proc; 
  int nfft2_proc  = nfft_proc/2;
  int nkf1        = cp->cp_para_fft_pkg3d_sm.nkf1;
  int nkf2        = cp->cp_para_fft_pkg3d_sm.nkf2;
  int nkf3        = cp->cp_para_fft_pkg3d_sm.nkf3;

  int skc_fft_ka_proc      = cp->cp_para_fft_pkg3d_sm.skc_fft_ka_proc;
  int ekc_fft_ka_proc      = cp->cp_para_fft_pkg3d_sm.ekc_fft_ka_proc;
  int skb_fft_ka_proc      = cp->cp_para_fft_pkg3d_sm.skb_fft_ka_proc;
  int ekb_fft_ka_proc      = cp->cp_para_fft_pkg3d_sm.ekb_fft_ka_proc;

  int idiv,irem,ka,kb,kc,i,ioff,icount;
  double dxx,dyy,dzz,da,db,dc,x,y,z,sa,sb,sc;
  double vol,scale_fact;

  int *irecv,*idispl;

  FILE *fp_wan;

/*------------------------------------------------------------------------------- */
/* 0) Check the form of the coefficients  */

  if(icoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in orthogonal form    \n");
    printf("on state processor %d in cp_rho_calc_dvr_all_space \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1 && icoef_form == 0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed (not normal) \n");
    printf("form on state processor %d in cp_rho_calc_dvr_all_sace \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/


/*---------------------------------------------------------------------------*/
/* I) initialization */

  irecv  = (int *) malloc((np_states)*sizeof(int))-1;
  idispl = (int *) malloc((np_states)*sizeof(int))-1;

  dxx = hmat[1]/(2.0*(double) nkf1);
  dyy = hmat[5]/(2.0*(double) nkf2);
  dzz = hmat[9]/(2.0*(double) nkf3);

  da = 1.0/((double) nkf1);
  db = 1.0/((double) nkf2);
  dc = 1.0/((double) nkf3);

  for(i=1; i<= nfft2_proc ; i++){
    dvrc_all[i] = 0.0;
    dvrc_tmp[i] = 0.0;
  }

  vol = getdeth(hmat);
  scale_fact = sqrt(((double)(nkf1*nkf2*nkf3))/vol);

/*--------------------------------------------------------------------------*/
/* II) Select the state we want */

  ioff  = (iwrite_st-1)*nfft2_proc;
  for(i=1; i<= nfft2_proc ; i++){
    dvrc_tmp[i] = dvrc[(ioff+i)];
  }

/*--------------------------------------------------------------------------*/
/* III) Gather orbital on Proc 0: */

  if( np_states > 1){

    idiv    = (nkf2*nkf3)/np_states;
    irem    = (nkf2*nkf3)%np_states;

    for(i=0; i <irem; i++){ irecv[i+1] = (idiv+1)*nkf1;  }
    for(i=irem ; i < np_states; i++){ irecv[i+1] = idiv*nkf1;  }

    idispl[1] = 0;
    for(i=2; i <= np_states; i++){
      idispl[i] = idispl[i-1] + irecv[i-1];
    }

    Barrier(comm_states);
    Gatherv(&dvrc_tmp[1],nfft2_proc,MPI_DOUBLE,&dvrc_all[1],
            &irecv[1],&idispl[1],MPI_DOUBLE,0,comm_states);
  }else{
    for(i=1; i <= nfft2; i++){ dvrc_all[i] = dvrc_tmp[i]; }
  }

/*------------------------------------------------------------------------*/
/* IV) File open */ 

  if(myid_state==0){
    fp_wan = cfopen(general_data->filenames.worbname,"w");
  }


/*------------------------------------------------------------------------*/
/* V) print out the orbital */

  if(myid_state==0){
    icount=0;
    fprintf(fp_wan, "%d  %d  %d\n", nkf3,nkf2,nkf1);
    for(kc=1;kc<=nkf3;kc++){
      sc = dc*((double)(kc-1)) - 0.5;
      for(kb=1;kb<=nkf2;kb++){
        sb = db*((double)(kb-1)) - 0.5;
        for(ka=1;ka<=nkf1;ka++){
          sa = da*((double)(ka-1)) - 0.5;

          x  = sa*hmat[1]+sb*hmat[4]+sc*hmat[7];
          y  = sa*hmat[2]+sb*hmat[5]+sc*hmat[8];
          z  = sa*hmat[3]+sb*hmat[6]+sc*hmat[9];
          x += dxx;
          y += dyy;
          z += dzz;

          icount++;
          fprintf(fp_wan,"%.8g  %.8g  %.8g  %.10g\n",z,y,x,dvrc_all[icount]*scale_fact);
        }
      }
    }
  }

  if(np_states>1){Barrier(comm_states);}

  free(&irecv[1]);
  free(&idispl[1]);

/*===========================================================================*/
} /*end routine */
/*============================================================================*/

