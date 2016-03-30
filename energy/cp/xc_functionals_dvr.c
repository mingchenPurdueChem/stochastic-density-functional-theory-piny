/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: xc_functionals.c                               */
/* Contains all the exchange and correlation functionals                    */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_energy_cpcon_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_exc_real_space_dvr (CELL *cell,EWALD *ewald,PSEUDO *pseudo,
                            CPOPTS *cpopts,CPCOEFFS_INFO *cpcoeffs_info,
                            CPSCR_WAVE *cpscr_wave,PARA_FFT_PKG3D *cp_para_fft_pkg, 
                            CPCOEFFS_POS_DVR *cpcoeffs_pos_dvr,
                            COMMUNICATE *communicate, CPSCR *cpscr,
                            DVR_MATRIX *dvr_matrix, double *zfft, 
                            double *zfft_tmp,
                            CP_COMM_STATE_PKG *cp_comm_state_pkg_dvr_up,
                            CP_COMM_STATE_PKG *cp_comm_state_pkg_dvr_dn,
                            double *cp_exc_ret, double *cp_muxc_ret) 

/*=======================================================================*/
{/*Begin Routine*/
/*=======================================================================*/

#include "../typ_defs/typ_mask.h"

  /* Local Pointers */
  int static ifirst = 1;
 
  int cp_lsda           = cpopts->cp_lsda;
  int cp_gga            = cpopts->cp_gga;
  int cp_b3             = cpopts->cp_b3;
  int cp_lyp            = cpopts->cp_lyp;
  int cp_lypm1          = cpopts->cp_lypm1;

  int nstate_up_proc = cpcoeffs_info->nstate_up_proc;
  int nstate_up      = cpcoeffs_info->nstate_up;
  int nstate_dn_proc = cpcoeffs_info->nstate_dn_proc;
  int nstate_dn      = cpcoeffs_info->nstate_dn;

  int np_states         = communicate->np_states;
  int myid_state        = communicate->myid_state;
  MPI_Comm  comm_states = communicate->comm_states;

  int nfft              = cp_para_fft_pkg->nfft;
  int nfft_proc         = cp_para_fft_pkg->nfft_proc;
  double erfcg_corr     = pseudo->erfcg_corr;

  double *dvrc_up       = cpcoeffs_pos_dvr->dvrc_up;
  double *dvrfc_up      = cpcoeffs_pos_dvr->dvrfc_up;
  double *dvrc_dn       = cpcoeffs_pos_dvr->dvrc_dn;
  double *dvrfc_dn      = cpcoeffs_pos_dvr->dvrfc_dn;
  double *rho_up        = cpscr->cpscr_rho.rho_up;
  double *rho_dn        = cpscr->cpscr_rho.rho_dn;

  double *cre_up        = cpscr_wave->cre_up;
  double *cre_dn        = cpscr_wave->cre_dn;
 
  double vol_cp         = cell->vol_cp;
  double gga_cut        = pseudo->gga_cut;

/*==========================================================================*/
/* I) First time only calculation for exact exchange*/

  if(ifirst ==1 && cp_b3==1) {

    if( nstate_up_proc < 2 || nstate_dn_proc < 2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf(" No. of state per processor is smaller than 2.    \n");
      printf(" This will cause a problem in calc_singular_corr  \n");
      printf(" STOP!                                            \n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@@\n");
    }

    /* calc_singular_corr(cell,ewald,pseudo, cpcoeffs_info->dvr_grid_dens,
                       zfft,zfft_tmp,cpscr_wave->fft_map, cpscr_wave->fft_map_rev,
                       cp_para_fft_pkg3d,myid_state); */
  }


/*========================================================================*/
/* II) LDA part : Perdew-Zunger only for now */

   if(cp_lsda ==0){
     excpot_pz_lda_dvr(dvrfc_up,dvrc_up,rho_up, cp_exc_ret, cp_muxc_ret,
                       cell->vol_cp,nstate_up,nfft,nfft_proc,
                       cp_lyp,cp_lypm1,cp_b3,erfcg_corr);

     // printf("myid_state  %d,cp_exc_ret %f \n", myid_state,*cp_exc_ret);

   }else if(cp_lsda==1){

     /* excpot_pz_lsda_dvr(dvrfc_up,dvrc_up,rho_up, cp_exc_ret, cp_muxc_ret,
                       cell->vol_cp,nstate_up,nfft,nfft_proc,
                       cp_lyp,cp_lypm1,cp_b3,erfcg_corr);*/

  }

/*=======================================================================*/
/* III) GGA + Hybrid */

  if(cp_gga==1){

     /* exact exchange for hybrid */
    if(cp_b3==1){

       /* calc_exact_exchange(cp_exc_ret,cell,pseudo,cpscr, cpcoeffs_info,cp_lsda,
                           cpcoeffs_pos_dvr, cp_para_fft_pkg3d,communicate); */

    }


    /* NEED TO Clean Up these routines:
      i) Is is possible to reduce the number of transpose_bck,fwd..
      ii) Do we need seperate routines for cp_lsda or call these routine again
          for down states
    */

    grad_dens_dvr_par(dvrc_up,dvr_matrix, &(cpscr->cpscr_grho),
                      cpcoeffs_info, cre_up, cp_comm_state_pkg_dvr_up,
                      zfft,zfft_tmp, communicate, cp_para_fft_pkg);


    grad_corr_lda_dvr_par(cpscr,cpcoeffs_info,dvrc_up,dvrfc_up, 
                          ewald,cell,cp_exc_ret,cp_muxc_ret, vol_cp,
                          cpopts,pseudo->gga_cut,zfft,zfft_tmp,dvr_matrix,
                          communicate, cp_comm_state_pkg_dvr_up, 
                          cp_para_fft_pkg,nstate_up);
    
    if(cp_lsda==1){
      if(myid_state==0){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        printf("             IMPLEMENT CP_LSDA for GGA \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      }
      fflush(stdout);
      exit(1);
    }/*endifr*/
  
  }/*endif gga*/

/*==========================================================================*/
}/*end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void  excpot_pz_lda_dvr(double *dvrfc,double *dvrc,double *rho,
                        double *exc_ret,double *muxc_ret,double vol,
                        int nstate,int nfft,int nfft_proc,
                        int cp_lyp,int cp_lypm1, int cp_b3_opt, double erfcg_corr)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*==========================================================================*/

/* Static variables               */
/* lda correlation parameters     */
  static double gamma = -0.14230;
  static double beta1 =  1.05290;
  static double beta2 =  0.33340;
  static double au    =  0.03110;
  static double bu    = -0.0480;
  static double cu    =  0.00200;
  static double du    = -0.01160;

/* Local variables */

  int nfft2      = nfft/2;
  int nfft2_proc = nfft_proc/2;
  int i,ioff;
  double ftemp,fcons;
  double ex_scale;

  double vscale;
  double pi,rho_r,power,xfact,cfact,rs,fpi,fpin;
  double ex,ec,mufact,rat1,rat2,cf1,cf2,sqtrs,vxc;
  double ffact;
  double lnrs;
  double lyp_fact;
  int    is,ind;

  double b3_corr,ex_factor;
  double ex_scale_g = 0.2;

/*=================================================================*/
/* I) Assign use constants                                         */

   pi       =  M_PI;
   fpi      =  4.0*pi;
   rat1     = (7.0/6.0)*beta1;
   rat2     = (4.0/3.0)*beta2;
   power    = 1.0/3.0;
   lyp_fact = (double)(1-(cp_lyp || cp_lypm1));
   ex_scale=1.0000000000;

   if(cp_b3_opt==1) {
     ex_scale = 0.8000000000;
     lyp_fact=1.0000000000;
   }

   ex       = 0.0;
   ec       = 0.0;
   vxc      = 0.0;

   vscale = vol/((double)(nfft2));
   ex_factor = vscale*vscale*ex_scale_g*erfcg_corr;

/*=================================================================*/
/* II) compute the LDA energy and force contribution */

   b3_corr=0.0;

   for(i=1 ; i<= nfft2_proc; i++){
     rho_r = rho[i];
     fpin  = fpi*rho_r;
     rs    = pow((3.0/fpin),power);
     sqtrs = sqrt(rs);

     xfact = pow(((3.0*rho_r/pi)),power);
     xfact = -xfact;

     b3_corr += rho_r*rho_r;

     if(rs >= 1.0){
       cf1 = 1.0 + beta1*sqtrs + beta2*rs;
       cf2 = 1.0 + rat1*sqtrs + rat2*rs;
       cfact = gamma/cf1;
       mufact = cfact*(cf2/cf1);   /* derivative of Ec = ec*rho wrt to rho */
     }else{
       lnrs   = log(rs);
       cfact  = au*lnrs + bu + cu*rs*lnrs + du*rs;
       mufact = au*lnrs + (bu - power*au) + 2.0*power*cu*rs*lnrs
               + power*(2.0*du - cu)*rs;   /* derivative of Ec = ec*rho wrt to rho */
     }/*endif*/

     /* For B3LYP, 0.8Ex[slater] */
     xfact=xfact*ex_scale;

     /* For B3LYP, Ec[PZ] */
     cfact = cfact*lyp_fact;  
     ex    = ex + xfact*rho_r;
     ec    = ec + cfact*rho_r;

     /* lyp_fact = 0 if gga on */
     mufact = mufact*lyp_fact;  
     vxc   += (xfact + mufact)*rho_r;

     /* forces on coefficients */
     ffact = -2.0*(xfact + mufact);
     for(is=1; is <= nstate; is++){
       ind = (is-1)*nfft2_proc + i;
       //test dvrfc
       //  if((i==44707)&&is==2){
       //printf("nfft2_proc i %d %d ind %d dvrfc[ind] %f ffact %f dvrc[ind] %f vscale %f 1\n", nfft2_proc,i,ind,dvrfc[ind],ffact,dvrc[ind],vscale);}
       dvrfc[ind] += ffact*dvrc[ind]*vscale;
       //test dvrfc
       //if((i==44707)&&is==2){
       //printf("nfft2_proc i %d %d ind %d dvrfc[ind] %f ffact %f dvrc[ind] %f vscale %f 2\n", nfft2_proc,i,ind,dvrfc[ind],ffact,dvrc[ind],vscale);}
       if(cp_b3_opt==1) {
         dvrfc[ind] += 2.0*dvrc[ind]*rho_r*ex_factor;
       }
     }
  }/*endfor my grid points*/

  ex *= 0.750;

  *exc_ret += (ex+ec)*vscale;

  if(cp_b3_opt==1) {
    *exc_ret += -0.5*b3_corr*ex_factor;
    /* muxc_ret = ??? */
   }

   *muxc_ret += vxc*vscale;

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/

void grad_dens_dvr_par(double *dvrc,DVR_MATRIX *dvr_matrix,
                      CPSCR_GRHO *cpscr_grho,CPCOEFFS_INFO *cpcoeffs_info,
                      double *cre_up, CP_COMM_STATE_PKG *cp_comm_state_pkg,
                      double *zfft,double *zfft_tmp,
                      COMMUNICATE *communicate,PARA_FFT_PKG3D *cp_para_fft_pkg)

/*==========================================================================*/
 {/*begin routine*/
/*==========================================================================*/


#include "../typ_defs/typ_mask.h"

  /*   Local Pointer declarations   */

  double *grad_rho_tmp     = zfft;
  double *grad_rho_tmp_all = zfft_tmp;

  double *dGx   = dvr_matrix->dGx;
  double *dGy   = dvr_matrix->dGy;
  double *dGz   = dvr_matrix->dGz;

  int  nkf1 = cp_para_fft_pkg->nkf1;
  int  nkf2 = cp_para_fft_pkg->nkf2;
  int  nkf3 = cp_para_fft_pkg->nkf3;

  int myid                  = communicate->myid_state;
  int nproc                 = communicate->np_states;
  MPI_Comm comm             = communicate->comm_states;


  /* Gradient of rho (del rho) */
  double *d_rhox            = cpscr_grho->d_rhox_up;
  double *d_rhoy            = cpscr_grho->d_rhoy_up;
  double *d_rhoz            = cpscr_grho->d_rhoz_up;

  /* derivative of del rho wrt coef C(is,ijk) */
  double *d_rhox_c          = cpscr_grho->d_rhox_up_c;
  double *d_rhoy_c          = cpscr_grho->d_rhoy_up_c;
  double *d_rhoz_c          = cpscr_grho->d_rhoz_up_c;

  int ncoef            = cp_para_fft_pkg->nfft/2;
  int nstate_proc      = cpcoeffs_info->nstate_up_proc;
  int ncoef_proc       = cp_para_fft_pkg->nfft_proc/2;
  int *ioff_up_st      = cpcoeffs_info->ioff_up;
  int icoef_form_tmp ;


  int ia,ib,ic,iap,ibp,icp,i,is;
  int index1,index2,index3,ioff;
  int icount;
  double sum;


/*--------------------------------------------------------------------------*/
/* I) Zero Arrays */

  for(i=1; i<= ncoef; i++){
    d_rhox[i] = d_rhoy[i] = d_rhoz[i] = 0.00;
    grad_rho_tmp[i] = 0.0;
    grad_rho_tmp_all[i]=0.0;
  }

  for(is=1; is <= nstate_proc; is++){
    for(i=1; i<= ncoef; i++){
      index1 = ioff_up_st[is] + i;
      d_rhox_c[index1] = d_rhoy_c[index1] = d_rhoz_c[index1] = 0.00;
    }
  }

/*--------------------------------------------------------------------------*/
/* II) convert into nstate_proc form                                            */

  if(nproc > 1){
     Barrier(comm);
     icoef_form_tmp=1;
     cp_transpose_bck_dvr(dvrc,&(icoef_form_tmp), cre_up,
                          cp_comm_state_pkg, nkf1,nkf2,nkf3);
  }

/*--------------------------------------------------------------------------*/
/* III) density gradient with respect to X                                       */

  icount=0;
  for(is=1; is <= nstate_proc; is++){

    for(ic=1; ic <= nkf3; ic++){
      for(ib=1; ib <= nkf2; ib++){
        for(ia=1;  ia  <= nkf1; ia++){
          index1 = (ic-1)*nkf1*nkf2 + (ib-1)*nkf1 + ia;
          sum=0.0;
          for(iap=1; iap <= nkf1; iap++){
            index2 = (iap-1)*nkf1 + ia;
            index3 = ioff_up_st[is]
                     + (ic-1)*nkf1*nkf2 + (ib-1)*nkf1 + iap;
            sum+= dvrc[index3]*dGx[index2];
          }/* iap */
          grad_rho_tmp[index1]+=sum*dvrc[index1+ioff_up_st[is]];
          icount+=1;
          d_rhox_c[icount]=2.0*sum;
        }/* ia */
      }/* ib */
    }/* ic */
  }/* is */

  /*++++++++MA++++++Test++++++*/

  /* printf("myid %d, nkf1 %d, nkf2 %d, nkf3 %d, nstate_proc %d, nproc %d\n", myid,nkf1,nkf2,nkf3,nstate_proc,nproc);
  for(i=1;i<=nstate_proc;i++){
    printf("myid %d is %d ioff_up_st %d %f\n",myid, i,ioff_up_st[i],dvrc[nkf1*nkf2+nkf1+20]);
  }
  
  for(i=1;i<=ncoef;i+=20){
    printf("myid %d i %d grad_rho_tmp %f\n",myid, i,grad_rho_tmp[i]);
    }*/

  /* sum over state*/
  if(nproc > 1){
    Allreduce(&(grad_rho_tmp[1]),&(grad_rho_tmp_all[1]),ncoef,
             MPI_DOUBLE,MPI_SUM,0,comm);
    for(i=1; i <= ncoef; i++){d_rhox[i] = 2.0*grad_rho_tmp_all[i];}
  }else{
    for(i=1; i <= ncoef; i++){d_rhox[i] = 2.0*grad_rho_tmp[i];}
  }
  /* for(i=1;i<=ncoef;i+=20){
    printf("myid %d i %d grad_rho_tmp_all %f\n",myid, i, grad_rho_tmp_all[i]);
    }*/
/*--------------------------------------------------------------------------*/
/* IV) density gradient with respect to Y                                       */

  for(i=1; i<= ncoef; i++){ grad_rho_tmp_all[i] = 0.0;}

  for(is=1; is <= nstate_proc; is++){

    for(i=1; i<= ncoef; i++){ grad_rho_tmp[i] = 0.0;}

    /* reorder the coefficients */
    ioff=0;
    for (ic=1; ic<= nkf3;ic++){
      for(ia=1; ia <= nkf1; ia++){
        for(ib=1; ib <= nkf2; ib++){
          ioff++;
          index1=ioff_up_st[is]+(ic-1)*nkf1*nkf2+(ib-1)*nkf1+ia;
          cre_up[ioff]=dvrc[index1];
        }
      }
    }

    ioff=0;
    for(ic=1; ic <= nkf3; ic++){
      for(ib=1; ib <= nkf2; ib++){
        for(ia=1; ia <= nkf1; ia++){
          ioff++;
          for(ibp = 1; ibp <=nkf2; ibp++){
            index2 = (ibp-1)*nkf2 + ib;
            index3 = (ic-1)*nkf1*nkf2 + (ia-1)*nkf2 + ibp;
            grad_rho_tmp[ioff] += cre_up[index3]*dGy[index2];
          }/* ibp */
        }/* ia */
      }/* ib  */
    }/* ic */

    icount=0;
    for(ic=1; ic <= nkf3; ic++){
      for(ib=1; ib <= nkf2; ib++){
        for(ia=1; ia <= nkf1; ia++){
          icount++;
          d_rhoy_c[(is-1)*ncoef+icount]=2.0*grad_rho_tmp[icount];
          grad_rho_tmp_all[icount] += dvrc[ioff_up_st[is]+icount]
                                *grad_rho_tmp[icount];
        }/* ia */
      }/* ib */
    }/* ic */

  }/* is */

  /* sum over states*/
  if(nproc > 1){

    icount=0;
    for(ic=1; ic <= nkf3; ic++){
      for(ib=1; ib <= nkf2; ib++){
        for(ia=1; ia <= nkf1; ia++){
          icount++;
          grad_rho_tmp[icount] = grad_rho_tmp_all[icount];
        }
      }
    }

    Allreduce(&(grad_rho_tmp[1]),&(grad_rho_tmp_all[1]),ncoef,
             MPI_DOUBLE,MPI_SUM,0,comm);
  }

  for(i=1; i <= ncoef; i++){d_rhoy[i] = 2.0*grad_rho_tmp_all[i];}


/*--------------------------------------------------------------------------*/
/* V) density gradient with respect to Z                                       */

  for(i=1; i<= ncoef; i++){ grad_rho_tmp_all[i] = 0.0;}

  for(is=1; is <= nstate_proc; is++){

    for(i=1; i<= ncoef; i++){ grad_rho_tmp[i] = 0.0;}

    /* reorder */
    ioff=0;
    for (ib=1;ib<=nkf2;ib++){
      for (ia=1;ia<=nkf1;ia++){
        for (ic=1;ic<=nkf3;ic++){
          ioff++;
          index1=ioff_up_st[is]+(ic-1)*nkf1*nkf2+(ib-1)*nkf1+ia;
          cre_up[ioff]=dvrc[index1];
        }
      }
    }

    ioff=0;
    for(ic=1; ic <= nkf3; ic++){
      for(ib=1; ib <= nkf2; ib++){
        for(ia=1; ia <= nkf1; ia++){
          ioff++;
          for(icp=1;icp<=nkf3;icp++){
            index2=(icp-1)*nkf3+ic;
            index3 = (ib-1)*nkf1*nkf3+(ia-1)*nkf3+icp;
            grad_rho_tmp[ioff] += cre_up[index3]*dGz[index2];
          }/* icp */
        }/* ia */
      }/* ib  */
    }/* ic */

    icount=0;
    for(ic=1; ic <= nkf3; ic++){
      for(ib=1; ib <= nkf2; ib++){
        for(ia=1; ia <= nkf1; ia++){
          icount++;
          d_rhoz_c[(is-1)*ncoef+icount]=2.0*grad_rho_tmp[icount];
          grad_rho_tmp_all[icount] += dvrc[ioff_up_st[is]+icount]
                                *grad_rho_tmp[icount];
        }/* ia */
      }/* ib */
    }/* ic */

  }/* is */

  /* sum over states*/
  if(nproc > 1){
    icount=0;
    for(ic=1; ic <= nkf3; ic++){
      for(ib=1; ib <= nkf2; ib++){
        for(ia=1; ia <= nkf1; ia++){
          icount++;
          grad_rho_tmp[icount] = grad_rho_tmp_all[icount];
        }
      }
    }

    Allreduce(&(grad_rho_tmp[1]),&(grad_rho_tmp_all[1]),ncoef,
             MPI_DOUBLE,MPI_SUM,0,comm);
  }

  for(i=1; i <= ncoef; i++){d_rhoz[i] = 2.0*grad_rho_tmp_all[i];}

/*--------------------------------------------------------------------------*/
/* VI) convert back to ncoef_proc form                                          */

  if(nproc > 1){
     Barrier(comm);
     icoef_form_tmp=0;
     cp_transpose_fwd_dvr(dvrc,&(icoef_form_tmp), cre_up,
                         cp_comm_state_pkg, nkf1,nkf2,nkf3);
   }


/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

void grad_corr_lda_dvr_par(CPSCR *cpscr, CPCOEFFS_INFO *cpcoeffs_info,
                          double *dvrc,double *dvrfc,EWALD *ewald, CELL *cell,
                          double *exc,double *muxc,double vol,CPOPTS *cpopts,
                          double gc_cut,double *dvrfc_tmp,double *dvrfc_tmp_all,
                          DVR_MATRIX *dvr_matrix,COMMUNICATE *communicate,
                          CP_COMM_STATE_PKG *cp_comm_state_pkg,
                          PARA_FFT_PKG3D *cp_para_fft_pkg3d_gen,int nstate)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*                 Local variable declarations                           */

#include "../typ_defs/typ_mask.h"
/*------------------------------*/
/* density and gradient */

   double rho_p,g_rho2,g_rhoi,g_rho;
   double unit_gx,unit_gy,unit_gz;

/*------------------------------*/
/*exchange and correlation energy */

   double  ex_gc,ec_gc,muxc_loc;
/*------------------------------*/
/*  cg functional (f_xc), df_xc/dn, and df_xc/d|grad n|  */

   double fxc_x,dfxc_dn_x,dfxc_dgn_x;
   double fxc_c,dfxc_dn_c,dfxc_dgn_c;

   double dfxc_dgn_x_c;

/*------------------------------*/
/* Switching function */

   double srho,dsrho;
   double rho_heali,rho_cut,rsw;

/*------------------------------*/
/* Forces */
   double force_dn_xc,force_dgn_xc;
   double force_dgn_xc_x,force_dgn_xc_y,force_dgn_xc_z;

/*------------------------------*/
/* loop variables */
   int i,m,k,j,is,index;
   int iup=1;

/*----------------------------------------------------------*/
/* Becke beta parameter (can vary depending on functional)  */
   double beta = 0.0042;

/* Assign local pointers                                    */
    double *cre_up     = cpscr->cpscr_wave.cre_up;
    double *del_rho_up = cpscr->cpscr_wave.cim_up;

    double *d_rho_x_c  = cpscr->cpscr_grho.d_rhox_up_c;
    double *d_rho_y_c  = cpscr->cpscr_grho.d_rhoy_up_c;
    double *d_rho_z_c  = cpscr->cpscr_grho.d_rhoz_up_c;

    double    *rho_up          = cpscr->cpscr_rho.rho_up;
    double    *del_rho_up_x    = cpscr->cpscr_grho.d_rhox_up;
    double    *del_rho_up_y    = cpscr->cpscr_grho.d_rhoy_up;
    double    *del_rho_up_z    = cpscr->cpscr_grho.d_rhoz_up;


    int cp_becke               = cpopts->cp_becke;
    int cp_b3                  = cpopts->cp_b3;
    int cp_fila_1x             = cpopts->cp_fila_1x;
    int cp_fila_2x             = cpopts->cp_fila_2x;
    int cp_pw91x               = cpopts->cp_pw91x;
    int cp_pw91c               = cpopts->cp_pw91c;
    int cp_lyp                 = cpopts->cp_lyp;
    int cp_lypm1               = cpopts->cp_lypm1;

    int   nfft_proc            = cp_para_fft_pkg3d_gen->nfft_proc;
    int   nfft                 = cp_para_fft_pkg3d_gen->nfft;
    int   nfft2                = nfft/2;
    int   nfft2_proc           = nfft_proc/2;

    int  nkf1 = cp_para_fft_pkg3d_gen->nkf1;
    int  nkf2 = cp_para_fft_pkg3d_gen->nkf2;
    int  nkf3 = cp_para_fft_pkg3d_gen->nkf3;
    int  ia,ib,ic,iap,ibp,icp,index1,index2,ix,iy,iz;

    int myid                  = communicate->myid_state;
    int nproc                 = communicate->np_states;
    MPI_Comm comm             = communicate->comm_states;

 int *skb_fft_ka_proc_all  = cp_para_fft_pkg3d_gen->skb_fft_ka_proc_all;
 int *skc_fft_ka_proc_all  = cp_para_fft_pkg3d_gen->skc_fft_ka_proc_all;
 int *ekb_fft_ka_proc_all  = cp_para_fft_pkg3d_gen->ekb_fft_ka_proc_all;
 int *ekc_fft_ka_proc_all  = cp_para_fft_pkg3d_gen->ekc_fft_ka_proc_all;

 int  skc_fft_ka_proc      = cp_para_fft_pkg3d_gen->skc_fft_ka_proc;
 int  ekc_fft_ka_proc      = cp_para_fft_pkg3d_gen->ekc_fft_ka_proc;
 int  skb_fft_ka_proc      = cp_para_fft_pkg3d_gen->skb_fft_ka_proc;
 int  ekb_fft_ka_proc      = cp_para_fft_pkg3d_gen->ekb_fft_ka_proc;

 double *dGx   = dvr_matrix->dGx;
 double *dGy   = dvr_matrix->dGy;
 double *dGz   = dvr_matrix->dGz;
 double volscale = cpcoeffs_info->scale_fact;
 int kb_str,kb_end;
 int idiv,irem,icount,itot,iproc,ioff;

 int icoef_start_upm       = cpcoeffs_info->icoef_start_up - 1;
 int ncoef_proc            = cp_para_fft_pkg3d_gen->nfft_proc/2;

 int nstate_proc   = cpcoeffs_info->nstate_up_proc;
 int *ioff_up_st      = cpcoeffs_info->ioff_up;
 int icoef_form_tmp;

 double *del_rho_up_tmp_x = dvrfc_tmp;
 double *del_rho_up_tmp_y = dvrfc_tmp_all;
 double *del_rho_up_tmp_z = &dvrfc_tmp_all[nfft2];
 double *force_tmp = &dvrfc_tmp[nfft2];

 /* for B3LYP */
 double ex_scale,ec_scale;

/*======================================================================*/
/* I) Initialize and malloc variables                                      */

   rho_heali=1.0/(3.9*gc_cut);
   rho_cut = 0.1*gc_cut;

   ex_gc         = 0.0;
   ec_gc         = 0.0;
   muxc_loc      = 0.0;

   fxc_x         = 0.0;
   dfxc_dn_x     = 0.0;
   dfxc_dgn_x    = 0.0;

   fxc_c         = 0.0;
   dfxc_dn_c     = 0.0;
   dfxc_dgn_c    = 0.0;

   ex_scale=1.0000000;
   ec_scale=1.0000000;

   if(cp_b3==1){
     ex_scale=0.72000000;
     ec_scale=0.81000000;
   }

   volscale = vol/((double)(nfft2));

  for(icount=1;icount<=nfft2;icount++){
    del_rho_up_tmp_x[icount]=0.0;
    del_rho_up_tmp_y[icount]=0.0;
    del_rho_up_tmp_z[icount]=0.0;
  }

/*=======================================================================*/
/* ncoef_proc operation for Exc                                          */
/* Loop over real space grid                                             */

  icount = 0;
  for(ic = skc_fft_ka_proc; ic <= ekc_fft_ka_proc; ic++){
  kb_str = (ic == skc_fft_ka_proc ? skb_fft_ka_proc : 1);
  kb_end = (ic == ekc_fft_ka_proc ? ekb_fft_ka_proc : nkf2);
  for(ib=kb_str; ib <= kb_end; ib++){
  for(ia = 1; ia <= nkf1; ia++){

    icount++;
    index=(ic-1)*nkf2*nkf1+(ib-1)*nkf1+ia;
    rho_p = rho_up[icount];

    if(rho_p > rho_cut) {
      g_rho2 = del_rho_up_x[index]*del_rho_up_x[index]
             + del_rho_up_y[index]*del_rho_up_y[index]
             + del_rho_up_z[index]*del_rho_up_z[index];


/*------------------------------------------------------------------------*/
/* iv) Switching functional                                               */

      rsw   = (rho_p - rho_cut)*rho_heali;
      rsw   = MIN(rsw,1.0);
      rsw   = MAX(rsw,0.0);
      srho  = rsw*rsw*(3.0-2.0*rsw);
      dsrho = 6.0*rsw*(1.0-rsw)*rho_heali;

/*------------------------------------------------------------------------*/
/* v) Exchange functional                                                 */

      if(cp_becke==1 || cp_b3==1) becke_gcx_lda(rho_p,g_rho2,&fxc_x,&dfxc_dn_x,&dfxc_dgn_x,beta);
      if(cp_fila_1x==1) fila_1x_lda(rho_p,g_rho2,&fxc_x,&dfxc_dn_x,&dfxc_dgn_x);
      if(cp_fila_2x==1) fila_2x_lda(rho_p,g_rho2,&fxc_x,&dfxc_dn_x,&dfxc_dgn_x);
      if(cp_pw91x==1) pw91_gcx(rho_p,g_rho2,&fxc_x,&dfxc_dn_x,&dfxc_dgn_x);

      fxc_x *= ex_scale;
      dfxc_dn_x *= ex_scale;
      dfxc_dgn_x *= ex_scale;

      dfxc_dn_x    = (dfxc_dn_x*srho + fxc_x*dsrho);
      dfxc_dgn_x  *= srho;
      fxc_x       *= srho;

/*-------------------------------------------------------------------------*/
/* vi) Correlation functional                                              */

      if(cp_pw91c==1)  pw91_gcc(rho_p,g_rho2,&fxc_c,&dfxc_dn_c,&dfxc_dgn_c);
      if(cp_lyp==1) {
        if(cp_b3==1){
           b3_lyp_gcc(rho_p,g_rho2,&fxc_c,&dfxc_dn_c,&dfxc_dgn_c);
        }else{
           lyp_gcc(rho_p,g_rho2,&fxc_c,&dfxc_dn_c,&dfxc_dgn_c);
        }
      }

      fxc_c *= ec_scale;
      dfxc_dn_c *= ec_scale;
      dfxc_dgn_c *= ec_scale;

      dfxc_dn_c   = (dfxc_dn_c*srho + fxc_c*dsrho);
      dfxc_dgn_c *= srho;
      fxc_c      *= srho;

/*-------------------------------------------------------------------------*/
/* vii) Add contributions to exchange and correlation energy               */

      ex_gc +=  fxc_x;
      ec_gc +=  fxc_c;

/*-------------------------------------------------------------------------*/
/* viii) Get non-gradiant contribution to forces                           */

      muxc_loc   += (dfxc_dn_x + dfxc_dn_c)*rho_p;

/*-------------------------------------------------------------------------*/
/* x) Get unit vector del rho/|del rho| and multiply it                    */
/*                                                 times fxn derivative    */

      g_rho = sqrt(g_rho2);
      g_rhoi = 1.0/g_rho;
      unit_gx = del_rho_up_x[index]*g_rhoi;
      unit_gy = del_rho_up_y[index]*g_rhoi;
      unit_gz = del_rho_up_z[index]*g_rhoi;

      dfxc_dgn_x_c = (dfxc_dgn_x + dfxc_dgn_c);
      if(nproc > 1){
        del_rho_up_tmp_x[index] = unit_gx*dfxc_dgn_x_c;
        del_rho_up_tmp_y[index] = unit_gy*dfxc_dgn_x_c;
        del_rho_up_tmp_z[index] = unit_gz*dfxc_dgn_x_c;
      }else{
        del_rho_up_x[index] = unit_gx*dfxc_dgn_x_c;
        del_rho_up_y[index] = unit_gy*dfxc_dgn_x_c;
        del_rho_up_z[index] = unit_gz*dfxc_dgn_x_c;
      }

/*-------------------------------------------------------------------------*/
/* Calculate (partial Exc/partial rho ) contribution to the coef force */

      for(is=1; is <= nstate; is++){
        index2 = (is-1)*nfft2_proc + icount;
        force_dn_xc  = (dfxc_dn_x + dfxc_dn_c)*2.0*dvrc[index2];

        muxc_loc -= force_dn_xc*rho_p;
        dvrfc[index2] -= (force_dn_xc)*volscale;
      }
/*-------------------------------------------------------------------------*/
    }else{
      if(nproc > 1){
        del_rho_up_tmp_x[index] = 0.0;
        del_rho_up_tmp_y[index] = 0.0;
        del_rho_up_tmp_z[index] = 0.0;
      }else{
        del_rho_up_x[icount] = 0.0;
        del_rho_up_y[icount] = 0.0;
        del_rho_up_z[icount] = 0.0;
      }
    }/* endif rho_p>gc_cut*/

  }/* ia */
  }/* ib */
  }/* ic */

/*==========================================================================*/
/* Parallel communications                                                  */

  if(nproc > 1){
    Allreduce(&(del_rho_up_tmp_x[1]),&(del_rho_up_x[1]),nfft2,
             MPI_DOUBLE,MPI_SUM,0,comm);
    Allreduce(&(del_rho_up_tmp_y[1]),&(del_rho_up_y[1]),nfft2,
             MPI_DOUBLE,MPI_SUM,0,comm);
    Allreduce(&(del_rho_up_tmp_z[1]),&(del_rho_up_z[1]),nfft2,
             MPI_DOUBLE,MPI_SUM,0,comm);

  }

/*==========================================================================*/
/* change the force to nstate_proc form                                     */

  if(nproc > 1){
     Barrier(comm);
     icoef_form_tmp=1;
     cp_transpose_bck_dvr(dvrc,&(icoef_form_tmp), cre_up,
                          cp_comm_state_pkg, nkf1,nkf2,nkf3);
     icoef_form_tmp=1;
     cp_transpose_bck_dvr(dvrfc,&(icoef_form_tmp), cre_up,
                          cp_comm_state_pkg, nkf1,nkf2,nkf3);
   }

/*-------------------------------------------------------------------------*/
/* Calculate (partial Exc/partial |del rho| ) contribution to  coef force */

  /*X component*/

  for(is=1; is<= nstate_proc; is++){
    for(ic=1; ic <= nkf3; ic++){
      for(ib=1; ib <= nkf2; ib++){
        icount=(ic-1)*nkf1*nkf2 + (ib-1)*nkf1;
        for(ia=1; ia<= nkf1; ia++){
          index1 = ioff_up_st[is] + icount + ia;
          force_dgn_xc_x=0.0;
          for(iap=1; iap <= nkf1; iap++){
            index2 = ioff_up_st[is] + icount + iap;
            ix = (ia-1)*nkf1 + iap;
            index = icount + iap;
            force_dgn_xc_x += (del_rho_up_x[index]*(2.0*dvrc[index2]*dGx[ix]));
          }
          dvrfc[index1] -= (force_dgn_xc_x)*volscale;
    }}}
  }/*endfor is*/

  ioff=0;
  for (ic=1; ic<= nkf3;ic++){
    for(ia=1; ia <= nkf1; ia++){
      for(ib=1; ib <= nkf2; ib++){
        ioff++;
        index = (ic-1)*nkf1*nkf2+(ib-1)*nkf1+ia;
        del_rho_up[ioff] = del_rho_up_y[index];
      }
    }

  }


   /* Y COMPONENT */
  for(is=1; is <= nstate_proc; is++){

    ioff=0;
    for (ic=1; ic<= nkf3;ic++){
      for(ia=1; ia <= nkf1; ia++){
        for(ib=1; ib <= nkf2; ib++){
          ioff++;
          index1=ioff_up_st[is]+(ic-1)*nkf1*nkf2+(ib-1)*nkf1+ia;
          cre_up[ioff]=dvrc[index1];
        }
      }
    }


    ioff=0;
    for(ic = 1; ic <= nkf3; ic++){
      for(ib=1; ib <= nkf2; ib++){
        for(ia=1; ia <= nkf1; ia++){
          ioff++;
          force_dgn_xc_y=0.0;
          for(ibp=1; ibp <= nkf2; ibp++){
            index = (ic-1)*nkf1*nkf2+(ia-1)*nkf2+ibp;
            iy = (ib-1)*nkf2 + ibp;
            force_dgn_xc_y += (del_rho_up[index]*(2.0*cre_up[index]*dGy[iy]));
          }
          dvrfc[ioff_up_st[is]+ioff] -= (force_dgn_xc_y)*volscale;
    }}}

  }/* endfor is*/

  ioff=0;
  for (ib=1;ib<=nkf2;ib++){
    for (ia=1;ia<=nkf1;ia++){
      for (ic=1;ic<=nkf3;ic++){
        ioff++;
        index = (ic-1)*nkf1*nkf2+(ib-1)*nkf1+ia;
        del_rho_up[ioff] = del_rho_up_z[index];
      }
    }
  }


   /* Z COMPONENT */
  for(is=1; is<= nstate_proc; is++){

    ioff=0;
    for (ib=1;ib<=nkf2;ib++){
      for (ia=1;ia<=nkf1;ia++){
        for (ic=1;ic<=nkf3;ic++){
          ioff++;
          index1=ioff_up_st[is]+(ic-1)*nkf1*nkf2+(ib-1)*nkf1+ia;
          cre_up[ioff]=dvrc[index1];
        }
      }
    }

    ioff=0;
    for(ic=1; ic <= nkf3; ic++){
      for(ib= 1; ib <= nkf2; ib++){
        for(ia=1; ia<= nkf1; ia++){
          ioff++;
          force_dgn_xc_z=0.0;
          for(icp= 1; icp <= nkf3; icp++){
            iz = (ic-1)*nkf3 + icp;
            index = (ib-1)*nkf1*nkf3+(ia-1)*nkf3+icp;
            force_dgn_xc_z += (del_rho_up[index]*(2.0*cre_up[index]*dGz[iz]));
          }
          dvrfc[ioff_up_st[is]+ioff] -= (force_dgn_xc_z)*volscale;

    }}}
  }/*endfor is*/

   /* DIAGONAL TERM OTHER PIECE */
  for(is=1; is<= nstate_proc; is++){
    for(ic=1; ic <= nkf3; ic++){
      for(ib=1; ib <= nkf2; ib++){
        for(ia=1; ia<= nkf1; ia++){
          index1 = ioff_up_st[is] + (ic-1)*nkf1*nkf2 + (ib-1)*nkf1 + ia;
          index = (ic-1)*nkf1*nkf2 + (ib-1)*nkf1 + ia ;

          ix = (ia-1)*nkf1 + ia;
          force_dgn_xc_x = (del_rho_up_x[index]*(d_rho_x_c[index1]));

          iy = (ib-1)*nkf2 + ib;
          force_dgn_xc_y = (del_rho_up_y[index]*(d_rho_y_c[index1]));

          iz = (ic-1)*nkf3 + ic;
          force_dgn_xc_z = (del_rho_up_z[index]*(d_rho_z_c[index1]));

          dvrfc[index1] -= (force_dgn_xc_x
                         +  force_dgn_xc_y
                         +  force_dgn_xc_z)*volscale;

  }}}}

/*===========================================================================*/
/*IV) add gradient corrections to energy and kohn-sham potential             */

    ex_gc *= volscale;
    ec_gc *= volscale;
    muxc_loc *= volscale;
    *exc +=  (ex_gc+ec_gc);
    *muxc += muxc_loc;

/*==========================================================================*/
 }/* end routine */
/*==========================================================================*/

void  b3_lyp_gcc(double rho,double g_rho2,double *fn,double *df_dn,double *df_dgn)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/
/*  local variables   */
   double rho2,drho;
   double denom,denom2,cexp,omega,domega,delta,ddelta;
   double arho,darho;
   double dprho;
   double rho13,rho113,rho43;
   double term1,brack1,brack2,brack3,term2;
   double dterm1,dterm2;

/* static variables */
   static double a = 0.049180;
   static double b = 0.1320;
   static double c = 0.25330;
   static double d = 0.3490;
   static double pi = M_PI;
   static double pow13;
   static double pow113;
   static double pow43;
   static double pow53;
   static double pow83;
   static double cf;

/*=======================================================================*/
/* I) Get square root density coefficient                                */

   pow13 = -1.0/3.0;
   pow113 = 11.0*pow13;
   pow43 = 4.0*pow13;
   pow53 = -5.0*pow13;
   pow83 = -8.0*pow13;
   cf = 0.30*pow((3.0*pi*pi),(2.0/3.0));
   drho = sqrt(g_rho2);

/*=======================================================================*/
/*  II) evaluate some grouped quantities                                 */

   rho13  = pow(rho,pow13);
   rho113 = pow(rho,pow113);
   rho43  = rho13/rho;
   denom  = 1.0 + d*rho13;
   denom2 = denom*denom;
   cexp   = exp(-c*rho13);

   omega  = (cexp*rho113)/denom;
   brack1 = (d*rho43 + denom*(c*rho43 - 11.0/rho))/3.0;
   domega = omega*brack1/denom;

   delta  = c*rho13 + d*rho13/denom;
   brack1 = (d*rho13 - denom)*d*rho43/denom2;
   ddelta = (-c*rho43 + brack1)/3.0;

/*=======================================================================*/
/* III) get all the brackets                                             */

   rho2   = rho*rho;
   term1 = 0.0;
   brack1 = ((5.0 - 7.0/6.0*delta)/3.0)*g_rho2;
   brack2 = brack1;
   brack3 = (0.250*brack2 - 11.0*g_rho2/24.0)*rho2;
   term2  = -a*b*omega*brack3;

   arho   = 0.250*brack2 - 11.0*g_rho2/24.0;
   brack1 = -7.0/18.0*ddelta*g_rho2;
   darho  = 0.250*brack1;

   dprho  = 2.0*arho*rho + darho*rho2;
   dterm1 = 0.0;
   dterm2 = -a*b*(domega*brack3 + omega*dprho);

/*=======================================================================*/
/* IV) evaluate the function                                             */

   *fn = term1 + term2;
/*=======================================================================*/
/*  V) evaluate the derivatives of the funcion                           */

    *df_dn = dterm1 + dterm2;
    brack1 = 0.50*((5.0 - 7.0/6.0*delta)/3.0) - 11.0/12.0;
    *df_dgn = -a*b*omega*brack1*rho2*drho;

/*==============================================================================*/
}/*end routine*/
/*==============================================================================*/

