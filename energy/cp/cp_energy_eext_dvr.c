/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: cp_energy_pot.c                                */
/*                                                                          */
/* This routine calls the required force and PE routines                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/



#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void eext_loc_real_space_dvr(CLATOMS_POS *clatoms_pos,ATOMMAPS *atommaps,
                             PSEUDO *pseudo,EWD_SCR *ewd_scr,CELL *cell,
                             COMMUNICATE *communicate, double *rho,
                             double *dvrc,double *dvrfc,int natm,int nstate,
                             double *eext_ret,PARA_FFT_PKG3D *cp_para_fft_pkg)

/*======================================================================*/
{/*Begin Routine*/
/*======================================================================*/
/*         Local Variable declarations                                  */

#include "../typ_defs/typ_mask.h"

/*--------------------------------------------*/
/*         Local Pointer declarations         */
  int natm_typ      = atommaps->natm_typ;
  int *iatm_atm_typ = atommaps->iatm_atm_typ;

  int    n_ang_max  = pseudo->n_ang_max;
  int    *loc_opt   = pseudo->loc_opt;

  double *vps0      = pseudo->vps0;
  double *vps1      = pseudo->vps1;
  double *vps2      = pseudo->vps2;
  double *vps3      = pseudo->vps3;

  double *dvps0     = pseudo->dvps0;
  double *dvps1     = pseudo->dvps1;
  double *dvps2     = pseudo->dvps2;
  double *dvps3     = pseudo->dvps3;

  int    *nsplin_r  = pseudo->nsplin_r;
  double *rmin_spl  = pseudo->rmin;
  double *dr_spl    = pseudo->dr_spl;
  double *x         = clatoms_pos->x;
  double *y         = clatoms_pos->y;
  double *z         = clatoms_pos->z;

  double *fx        = clatoms_pos->fx;
  double *fy        = clatoms_pos->fy;
  double *fz        = clatoms_pos->fz;

  double *fx_tmp    = ewd_scr->fx;
  double *fy_tmp    = ewd_scr->fy;
  double *fz_tmp    = ewd_scr->fz;

  double *fx_all    = ewd_scr->fx2;
  double *fy_all    = ewd_scr->fy2;
  double *fz_all    = ewd_scr->fz2;

  double *hmat      = cell->hmat_cp;
  double vol        = cell->vol_cp;

  int nkf1            = cp_para_fft_pkg->nkf1;
  int nkf2            = cp_para_fft_pkg->nkf2;
  int nkf3            = cp_para_fft_pkg->nkf3;
  int skc_fft_ka_proc = cp_para_fft_pkg->skc_fft_ka_proc;
  int ekc_fft_ka_proc = cp_para_fft_pkg->ekc_fft_ka_proc;
  int skb_fft_ka_proc = cp_para_fft_pkg->skb_fft_ka_proc;
  int ekb_fft_ka_proc = cp_para_fft_pkg->ekb_fft_ka_proc;
  int nfft_proc       = cp_para_fft_pkg->nfft_proc;
  int nfft2           = cp_para_fft_pkg->nfft/2;
  int nfft2_proc      = nfft_proc/2;

  int np_states          = communicate->np_states;
  MPI_Comm comm_states   = communicate->comm_states;

  /*long range contribution 
  double *z_1          = pseudo->z_1;
  double *z_2          = pseudo->z_2;
  double *alp_1        = pseudo->alp_1;
  double *alp_2        = pseudo->alp_2; */

/*--------------------------------------------*/
/* Local Variable Declarations */

  int i,iii,ka,kb,kc,kb_str,kb_end;
  int ind_rho,index_atm;
  int loc_opt_now;
  int is, iis,iatm,index,icount,iatm_typ_now;

  double sa,sb,sc,da,db,dc,grid_x,grid_y,grid_z;
  double dx,dy,dz,dxx,dyy,dzz,pi,r,vloc,vtemp,dvtemp,eext_loc;
  double dfact,tdfact;

/* ------------------------------------------------------------------ */
/*  0) Assign Local Variables */

  pi = M_PI;

  da = 1.0/((double) nkf1);
  db = 1.0/((double) nkf2);
  dc = 1.0/((double) nkf3);

  dxx = hmat[1]/(2.0*(double) nkf1);
  dyy = hmat[5]/(2.0*(double) nkf2);
  dzz = hmat[9]/(2.0*(double) nkf3);

  dfact   = vol/((double)nfft2);
  tdfact  = 2.0*dfact;
  eext_loc = 0.00;

  for(iatm = 1; iatm <= natm; iatm++){
    fx_tmp[iatm] = 0.00;
    fy_tmp[iatm] = 0.00;
    fz_tmp[iatm] = 0.00;
  }

/*------------------------------------------------------------------------ */
/* 1) Evaluate the local external energy and force at each grid point */

  icount = 0;

  for(kc=skc_fft_ka_proc;kc<=ekc_fft_ka_proc;kc++){
    kb_str = (kc==skc_fft_ka_proc ? skb_fft_ka_proc : 1);
    kb_end = (kc==ekc_fft_ka_proc ? ekb_fft_ka_proc : nkf2);
    for(kb=kb_str;kb<=kb_end;kb++){
      for(ka=1;ka<=nkf1;ka++){

        /* Obtain coordinate of real space grid point */
        sa = da*((double)(ka-1)) - 0.5;
        sb = db*((double)(kb-1)) - 0.5;
        sc = dc*((double)(kc-1)) - 0.5;

        grid_x  = sa*hmat[1]+sb*hmat[4]+sc*hmat[7];
        grid_y  = sa*hmat[2]+sb*hmat[5]+sc*hmat[8];
        grid_z  = sa*hmat[3]+sb*hmat[6]+sc*hmat[9];

        grid_x += dxx;
        grid_y += dyy;
        grid_z += dzz;

        ind_rho = icount + ka;

        vloc  = 0.00000;

        for(iatm = 1; iatm <= natm; iatm++){
          /* Calculate distance between atom and grid point */
          dx = grid_x - x[iatm];
          dy = grid_y - y[iatm];
          dz = grid_z - z[iatm];
          r  = sqrt(dx*dx + dy*dy + dz*dz);

          /* calculate the offset for this atom in pseudopotential array */
          iatm_typ_now = iatm_atm_typ[iatm];
          loc_opt_now  = loc_opt[iatm_typ_now];

          index_atm = 0;
          for(i=1; i< iatm_typ_now; i++){
            index_atm += (n_ang_max + 1)*nsplin_r[i];
          }
          index_atm += (loc_opt_now)*nsplin_r[iatm_typ_now];

          /* compute local pp and its derivative */
          get_vpsnow_dvr(index_atm,nsplin_r[iatm_typ_now],
                         rmin_spl[iatm_typ_now],dr_spl[iatm_typ_now],r,
                         vps0,vps1,vps2,vps3,&vtemp);

          get_dvpsnow_dvr(index_atm,nsplin_r[iatm_typ_now],
                          rmin_spl[iatm_typ_now],dr_spl[iatm_typ_now],r,
                          vps0,vps1,vps2,vps3,&dvtemp);

          /* add atom force contribution */
          fx_tmp[iatm] += (rho[icount + ka]*dvtemp*dx/r);
          fy_tmp[iatm] += (rho[icount + ka]*dvtemp*dy/r);
          fz_tmp[iatm] += (rho[icount + ka]*dvtemp*dz/r);

          /* accumulate the potential over all atoms*/
          vloc += vtemp;

        }/*endfor atoms*/

        /* Calculate the force on the coefficient */
        for(is=1; is <= nstate; is++){
          iis = (is-1)*nfft2_proc + (icount + ka);
          dvrfc[iis] -= (tdfact*dvrc[iis]*vloc);
        }

        /* accumulate the local energy over grid */
        eext_loc += rho[(icount + ka)]*vloc;

      }/*endfor grid points x*/
      icount += nkf1;
  }}/*endfor grid points y z*/

/*------------------------------------------------------------------------------*/
/* 2) Update the local energy and atom forces  */

  if(np_states > 1){
    Allreduce(&(fx_tmp[1]),&(fx_all[1]),natm,MPI_DOUBLE,MPI_SUM,0,comm_states);
    Allreduce(&(fy_tmp[1]),&(fy_all[1]),natm,MPI_DOUBLE,MPI_SUM,0,comm_states);
    Allreduce(&(fz_tmp[1]),&(fz_all[1]),natm,MPI_DOUBLE,MPI_SUM,0,comm_states);
  }else{
    for(iatm=1; iatm <= natm; iatm++){
      fx_all[iatm] = fx_tmp[iatm];
      fy_all[iatm] = fy_tmp[iatm];
      fz_all[iatm] = fz_tmp[iatm];
    }
  }

  for(iatm=1; iatm <= natm; iatm++){
    fx[iatm] += fx_all[iatm]*dfact;
    fy[iatm] += fy_all[iatm]*dfact;
    fz[iatm] += fz_all[iatm]*dfact;
  }

  *eext_ret += eext_loc*dfact;

/*======================================================================*/
}/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void eext_loc_recip_space_dvr(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                              CELL *cell,EWALD *ewald,CPEWALD *cpewald,
                              CPSCR *cpscr,PSEUDO *pseudo,EWD_SCR *ewd_scr,
                              ATOMMAPS *atommaps, double *fcre,double *fcim,
                              double *vrecip_ret,double *veext_ret,
                              COMMUNICATE *communicate,FOR_SCR *for_scr,
                              int nkf1,int nkf2,int nkf3)

/*=======================================================================*/
/*         Begin Routine                                                 */
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"

/*--------------------------------------------*/
/*         Local Pointer declarations         */
  int natm_tot      = clatoms_info->natm_tot;
  int *iatm_typ     = atommaps->iatm_atm_typ;;
  double *x         = clatoms_pos->x;
  double *y         = clatoms_pos->y;
  double *z         = clatoms_pos->z;
  double *fx        = clatoms_pos->fx;
  double *fy        = clatoms_pos->fy;
  double *fz        = clatoms_pos->fz;

  double *fx_tmp    = ewd_scr->fx;
  double *fy_tmp    = ewd_scr->fy;
  double *fz_tmp    = ewd_scr->fz;

  double *fx_all    = ewd_scr->fx2;
  double *fy_all    = ewd_scr->fy2;
  double *fz_all    = ewd_scr->fz2;

  double *vextr     = cpscr->cpscr_loc.vextr;
  double *vexti     = cpscr->cpscr_loc.vexti;

  double *q         = clatoms_info->q;
  int natm_typ      = atommaps->natm_typ;
  int *index_atm    = for_scr->index_atm;

  double *vps0      = pseudo->vps0;
  double *vps1      = pseudo->vps1;
  double *vps2      = pseudo->vps2;
  double *vps3      = pseudo->vps3;
  double *gzvps     = pseudo->gzvps;
  double *gzvps0    = pseudo->gzvps0;

  int    *nsplin_r  = pseudo->nsplin_r;
  double *rmin_spl  = pseudo->rmin;
  double *dr_spl    = pseudo->dr_spl;

  int     nsplin_g  = pseudo->nsplin_g;
  double  gmin_spl  = pseudo->gmin_spl;
  double  dg_spl    = pseudo->dg_spl;
  int n_ang_max     = pseudo->n_ang_max;
  int *loc_opt      = pseudo->loc_opt;

  int ioff_nr;
  int iatm_typ_now;

  /*----------------------*/
  /* G-vector information */
  int *kastore   = ewald->kastr;
  int *kbstore   = ewald->kbstr;
  int *kcstore   = ewald->kcstr;
  int *ibreak1   = ewald->ibrk1;
  int *ibreak2   = ewald->ibrk2;
  double *rhocr  = cpscr->cpscr_rho.rhocr_up;
  double *rhoci  = cpscr->cpscr_rho.rhoci_up;
  double *ak2    = cpewald->ak2;
  int nktot      = ewald->nktot;
  double *hmat   = cell->hmat;
  double *hmati  = cell->hmati;

  /* long range stuff 
  double *z_1          = pseudo->z_1;
  double *z_2          = pseudo->z_2;
  double *alp_1        = pseudo->alp_1;
  double *alp_2        = pseudo->alp_2;
  double alp1_2,alp2_2,falp1_2,falp2_2; */

  /*---------------------------------*/
  /* Ewald and ewald scr information */

  double alp_ewald  = ewald->alp_ewd;
  double *cossc     = ewd_scr->cossc;
  double *sinsc     = ewd_scr->sinsc;
  double *helr      = ewd_scr->helr;
  double *heli      = ewd_scr->heli;
  double *vtemp     = ewd_scr->temp;

  double *ewd_scr_x = ewd_scr->x;
  double *ewd_scr_y = ewd_scr->y;
  double *ewd_scr_z = ewd_scr->z;

  /*---------------------------*/
  /* Communication information */

  int myid_state       = communicate->myid_state;
  int np_states        = communicate->np_states;
  MPI_Comm comm_states = communicate->comm_states;

  /*---------------------------*/
  /* Local Variables           */

  int istart,ngo,irem,idiv;
  int ipart,jpart,iii,itype,i;
  int icount,koff;

  double falp2,vol,rvol,pivol,fpi,arg,q_sum1;
  double aka,akb,akc,xk,yk,zk,atemp,btemp,ctemp;
  double xtemp,ytemp,ztemp;
  double sumr,sumi,g,g2,preg,prep,tpi,pi;
  double sumr_h,sumi_h;
  double srx,sry,srz,six,siy,siz,temp,smag;
  double vrecip,veext;
  int nfft2_proc;

/*======================================================================*/
/* I) Get some useful constants                                         */

  pi    = M_PI;
  tpi   = 2.0*pi;
  fpi   = 4.0*pi;

  vol   = getdeth(hmat);

  rvol  = 1.0/vol;
  pivol = vol/2.0/pi; /* 4.0 changed to 2.0 DY FULL SPACE NOT HALF SPACE */
  falp2 = 4.0*alp_ewald*alp_ewald;

/*======================================================================*/
/* II) Find cos and sin of exp(igR_I) components of the particles       */
/*    ( hmati rvec = svec   r=(x,y,z) s=(a,b,c) )                       */

  for(ipart=1;ipart<=natm_tot;ipart++){
    fx_tmp[ipart] = 0.0;
    fy_tmp[ipart] = 0.0;
    fz_tmp[ipart] = 0.0;
  }/*endfor*/

  for(ipart=1;ipart<=natm_tot;ipart++){
    xtemp = x[ipart];
    ytemp = y[ipart];
    ztemp = z[ipart];

    ewd_scr_x[ipart] = xtemp*hmati[1] + ytemp*hmati[4] + ztemp*hmati[7];
    ewd_scr_y[ipart] = xtemp*hmati[2] + ytemp*hmati[5] + ztemp*hmati[8];
    ewd_scr_z[ipart] = xtemp*hmati[3] + ytemp*hmati[6] + ztemp*hmati[9];

    ctemp            = ewd_scr_z[ipart]*tpi;
    cossc[ipart]     = cos(ctemp);
    sinsc[ipart]     = sin(ctemp);
  }/*endfor*/

/*======================================================================*/
/* II) Perform the ewald sum/ CP-potential calculation */

  vrecip  = 0.0;
  veext   = 0.0;

  idiv    = (nkf2*nkf3)/np_states;
  irem    = (nkf2*nkf3)%np_states;

  nfft2_proc = ( myid_state < irem ? (idiv+1)*nkf1 : idiv*nkf1);

  ngo = nfft2_proc;
  istart = (myid_state <= irem ? myid_state*(idiv+1)*nkf1 + 1
                               : irem*(idiv+1)*nkf1+(myid_state-irem)*idiv*nkf1+1);

  koff = istart-1;

  if(np_states==myid_state+1){ngo--;} /*for PE with g=0 component*/

  /* GO OVER g-VECTORS */
  for(icount=1;icount<= ngo;icount++){

   /* (a) Get the k vectors */
    aka = (double)(kastore[(icount+koff)]);
    akb = (double)(kbstore[(icount+koff)]);
    akc = (double)(kcstore[(icount+koff)]);

    xk = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
    yk = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
    zk = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;

    g2 = xk*xk+yk*yk+zk*zk;
    g  = sqrt(g2);
    ak2[icount] = g2;

    /* (b) compute the useful vectors if break point number one */
    /*     or you are just starting out*/

    for(ipart=1;ipart<=natm_tot;ipart++){
      atemp = ewd_scr_x[ipart];
      btemp = ewd_scr_y[ipart];
      ctemp = ewd_scr_z[ipart];

      arg = (aka*atemp + akb*btemp + akc*ctemp)*tpi;
      helr[ipart] = cos(arg);
      heli[ipart] = sin(arg);
    }

    /* (c) Calculate the offsets in PP array for all atoms*/

    for(itype=1; itype <= natm_typ; itype++){
      ioff_nr = 0;
      for(i=1; i< itype; i++){ 
        ioff_nr += (n_ang_max+1)*nsplin_r[i];
      }
      ioff_nr += loc_opt[itype]*nsplin_r[itype];
      index_atm[itype] = ioff_nr;
    }/*endfor*/

    /* (d) Calculate the local external potential at this g
           for all atoms */

    get_vpsnow(index_atm,nsplin_g,gmin_spl,dg_spl,g,
               vps0,vps1,vps2,vps3,vtemp,iatm_typ,natm_typ,natm_tot,1);

    for(ipart=1;ipart<=natm_tot;ipart++){
      vtemp[ipart] *= rvol;
    }

    /* (e) This is for force: */
    vextr[icount]  =  ddot1(natm_tot,helr,1,vtemp,1);
    vexti[icount]  = -ddot1(natm_tot,heli,1,vtemp,1);
    fcre[icount]   = vextr[icount];
    fcim[icount]   = vexti[icount];

    /* (f) accumulate the energy */
    veext += (vextr[icount]*rhocr[icount] + vexti[icount]*rhoci[icount]);

    /* (g) Get the real and imag parts of the structure factor  */
    sumr = ddot1(natm_tot,helr,1,q,1);
    sumi = ddot1(natm_tot,heli,1,q,1);
    smag = sumr*sumr+sumi*sumi;

    /* (h) accumulate ewald energy */
    preg    = exp(-g2/falp2)/(g2*pivol);
    vrecip  = vrecip + smag*preg;

    /* (i) get the atom force */
    sumr_h = sumr;
    sumi_h = sumi;
    sumr = sumr*preg*2.0;
    sumi = sumi*preg*2.0;
    for(ipart=1;ipart<=natm_tot;ipart++){
      srx = xk*(sumr*q[ipart] + rhocr[icount]*vtemp[ipart]);
      sry = yk*(sumr*q[ipart] + rhocr[icount]*vtemp[ipart]);
      srz = zk*(sumr*q[ipart] + rhocr[icount]*vtemp[ipart]);
      six = xk*(sumi*q[ipart] - rhoci[icount]*vtemp[ipart]);
      siy = yk*(sumi*q[ipart] - rhoci[icount]*vtemp[ipart]);
      siz = zk*(sumi*q[ipart] - rhoci[icount]*vtemp[ipart]);

      fx_tmp[ipart] += (srx*heli[ipart]  - six*helr[ipart]);
      fy_tmp[ipart] += (sry*heli[ipart]  - siy*helr[ipart]);
      fz_tmp[ipart] += (srz*heli[ipart]  - siz*helr[ipart]);
    }/*endfor*/

  }/*endfor:icount loop over k vectors */

/*======================================================================*/
/* III) Taking care of g=0 term (local pseudopotential)*/

  if((myid_state+1)==np_states){
    icount = ngo+1;

    for(ipart=1;ipart<=natm_tot;ipart++){
      iatm_typ_now  = iatm_typ[ipart];
      vtemp[ipart]  = gzvps[iatm_typ_now];
    }/*endfor*/

    ak2[(ngo+1)]   = 0.0;
    vextr[(ngo+1)] = dsum1(natm_tot,vtemp,1)*rvol;
    vexti[(ngo+1)] = 0.0;

    fcre[icount]   = vextr[(ngo+1)];
    fcim[icount]   = vexti[(ngo+1)];

    veext +=  vextr[(ngo+1)]*rhocr[icount];
  }/*endif*/

/*======================================================================*/
/* IV) Collect the forces */

  if(np_states > 1){
   Allreduce(&(fx_tmp[1]),&(fx_all[1]),natm_tot,MPI_DOUBLE,
             MPI_SUM,0,comm_states);
   Allreduce(&(fy_tmp[1]),&(fy_all[1]),natm_tot,MPI_DOUBLE,
             MPI_SUM,0,comm_states);
   Allreduce(&(fz_tmp[1]),&(fz_all[1]),natm_tot,MPI_DOUBLE,
             MPI_SUM,0,comm_states);
  }else{
    for(ipart=1; ipart <= natm_tot; ipart++){
      fx_all[ipart] = fx_tmp[ipart];
      fy_all[ipart] = fy_tmp[ipart];
      fz_all[ipart] = fz_tmp[ipart];
    }
  }

  for(ipart=1;ipart<=natm_tot;ipart++){
    fx[ipart] += fx_all[ipart];
    fy[ipart] += fy_all[ipart];
    fz[ipart] += fz_all[ipart];
  }

  if(np_states > 1) Barrier(comm_states);

/*======================================================================*/
/* V) Finally, store the final value of vrecip */

   *vrecip_ret  = vrecip;
   *veext_ret  += veext;

/*======================================================================*/
}/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vps_atm_list_dvr(PSEUDO *pseudo,CELL *cell, PARA_FFT_PKG3D *cp_para_fft_pkg,
                      CPSCR_NONLOC *cpscr_nonloc,ATOMMAPS *atommaps,int natm_tot,
                      int nstate)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/

/*         Local Variable declarations                                   */

  int iatm,iang,ishift,loc_now,iii,i,j;
  int icount,ioff,isame;
  int ishift2;
  int np_nonloc_cp_box;
  int np_nonloc_cp_box_kb;
  int count,count_np_nonloc_cp_box_kb;
  int *inonloc_index;
  int iatm_typ_now,natm_typ_nl;

 /* Local pointers */
 /*----------------*/

  int *iatm_typ         = atommaps->iatm_atm_typ;
  int *ip_nl            = pseudo->ip_nl;
  int *np_nl            = pseudo->np_nl;/* Number of particles in l channel*/
  int *ip_nl_rev        = pseudo->ip_nl_rev;     /*length natm*nang_max*/

  /* list: index from full atom array into nonlocal cp atom array*/
  int *ivps_label       = pseudo->ivps_label;
  int *n_ang            = pseudo->n_ang;
  int *loc_opt          = pseudo->loc_opt;
  int n_ang_max         = pseudo->n_ang_max;

  int *iatm_typ_nl,*iatm_typ_nl_rev,*imap_atm_typ_nl;

 /*for nonlocal truncation */

  double dxx,dyy,dzz;
  double *hmat             = cell->hmat;
  double *hmat_cp          = cell->hmat_cp;

  int myid                 = cp_para_fft_pkg->myid;
  int nkf1                 = cp_para_fft_pkg->nkf1;
  int nkf2                 = cp_para_fft_pkg->nkf2;
  int nkf3                 = cp_para_fft_pkg->nkf3;
  int skc_fft_ka_proc      = cp_para_fft_pkg->skc_fft_ka_proc;
  int ekc_fft_ka_proc      = cp_para_fft_pkg->ekc_fft_ka_proc;
  int skb_fft_ka_proc      = cp_para_fft_pkg->skb_fft_ka_proc;
  int ekb_fft_ka_proc      = cp_para_fft_pkg->ekb_fft_ka_proc;

  int nmax_wan_orb  = cpscr_nonloc->nmax_wan_orb;
  int nloc_trunc_on = cpscr_nonloc->nloc_trunc_on;
  int nloc_wan_on   = cpscr_nonloc->nloc_wan_on;
  int nrmax_nl,ngrid;
  double *rcut_nl = pseudo->rcut_nl;
  double rcut_max,memory_nl;

/*=======================================================================*/
/* I) Make a list of the local/non-local atoms  (1 GRID ONLY)            */

  np_nonloc_cp_box    = 0;
  np_nonloc_cp_box_kb = 0;

  inonloc_index = (int *) cmalloc(natm_tot*sizeof(int))-1;

  /* Count up the number of atoms with KB non-local pseudo potentials */

  icount = 0;
  for(iatm=1;iatm<=natm_tot;iatm++){
    iatm_typ_now = iatm_typ[iatm];
    if(ivps_label[iatm_typ_now]==1){  /* Kleinman-Bylander ONLY */
      loc_now = loc_opt[iatm_typ_now];
      if(loc_now > 0){
        np_nonloc_cp_box_kb++;
        icount++;
        inonloc_index[icount] = iatm;
      }
    }
  }

  /* only KB-type nonlocal is available for DVR */

  pseudo->np_nonloc_cp_box    = np_nonloc_cp_box_kb;
  pseudo->np_nonloc_cp_box_kb = np_nonloc_cp_box_kb;
  pseudo->np_nonloc_cp_box_gh = 0;

/*==================================================================*/
/* II) Assign number of KB atoms in each angular momentum channel       */

  for(iang=1;iang<=(n_ang_max+1);iang++){
    np_nl[iang] = 0;
  }

  count_np_nonloc_cp_box_kb = 0;

  for(iatm=1;iatm<=natm_tot;iatm++){
    iatm_typ_now = iatm_typ[iatm];

    if(ivps_label[iatm_typ_now]==1){ /* KLEINMAN-BYLANDER */

      loc_now = loc_opt[iatm_typ_now];
      if(loc_now > 0 ){count_np_nonloc_cp_box_kb++;}

      for(iang=1;iang<=loc_now;iang++){
        np_nl[iang]  += 1;
        ishift        = (iang - 1)*natm_tot + np_nl[iang];
        ip_nl[ishift] = iatm;
        ip_nl_rev[ishift] = count_np_nonloc_cp_box_kb;
      }

      for(iang=loc_now+2;iang<=(n_ang[iatm_typ_now]+1);iang++){
        np_nl[iang]  += 1;
        ishift        = (iang - 1)*natm_tot + np_nl[iang];
        ip_nl[ishift] = iatm;
	ip_nl_rev[ishift] = count_np_nonloc_cp_box_kb;
      }

    }/*endif: atm has a KB pseudopotential*/
  }/*endfor*/


/*--------------------------------------------------------------------------*/
/* III) assign global and NL atom type */

  atommaps->iatm_atm_typ_nl = (int *) cmalloc(np_nonloc_cp_box_kb*sizeof(int))-1;
  atommaps->iatm_atm_typ_nl_rev = (int *) cmalloc(np_nonloc_cp_box_kb*sizeof(int))-1;

  iatm_typ_nl     = atommaps->iatm_atm_typ_nl;
  iatm_typ_nl_rev = atommaps->iatm_atm_typ_nl_rev;

  /*Nonlocal atom's  global atom type  */
  for(i=1; i<= np_nonloc_cp_box_kb; i++){
    iatm = inonloc_index[i];
    iatm_typ_nl[i] = iatm_typ[iatm];
  }/*endfor*/

  /*Count the number of non-local atom types*/
  natm_typ_nl = 0;

  for(i=1; i<= np_nonloc_cp_box_kb; i++){
    isame = 0;
    for(j=(i-1); j >= 1; j--){
      if(iatm_typ_nl[j] == iatm_typ_nl[i]){
        isame = 1;
        break;
      }/*endif*/
    }/*endfor*/
    if(isame == 0){natm_typ_nl++;}/*endif*/
  }/*endfor*/

  pseudo->natm_typ_nl = natm_typ_nl;

  /* Assign atom's nonlocal atom type */
  if( np_nonloc_cp_box_kb > 0){
    iatm_typ_nl_rev[1] = 1;
    icount = 1;
  }/*endif*/
  for(i=2; i<= np_nonloc_cp_box_kb; i++){
    for(j=(i-1); j>= 1;j--){
      isame = 0;
      if(iatm_typ_nl[i] == iatm_typ_nl[j]){
        isame = 1;
        break;
      }/*endif*/
    }/*endfor*/
    if(isame==0){/* not the same*/
      icount++;
      iatm_typ_nl_rev[i] = icount;
    }else{/*the same*/
      iatm_typ_nl_rev[i] = iatm_typ_nl_rev[j];
    }
  }/*endfor*/

/*-------------------------------------------------------------------------*/
/* IV) Create a map: nonlocal atom type j is of global atom type k */

  atommaps->imap_atm_typ_nl = (int *) cmalloc(natm_typ_nl*sizeof(int))-1;
  imap_atm_typ_nl  = atommaps->imap_atm_typ_nl;

  if( np_nonloc_cp_box_kb > 0){
    icount = 1;
    imap_atm_typ_nl[1] = iatm_typ_nl[1];
  }

  for(i=2; i<= np_nonloc_cp_box_kb; i++){
    isame = 0;
    for(j=(i-1); j >= 1; j--){
      if(iatm_typ_nl_rev[i] == iatm_typ_nl_rev[j]){
        isame = 1;
        break;
      }/*endif*/
    }/*endfor*/
    if(isame == 0){
      icount++;
      imap_atm_typ_nl[icount] = iatm_typ_nl[i];
    }/*endif*/
  }/*endfor*/

/*-----------------------------------------------------------------------------*/
/* V) allocate more memory for the nonlocal truncation                         */

  if(nloc_trunc_on==1 || nloc_wan_on==1){

    ngrid = (ekc_fft_ka_proc - skc_fft_ka_proc - 1)*nkf1*nkf2
            + (nkf2 - skb_fft_ka_proc + 1)*nkf1 + (ekb_fft_ka_proc)*nkf1;

    rcut_max=0.0;
    for(i=1;i<= np_nonloc_cp_box_kb;i++){
      iatm = ip_nl[i];
      iatm_typ_now = iatm_typ[iatm];
      if(rcut_max < rcut_nl[iatm_typ_now]) {
        rcut_max = rcut_nl[iatm_typ_now];
      }
    }

    dxx = hmat[1]/(2.0*(double) nkf1);
    dyy = hmat[5]/(2.0*(double) nkf2);
    dzz = hmat[9]/(2.0*(double) nkf3);

    cpscr_nonloc->nrmax_nl = (int)((rcut_max/dxx)*(rcut_max/dyy)*(rcut_max/dzz)*0.53);

    nrmax_nl = cpscr_nonloc->nrmax_nl;
    if(nloc_trunc_on==0 || nrmax_nl > ngrid || nrmax_nl < 0) {
      cpscr_nonloc->nrmax_nl = ngrid;
    }
    nrmax_nl = cpscr_nonloc->nrmax_nl;

    cpscr_nonloc->dvrc_wan = (double *) cmalloc(nmax_wan_orb*nrmax_nl*sizeof(double))-1;
    cpscr_nonloc->igrid_map = (int *) cmalloc(nrmax_nl*sizeof(int))-1;
    if(nloc_wan_on==1){
      cpscr_nonloc->my_wan_orb = (int *) cmalloc(nstate*np_nonloc_cp_box_kb*sizeof(int))-1;
      cpscr_nonloc->num_wan_orb = (int *) cmalloc(np_nonloc_cp_box_kb*sizeof(int))-1;
      cpscr_nonloc->cnum_wan_orb = (int *) cmalloc((np_nonloc_cp_box_kb+1)*sizeof(int))-1;
    }

    memory_nl = nmax_wan_orb*nrmax_nl*sizeof(double)*(1.0e-6);
    printf("ADDITIONAL MEMORY ALLOCATION FOR NON-LOCAL: %g MB \n",memory_nl);
  }/* endif nloc_trunc or nloc_wan */

  cfree(&(inonloc_index[1]));

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void get_vpsnow_dvr(int index_atm,int nsplin_r, double rmin_spl,
                    double dr_spl,double r,double *vps0,double *vps1,
                    double *vps2,double *vps3,double *vtemp)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*        Local variable declarations            */

#include "../typ_defs/typ_mask.h"

  int ipart,itype,iii;
  int index_now;
  double vtemp_atyp[200];
  double h,h0;
  double partem1,partem2,partem3,partem4;

/*==========================================================================*/
/* Loop over atom types to calculate pseudopotential                        */

    iii = (r-rmin_spl)/dr_spl + 1 ;
    iii = MIN(iii,nsplin_r);
    iii = MAX(iii,1);

    h0  = (double)(iii-1)*dr_spl + rmin_spl;
    h   = r-h0;
    index_now = index_atm + iii;

    partem1 = vps0[index_now];
    partem2 = vps1[index_now];
    partem3 = vps2[index_now];
    partem4 = vps3[index_now];

    *vtemp = ((partem4*h+partem3)*h+partem2)*h+ partem1;

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

void get_dvpsnow_dvr(int index_atm,int nsplin_r,double rmin_spl,
                     double dr_spl,double r,double *vps0,double *vps1,
                     double *vps2,double *vps3,double *dvtemp)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*        Local variable declarations            */

#include "../typ_defs/typ_mask.h"

  int ipart,itype,iii;
  int index_now;
  double vtemp_atyp[200];
  double h,h0;
  double partem1,partem2,partem3,partem4;

/*==========================================================================*/
/* Loop over atom types to calculate pseudopotential                        */

    iii = (r-rmin_spl)/dr_spl + 1 ;
    iii = MIN(iii,nsplin_r);
    iii = MAX(iii,1);

    h0  = (double)(iii-1)*dr_spl + rmin_spl;
    h   = r-h0;
    index_now = index_atm + iii;

    partem1 = vps0[index_now];
    partem2 = vps1[index_now];
    partem3 = vps2[index_now];
    partem4 = vps3[index_now];

    *dvtemp = (3.0*partem4*h + 2.0*partem3)*h+ partem2 ;

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_eext_nonloc_dvr(ATOMMAPS *atommaps,PSEUDO *pseudo,
                             CLATOMS_POS *clatoms_pos, CPSCR_NONLOC *cpscr_nonloc,
                             CPSCR_WAVE *cpscr_wave, EWD_SCR *ewd_scr,
                             CPOPTS *cpopts, CPCOEFFS_INFO *cpcoeffs_info,
                             CPCOEFFS_POS_DVR *cpcoeffs_pos_dvr,CELL *cell,
                             COMMUNICATE *communicate,
                             PARA_FFT_PKG3D *cp_para_fft_pkg, int natm_tot,
                             double *enl,int itime_nl, double **wan_cent,
                             double r_cut_wan)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*        Local Pointer declarations             */

  int cp_lsda     = cpopts->cp_lsda;
  int nstate_up   = cpcoeffs_info->nstate_up;
  int nstate_dn   = cpcoeffs_info->nstate_dn;

  int nloc_trunc_on =   cpscr_nonloc->nloc_trunc_on;
  int nloc_wan_on   =   cpscr_nonloc->nloc_wan_on;

  int *np_nl      = pseudo->np_nl;
  int n_ang_max   = pseudo->n_ang_max_kb;

  double *dvrc_up  = cpcoeffs_pos_dvr->dvrc_up;
  double *dvrc_dn  = cpcoeffs_pos_dvr->dvrc_dn;
  double *dvrfc_up = cpcoeffs_pos_dvr->dvrfc_up;
  double *dvrfc_dn = cpcoeffs_pos_dvr->dvrfc_dn;

  int myid_state       = communicate->myid_state;

  /* Local Variable declarations */

  int i;
  int nl_max,np_nlmax_kb;
  double cp_enl,cp_enl_ret;

/*======================================================================*/
/* I) determin some maximum values */

  /*  maximum open non-local angular momentum channel*/

  nl_max = -1;
  for(i=1;i<=(n_ang_max +1);i++){
    if(np_nl[i]>0){nl_max=i-1;}
  }/*endfor*/

  /* maximum number of atoms in any open angular momentum channel  */

  np_nlmax_kb = 1;
  for(i = 1;i<=(nl_max + 1);i++){
   np_nlmax_kb = MAX(np_nlmax_kb,np_nl[i]);
  }/*endfor*/


/*----------------------------------------------------------------------*/
/* II) compute the nonlocal energy/force */

  cp_enl = cp_enl_ret = 0.00000;

/*
  if(nloc_trunc_on==1 && nloc_wan_on==1){
    eext_nonloc_dvr_wannier(atommaps,pseudo,clatoms_pos,cpscr_nonloc,
                            cpscr_wave,ewd_scr,cell,cp_para_fft_pkg,
                            communicate,dvrc_up,dvrfc_up,
                            nstate_up,natm_tot,
                            nl_max,np_nlmax_kb,&cp_enl_ret,wan_cent,r_cut_wan);

  }else if (nloc_trunc_on==1 && nloc_wan_on==0){
    eext_nonloc_dvr_trunc(atommaps,pseudo,clatoms_pos,cpscr_nonloc,
                          cpscr_wave,ewd_scr,cell,cp_para_fft_pkg,
                          communicate,dvrc_up,dvrfc_up,
                          nstate_up,natm_tot,
                          nl_max,np_nlmax_kb,&cp_enl_ret);
  }else{
*/
    eext_nonloc_dvr_mm(atommaps,pseudo,clatoms_pos,cpscr_nonloc,
                       cpscr_wave,ewd_scr,cell,cp_para_fft_pkg,
                       communicate,dvrc_up,dvrfc_up,nstate_up,
                       natm_tot,nl_max,np_nlmax_kb,&cp_enl_ret,
                       (cpscr_nonloc->vnlre_up), (cpscr_nonloc->vnlim_up),
                       (cpscr_nonloc->dvnlre_x_up), (cpscr_nonloc->dvnlim_x_up),
                       (cpscr_nonloc->dvnlre_y_up), (cpscr_nonloc->dvnlim_y_up),
                       (cpscr_nonloc->dvnlre_z_up), (cpscr_nonloc->dvnlim_z_up));

    cp_enl += cp_enl_ret;

    if(cp_lsda==1){
      eext_nonloc_dvr_mm(atommaps,pseudo,clatoms_pos,cpscr_nonloc,
                         cpscr_wave,ewd_scr,cell,cp_para_fft_pkg,
                         communicate,dvrc_dn,dvrfc_dn,nstate_dn,
                         natm_tot,nl_max,np_nlmax_kb,&cp_enl_ret,
                         (cpscr_nonloc->vnlre_dn), (cpscr_nonloc->vnlim_dn),
                         (cpscr_nonloc->dvnlre_x_dn), (cpscr_nonloc->dvnlim_x_dn),
                         (cpscr_nonloc->dvnlre_y_dn), (cpscr_nonloc->dvnlim_y_dn),
                         (cpscr_nonloc->dvnlre_z_dn), (cpscr_nonloc->dvnlim_z_dn));

      cp_enl += cp_enl_ret;
    }

/*
  }
*/

  *enl = cp_enl/(double)(communicate->np_states);

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void eext_nonloc_dvr_mm(ATOMMAPS *atommaps,PSEUDO *pseudo, 
                        CLATOMS_POS *clatoms_pos, CPSCR_NONLOC *cpscr_nonloc,
                        CPSCR_WAVE *cpscr_wave, EWD_SCR *ewd_scr,
                        CELL *cell, PARA_FFT_PKG3D *cp_para_fft_pkg,
                        COMMUNICATE *communicate, double *dvrc, double *dvrfc,
                        int nstate, int natm_tot, int nl_max, int np_nlmax_kb,
                        double *pcp_enl, double *vnlreal, double *vnlimag,
                        double *dvnlreal_x, double *dvnlimag_x, double *dvnlreal_y,
                        double *dvnlimag_y, double *dvnlreal_z, double *dvnlimag_z)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*        Local Pointer declarations             */

#include "../typ_defs/typ_mask.h"

  int *iatm_typ            = atommaps->iatm_atm_typ;
  int *iatm_typ_nl         = atommaps->iatm_atm_typ_nl;
  int natm_typ             = atommaps->natm_typ;

  int np_nonloc            = pseudo->np_nonloc_cp_box_kb;
  int *np_nl               = pseudo->np_nl;
  int *ip_nl               = pseudo->ip_nl;
  int *ip_nl_rev           = pseudo->ip_nl_rev;
  int n_ang_max            = pseudo->n_ang_max_kb;
  double *vpsnorm          = pseudo->vpsnorm;
  int np_nonloc_cp_box_kb  = pseudo->np_nonloc_cp_box_kb;

  int *nsplin_r            = pseudo->nsplin_r;
  double *rmin_spl         = pseudo->rmin;
  double *dr_spl           = pseudo->dr_spl;
  double *vps0             = pseudo->vps0;
  double *vps1             = pseudo->vps1;
  double *vps2             = pseudo->vps2;
  double *vps3             = pseudo->vps3;
  double *dvps0            = pseudo->dvps0;
  double *dvps1            = pseudo->dvps1;
  double *dvps2            = pseudo->dvps2;
  double *dvps3            = pseudo->dvps3;

  double *x                = clatoms_pos->x;
  double *y                = clatoms_pos->y;
  double *z                = clatoms_pos->z;
  double *fx               = clatoms_pos->fx;
  double *fy               = clatoms_pos->fy;
  double *fz               = clatoms_pos->fz;

  double *vscr             = ewd_scr->fx2;
  double *vnorm            = ewd_scr->fy2;
  double *vnorm_now        = ewd_scr->fz2;

  double *fx_tmp    = ewd_scr->fx;
  double *fy_tmp    = ewd_scr->fy;
  double *fz_tmp    = ewd_scr->fz;

  /* need to define and allocate vnlreal_tmp, dvnlreal_x_tmp etc.. */

  double *vnlreal_tmp      = cpscr_nonloc->vnlre_tmp;
  double *vnlimag_tmp      = cpscr_nonloc->vnlim_tmp;

  double *dvnlreal_x_tmp   = cpscr_nonloc->dvnlre_x_tmp;
  double *dvnlreal_y_tmp   = cpscr_nonloc->dvnlre_y_tmp;
  double *dvnlreal_z_tmp   = cpscr_nonloc->dvnlre_z_tmp;
  double *dvnlimag_x_tmp   = cpscr_nonloc->dvnlim_x_tmp;
  double *dvnlimag_y_tmp   = cpscr_nonloc->dvnlim_y_tmp;
  double *dvnlimag_z_tmp   = cpscr_nonloc->dvnlim_z_tmp;

  double *hmat             = cell->hmat;
  double *hmat_cp          = cell->hmat_cp;
  int    iperd             = cell->iperd;

  int nkf1                 = cp_para_fft_pkg->nkf1;
  int nkf2                 = cp_para_fft_pkg->nkf2;
  int nkf3                 = cp_para_fft_pkg->nkf3;
  int skc_fft_ka_proc      = cp_para_fft_pkg->skc_fft_ka_proc;
  int ekc_fft_ka_proc      = cp_para_fft_pkg->ekc_fft_ka_proc;
  int skb_fft_ka_proc      = cp_para_fft_pkg->skb_fft_ka_proc;
  int ekb_fft_ka_proc      = cp_para_fft_pkg->ekb_fft_ka_proc;
  int nfft_proc            = cp_para_fft_pkg->nfft_proc;
  int nfft2                = cp_para_fft_pkg->nfft/2;
  int nfft2_proc           = nfft_proc/2;

  int np_states            = communicate->np_states;
  int myid_state           = communicate->myid_state; 
  MPI_Comm comm_states     = communicate->comm_states;

  int ioffset  = (nfft2_proc)*(nstate)/2;

  double *vnlreal_tmp_all     = cpscr_wave->zfft;
  double *vnlimag_tmp_all     = cpscr_wave->zfft+ioffset;

  double *dvnlreal_x_tmp_all  = cpscr_wave->zfft_tmp;
  double *dvnlimag_x_tmp_all  = cpscr_wave->zfft_tmp + ioffset;
  double *dvnlreal_y_tmp_all  = cpscr_wave->cre_up;
  double *dvnlimag_y_tmp_all  = cpscr_wave->cre_up+ioffset;
  double *dvnlreal_z_tmp_all  = cpscr_wave->cim_up;
  double *dvnlimag_z_tmp_all  = cpscr_wave->cim_up+ioffset; 

  /* This is not an error. Just for convenience */
  double *freal_tmp         = cpscr_wave->cre_up;
  double *fimag_tmp         = cpscr_wave->cim_up;

/*        Local Variable declarations             */

  int ioff_nr,ind_nr,ioff_ang,ioff;
  int i,iii,is,m,iatm,iang,im,ind_atm,iatm_typ_now,ntot,nlm_tot,icount,ind_lm;
  int lp1,i_shift;
  int ka,kb,kc,kb_str,kb_end;

  int index;

  double sa,sb,sc,da,db,dc, grid_x,grid_y,grid_z,dxx,dyy,dzz;
  double dx,dy,dz,r, vtemp,dvtemp,cp_enl;
  double rvol_cp,vol_cp,tdfact_sq,dfact;

  YLM_CONS ylm_cons;
  double fpi;
  double ylmr[21], ylmi[21];
  double dylmr_x[21], dylmi_x[21];
  double dylmr_y[21], dylmi_y[21];
  double dylmr_z[21], dylmi_z[21];

  double rinv,ylmr_loc,ylmi_loc;

  int itransp  = 0;
  int inormal  = 1;
  double beta  = 1.0;
  double alpha = 1.0;
  int ldx,ldy,ldz,nrz,ncz,nxy,lmi,ngrid;

/*======================================================================*/
/* I) Set constants for spherical harmonics                             */

  fpi                 = 4.0*M_PI;
  ylm_cons.rt_fpi     = 1.0/sqrt(fpi);
  ylm_cons.rt_thrfpi  = sqrt(3.0/fpi);
  ylm_cons.rt_threpi  = sqrt(1.50/fpi);
  ylm_cons.hrt_fivfpi = 0.50*sqrt(5.0/fpi);
  ylm_cons.rt_fiftepi = sqrt(7.50/fpi);
  ylm_cons.hrt_sevfpi = 0.50*sqrt(7.0/fpi);
  ylm_cons.hrt_toepi  = 0.50*sqrt(10.50/fpi)/sqrt(2.0);
  ylm_cons.hrt_ohffpi = 0.50*sqrt(105.0/fpi)/sqrt(2.0);
  ylm_cons.hrt_tfepi  = 0.50*sqrt(17.50/fpi)/sqrt(2.0);

  ylmr[1] =  ylm_cons.rt_fpi;
  ylmi[1] =  0.000000000000000;

  dylmr_x[1] = 0.00;
  dylmi_x[1] = 0.00;
  dylmr_y[1] = 0.00;
  dylmi_y[1] = 0.00;
  dylmr_z[1] = 0.00;
  dylmi_z[1] = 0.00;

/*=====================================================================*/
/* II) Cell/grid information */

  vol_cp  = getdeth(hmat_cp);
  rvol_cp = 1.0/vol_cp;

  dfact  = vol_cp/((double)nfft2);
  tdfact_sq = 2.0*dfact*dfact;

  dxx = hmat[1]/(2.0*(double) nkf1);
  dyy = hmat[5]/(2.0*(double) nkf2);
  dzz = hmat[9]/(2.0*(double) nkf3);

  da = 1.0/((double) nkf1);
  db = 1.0/((double) nkf2);
  dc = 1.0/((double) nkf3);

/*======================================================================*/
/* III) Zero the non-local tensors                                       */

  nlm_tot = (n_ang_max + 1)*(n_ang_max + 1);

  ntot = nlm_tot*nstate*np_nlmax_kb;

  for(i=1;i<=ntot;i++){
    vnlreal_tmp[i]     = 0.0;
    vnlimag_tmp[i]     = 0.0;
    vnlreal[i]         = 0.0;
    vnlimag[i]         = 0.0;
    dvnlreal_x[i]      = 0.0;
    dvnlreal_y[i]      = 0.0;
    dvnlreal_z[i]      = 0.0;
    dvnlimag_x[i]      = 0.0;
    dvnlimag_y[i]      = 0.0;
    dvnlimag_z[i]      = 0.0;
    dvnlreal_x_tmp[i]  = 0.0;
    dvnlreal_y_tmp[i]  = 0.0;
    dvnlreal_z_tmp[i]  = 0.0;
    dvnlimag_x_tmp[i]  = 0.0;
    dvnlimag_y_tmp[i]  = 0.0;
    dvnlimag_z_tmp[i]  = 0.0;
  }/*endfor*/

/************************************************************************/
/* Create Zilmn                                                         */
/* Zilmn = sum_atoms sum_grid_points Psi_n*Vkb_l*Phi_kb_l*Ylm(omega_r)  */
/************************************************************************/

/*=======================================================================*/
/* IV) Compute pseudopotential and its derivative */

  lmi = 0;
  icount = 1;
  for(iang=0; iang <= nl_max ; iang++){
    ioff     = natm_tot*iang;

    for(im=0; im < (2*iang + 1) ; im++){
      ind_lm = (iang*iang) + im + 1;

      /* loop over atoms with non-local pseudo component in l channel */
      for(iatm = 1; iatm <= np_nl[(iang+1)]; iatm++){
        lmi++;
        ind_atm  = ip_nl[ioff + iatm];
        iatm_typ_now = iatm_typ[ind_atm];

        /* loop over grid points */
        for(kc=skc_fft_ka_proc;kc<=ekc_fft_ka_proc;kc++){
          kb_str = (kc==skc_fft_ka_proc ? skb_fft_ka_proc : 1);
          kb_end = (kc==ekc_fft_ka_proc ? ekb_fft_ka_proc : nkf2);
          sc = dc*((double)(kc-1)) - 0.5;

          for(kb=kb_str;kb<=kb_end;kb++){
            sb = db*((double)(kb-1)) - 0.5;

            for(ka=1;ka<=nkf1;ka++){
              sa = da*((double)(ka-1)) - 0.5;

              grid_x  = sa*hmat[1]+sb*hmat[4]+sc*hmat[7];
              grid_y  = sa*hmat[2]+sb*hmat[5]+sc*hmat[8];
              grid_z  = sa*hmat[3]+sb*hmat[6]+sc*hmat[9];

              grid_x += dxx;
              grid_y += dyy;
              grid_z += dzz;

              dx = grid_x - x[ind_atm];
              dy = grid_y - y[ind_atm];
              dz = grid_z - z[ind_atm];

              if(iperd==3) {period_one(1,&dx,&dy,&dz,cell);}

              r     = sqrt(dx*dx + dy*dy + dz*dz);
              rinv  = 1.0/r;

              /* calculate offset into pseudo array for spline lookup */
              ioff_nr = 0;
              for(i=1; i< iatm_typ_now; i++){ ioff_nr += nsplin_r[i];}

              ioff_nr *= (n_ang_max + 1);
              ind_nr   = ioff_nr + iang*nsplin_r[iatm_typ_now];

              /* Look up KB potential Dv_l*Phi_l at grid point r */
              get_vpsnow_dvr(ind_nr,nsplin_r[iatm_typ_now],
                             rmin_spl[iatm_typ_now],dr_spl[iatm_typ_now],r,
                             vps0,vps1,vps2,vps3,&vtemp);

              get_dvpsnow_dvr(ind_nr,nsplin_r[iatm_typ_now],
                              rmin_spl[iatm_typ_now],dr_spl[iatm_typ_now],r,
                              vps0,vps1,vps2,vps3,&dvtemp);

              /*  Obtain spherical harmonics at grid point r */
              if( iang > 0){
                get_ylm(dx,dy,dz,r,ylmr,ylmi,dylmr_x,dylmi_x,
                        dylmr_y,dylmi_y,dylmr_z,dylmi_z,&ylm_cons);
                ioff_ang = iang + im + 1;
                ylmr_loc = ylmr[ioff_ang];
                ylmi_loc = ylmi[ioff_ang];

                vnlreal_tmp_all[icount] = (vtemp*ylmr_loc);
                vnlimag_tmp_all[icount] = (vtemp*ylmi_loc);

                dvnlreal_x_tmp_all[icount] = (dvtemp*(dx*rinv)*ylmr_loc
                                           -  vtemp*dylmr_x[ioff_ang]);
                dvnlimag_x_tmp_all[icount] = (dvtemp*(dx*rinv)*ylmi_loc
                                           -  vtemp*dylmi_x[ioff_ang]);
                dvnlreal_y_tmp_all[icount] = (dvtemp*(dy*rinv)*ylmr_loc
                                           -  vtemp*dylmr_y[ioff_ang]);
                dvnlimag_y_tmp_all[icount] = (dvtemp*(dy*rinv)*ylmi_loc
                                           -  vtemp*dylmi_y[ioff_ang]);
                dvnlreal_z_tmp_all[icount] = (dvtemp*(dz*rinv)*ylmr_loc
                                           -  vtemp*dylmr_z[ioff_ang]);
                dvnlimag_z_tmp_all[icount] = (dvtemp*(dz*rinv)*ylmi_loc
                                           -  vtemp*dylmi_z[ioff_ang]);

                icount++;

              }else{  /* L=0 TERM */

                ioff_ang = iang + im + 1;
                ylmr_loc = ylmr[ioff_ang];
                ylmi_loc = 0.000000;

                vnlreal_tmp_all[icount] = (vtemp*ylmr_loc);
                vnlimag_tmp_all[icount] = 0.000000000;
                dvnlreal_x_tmp_all[icount] = (dvtemp*(dx*rinv)*ylmr_loc
                                           -  vtemp*dylmr_x[ioff_ang]);
                dvnlimag_x_tmp_all[icount] = 0.000000000;
                dvnlreal_y_tmp_all[icount] = (dvtemp*(dy*rinv)*ylmr_loc
                                           -  vtemp*dylmr_y[ioff_ang]);
                dvnlimag_y_tmp_all[icount] = 0.000000000;
                dvnlreal_z_tmp_all[icount] = (dvtemp*(dz*rinv)*ylmr_loc
                                           -  vtemp*dylmr_z[ioff_ang]);
                dvnlimag_z_tmp_all[icount] = 0.000000000;
                icount++;
              }/*endif ang > 1*/

        } } }/* end loop over grid points ka kb kc */

      }/* end loop over non-local atoms in channel l*/
    }/* end loop over m components of l channel */
  }/* end loop over angular momentum channels */

                /*****************************************/
                /* Time for Matrix multiplication..      */
                /*****************************************/

/*========================================================================*/
/* V) compute Z_iIlm */

  ngrid = (ekc_fft_ka_proc - skc_fft_ka_proc - 1)*nkf1*nkf2
           + (nkf2 - skb_fft_ka_proc + 1)*nkf1
           + (ekb_fft_ka_proc)*nkf1;

  ldx = ngrid;
  ldy = ngrid;
  ldz = lmi;
  nrz = lmi;
  ncz = nstate;
  nxy = ngrid;

  GEN_MATMUL(&(vnlreal_tmp_all[1]),&ldx,&itransp,
                &(dvrc[1]),&ldy,&inormal,
                &(vnlreal_tmp[1]),&ldz,
                &nrz,&ncz,&nxy,&alpha,&beta);

  GEN_MATMUL(&(vnlimag_tmp_all[1]),&ldx,&itransp,
                &(dvrc[1]),&ldy,&inormal,
                &(vnlimag_tmp[1]),&ldz,
                &nrz,&ncz,&nxy,&alpha,&beta);

  if(np_states > 1){
    Allreduce(&(vnlreal_tmp[1]),&(vnlreal[1]),(lmi*nstate),MPI_DOUBLE,
                MPI_SUM,0,comm_states);
    Allreduce(&(vnlimag_tmp[1]),&(vnlimag[1]),(lmi*nstate),MPI_DOUBLE,
                MPI_SUM,0,comm_states);
  }else{
    for(i=1; i <= (lmi*nstate) ; i++){
      vnlreal[i] = vnlreal_tmp[i];
      vnlimag[i] = vnlimag_tmp[i];
    }
  }

/*==============================================================================*/
/* VI) Forces on Ions  */

  GEN_MATMUL(&(dvnlreal_x_tmp_all[1]),&ldx,&itransp,
                &(dvrc[1]),&ldy,&inormal,
                &(dvnlreal_x_tmp[1]),&ldz,
                &nrz,&ncz,&nxy,&alpha,&beta);

  GEN_MATMUL(&(dvnlimag_x_tmp_all[1]),&ldx,&itransp,
                &(dvrc[1]),&ldy,&inormal,
                &(dvnlimag_x_tmp[1]),&ldz,
                &nrz,&ncz,&nxy,&alpha,&beta);

  GEN_MATMUL(&(dvnlreal_y_tmp_all[1]),&ldx,&itransp,
                &(dvrc[1]),&ldy,&inormal,
                &(dvnlreal_y_tmp[1]),&ldz,
                &nrz,&ncz,&nxy,&alpha,&beta);

  GEN_MATMUL(&(dvnlimag_y_tmp_all[1]),&ldx,&itransp,
                &(dvrc[1]),&ldy,&inormal,
                &(dvnlimag_y_tmp[1]),&ldz,
                &nrz,&ncz,&nxy,&alpha,&beta);

  GEN_MATMUL(&(dvnlreal_z_tmp_all[1]),&ldx,&itransp,
                &(dvrc[1]),&ldy,&inormal,
                &(dvnlreal_z_tmp[1]),&ldz,
                &nrz,&ncz,&nxy,&alpha,&beta);

  GEN_MATMUL(&(dvnlimag_z_tmp_all[1]),&ldx,&itransp,
                &(dvrc[1]),&ldy,&inormal,
                &(dvnlimag_z_tmp[1]),&ldz,
                &nrz,&ncz,&nxy,&alpha,&beta);

  if(np_states > 1){
    Allreduce(&(dvnlreal_x_tmp[1]),&(dvnlreal_x[1]),(lmi*nstate),MPI_DOUBLE,
              MPI_SUM,0,comm_states);
    Allreduce(&(dvnlimag_x_tmp[1]),&(dvnlimag_x[1]),(lmi*nstate),MPI_DOUBLE,
              MPI_SUM,0,comm_states);

    Allreduce(&(dvnlreal_y_tmp[1]),&(dvnlreal_y[1]),(lmi*nstate),MPI_DOUBLE,
              MPI_SUM,0,comm_states);
    Allreduce(&(dvnlimag_y_tmp[1]),&(dvnlimag_y[1]),(lmi*nstate),MPI_DOUBLE,
              MPI_SUM,0,comm_states);

    Allreduce(&(dvnlreal_z_tmp[1]),&(dvnlreal_z[1]),(lmi*nstate),MPI_DOUBLE,
              MPI_SUM,0,comm_states);
    Allreduce(&(dvnlimag_z_tmp[1]),&(dvnlimag_z[1]),(lmi*nstate),MPI_DOUBLE,
              MPI_SUM,0,comm_states);

  }else{
    for(i=1; i <= (lmi*nstate) ; i++){
      dvnlreal_x[i] = dvnlreal_x_tmp[i];
      dvnlimag_x[i] = dvnlimag_x_tmp[i];

      dvnlreal_y[i] = dvnlreal_y_tmp[i];
      dvnlimag_y[i] = dvnlimag_y_tmp[i];

      dvnlreal_z[i] = dvnlreal_z_tmp[i];
      dvnlimag_z[i] = dvnlimag_z_tmp[i];
    }
  }

/*========================================================================*/
/* VII) Calculate Non-local energy and ion forces */

  cp_enl = 0.000;

  icount = 1;

  /* sum over i (state) */
  for(is=1;is<=nstate;is++){
    /* sum over l (angular momentum) */
    for(iang=0;iang <= nl_max;iang++){
      lp1  = iang+1;
      ioff = natm_tot*iang;

      /* get the normalization factor */
      get_vpsnorm(vscr,vpsnorm,vnorm,iatm_typ_nl,natm_typ,np_nonloc_cp_box_kb,
                  iang,n_ang_max,1,1,1);

      i_shift = iang*natm_tot;

      for(iatm=1 ;iatm<=np_nl[lp1];iatm++){
        vnorm_now[iatm] = vnorm[ip_nl_rev[(iatm+i_shift)]];
      }

      /* sum over m (ang. projection) */
      for(im=0; im < (2*iang+1);im++){
        /* sum over I (ions) */
        for(iatm=1;iatm <= np_nl[(iang+1)];iatm++){

          
          cp_enl += (vnlreal[icount]*vnlreal[icount]
                    +vnlimag[icount]*vnlimag[icount])*vnorm_now[iatm];

          ind_atm  = ip_nl[ioff + iatm];

          fx[ind_atm] += (dvnlreal_x[icount]*vnlreal[icount]
                       + dvnlimag_x[icount]*vnlimag[icount])*vnorm_now[iatm]*tdfact_sq;
          fy[ind_atm] += (dvnlreal_y[icount]*vnlreal[icount]
                       + dvnlimag_y[icount]*vnlimag[icount])*vnorm_now[iatm]*tdfact_sq;
          fz[ind_atm] += (dvnlreal_z[icount]*vnlreal[icount]
                       + dvnlimag_z[icount]*vnlimag[icount])*vnorm_now[iatm]*tdfact_sq;

          vnlreal[icount] *= vnorm_now[iatm];
          vnlimag[icount] *= vnorm_now[iatm];

          icount++;

        }/*end loop over I atoms */
      }/*end loop over m channels */
    }/*end loop over l channels */
  }/*end loop over i states */

  *pcp_enl = cp_enl*dfact*dfact;

/*==============================================================================*/
/* VIII) Forces on coefficients */

  for(i=1; i <= nstate*ngrid; i++){
    freal_tmp[i] =  0.00;
    fimag_tmp[i] =  0.00;
  }

  ldx = ngrid;
  ldy = lmi;
  ldz = ngrid;
  nrz = ngrid;
  ncz = nstate;
  nxy = lmi;

  GEN_MATMUL(&(vnlreal_tmp_all[1]),&ldx,&inormal,
                &(vnlreal[1]),&ldy,&inormal,
                &(freal_tmp[1]),&ldz,
                &nrz,&ncz,&nxy,&alpha,&beta);

  icount = 1;
  for(is=1; is <= nstate; is++){
    for(i=1; i <= ngrid; i++){
      dvrfc[icount] -= tdfact_sq*freal_tmp[icount];
      icount++;
    }
  }

  GEN_MATMUL(&(vnlimag_tmp_all[1]),&ldx,&inormal,
               &(vnlimag[1]),&ldy,&inormal,
                &(fimag_tmp[1]),&ldz,
                &nrz,&ncz,&nxy,&alpha,&beta);

  icount = 1;
  for(is=1; is <= nstate; is++){
    for(i=1; i <= ngrid; i++){
      dvrfc[icount] -= tdfact_sq*fimag_tmp[icount];
      icount++;
    }
  }

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/
