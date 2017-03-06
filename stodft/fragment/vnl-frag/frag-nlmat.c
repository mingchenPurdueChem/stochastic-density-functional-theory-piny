/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: frag-nlmat.c                                 */
/*                                                                          */
/* This routine calculate the kinectic energy, non-local pseudo-potential   */
/* energy and nuclei force correction.                                      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_frag_local.h"
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
calcNonLocalMatrix(CP *cp, CP *cpMini)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
/**************************************************************************/
/* This function calculate non-local pseudopotential matrix w.r.t. frag-  */
/* -ment MO, as well as the force component.				  */
/**************************************************************************/

/*======================================================================*/
/* 0) Check the forms                                                   */

  if(icoef_orth_up!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The Up coefficients must be in orthogonal form    \n");
    printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/
  if(cp_lsda==1 && nstate_dn != 0){
    if(icoef_orth_dn!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The dn coefficients must be in orthogonal form    \n");
      printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
  }/*endif*/

  if(np_states>1){
    if((icoef_form_up+ifcoef_form_up)!=0){
      //printf("icoef_form_up %i ifcoef_form_up %i\n",icoef_form_up,ifcoef_form_up);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The up coefs and coef forces must not be in transposed form\n");
      printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
    if(cp_lsda==1 && nstate_dn != 0){
      if((icoef_form_dn+ifcoef_form_dn)!=0){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("The dn coefs and coef forces must not be in transposed form\n");
        printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);exit(1);
      }/*endif*/
    }/*endif*/
  }/*endif*/
/*======================================================================*/
/* I) Get some useful constants                                         */
  tpi                 = 2.0*M_PI;
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
/*======================================================================*/
/* II) Determine the maximum open non-local angular momentum channel     */
  nl_max = -1;
  for(i=1;i<=(n_ang_max_kb +1);i++){
    if(np_nl[i]>0){nl_max=i-1;}
  }/*endfor*/
  nl_chan_max = (nl_max + 1)*(nl_max + 1);
/*======================================================================*/
/* III) Determine the maximum open non-local angular momentum channel   */
/*      for Kleinman-Bylander and Goedecker pseudo potentials           */
  nl_max_kb = -1;
  for(i=1;i<=(n_ang_max_kb+1);i++){
    if(np_nl[i]>0){nl_max_kb=i-1;}
  }/*endfor*/
  nl_max_gh = -1;
  for(i=1;i<=(n_ang_max_gh+1);i++){
    if(np_nl_gh[i]>0){nl_max_gh=i-1;}
  }/*endfor*/
/*======================================================================*/
/* IV) Determine the maximum number of atoms in any                     */
/*       open angular momentum channel                                  */
  np_nlmax_kb = 1;
  for(i = 1;i<=(nl_max_kb+1);i++){
    np_nlmax_kb = MAX(np_nlmax_kb,np_nl[i]);
  }/*endfor*/
  np_nlmax_gh = 1;
  for(i = 1;i<=(nl_max_gh+1);i++){
    np_nlmax_gh = MAX(np_nlmax_gh,np_nl_gh[i]);
  }/*endfor*/
  np_nlmax_all = (np_nlmax_gh > np_nlmax_kb ? np_nlmax_gh : np_nlmax_kb);
/*======================================================================*/
/* IV) Find cos and sin of sc components of the particles               */
/*    ( hmati rvec = svec   r=(x,y,z) s=(a,b,c) )                       */
  for(ipart=1;ipart<= np_nonloc_cp_box_kb;ipart++){
    iatm = ip_nl[ipart];
    dx  = x[iatm] - cp_box_center[1];
    dy  = y[iatm] - cp_box_center[2];
    dz  = z[iatm] - cp_box_center[3];
    asx = dx*hmati_big[1]+dy*hmati_big[4]+dz*hmati_big[7];
    asy = dx*hmati_big[2]+dy*hmati_big[5]+dz*hmati_big[8];
    asz = dx*hmati_big[3]+dy*hmati_big[6]+dz*hmati_big[9];
    sx  = asx - NINT(asx);
    sy  = asy - NINT(asy);
    sz  = asz - NINT(asz);
    dx  = sx*hmat_big[1]+sy*hmat_big[4]+sz*hmat_big[7];
    dy  = sx*hmat_big[2]+sy*hmat_big[5]+sz*hmat_big[8];
    dz  = sx*hmat_big[3]+sy*hmat_big[6]+sz*hmat_big[9];
    xtemp = dx + cp_box_center_rel[1];
    ytemp = dy + cp_box_center_rel[2];
    ztemp = dz + cp_box_center_rel[3];
    ewd_scr_x[ipart] = xtemp*hmati_cp[1]
                     + ytemp*hmati_cp[4]
                     + ztemp*hmati_cp[7];
    ewd_scr_y[ipart] = xtemp*hmati_cp[2]
                     + ytemp*hmati_cp[5]
                     + ztemp*hmati_cp[8];
    ewd_scr_z[ipart] = xtemp*hmati_cp[3]
                     + ytemp*hmati_cp[6]
                     + ztemp*hmati_cp[9];
    ctemp = ewd_scr_z[ipart]*tpi;
    cossc[ipart] = cos(ctemp);
    sinsc[ipart] = sin(ctemp);
  }/*endfor*/
/*======================================================================*/
/* V) Perform the ewald sum/ non-local potential calculation            */
  for(icount=1;icount<=nktot_sm;icount++){
/*----------------------------------------------------------------------*/
/* i) Get the k vectors                                                 */
    aka = (double)(kastore_sm[icount]);
    akb = (double)(kbstore_sm[icount]);
    akc = (double)(kcstore_sm[icount]);
    xk = (aka*hmati_cp[1]+akb*hmati_cp[2]+akc*hmati_cp[3])*tpi;
    yk = (aka*hmati_cp[4]+akb*hmati_cp[5]+akc*hmati_cp[6])*tpi;
    zk = (aka*hmati_cp[7]+akb*hmati_cp[8]+akc*hmati_cp[9])*tpi;
    g2 = xk*xk+yk*yk+zk*zk;
    g  = sqrt(g2);
    ak2_sm[icount] = g2;
/*----------------------------------------------------------------------*/
/* ii) If break point number one calculate the helpful vectors          */
    if(ibreak1_sm[icount]==1){
     for(ipart=1;ipart<=np_nonloc_cp_box_kb;ipart++){
       atemp = ewd_scr_x[ipart];
       btemp = ewd_scr_y[ipart];
       ctemp = ewd_scr_z[ipart];
       arg = (aka*atemp + akb*btemp + akc*ctemp)*tpi;
       helr[ipart] = cos(arg);
       heli[ipart] = sin(arg);
     }/*endfor*/
   }/*endif*/
/*----------------------------------------------------------------------*/
/* iii) nonlocal matrix                                                 */
/*     use structure factor and wavefunction to compute                 */
/*     the array znl and its derivative dznl for the                    */
/*     calculation of the nonlocal forces and energy                    */
/*     (done separately for each angular momentum component)            */
    controlNlmatFrag(clatoms_info,cpcoeffs_info,cpcoeffs_pos,
                  cpscr,cpopts,pseudo,ewd_scr,atommaps,
                  np_nlmax,nl_max,index_atm,g,xk,yk,zk,
                  icount,&ylm_cons);

/*----------------------------------------------------------------------*/
/* iv) If break point two, increment the helpful vectors                */

    if(ibreak2_sm[icount]==1){
      for(ipart=1;ipart<=np_nonloc_cp_box_kb;ipart++){
	temp = helr[ipart];
	helr[ipart] = helr[ipart]*cossc[ipart] - heli[ipart]*sinsc[ipart];
	heli[ipart] = heli[ipart]*cossc[ipart] + temp*sinsc[ipart];
      }/*endfor*/
    }/*endif*/
  }/*endfor:icount loop over k vectors */
/*======================================================================*/
/* VI) g=0 term                                                         */

  ak2_sm[ncoef] = 0.0;
  if(np_nl[1]>0){

    ylmr[1] = 1.0/sqrt(fpi);
    for(irad=1;irad<=nrad_max_l[1];irad++){

      for(ipart=1;ipart<=np_nonloc_cp_box_kb;ipart++){
        iii = n_rad_max*(iatm_typ_nl[ipart]-1) + irad;
        vtemp[ipart] = gzvps0[iii];
      }/*endfor*/

      for(ipart=np_nl_rad_str[1][irad];ipart<=np_nl[1];ipart++){
        vtemp_now[ipart] = vtemp[ipart]*ylmr[1];
      }/*endfor*/

      get_nlmat0(ncoef,ncoef,nstate_up,np_nlmax,
                 np_nl_rad_str[1][irad],np_nl[1],nl_chan_max,irad,
                 creal_up,cimag_up,vtemp_now,vnlreal_up,vnlimag_up,ylmr[1]);
      if(cp_lsda==1){
        get_nlmat0(ncoef,ncoef,nstate_dn,np_nlmax,
                  np_nl_rad_str[1][irad],np_nl[1],nl_chan_max,irad,
                  creal_dn,cimag_dn,vtemp_now,vnlreal_dn,vnlimag_dn,ylmr[1]);
      }/*endif*/

    }/*endfor*/

  }/*endif: l=0 nonlocal*/

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void controlNlmatFrag(CLATOMS_INFO *clatoms_info,
                   CPCOEFFS_INFO *cpcoeffs_info,
                   CPCOEFFS_POS *cpcoeffs_pos,
                   CPSCR *cpscr,CPOPTS *cpopts,PSEUDO *pseudo,
                   EWD_SCR *ewd_scr,ATOMMAPS *atommaps,FRAGINFO *fragInfo,
                   int np_nlmax,int nl_max,int *index_atm,
                   double g,double xk,double yk,double zk,
                   int ismcount,YLM_CONS *ylm_cons)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/

/*=======================================================================*/
/*         Local Variable declarations                                   */

  int ind_lm,ipart,lp1,irad,jrad,ipart_nl;
  double sgn_l;
  double ylmr[21], ylmi[21];
  double dylmr_gx[21], dylmi_gx[21];
  double dylmr_gy[21], dylmi_gy[21];
  double dylmr_gz[21], dylmi_gz[21];
  int itype,itype_nl,l,m,i_shift,iii;
  int ktemp,ltemp;
  int nl_chan_max;
/* Local pointers */
  int npart                = clatoms_info->natm_tot;
  int natm_typ             = atommaps->natm_typ;
  int *iatm_typ            = atommaps->iatm_atm_typ;
  int *iatm_typ_nl         = atommaps->iatm_atm_typ_nl;
  int *iatm_typ_nl_rev     = atommaps->iatm_atm_typ_nl_rev;
  int *imap_atm_typ_nl     = atommaps->imap_atm_typ_nl;
  int  natm_typ_nl         = pseudo->natm_typ_nl;

  int nstate_up            = cpcoeffs_info->nstate_up_proc;
  int nstate_dn            = cpcoeffs_info->nstate_dn_proc;
  int ncoef                = cpcoeffs_info->ncoef;
  double *creal_up         = cpcoeffs_pos->cre_up;
  double *cimag_up         = cpcoeffs_pos->cim_up;
  double *creal_dn         = cpcoeffs_pos->cre_dn;
  double *cimag_dn         = cpcoeffs_pos->cim_dn;
  int cp_lsda              = cpopts->cp_lsda;
  int cp_hess_calc         = cpopts->cp_hess_calc;
  int cp_ptens             = cpopts->cp_ptens_calc;
  int atm_hess_calc        = clatoms_info->hess_calc;

  int *ip_nl               = pseudo->ip_nl;
  int *ip_nl_rev           = pseudo->ip_nl_rev;
  int *np_nl               = pseudo->np_nl;
  int nsplin_g             = pseudo->nsplin_g;
  double dg_spl            = pseudo->dg_spl;
  double gmin_spl          = pseudo->gmin_spl;
  double *vps0             = pseudo->vps0;
  double *vps1             = pseudo->vps1;
  double *vps2             = pseudo->vps2;
  double *vps3             = pseudo->vps3;
  double *dvps0            = pseudo->dvps0;
  double *dvps1            = pseudo->dvps1;
  double *dvps2            = pseudo->dvps2;
  double *dvps3            = pseudo->dvps3;
  int n_ang_max            = pseudo->n_ang_max;
  int n_ang_max_kb         = pseudo->n_ang_max_kb;
  int n_rad_max            = pseudo->n_rad_max;
  int *nrad_max_l          = pseudo->nrad_max_l;
  int **np_nl_rad_str      = pseudo->np_nl_rad_str;

  int np_nonloc_cp_box_kb  = pseudo->np_nonloc_cp_box_kb;

  int np_nonloc_cp_box_kb  = pseudo->np_nonloc_cp_box_kb;

  double *vnlreal_up       = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up       = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn       = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn       = cpscr->cpscr_nonloc.vnlim_dn;
  double *dvnlreal_x_up    = cpscr->cpscr_nonloc.dvnlre_x_up;
  double *dvnlreal_y_up    = cpscr->cpscr_nonloc.dvnlre_y_up;
  double *dvnlreal_z_up    = cpscr->cpscr_nonloc.dvnlre_z_up;
  double *dvnlimag_x_up    = cpscr->cpscr_nonloc.dvnlim_x_up;
  double *dvnlimag_y_up    = cpscr->cpscr_nonloc.dvnlim_y_up;
  double *dvnlimag_z_up    = cpscr->cpscr_nonloc.dvnlim_z_up;
  double *dvnlreal_x_dn    = cpscr->cpscr_nonloc.dvnlre_x_dn;
  double *dvnlreal_y_dn    = cpscr->cpscr_nonloc.dvnlre_y_dn;
  double *dvnlreal_z_dn    = cpscr->cpscr_nonloc.dvnlre_z_dn;
  double *dvnlimag_x_dn    = cpscr->cpscr_nonloc.dvnlim_x_dn;
  double *dvnlimag_y_dn    = cpscr->cpscr_nonloc.dvnlim_y_dn;
  double *dvnlimag_z_dn    = cpscr->cpscr_nonloc.dvnlim_z_dn;
  double *dvnlreal_gxgx_up = cpscr->cpscr_nonloc.dvnlre_gxgx_up;
  double *dvnlimag_gxgx_up = cpscr->cpscr_nonloc.dvnlim_gxgx_up;
  double *dvnlreal_gygy_up = cpscr->cpscr_nonloc.dvnlre_gygy_up;
  double *dvnlimag_gygy_up = cpscr->cpscr_nonloc.dvnlim_gygy_up;
  double *dvnlreal_gzgz_up = cpscr->cpscr_nonloc.dvnlre_gzgz_up;
  double *dvnlimag_gzgz_up = cpscr->cpscr_nonloc.dvnlim_gzgz_up;
  double *dvnlreal_gxgy_up = cpscr->cpscr_nonloc.dvnlre_gxgy_up;
  double *dvnlimag_gxgy_up = cpscr->cpscr_nonloc.dvnlim_gxgy_up;
  double *dvnlreal_gygz_up = cpscr->cpscr_nonloc.dvnlre_gygz_up;
  double *dvnlimag_gygz_up = cpscr->cpscr_nonloc.dvnlim_gygz_up;
  double *dvnlreal_gxgz_up = cpscr->cpscr_nonloc.dvnlre_gxgz_up;
  double *dvnlimag_gxgz_up = cpscr->cpscr_nonloc.dvnlim_gxgz_up;
  double *dvnlreal_gxgx_dn = cpscr->cpscr_nonloc.dvnlre_gxgx_dn;
  double *dvnlimag_gxgx_dn = cpscr->cpscr_nonloc.dvnlim_gxgx_dn;
  double *dvnlreal_gygy_dn = cpscr->cpscr_nonloc.dvnlre_gygy_dn;
  double *dvnlimag_gygy_dn = cpscr->cpscr_nonloc.dvnlim_gygy_dn;
  double *dvnlreal_gzgz_dn = cpscr->cpscr_nonloc.dvnlre_gzgz_dn;
  double *dvnlimag_gzgz_dn = cpscr->cpscr_nonloc.dvnlim_gzgz_dn;
  double *dvnlreal_gxgy_dn = cpscr->cpscr_nonloc.dvnlre_gxgy_dn;
  double *dvnlimag_gxgy_dn = cpscr->cpscr_nonloc.dvnlim_gxgy_dn;
  double *dvnlreal_gygz_dn = cpscr->cpscr_nonloc.dvnlre_gygz_dn;
  double *dvnlimag_gygz_dn = cpscr->cpscr_nonloc.dvnlim_gygz_dn;
  double *dvnlreal_gxgz_dn = cpscr->cpscr_nonloc.dvnlre_gxgz_dn;
  double *dvnlimag_gxgz_dn = cpscr->cpscr_nonloc.dvnlim_gxgz_dn;

  double *helr             = ewd_scr->helr;
  double *heli             = ewd_scr->heli;
  double *helr_now         = ewd_scr->helr_now;
  double *heli_now         = ewd_scr->heli_now;
  double *dheli_now        = ewd_scr->temp;
  double *dhelr_now        = ewd_scr->vtemp_now;
  double *vtemp            = ewd_scr->fx2;
  double *dvtemp           = ewd_scr->q;
  double *scr1             = ewd_scr->fx2; /* yes the same as vtemp, its ok */
  double *scr2             = ewd_scr->q;   /* yes the same as dvtemp,its ok */
  double *scr3             = ewd_scr->fz2;

  //fragment related
  int iFrag = fragInfo->iFrag;
  double *vnlMatrixUp;
  double *vnlMatrixDn;
  double *vnlForceMatrixUp;
  double *vnlForceMatrixDn;


  vnlMatrixUp = fragInfo->vnlMatrixUp[iFrag];
  vnlForceMatrixUp = fragInfo->vnlForceMatrixUp[iFrag];
  if(cp_lsda==1){
    vnlMatrixDn = fragInfo->vnlMatrixDn[iFrag];
    vnlForceMatrixDn = fragInfo->vnlForceMatrixDn[iFrag];
  }

/*======================================================================*/
/* I) Get the ylm(g)                                                   */

  get_ylm(xk,yk,zk,g,ylmr,ylmi,dylmr_gx,dylmi_gx,dylmr_gy,dylmi_gy,
          dylmr_gz,dylmi_gz,ylm_cons);
  nl_chan_max = (nl_max+1)*(nl_max+1);
/*======================================================================*/
/* II) Calculate the nl-pseudoponential matrix elements by looping over  */
/* the channels, l, and then the 2l+1 components of the channel         */

  for(l=0;l<=nl_max;l++){
    lp1 = l+1;
    if(np_nl[lp1]>0){
      for(irad=1;irad<=nrad_max_l[lp1];irad++){
/*---------------------------------------------------------------------*/
/* i) Get the bessel transform of the pseudopotential at this g vector */
/*    and scale the structure factor appropriately                     */
        for(itype=1;itype<=natm_typ_nl;itype++){
          itype_nl = imap_atm_typ_nl[itype];
          index_atm[itype] = (itype_nl-1)*nsplin_g*(n_ang_max+1)*n_rad_max
                             +l*nsplin_g*n_rad_max+(irad-1)*nsplin_g;
        }/*endfor*/
        getVpsNowFrag(index_atm,nsplin_g,gmin_spl,dg_spl,g,
                      vps0,vps1,vps2,vps3,vtemp,iatm_typ_nl_rev,
                      natm_typ_nl,np_nonloc_cp_box_kb,
                      np_nl_rad_str[lp1][irad]);
	i_shift = l*npart;
	if(cp_ptens==0){
	  for(ipart=np_nl_rad_str[lp1][irad];ipart<=np_nl[lp1];ipart++){
	    ktemp             =  ipart+i_shift;
	    ltemp             =  ip_nl_rev[ktemp];
	    helr_now[ipart]   =  helr[ltemp]*vtemp[ltemp];
	    heli_now[ipart]   =  heli[ltemp]*vtemp[ltemp];
	  }/*endfor*/
	}else{
	  /*PRESSURE TENSOR CALC*/
	  get_vpsnow(index_atm,nsplin_g,gmin_spl,dg_spl,g,
		     dvps0,dvps1,dvps2,dvps3,dvtemp,iatm_typ,
		     natm_typ,npart,np_nl_rad_str[lp1][irad]);
	  for(ipart=np_nl_rad_str[lp1][irad];ipart<=np_nl[lp1];ipart++){
	    ktemp              = ipart+i_shift;
	    ltemp              =  ip_nl_rev[ktemp];
	    helr_now[ipart]    = helr[ltemp]*vtemp[ltemp];
	    heli_now[ipart]    = heli[ltemp]*vtemp[ltemp];
	    dhelr_now[ipart]   = helr[ltemp]*dvtemp[ltemp];
	    dheli_now[ipart]   = heli[ltemp]*dvtemp[ltemp];
	  }/*endfor*/
	}/*endif*/
/*---------------------------------------------------------------------*/
/* ii) Loop over m components of the channel and get the matrix elements */
/*     vtemp,dvtemp, */
	sgn_l = 1.0;
	if((l%2)==1)sgn_l = -1.0;
	for(m=1;m<=(2*l+1);m++){
	  ind_lm = m + l*l;
	  /* // I dont have cp_ptens and atm_hess
	  if(cp_ptens==1) {
	    get_nlmat_pv(ncoef,ismcount,nstate_up,ind_lm,irad,
	      np_nlmax,np_nl[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
	      helr_now,heli_now,
	      creal_up,cimag_up,dhelr_now,dheli_now,
	      vnlreal_up,vnlimag_up,dvnlreal_x_up,dvnlreal_y_up,
	      dvnlreal_z_up,dvnlimag_x_up,dvnlimag_y_up,
	      dvnlimag_z_up,dvnlreal_gxgx_up,dvnlreal_gxgy_up,
	      dvnlreal_gxgz_up,dvnlreal_gygy_up,dvnlreal_gygz_up,
	      dvnlreal_gzgz_up,
	      dvnlimag_gxgx_up,dvnlimag_gxgy_up,dvnlimag_gxgz_up,
	      dvnlimag_gygy_up,dvnlimag_gygz_up,dvnlimag_gzgz_up,
	      xk,yk,zk,
	      ylmr[ind_lm],ylmi[ind_lm],
	      dylmr_gx[ind_lm],dylmi_gx[ind_lm],
	      dylmr_gy[ind_lm],dylmi_gy[ind_lm],
	      dylmr_gz[ind_lm],dylmi_gz[ind_lm],
	      sgn_l,scr1,scr2,scr3);
	  }// endif
	  if(atm_hess_calc == 3){
	    get_nlmat_hess(ncoef,ismcount,nstate_up,ind_lm,irad,
	      np_nlmax,np_nl[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
	      helr_now,heli_now,
	      creal_up,cimag_up,dhelr_now,dheli_now,
	      vnlreal_up,vnlimag_up,dvnlreal_x_up,dvnlreal_y_up,
	      dvnlreal_z_up,dvnlimag_x_up,dvnlimag_y_up,
	      dvnlimag_z_up,dvnlreal_gxgx_up,dvnlreal_gxgy_up,
	      dvnlreal_gxgz_up,dvnlreal_gygy_up,dvnlreal_gygz_up,
	      dvnlreal_gzgz_up,
	      dvnlimag_gxgx_up,dvnlimag_gxgy_up,dvnlimag_gxgz_up,
	      dvnlimag_gygy_up,dvnlimag_gygz_up,dvnlimag_gzgz_up,
	      xk,yk,zk,
	      ylmr[ind_lm],ylmi[ind_lm],
	      sgn_l,scr1,scr2);
	  }// endif
	  */
	  if(atm_hess_calc != 3 && cp_ptens == 0){
	    get_nlmat(ncoef,ismcount,nstate_up,ind_lm,irad,
	      np_nlmax,np_nl[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
	      helr_now,heli_now,
	      creal_up,cimag_up,
	      vnlreal_up,vnlimag_up,
	      dvnlreal_x_up,dvnlreal_y_up,dvnlreal_z_up,
	      dvnlimag_x_up,dvnlimag_y_up,dvnlimag_z_up,
	      xk,yk,zk,
	      ylmr[ind_lm],ylmi[ind_lm],sgn_l,scr1,scr2);
	  }/*endif*/


	  if(cp_lsda==1) {
	    /*
	    if(cp_ptens==1) {
	      get_nlmat_pv(ncoef,ismcount,nstate_dn,ind_lm,irad,
		np_nlmax,np_nl[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
		helr_now,heli_now,
		creal_dn,cimag_dn,dhelr_now,dheli_now,
		vnlreal_dn,vnlimag_dn,
		dvnlreal_x_dn,dvnlreal_y_dn,dvnlreal_z_dn,
		dvnlimag_x_dn,dvnlimag_y_dn,dvnlimag_z_dn,
		dvnlreal_gxgx_dn,dvnlreal_gxgy_dn,dvnlreal_gxgz_dn,
		dvnlreal_gygy_dn,dvnlreal_gygz_dn,dvnlreal_gzgz_dn,
		dvnlimag_gxgx_dn,dvnlimag_gxgy_dn,dvnlimag_gxgz_dn,
		dvnlimag_gygy_dn,dvnlimag_gygz_dn,dvnlimag_gzgz_dn,
		xk,yk,zk,
		ylmr[ind_lm],ylmi[ind_lm],
		dylmr_gx[ind_lm],dylmi_gx[ind_lm],
		dylmr_gy[ind_lm],dylmi_gy[ind_lm],
		dylmr_gz[ind_lm],dylmi_gz[ind_lm],sgn_l,scr1,scr2,scr3);
	    }// endif 
	    if(atm_hess_calc == 3){
	      get_nlmat_hess(ncoef,ismcount,nstate_dn,ind_lm,irad,
		np_nlmax,np_nl[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
		helr_now,heli_now,
		creal_dn,cimag_dn,dhelr_now,dheli_now,
		vnlreal_dn,vnlimag_dn,
		dvnlreal_x_dn,dvnlreal_y_dn,dvnlreal_z_dn,
		dvnlimag_x_dn,dvnlimag_y_dn,dvnlimag_z_dn,
		dvnlreal_gxgx_dn,dvnlreal_gxgy_dn,dvnlreal_gxgz_dn,
		dvnlreal_gygy_dn,dvnlreal_gygz_dn,dvnlreal_gzgz_dn,
		dvnlimag_gxgx_dn,dvnlimag_gxgy_dn,dvnlimag_gxgz_dn,
		dvnlimag_gygy_dn,dvnlimag_gygz_dn,dvnlimag_gzgz_dn,
		xk,yk,zk,
		ylmr[ind_lm],ylmi[ind_lm],
		sgn_l,scr1,scr2);
	    }// endif
	    */
	    if(atm_hess_calc != 3 && cp_ptens == 0){
	       get_nlmat(ncoef,ismcount,nstate_dn,ind_lm,irad,
		 np_nlmax,np_nl[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
		 helr_now,heli_now,
		 creal_dn,cimag_dn,
		 vnlreal_dn,vnlimag_dn,
		 dvnlreal_x_dn,dvnlreal_y_dn,dvnlreal_z_dn,
		 dvnlimag_x_dn,dvnlimag_y_dn,dvnlimag_z_dn,xk,yk,zk,
		 ylmr[ind_lm],ylmi[ind_lm],sgn_l,scr1,scr2);
	    }//endif:pvten
	  }//endif:lsda
        }//endfor:m quantum number loop
      }//endfor: irad quantum number loop 
    }//endif: this channel has atoms in it
  }//endfor:l quantum number loop

/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void getnl_pot_pv_fatm(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                  CELL *cell,CPCOEFFS_INFO *cpcoeffs_info,CPSCR *cpscr,
                  EWD_SCR *ewd_scr,CPOPTS *cpopts,
                  PSEUDO *pseudo,ATOMMAPS *atommaps,
                  double *cp_enl_ret, int np_nlmax,double *pvten)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
  int i,l,m,i_shift,ipart,is,iii,ioff,lp1;
  int iatm;
  int ind_loc,nl_max,irad,jrad;
  int nl_chan_max;
  double rvol_cp,vol_cp,cp_enl;
  double p11,p22,p33,p12,p13,p23;

/* Local pointers */
  int npart                = clatoms_info->natm_tot;
  double *fx               = clatoms_pos->fx;
  double *fy               = clatoms_pos->fy;
  double *fz               = clatoms_pos->fz;
  double *hess_xx          = clatoms_pos->hess_xx;
  double *hess_xy          = clatoms_pos->hess_xy;
  double *hess_xz          = clatoms_pos->hess_xz;
  double *hess_yy          = clatoms_pos->hess_yy;
  double *hess_yz          = clatoms_pos->hess_yz;
  double *hess_zz          = clatoms_pos->hess_zz;
  int natm_typ             = atommaps->natm_typ;
  int *iatm_typ            = atommaps->iatm_atm_typ;
  int *iatm_typ_nl         = atommaps->iatm_atm_typ_nl;
  double *hmat_cp          = cell->hmat_cp;

  int cp_ptens             = cpopts->cp_ptens_calc;
  int cp_lsda              = cpopts->cp_lsda;
  int atm_hess_calc        = clatoms_info->hess_calc;
  int nstate_up            = cpcoeffs_info->nstate_up_proc;
  int nstate_dn            = cpcoeffs_info->nstate_dn_proc;

  int np_nonloc_cp_box_kb  = pseudo->np_nonloc_cp_box_kb;
  double *vpsnorm          = pseudo->vpsnorm;
  int n_ang_max            = pseudo->n_ang_max;
  int n_ang_max_kb         = pseudo->n_ang_max_kb;
  int n_rad_max            = pseudo->n_rad_max;
  int *loc_opt             = pseudo->loc_opt;
  int *ip_nl               = pseudo->ip_nl;
  int *ip_nl_rev           = pseudo->ip_nl_rev;
  int *np_nl               = pseudo->np_nl;
  int *nrad_max_l          = pseudo->nrad_max_l;
  int **np_nl_rad_str      = pseudo->np_nl_rad_str;

  double *vnlreal_up       = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up       = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn       = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn       = cpscr->cpscr_nonloc.vnlim_dn;
  double *dvnlreal_x_up    = cpscr->cpscr_nonloc.dvnlre_x_up;
  double *dvnlreal_y_up    = cpscr->cpscr_nonloc.dvnlre_y_up;
  double *dvnlreal_z_up    = cpscr->cpscr_nonloc.dvnlre_z_up;
  double *dvnlimag_x_up    = cpscr->cpscr_nonloc.dvnlim_x_up;
  double *dvnlimag_y_up    = cpscr->cpscr_nonloc.dvnlim_y_up;
  double *dvnlimag_z_up    = cpscr->cpscr_nonloc.dvnlim_z_up;
  double *dvnlreal_x_dn    = cpscr->cpscr_nonloc.dvnlre_x_dn;
  double *dvnlreal_y_dn    = cpscr->cpscr_nonloc.dvnlre_y_dn;
  double *dvnlreal_z_dn    = cpscr->cpscr_nonloc.dvnlre_z_dn;
  double *dvnlimag_x_dn    = cpscr->cpscr_nonloc.dvnlim_x_dn;
  double *dvnlimag_y_dn    = cpscr->cpscr_nonloc.dvnlim_y_dn;
  double *dvnlimag_z_dn    = cpscr->cpscr_nonloc.dvnlim_z_dn;
  double *dvnlreal_gxgx_up = cpscr->cpscr_nonloc.dvnlre_gxgx_up;
  double *dvnlimag_gxgx_up = cpscr->cpscr_nonloc.dvnlim_gxgx_up;
  double *dvnlreal_gygy_up = cpscr->cpscr_nonloc.dvnlre_gygy_up;
  double *dvnlimag_gygy_up = cpscr->cpscr_nonloc.dvnlim_gygy_up;
  double *dvnlreal_gzgz_up = cpscr->cpscr_nonloc.dvnlre_gzgz_up;
  double *dvnlimag_gzgz_up = cpscr->cpscr_nonloc.dvnlim_gzgz_up;
  double *dvnlreal_gxgy_up = cpscr->cpscr_nonloc.dvnlre_gxgy_up;
  double *dvnlimag_gxgy_up = cpscr->cpscr_nonloc.dvnlim_gxgy_up;
  double *dvnlreal_gygz_up = cpscr->cpscr_nonloc.dvnlre_gygz_up;
  double *dvnlimag_gygz_up = cpscr->cpscr_nonloc.dvnlim_gygz_up;
  double *dvnlreal_gxgz_up = cpscr->cpscr_nonloc.dvnlre_gxgz_up;
  double *dvnlimag_gxgz_up = cpscr->cpscr_nonloc.dvnlim_gxgz_up;
  double *dvnlreal_gxgx_dn = cpscr->cpscr_nonloc.dvnlre_gxgx_dn;
  double *dvnlimag_gxgx_dn = cpscr->cpscr_nonloc.dvnlim_gxgx_dn;
  double *dvnlreal_gygy_dn = cpscr->cpscr_nonloc.dvnlre_gygy_dn;
  double *dvnlimag_gygy_dn = cpscr->cpscr_nonloc.dvnlim_gygy_dn;
  double *dvnlreal_gzgz_dn = cpscr->cpscr_nonloc.dvnlre_gzgz_dn;
  double *dvnlimag_gzgz_dn = cpscr->cpscr_nonloc.dvnlim_gzgz_dn;
  double *dvnlreal_gxgy_dn = cpscr->cpscr_nonloc.dvnlre_gxgy_dn;
  double *dvnlimag_gxgy_dn = cpscr->cpscr_nonloc.dvnlim_gxgy_dn;
  double *dvnlreal_gygz_dn = cpscr->cpscr_nonloc.dvnlre_gygz_dn;
  double *dvnlimag_gygz_dn = cpscr->cpscr_nonloc.dvnlim_gygz_dn;
  double *dvnlreal_gxgz_dn = cpscr->cpscr_nonloc.dvnlre_gxgz_dn;
  double *dvnlimag_gxgz_dn = cpscr->cpscr_nonloc.dvnlim_gxgz_dn;

  double *fxtemp           = ewd_scr->fx;
  double *fytemp           = ewd_scr->fy;
  double *fztemp           = ewd_scr->fz;
  double *vscr             = ewd_scr->fx2;
  double *vnorm            = ewd_scr->fy2;
  double *vnorm_now        = ewd_scr->fz2;

/*======================================================================*/
/* I) Useful constants                                                  */

  vol_cp   = getdeth(hmat_cp);
  rvol_cp  = 1.0/vol_cp;
  cp_enl   = 0.0;
  nl_max   = -1;
  for(i=1;i<=(n_ang_max_kb+1);i++){
    if(np_nl[i]>0){nl_max=i-1;}
  }/*endfor*/
  nl_chan_max = (nl_max+1)*(nl_max+1);

/*======================================================================*/
/* II) Loop over the open channels, the states and get the nl potent,   */
/*     pvten and  particle forces                                       */

  for(l=0;l<=nl_max;l++){
    lp1 = l+1;
    if(np_nl[lp1]>0){
      for(irad=1;irad<=nrad_max_l[lp1];irad++){
        for(jrad=irad;jrad<=nrad_max_l[lp1];jrad++){
/*-----------------------------------------------------------------------*/
/* i) Get the normalization scaled by the volume                        */
	  get_vpsnorm(vscr,vpsnorm,vnorm,iatm_typ_nl,natm_typ,np_nonloc_cp_box_kb,l,
		      n_ang_max,irad,jrad,n_rad_max);
	  i_shift = l*npart;
	  for(ipart=np_nl_rad_str[lp1][jrad];ipart<=np_nl[lp1];ipart++){
	    vnorm_now[ipart] = vnorm[ip_nl_rev[(ipart+i_shift)]]*rvol_cp;
	  }/*endfor*/
/*-----------------------------------------------------------------------*/
/* ii) Sum the contributions over the 2l+1 directions and the states    */
	   sumnl_pot_pv_fatm_hess(npart,nstate_up,np_nlmax,nl_chan_max,
				  np_nl[lp1],l,np_nl_rad_str[lp1][jrad],
				  irad,jrad,
				  ip_nl,vnorm_now,vnlreal_up,vnlimag_up,
				  dvnlreal_gxgx_up,dvnlimag_gxgx_up,
				  dvnlreal_gygy_up,dvnlimag_gygy_up,
				  dvnlreal_gzgz_up,dvnlimag_gzgz_up,
				  dvnlreal_gxgy_up,dvnlimag_gxgy_up,
				  dvnlreal_gxgz_up,dvnlimag_gxgz_up,
				  dvnlreal_gygz_up,dvnlimag_gygz_up,
				  dvnlreal_x_up,dvnlimag_x_up,
				  dvnlreal_y_up,dvnlimag_y_up,
				  dvnlreal_z_up,dvnlimag_z_up,
				  fx,fy,fz,fxtemp,fytemp,fztemp,
				  hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
				  atm_hess_calc,cp_ptens,pvten,&cp_enl);
	   if(cp_lsda==1){
	     sumnl_pot_pv_fatm_hess(npart,nstate_dn,np_nlmax,nl_chan_max,
				    np_nl[lp1],l,np_nl_rad_str[lp1][jrad],
				    irad,jrad,
				    ip_nl,vnorm_now,vnlreal_dn,vnlimag_dn,
				    dvnlreal_gxgx_dn,dvnlimag_gxgx_dn,
				    dvnlreal_gygy_dn,dvnlimag_gygy_dn,
				    dvnlreal_gzgz_dn,dvnlimag_gzgz_dn,
				    dvnlreal_gxgy_dn,dvnlimag_gxgy_dn,
				    dvnlreal_gxgz_dn,dvnlimag_gxgz_dn,
				    dvnlreal_gygz_dn,dvnlimag_gygz_dn,
				    dvnlreal_x_dn,dvnlimag_x_dn,
				    dvnlreal_y_dn,dvnlimag_y_dn,
				    dvnlreal_z_dn,dvnlimag_z_dn,
				    fx,fy,fz,fxtemp,fytemp,fztemp,
				    hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
				    atm_hess_calc,cp_ptens,pvten,&cp_enl);
	   }/*endif*/
	}//endfor jrad
      }//endfor: radial channels irad
    }//endif: l channel open
  }//endfor: l channels

/*======================================================================*/
/* III) Assign the non-local energy  and add it to the pvten            */

  *cp_enl_ret = cp_enl;
  if(cp_ptens==1){
    pvten[1] += cp_enl;
    pvten[5] += cp_enl;
    pvten[9] += cp_enl;
  }/*endif*/

/*======================================================================*/
  }/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_nlmat(int ncoef,int ismcount,int nstate,int ind_lm,int irad,
               int np_nlmax,int np_nl,int np_nl_rad_str,int nl_chan_max,
               double *helr,double *heli,
               double *creal,double *cimag,
               double *vnlreal,double *vnlimag,
               double *dvnlreal_x,double *dvnlreal_y,double *dvnlreal_z,
               double *dvnlimag_x,double *dvnlimag_y,double *dvnlimag_z,
               double xk,double yk,double zk,double ylmr,double ylmi,
               double sgn_l,double *cbyhelr,double *cbyheli)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*        Local variable declarations            */


  int is,i,ipart,iii;
  int ind_loc_c,ioff_v,ind_loc_v;
  int isgn_l;
  double cre_ind,cim_ind;
  double tylmr, tylmi;
  double tmp1,tmp2,tmp3,tmp4;
/*==========================================================================*/
/* I) Loop over states and particles:                                       */
/*     Here helr,heli have vtemp mutliplied in them                         */
/*     See control routine                                                  */
/* storage : atm,state,l,m,n */
  isgn_l = (int) (0.5*(sgn_l+1.0)+1.0);
  tylmr = 2.0*ylmr;
  tylmi = 2.0*ylmi;
  for(is=1;is<=nstate;is++){
    ind_loc_c = ismcount + ncoef*(is-1);
    cre_ind   = creal[ind_loc_c];
    cim_ind   = cimag[ind_loc_c];
    ioff_v    = (is-1)*np_nlmax+(ind_lm-1)*nstate*np_nlmax
               +(irad-1)*nl_chan_max*nstate*np_nlmax;
    switch(isgn_l){
      case 1:
        for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
          cbyhelr[ipart]     =  cre_ind*helr[ipart] - cim_ind*heli[ipart];
          cbyheli[ipart]     =  cre_ind*heli[ipart] + cim_ind*helr[ipart];
        }// endfor loop over particles
        for(i=np_nl_rad_str;i<=np_nl;i++){
          tmp1        = -cbyheli[i]*tylmi;
          tmp2        =  cbyheli[i]*tylmr;
          vnlreal[(i+ioff_v)]    += tmp1;
          vnlimag[(i+ioff_v)]    += tmp2;
        }//endfor
        for(i=np_nl_rad_str;i<=np_nl;i++){
          tmp3        = -cbyhelr[i]*tylmi;
          dvnlreal_x[(i+ioff_v)] += (tmp3*xk);
          dvnlreal_y[(i+ioff_v)] += (tmp3*yk);
          dvnlreal_z[(i+ioff_v)] += (tmp3*zk);
        }//endfor
        for(i=np_nl_rad_str;i<=np_nl;i++){
          tmp4        =  cbyhelr[i]*tylmr;
          dvnlimag_x[(i+ioff_v)] += (tmp4*xk);
          dvnlimag_y[(i+ioff_v)] += (tmp4*yk);
          dvnlimag_z[(i+ioff_v)] += (tmp4*zk);
        }//endfor
        break;
      case 2:
        for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
          cbyhelr[ipart]     =  cre_ind*helr[ipart] - cim_ind*heli[ipart];
          cbyheli[ipart]     =  cre_ind*heli[ipart] + cim_ind*helr[ipart];
        }// endfor : loop over particles
        for(i=np_nl_rad_str;i<=np_nl;i++){
          tmp1        =  cbyhelr[i]*tylmr;
          tmp2        =  cbyhelr[i]*tylmi;
          vnlreal[(i+ioff_v)]    += tmp1;
          vnlimag[(i+ioff_v)]    += tmp2;
        }//endfor
        for(i=np_nl_rad_str;i<=np_nl;i++){
          tmp3        = -cbyheli[i]*tylmr;
          dvnlreal_x[(i+ioff_v)] += (tmp3*xk);
          dvnlreal_y[(i+ioff_v)] += (tmp3*yk);
          dvnlreal_z[(i+ioff_v)] += (tmp3*zk);
        }// endfor
        for(i=np_nl_rad_str;i<=np_nl;i++){
          tmp4        = -cbyheli[i]*tylmi;
          dvnlimag_x[(i+ioff_v)] += (tmp4*xk);
          dvnlimag_y[(i+ioff_v)] += (tmp4*yk);
          dvnlimag_z[(i+ioff_v)] += (tmp4*zk);
        }// endfor
	break;
    }//end switch: sgn of m in Ylm
  }// endfor : loop over states

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void getnl_pot_pv_fatm(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                  CELL *cell,CPCOEFFS_INFO *cpcoeffs_info,CPSCR *cpscr,
                  EWD_SCR *ewd_scr,CPOPTS *cpopts,
                  PSEUDO *pseudo,ATOMMAPS *atommaps,
                  double *cp_enl_ret, int np_nlmax,double *pvten,double *vnlMatrix)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int i,l,m,i_shift,ipart,is,iii,ioff,lp1;
  int iatm;
  int ind_loc,nl_max,irad,jrad;
  int nl_chan_max;
  double rvol_cp,vol_cp,cp_enl;
  double p11,p22,p33,p12,p13,p23;

/* Local pointers */
  int npart                = clatoms_info->natm_tot;
  double *fx               = clatoms_pos->fx;
  double *fy               = clatoms_pos->fy;
  double *fz               = clatoms_pos->fz;
  double *hess_xx          = clatoms_pos->hess_xx;
  double *hess_xy          = clatoms_pos->hess_xy;
  double *hess_xz          = clatoms_pos->hess_xz;
  double *hess_yy          = clatoms_pos->hess_yy;
  double *hess_yz          = clatoms_pos->hess_yz;
  double *hess_zz          = clatoms_pos->hess_zz;
  int natm_typ             = atommaps->natm_typ;
  int *iatm_typ            = atommaps->iatm_atm_typ;
  int *iatm_typ_nl         = atommaps->iatm_atm_typ_nl;
  double *hmat_cp          = cell->hmat_cp;

  int cp_ptens             = cpopts->cp_ptens_calc;
  int cp_lsda              = cpopts->cp_lsda;
  int atm_hess_calc        = clatoms_info->hess_calc;
  int nstate_up            = cpcoeffs_info->nstate_up_proc;
  int nstate_dn            = cpcoeffs_info->nstate_dn_proc;

  int np_nonloc_cp_box_kb  = pseudo->np_nonloc_cp_box_kb;
  double *vpsnorm          = pseudo->vpsnorm;
  int n_ang_max            = pseudo->n_ang_max;
  int n_ang_max_kb         = pseudo->n_ang_max_kb;
  int n_rad_max            = pseudo->n_rad_max;
  int *loc_opt             = pseudo->loc_opt;
  int *ip_nl               = pseudo->ip_nl;
  int *ip_nl_rev           = pseudo->ip_nl_rev;
  int *np_nl               = pseudo->np_nl;
  int *nrad_max_l          = pseudo->nrad_max_l;
  int **np_nl_rad_str      = pseudo->np_nl_rad_str;
  double *vnlreal_up       = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up       = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn       = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn       = cpscr->cpscr_nonloc.vnlim_dn;
  double *dvnlreal_x_up    = cpscr->cpscr_nonloc.dvnlre_x_up;
  double *dvnlreal_y_up    = cpscr->cpscr_nonloc.dvnlre_y_up;
  double *dvnlreal_z_up    = cpscr->cpscr_nonloc.dvnlre_z_up;
  double *dvnlimag_x_up    = cpscr->cpscr_nonloc.dvnlim_x_up;
  double *dvnlimag_y_up    = cpscr->cpscr_nonloc.dvnlim_y_up;
  double *dvnlimag_z_up    = cpscr->cpscr_nonloc.dvnlim_z_up;
  double *dvnlreal_x_dn    = cpscr->cpscr_nonloc.dvnlre_x_dn;
  double *dvnlreal_y_dn    = cpscr->cpscr_nonloc.dvnlre_y_dn;
  double *dvnlreal_z_dn    = cpscr->cpscr_nonloc.dvnlre_z_dn;
  double *dvnlimag_x_dn    = cpscr->cpscr_nonloc.dvnlim_x_dn;
  double *dvnlimag_y_dn    = cpscr->cpscr_nonloc.dvnlim_y_dn;
  double *dvnlimag_z_dn    = cpscr->cpscr_nonloc.dvnlim_z_dn;
  double *dvnlreal_gxgx_up = cpscr->cpscr_nonloc.dvnlre_gxgx_up;
  double *dvnlimag_gxgx_up = cpscr->cpscr_nonloc.dvnlim_gxgx_up;
  double *dvnlreal_gygy_up = cpscr->cpscr_nonloc.dvnlre_gygy_up;
  double *dvnlimag_gygy_up = cpscr->cpscr_nonloc.dvnlim_gygy_up;
  double *dvnlreal_gzgz_up = cpscr->cpscr_nonloc.dvnlre_gzgz_up;
  double *dvnlimag_gzgz_up = cpscr->cpscr_nonloc.dvnlim_gzgz_up;
  double *dvnlreal_gxgy_up = cpscr->cpscr_nonloc.dvnlre_gxgy_up;
  double *dvnlimag_gxgy_up = cpscr->cpscr_nonloc.dvnlim_gxgy_up;
  double *dvnlreal_gygz_up = cpscr->cpscr_nonloc.dvnlre_gygz_up;
  double *dvnlimag_gygz_up = cpscr->cpscr_nonloc.dvnlim_gygz_up;
  double *dvnlreal_gxgz_up = cpscr->cpscr_nonloc.dvnlre_gxgz_up;
  double *dvnlimag_gxgz_up = cpscr->cpscr_nonloc.dvnlim_gxgz_up;
  double *dvnlreal_gxgx_dn = cpscr->cpscr_nonloc.dvnlre_gxgx_dn;
  double *dvnlimag_gxgx_dn = cpscr->cpscr_nonloc.dvnlim_gxgx_dn;
  double *dvnlreal_gygy_dn = cpscr->cpscr_nonloc.dvnlre_gygy_dn;
  double *dvnlimag_gygy_dn = cpscr->cpscr_nonloc.dvnlim_gygy_dn;
  double *dvnlreal_gzgz_dn = cpscr->cpscr_nonloc.dvnlre_gzgz_dn;
  double *dvnlimag_gzgz_dn = cpscr->cpscr_nonloc.dvnlim_gzgz_dn;
  double *dvnlreal_gxgy_dn = cpscr->cpscr_nonloc.dvnlre_gxgy_dn;
  double *dvnlimag_gxgy_dn = cpscr->cpscr_nonloc.dvnlim_gxgy_dn;
  double *dvnlreal_gygz_dn = cpscr->cpscr_nonloc.dvnlre_gygz_dn;
  double *dvnlimag_gygz_dn = cpscr->cpscr_nonloc.dvnlim_gygz_dn;
  double *dvnlreal_gxgz_dn = cpscr->cpscr_nonloc.dvnlre_gxgz_dn;
  double *dvnlimag_gxgz_dn = cpscr->cpscr_nonloc.dvnlim_gxgz_dn;

  double *fxtemp           = ewd_scr->fx;
  double *fytemp           = ewd_scr->fy;
  double *fztemp           = ewd_scr->fz;
  double *vscr             = ewd_scr->fx2;
  double *vnorm            = ewd_scr->fy2;
  double *vnorm_now        = ewd_scr->fz2;

/*======================================================================*/
/* I) Useful constants                                                  */

  vol_cp   = getdeth(hmat_cp);
  rvol_cp  = 1.0/vol_cp;
  cp_enl   = 0.0;
  nl_max   = -1;
  for(i=1;i<=(n_ang_max_kb+1);i++){
   if(np_nl[i]>0){nl_max=i-1;}
  }//endfor
  nl_chan_max = (nl_max+1)*(nl_max+1);

/*======================================================================*/
/* II) Loop over the open channels, the states and get the nl potent,   */
/*     pvten and  particle forces                                       */

  for(l=0;l<=nl_max;l++){
    lp1 = l+1;
    if(np_nl[lp1]>0){
      for(irad=1;irad<=nrad_max_l[lp1];irad++){
        for(jrad=irad;jrad<=nrad_max_l[lp1];jrad++){
/*-----------------------------------------------------------------------*/
/* i) Get the normalization scaled by the volume                        */
	  get_vpsnorm(vscr,vpsnorm,vnorm,iatm_typ_nl,natm_typ,np_nonloc_cp_box_kb,l,
                     n_ang_max,irad,jrad,n_rad_max);
	  i_shift = l*npart;
	  for(ipart=np_nl_rad_str[lp1][jrad];ipart<=np_nl[lp1];ipart++){
            vnorm_now[ipart] = vnorm[ip_nl_rev[(ipart+i_shift)]]*rvol_cp;
          }//endfor
/*-----------------------------------------------------------------------*/
/* ii) Sum the contributions over the 2l+1 directions and the states    */
	  sumnl_pot_pv_fatm_hess_frag(npart,nstate_up,np_nlmax,nl_chan_max,
				    np_nl[lp1],l,np_nl_rad_str[lp1][jrad],
				    irad,jrad,
				    ip_nl,vnorm_now,vnlreal_up,vnlimag_up,
				    dvnlreal_gxgx_up,dvnlimag_gxgx_up,
				    dvnlreal_gygy_up,dvnlimag_gygy_up,
				    dvnlreal_gzgz_up,dvnlimag_gzgz_up,
				    dvnlreal_gxgy_up,dvnlimag_gxgy_up,
				    dvnlreal_gxgz_up,dvnlimag_gxgz_up,
				    dvnlreal_gygz_up,dvnlimag_gygz_up,
				    dvnlreal_x_up,dvnlimag_x_up,
				    dvnlreal_y_up,dvnlimag_y_up,
				    dvnlreal_z_up,dvnlimag_z_up,
				    fx,fy,fz,fxtemp,fytemp,fztemp,
				    hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
				    atm_hess_calc,cp_ptens,pvten,&cp_enl);
          if(cp_lsda==1){
            sumnl_pot_pv_fatm_hess_frag(npart,nstate_dn,np_nlmax,nl_chan_max,
				      np_nl[lp1],l,np_nl_rad_str[lp1][jrad],
				      irad,jrad,
				      ip_nl,vnorm_now,vnlreal_dn,vnlimag_dn,
				      dvnlreal_gxgx_dn,dvnlimag_gxgx_dn,
				      dvnlreal_gygy_dn,dvnlimag_gygy_dn,
				      dvnlreal_gzgz_dn,dvnlimag_gzgz_dn,
				      dvnlreal_gxgy_dn,dvnlimag_gxgy_dn,
				      dvnlreal_gxgz_dn,dvnlimag_gxgz_dn,
				      dvnlreal_gygz_dn,dvnlimag_gygz_dn,
				      dvnlreal_x_dn,dvnlimag_x_dn,
				      dvnlreal_y_dn,dvnlimag_y_dn,
				      dvnlreal_z_dn,dvnlimag_z_dn,
				      fx,fy,fz,fxtemp,fytemp,fztemp,
				      hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
				      atm_hess_calc,cp_ptens,pvten,&cp_enl);
          }//endif
	}//endfor radial channels jrad
      }//endfor: radial channels irad
    }//endif: l channel open 
  }//endfor: l channels

/*======================================================================*/
/* III) Assign the non-local energy  and add it to the pvten            */

  *cp_enl_ret = cp_enl;
  if(cp_ptens==1){
    pvten[1] += cp_enl;
    pvten[5] += cp_enl;
    pvten[9] += cp_enl;
  }/*endif*/

/*======================================================================*/
  }/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void sumnl_pot_pv_fatm_hess_frag(int npart,int nstate,int np_nlmax,
				int nl_chan_max,int np_nl,int l,int np_nl_rad_str,
				int irad, int jrad,
				int *ip_nl,double *vnorm_now,
				double *vnlreal,double *vnlimag,
				double *dvnlreal_gxgx,double *dvnlimag_gxgx,
				double *dvnlreal_gygy,double *dvnlimag_gygy,
				double *dvnlreal_gzgz,double *dvnlimag_gzgz,
				double *dvnlreal_gxgy,double *dvnlimag_gxgy,
				double *dvnlreal_gxgz,double *dvnlimag_gxgz,
				double *dvnlreal_gygz,double *dvnlimag_gygz,
				double *dvnlreal_x,double *dvnlimag_x,
				double *dvnlreal_y,double *dvnlimag_y,
				double *dvnlreal_z,double *dvnlimag_z,
				double *fx,double *fy,double *fz,
				double *fxtemp,double *fytemp,double *fztemp,
				double *hess_xx,double *hess_xy,double *hess_xz,
				double *hess_yy,double *hess_yz,double *hess_zz,
				int atm_hess_calc,
				int cp_ptens,double *pvten, double *cp_enl_ret,
				double *vnlFxMatrix,double *vnlFyMatrix,
				double *vnlFzMatrix,double *vnlMatrix)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
  int i,m,i_shift,ipart,is,js,iii,ioff,ltemp;
  int ind_loc,ind_lm;
  int ioff_i,ioff_j,ind_loc_i,ind_loc_j;
  int ioff_i_is,ioff_j_is,ioff_i_js,ioff_j_js;
  int ind_loc_i_is,ind_loc_i_js,ind_loc_j_is,ind_loc_j_js;
  int hess_ind;
  int ind_force_mat_1,ind_force_mat_2;
  double p11,p22,p33,p12,p13,p23,cp_enl_now;
  double cp_enl = *cp_enl_ret;
/*==========================================================================*/
/* I) Loop over the 2*l+1 directions and sum the nl contributions           */
    for(m = 1;m<=(2*l+1);m++){
      ind_lm = m + l*l;
      for(is=1;is<=nstate;is++){
	for(js=is;js<=nstate;js++){
/*----------------------------------------------------------------------*/
/*  i) Get the contrib to the non-local energy                          */
	  cp_enl_now = 0.0;
	  ioff_i_is = (is-1)*np_nlmax+(ind_lm-1)*nstate*np_nlmax
                   +(irad-1)*nl_chan_max*nstate*np_nlmax;
          ioff_j_is = (is-1)*np_nlmax+(ind_lm-1)*nstate*np_nlmax
                   +(jrad-1)*nl_chan_max*nstate*np_nlmax;
          ioff_i_js = (js-1)*np_nlmax+(ind_lm-1)*nstate*np_nlmax
                   +(irad-1)*nl_chan_max*nstate*np_nlmax;
          ioff_j_js = (js-1)*np_nlmax+(ind_lm-1)*nstate*np_nlmax
                   +(jrad-1)*nl_chan_max*nstate*np_nlmax;
	  if(irad==jrad){
            for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
              ind_loc_i_is = ipart + ioff_i_is;
	      ind_loc_i_js = ipart + ioff_i_js;
              cp_enl_now = (vnlreal[ind_loc_i_is]*vnlreal[ind_loc_i_js]
                           +vnlimag[ind_loc_i_is]*vnlimag[ind_loc_i_js])*vnorm_now[ipart];
	      vnlMatrix[(is-1)*nstate+js-1] += cp_enl_now;
	      vnlMatrix[(js-1)*nstate+is-1] += cp_enl_now;
              if(is==js)cp_enl += cp_enl_now;
            }//endfor
	  }
	  else{
            for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
              ind_loc_i_is = ipart + ioff_i_is;
              ind_loc_i_js = ipart + ioff_i_js;
              ind_loc_j_is = ipart + ioff_j_is;
              ind_loc_j_js = ipart + ioff_j_js;
              cp_enl_now = 2.0*(vnlreal[ind_loc_i_is]*vnlreal[ind_loc_j_js]
                            +vnlimag[ind_loc_i_is]*vnlimag[ind_loc_j_js])
                            *vnorm_now[ipart];
	      vnlMatrix[(is-1)*nstate+js-1] += cp_enl_now;
	      vnlMatrix[(js-1)*nstate+is-1] += cp_enl_now;
	      if(is==js)cp_enl += cp_enl_now;
            }//endfor
	  }
	  /*
	  cp_enl_now = 0.0;
	  ioff_i = (is-1)*np_nlmax+(ind_lm-1)*nstate*np_nlmax
		 +(irad-1)*nl_chan_max*nstate*np_nlmax;
	  ioff_j = (is-1)*np_nlmax+(ind_lm-1)*nstate*np_nlmax
		 +(jrad-1)*nl_chan_max*nstate*np_nlmax;
	  if(irad==jrad){
	    for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
	      ind_loc = ipart + ioff_i;
	      cp_enl_now = (vnlreal[ind_loc]*vnlreal[ind_loc]
			   +vnlimag[ind_loc]*vnlimag[ind_loc])*vnorm_now[ipart];
	      cp_enl += cp_enl_now;
	    }//endfor
	  }
	  else{
	    for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
	      ind_loc_i = ipart + ioff_i;
	      ind_loc_j = ipart + ioff_j;
	      cp_enl += 2.0*(vnlreal[ind_loc_i]*vnlreal[ind_loc_j]
			   +vnlimag[ind_loc_i]*vnlimag[ind_loc_j])
			   *vnorm_now[ipart];

	    }//endfor
	  }//endif
	  */

/*----------------------------------------------------------------------*/
/* ii) Get the contrib to non-local piece of the pressure tensor        */
	/* //Turn on ptens in the future
        if(cp_ptens==1){
         if(irad==jrad){
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc = ipart + ioff_i;
            p11 = 2.0*(dvnlreal_gxgx[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gxgx[ind_loc]*vnlimag[ind_loc])
                     *vnorm_now[ipart];
            p22 = 2.0*(dvnlreal_gygy[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gygy[ind_loc]*vnlimag[ind_loc])
                     *vnorm_now[ipart];
            p33 = 2.0*(dvnlreal_gzgz[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gzgz[ind_loc]*vnlimag[ind_loc])
                      *vnorm_now[ipart];
            p12 = 2.0*(dvnlreal_gxgy[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gxgy[ind_loc]*vnlimag[ind_loc])
                      *vnorm_now[ipart];
            p13 = 2.0*(dvnlreal_gxgz[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gxgz[ind_loc]*vnlimag[ind_loc])
                      *vnorm_now[ipart];
            p23 = 2.0*(dvnlreal_gygz[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gygz[ind_loc]*vnlimag[ind_loc])
                     *vnorm_now[ipart];
            pvten[1] += p11;
            pvten[5] += p22;
            pvten[9] += p33;
            pvten[4] += p12;
            pvten[2] += p12;
            pvten[7] += p13;
            pvten[3] += p13;
            pvten[8] += p23;
            pvten[6] += p23;
          }//endfor:ipart
         }else{
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc_i = ipart + ioff_i;
            ind_loc_j = ipart + ioff_j;
            p11 = 2.0*(dvnlreal_gxgx[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gxgx[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gxgx[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gxgx[ind_loc_j]*vnlimag[ind_loc_i])
                     *vnorm_now[ipart];
            p22 = 2.0*(dvnlreal_gygy[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gygy[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gygy[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gygy[ind_loc_j]*vnlimag[ind_loc_i])
                     *vnorm_now[ipart];
            p33 = 2.0*(dvnlreal_gzgz[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gzgz[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gzgz[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gzgz[ind_loc_j]*vnlimag[ind_loc_i])
                      *vnorm_now[ipart];
            p12 = 2.0*(dvnlreal_gxgy[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gxgy[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gxgy[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gxgy[ind_loc_j]*vnlimag[ind_loc_i])
                      *vnorm_now[ipart];
            p13 = 2.0*(dvnlreal_gxgz[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gxgz[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gxgz[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gxgz[ind_loc_j]*vnlimag[ind_loc_i])
                      *vnorm_now[ipart];
            p23 = 2.0*(dvnlreal_gygz[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gygz[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gygz[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gygz[ind_loc_j]*vnlimag[ind_loc_i])
                     *vnorm_now[ipart];
            pvten[1] += p11;
            pvten[5] += p22;
            pvten[9] += p33;
            pvten[4] += p12;
            pvten[2] += p12;
            pvten[7] += p13;
            pvten[3] += p13;
            pvten[8] += p23;
            pvten[6] += p23;
          }//endfor:ipart
         }//endif
        }//endif:cp_ptens on
	*/
/*----------------------------------------------------------------------*/
/* iii) Sum the non-local particle force                                  */

          if(irad==jrad){
            for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
              ind_loc_i_is = ipart + ioff_i_is;
              ind_loc_i_js = ipart + ioff_i_js;
              fxtemp[ipart] = -2.0*(dvnlreal_x[ind_loc_i_is]*vnlreal[ind_loc_i_js]
                               +dvnlimag_x[ind_loc_i_is]*vnlimag[ind_loc_i_js])
                               *vnorm_now[ipart];
              fytemp[ipart] = -2.0*(dvnlreal_y[ind_loc_i_is]*vnlreal[ind_loc_i_js]
                               +dvnlimag_y[ind_loc_i_is]*vnlimag[ind_loc_i_js])
                               *vnorm_now[ipart];
              fztemp[ipart] = -2.0*(dvnlreal_z[ind_loc_i_is]*vnlreal[ind_loc_i_js]
                               +dvnlimag_z[ind_loc_i_is]*vnlimag[ind_loc_i_js])
                               *vnorm_now[ipart];
            }//endfor
          }
          else{
            for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
              ind_loc_i_is = ipart + ioff_i_is;
              ind_loc_i_js = ipart + ioff_i_js;
              ind_loc_j_is = ipart + ioff_j_is;
              ind_loc_j_js = ipart + ioff_j_js;
              fxtemp[ipart] = -2.0*(dvnlreal_x[ind_loc_i_is]*vnlreal[ind_loc_j_js]
                                +dvnlreal_x[ind_loc_j_is]*vnlreal[ind_loc_i_js]
                                +dvnlimag_x[ind_loc_i_is]*vnlimag[ind_loc_j_js]
                                +dvnlimag_x[ind_loc_j_is]*vnlimag[ind_loc_i_js])
                                *vnorm_now[ipart];
              fytemp[ipart] = -2.0*(dvnlreal_y[ind_loc_i_is]*vnlreal[ind_loc_j_js]
                                +dvnlreal_y[ind_loc_j_is]*vnlreal[ind_loc_i_js]
                                +dvnlimag_y[ind_loc_i_is]*vnlimag[ind_loc_j_js]
                                +dvnlimag_y[ind_loc_j_is]*vnlimag[ind_loc_i_js])
                                *vnorm_now[ipart];
              fztemp[ipart] = -2.0*(dvnlreal_z[ind_loc_i_is]*vnlreal[ind_loc_j_js]
                                +dvnlreal_z[ind_loc_j_is]*vnlreal[ind_loc_i_js]
                                +dvnlimag_z[ind_loc_i_is]*vnlimag[ind_loc_j_js]
                                +dvnlimag_z[ind_loc_j_is]*vnlimag[ind_loc_i_js])
                               *vnorm_now[ipart];
            }//endfor
          }//endif
          i_shift = l*npart;
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ltemp = ip_nl[(ipart+i_shift)];
	    if(js==is){
              fx[ltemp] += fxtemp[ipart];
              fy[ltemp] += fytemp[ipart];
              fz[ltemp] += fztemp[ipart];
	    }
	    ind_force_mat_1 = ltemp*nstate*nstate+(is-1)*nstate+js-1;
	    ind_force_mat_2 = ltemp*nstate*nstate+(js-1)*nstate+is-1;
	    vnlFxMatrix[ind_force_mat_1] += fxtemp[ipart];
	    vnlFxMatrix[ind_force_mat_2] += fxtemp[ipart];
            vnlFyMatrix[ind_force_mat_1] += fytemp[ipart];
            vnlFyMatrix[ind_force_mat_2] += fytemp[ipart];
            vnlFzMatrix[ind_force_mat_1] += fztemp[ipart];
            vnlFzMatrix[ind_force_mat_2] += fztemp[ipart];
          }/*endfor:atomic forces*/

	  /*     
	  if(irad==jrad){
            for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
	      ind_loc = ipart + ioff_i;
	      fxtemp[ipart] = -2.0*(dvnlreal_x[ind_loc]*vnlreal[ind_loc]
                               +dvnlimag_x[ind_loc]*vnlimag[ind_loc])
                               *vnorm_now[ipart];
	      fytemp[ipart] = -2.0*(dvnlreal_y[ind_loc]*vnlreal[ind_loc]
                               +dvnlimag_y[ind_loc]*vnlimag[ind_loc])
                               *vnorm_now[ipart];
	      fztemp[ipart] = -2.0*(dvnlreal_z[ind_loc]*vnlreal[ind_loc]
                               +dvnlimag_z[ind_loc]*vnlimag[ind_loc])
                               *vnorm_now[ipart];
            }//endfor
	  }
	  else{
            for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
              ind_loc_i = ipart + ioff_i;
              ind_loc_j = ipart + ioff_j;
	      fxtemp[ipart] = -2.0*(dvnlreal_x[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_x[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_x[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_x[ind_loc_j]*vnlimag[ind_loc_i])
                                *vnorm_now[ipart];
	      fytemp[ipart] = -2.0*(dvnlreal_y[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_y[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_y[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_y[ind_loc_j]*vnlimag[ind_loc_i])
                                *vnorm_now[ipart];
	      fztemp[ipart] = -2.0*(dvnlreal_z[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_z[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_z[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_z[ind_loc_j]*vnlimag[ind_loc_i])
                               *vnorm_now[ipart];
            }//endfor
          }//endif
          i_shift = l*npart;
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ltemp = ip_nl[(ipart+i_shift)];
            fx[ltemp] += fxtemp[ipart];
            fy[ltemp] += fytemp[ipart];
            fz[ltemp] += fztemp[ipart];
          }//endfor:atomic forces
	  */

/*----------------------------------------------------------------------*/
/* iv) Sum the non-local atomic hessian                                  */

        /* // Hessian comes later
        if(atm_hess_calc == 3){
         if(irad==jrad){
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc = ipart + ioff_i;
            fxtemp[ipart] = 2.0*(dvnlreal_gxgx[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gxgx[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_x[ind_loc]*dvnlreal_x[ind_loc]
                                +dvnlimag_x[ind_loc]*dvnlimag_x[ind_loc])
                                *vnorm_now[ipart];
            fytemp[ipart] = 2.0*(dvnlreal_gxgy[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gxgy[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_x[ind_loc]*dvnlreal_y[ind_loc]
                                +dvnlimag_x[ind_loc]*dvnlimag_y[ind_loc])
                                *vnorm_now[ipart];
            fztemp[ipart] = 2.0*(dvnlreal_gxgz[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gxgz[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_x[ind_loc]*dvnlreal_z[ind_loc]
                                +dvnlimag_x[ind_loc]*dvnlimag_z[ind_loc])
                                *vnorm_now[ipart];
          }// endfor ipart
         } else {
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc_i = ipart + ioff_i;
            ind_loc_j = ipart + ioff_j;
            fxtemp[ipart] = 2.0*(dvnlreal_gxgx[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gxgx[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gxgx[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gxgx[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_x[ind_loc_i]*dvnlreal_x[ind_loc_j]
                                     +dvnlimag_x[ind_loc_i]*dvnlimag_x[ind_loc_j]))
                                *vnorm_now[ipart];
            fytemp[ipart] = 2.0*(dvnlreal_gxgy[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gxgy[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gxgy[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gxgy[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_x[ind_loc_i]*dvnlreal_y[ind_loc_j]
                                     +dvnlimag_x[ind_loc_i]*dvnlimag_y[ind_loc_j]))
                                *vnorm_now[ipart];
            fztemp[ipart] = 2.0*(dvnlreal_gxgz[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gxgz[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gxgz[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gxgz[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_x[ind_loc_i]*dvnlreal_z[ind_loc_j]
                                     +dvnlimag_x[ind_loc_i]*dvnlimag_z[ind_loc_j]))
                                *vnorm_now[ipart];
          }// endfor
         }// endif
         i_shift = l*npart;
         for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
           ltemp = ip_nl[(ipart+i_shift)];
           hess_ind = (ltemp-1)*npart + ltemp;
           hess_xx[hess_ind] += fxtemp[ipart];
           hess_xy[hess_ind] += fytemp[ipart];
           hess_xz[hess_ind] += fztemp[ipart];
         }//endfor:atomic forces

         if(irad==jrad){
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc = ipart + ioff_i;
            fxtemp[ipart] = 2.0*(dvnlreal_gygy[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gygy[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_y[ind_loc]*dvnlreal_y[ind_loc]
                                +dvnlimag_y[ind_loc]*dvnlimag_y[ind_loc])
                                *vnorm_now[ipart];
            fytemp[ipart] = 2.0*(dvnlreal_gygz[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gygz[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_y[ind_loc]*dvnlreal_z[ind_loc]
                                +dvnlimag_y[ind_loc]*dvnlimag_z[ind_loc])
                                *vnorm_now[ipart];
            fztemp[ipart] = 2.0*(dvnlreal_gzgz[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gzgz[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_z[ind_loc]*dvnlreal_z[ind_loc]
                                +dvnlimag_z[ind_loc]*dvnlimag_z[ind_loc])
                                *vnorm_now[ipart];
          }// endfor ipart
         } else {
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc_i = ipart + ioff_i;
            ind_loc_j = ipart + ioff_j;
            fxtemp[ipart] = 2.0*(dvnlreal_gygy[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gygy[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gygy[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gygy[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_y[ind_loc_i]*dvnlreal_y[ind_loc_j]
                                     +dvnlimag_y[ind_loc_i]*dvnlimag_y[ind_loc_j]))
                                *vnorm_now[ipart];
            fytemp[ipart] = 2.0*(dvnlreal_gygz[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gygz[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gygz[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gygz[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_y[ind_loc_i]*dvnlreal_z[ind_loc_j]
                                     +dvnlimag_y[ind_loc_i]*dvnlimag_z[ind_loc_j]))
                                *vnorm_now[ipart];
            fztemp[ipart] = 2.0*(dvnlreal_gzgz[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gzgz[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gzgz[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gzgz[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_z[ind_loc_i]*dvnlreal_z[ind_loc_j]
                                     +dvnlimag_z[ind_loc_i]*dvnlimag_z[ind_loc_j]))
                                *vnorm_now[ipart];
          }// endfor ipart
         }// endif
         i_shift = l*npart;
         for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
           ltemp = ip_nl[(ipart+i_shift)];
           hess_ind = (ltemp-1)*npart + ltemp;
           hess_yy[hess_ind] += fxtemp[ipart];
           hess_yz[hess_ind] += fytemp[ipart];
           hess_zz[hess_ind] += fztemp[ipart];
         }//endfor:atomic forces 

        }// endif: atm hess calc
	*/

      }/*endfor:loop over states*/
    }/*endfor: loop over the m channels */

/*==========================================================================*/
/* II) Set the return values                                               */

 *cp_enl_ret = cp_enl;

/*======================================================================*/
  }/*end routine*/
/*======================================================================*/

