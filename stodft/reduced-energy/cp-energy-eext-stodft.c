/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: cp-energy-eext-stodft.c                        */
/*                                                                          */
/* This routine wrapps all functions used within SCF. Nuclei forces are not */
/* calculated.                                                              */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_stodft_local.h"

#include "complex.h"
#define TIME_CP_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void controlEwdLocPreScf(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                         CELL *cell, PTENS *ptens, EWALD *ewald, CPEWALD *cpewald, 
                         CPSCR *cpscr, PSEUDO *pseudo, EWD_SCR *ewd_scr,  
                         CPOPTS *cpopts, ATOMMAPS *atommaps, double *vrecip_ret,
                         double *pseud_hess_loc,COMMUNICATE *communicate,
                         FOR_SCR *for_scr,int cp_dual_grid_opt,int idual_switch)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* KS potential from local pseudo pp can be calculated before SCF	 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
     
#include "../typ_defs/typ_mask.h"

  int idens_opt,ipseud_opt;
  int istart,ngo,irem,idiv;
  int ipart,jpart,iii,itype,i;
  int icount,koff,natm_use;
  int hess_ind;
  int realSparseOpt = cpewald->realSparseOpt;

  double falp2,falp_clus2,vol,rvol,pivol,fpi,arg,q_sum1,bgr;
  double aka,akb,akc,xk,yk,zk,atemp,btemp,ctemp;
  double xtemp,ytemp,ztemp;
  double sumr,sumi,g2,g4,preg,prep,tpi,pi,g;
  double sumr_h,sumi_h;
  double srx,sry,srz,six,siy,siz,temp,smag;
  double vrecip;
  double dx,dy,dz;
  double sx,sy,sz;
  double asx,asy,asz;
  double phase;
  double argp,fargp,argm,fargm,area;

/*--------------------------------------------*/
/*         Local Pointer declarations         */

  /*------------------*/
  /* Atom information */
  int natm_tot      = clatoms_info->natm_tot; 
  int hess_calc     = clatoms_info->hess_calc;
  double *x         = clatoms_pos->x;
  double *y         = clatoms_pos->y;
  double *z         = clatoms_pos->z;
  double *fx        = clatoms_pos->fx;     
  double *fy        = clatoms_pos->fy;
  double *fz        = clatoms_pos->fz;
  double *fx_tmp    = ewd_scr->fx2;     
  double *fy_tmp    = ewd_scr->fy2;     
  double *fz_tmp    = ewd_scr->fz2;     
  double *q         = clatoms_info->q;
  double *hess_xx   = clatoms_pos->hess_xx;
  double *hess_xy   = clatoms_pos->hess_xy;
  double *hess_xz   = clatoms_pos->hess_xz;
  double *hess_yy   = clatoms_pos->hess_yy;
  double *hess_yz   = clatoms_pos->hess_yz;
  double *hess_zz   = clatoms_pos->hess_zz;
  int natm_typ      = atommaps->natm_typ;
  int *index_atm    = for_scr->index_atm;
  int *iatm_typ;            /*Assigned below based on flags */
  int *iatm_typ_full;

  /*--------------------------------*/
  /* Cell and pressure information */
  int iperd                 = cell->iperd;
  int cp_ptens              = cpopts->cp_ptens_calc;
  double *pvten             = ptens->pvten_tmp;
  double *cp_box_center     = cell->cp_box_center;
  double *cp_box_center_rel = cell->cp_box_center_rel;
  double *hmat;             /* Assigned below based on flags */
  double *hmati;
  double *hmat_big          = cell->hmat;
  double *hmati_big         = cell->hmati;

  /*----------------------*/
  /* G-vector information */
  int *kastore;             /* Assigned below based on flags */
  int *kbstore;
  int *kcstore;
  int *ibreak1;
  int *ibreak2;
  double *vextr;
  double *vexti;
  double *vextr_loc;
  double *vexti_loc;
  double *dvextr;
  double *dvexti;
  double *rhocr;
  double *rhoci;
  double *ak2;
  int nktot;

  /*----------------------------------------------*/
  /* Pseudo-potential and Reduced Periodicity info*/
  int nsplin_g         = pseudo->nsplin_g;
  int n_rad_max        = pseudo->n_rad_max;
  double *clus_corr_r  = ewald->clus_corr_r;
  double *dclus_corr_r = ewald->dclus_corr_r;
  double alpha_conv_dual = pseudo->alpha_conv_dual;
  double dg_spl        = pseudo->dg_spl;
  double gmin_spl      = pseudo->gmin_spl;
  double *vps0         = pseudo->vps0;
  double *vps1         = pseudo->vps1;
  double *vps2         = pseudo->vps2;
  double *vps3         = pseudo->vps3;
  double *dvps0        = pseudo->dvps0;
  double *dvps1        = pseudo->dvps1;
  double *dvps2        = pseudo->dvps2;
  double *dvps3        = pseudo->dvps3;
  double *gzvps        = pseudo->gzvps;
  double *q_pseud      = pseudo->q_pseud;
  int n_ang_max        = pseudo->n_ang_max;
  int *loc_opt         = pseudo->loc_opt;
  int np_loc_cp_box    = pseudo->np_loc_cp_box;
  int *ip_loc_cp_box   = pseudo->ip_loc_cp_box;

  /*---------------------------------*/
  /* Ewald and ewald scr information */
  double alp_ewald  = ewald->alp_ewd;
  double alp_clus   = ewald->alp_clus;
  double *cossc     = ewd_scr->cossc;  
  double *sinsc     = ewd_scr->sinsc;
  double *helr      = ewd_scr->helr;   
  double *heli      = ewd_scr->heli;
  double *vtemp     = ewd_scr->temp;
  double *dvtemp    = ewd_scr->vtemp_now;
  double *ewd_scr_x = ewd_scr->x;
  double *ewd_scr_y = ewd_scr->y;
  double *ewd_scr_z = ewd_scr->z;
  double *q_tmp     = ewd_scr->q;

  /*---------------------------*/
  /* Communication information */
  int myid_state    = communicate->myid_state;
  int np_states     = communicate->np_states;
  MPI_Comm comm     = communicate->comm_states;
  
/*======================================================================*/
/* 0) Assign local pointers                                             */

  if(cp_dual_grid_opt < 2 || idual_switch == 0){
    /* large sparse grid when cp_dual_grid_opt == 2*/
    idens_opt = 0;
    ipseud_opt= (cp_dual_grid_opt==2 ? 0 : 1);

    //if(realSparseOpt==0){
    kastore   = ewald->kastr;
    kbstore   = ewald->kbstr;
    kcstore   = ewald->kcstr;
    ak2       = cpewald->ak2;
    nktot     = ewald->nktot;
    ibreak1   = ewald->ibrk1;
    ibreak2   = ewald->ibrk2;
    //printf("nktot %i\n",nktot);
    //}
    /*
    else{
      kastore = cpewald->kastr_sm;
      kbstore = cpewald->kbstr_sm;
      kcstore = cpewald->kcstr_sm;
      ak2 = cpewald->ak2_sm;
      nktot = cpewald->nktot_sm;
      ibreak1   = cpewald->ibrk1_sm;
      ibreak2   = cpewald->ibrk2_sm;
      printf("nktot %i\n",nktot);
    }
    */
    /*
    kastore   = ewald->kastr;
    kbstore   = ewald->kbstr;
    kcstore   = ewald->kcstr;
    ibreak1   = ewald->ibrk1;
    ibreak2   = ewald->ibrk2;
    */
    vextr     = cpscr->cpscr_loc.vextr;
    vexti     = cpscr->cpscr_loc.vexti;
    vextr_loc = cpscr->cpscr_loc.vextr_loc;
    vexti_loc = cpscr->cpscr_loc.vexti_loc;
    dvextr    = cpscr->cpscr_loc.dvextr;
    dvexti    = cpscr->cpscr_loc.dvexti;
    rhocr     = cpscr->cpscr_rho.rhocr_up;
    rhoci     = cpscr->cpscr_rho.rhoci_up;
    //ak2       = cpewald->ak2;
    //nktot     = ewald->nktot;
    hmat      = cell->hmat;
    hmati     = cell->hmati;
    natm_use  = natm_tot;
    iatm_typ  = atommaps->iatm_atm_typ;
  }else{
    /* small dense grid */
    idens_opt = 1;
    ipseud_opt= 1;
    kastore   = cpewald->kastr_dens_cp_box;
    kbstore   = cpewald->kbstr_dens_cp_box;
    kcstore   = cpewald->kcstr_dens_cp_box;
    ibreak1   = cpewald->ibrk1_dens_cp_box; /*DY edit*/
    ibreak2   = cpewald->ibrk2_dens_cp_box; /*DY edit*/
    vextr     = cpscr->cpscr_loc.vextr_dens_cp_box;
    vexti     = cpscr->cpscr_loc.vexti_dens_cp_box;
    vextr_loc = cpscr->cpscr_loc.vextr_dens_cp_box_loc;
    vexti_loc = cpscr->cpscr_loc.vexti_dens_cp_box_loc;
    rhocr     = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
    rhoci     = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
    ak2       = cpewald->ak2_dens_cp_box;
    nktot     = cpewald->nktot_dens_cp_box;
    hmat      = cell->hmat_cp;
    hmati     = cell->hmati_cp;
    natm_use  = np_loc_cp_box;
    iatm_typ  = for_scr->iexcl;
    iatm_typ_full = atommaps->iatm_atm_typ;
 }/*endif*/

/*======================================================================*/
/* I) Get more useful constants                                         */

  pi    = M_PI;
  tpi   = 2.0*pi;
  fpi   = 4.0*pi;

  vol   = getdeth(hmat);
  rvol  = 1.0/vol;
  pivol = vol/4.0/pi;
  falp2 = 4.0*alp_ewald*alp_ewald;
  falp_clus2 = 4.0*alp_clus*alp_clus;

/*======================================================================*/
/* II) Find cos and sin of sc components of the particles               */
/*    ( hmnati rvec = svec   r=(x,y,z) s=(a,b,c) )                       */

  if(idens_opt==0){
    for(ipart=1;ipart<=natm_use;ipart++){
      xtemp = x[ipart];
      ytemp = y[ipart];
      ztemp = z[ipart];
      q_tmp[ipart]     = q[ipart];
      ewd_scr_x[ipart] = xtemp*hmati[1]+ytemp*hmati[4]+ztemp*hmati[7];
      ewd_scr_y[ipart] = xtemp*hmati[2]+ytemp*hmati[5]+ztemp*hmati[8];
      ewd_scr_z[ipart] = xtemp*hmati[3]+ytemp*hmati[6]+ztemp*hmati[9];
      ctemp            = ewd_scr_z[ipart]*tpi;
      cossc[ipart]     = cos(ctemp);
      sinsc[ipart]     = sin(ctemp);
    }//endfor
  }
  else{
    for(ipart=1;ipart<=natm_use;ipart++){
      iatm_typ[ipart]  = iatm_typ_full[ip_loc_cp_box[ipart]];
      q_tmp[ipart]     = q[ip_loc_cp_box[ipart]];
      dx               = x[ip_loc_cp_box[ipart]]-cp_box_center[1];
      dy               = y[ip_loc_cp_box[ipart]]-cp_box_center[2];
      dz               = z[ip_loc_cp_box[ipart]]-cp_box_center[3];

      asx              = dx*hmati_big[1]+dy*hmati_big[4]+dz*hmati_big[7];
      asy              = dx*hmati_big[2]+dy*hmati_big[5]+dz*hmati_big[8];
      asz              = dx*hmati_big[3]+dy*hmati_big[6]+dz*hmati_big[9];
      sx               = asx-NINT(asx);
      sy               = asy-NINT(asy);
      sz               = asz-NINT(asz);
      dx               = sx*hmat_big[1]+sy*hmat_big[4]+sz*hmat_big[7];
      dy               = sx*hmat_big[2]+sy*hmat_big[5]+sz*hmat_big[8];
      dz               = sx*hmat_big[3]+sy*hmat_big[6]+sz*hmat_big[9];

      xtemp            = dx+cp_box_center_rel[1];
      ytemp            = dy+cp_box_center_rel[2];
      ztemp            = dz+cp_box_center_rel[3];
      ewd_scr_x[ipart] = xtemp*hmati[1]+ytemp*hmati[4]+ztemp*hmati[7];
      ewd_scr_y[ipart] = xtemp*hmati[2]+ytemp*hmati[5]+ztemp*hmati[8];
      ewd_scr_z[ipart] = xtemp*hmati[3]+ytemp*hmati[6]+ztemp*hmati[9];
      ctemp            = ewd_scr_z[ipart]*tpi;
      cossc[ipart]     = cos(ctemp);
      sinsc[ipart]     = sin(ctemp);
    }//endfor
  }//endif

/*======================================================================*/
/*======================================================================*/
/* Perform the ewald sum/ cp-potential calculation                      */

  idiv    = (nktot+1)/np_states;
  irem    = (nktot+1) % np_states;
  ngo     = (myid_state <  irem ? idiv+1 : idiv);
  istart  = (myid_state <= irem ? myid_state*(idiv+1)+1 : 
                                irem*(idiv+1)+1+(myid_state-irem)*idiv);
  koff    = istart-1;
  if(np_states==myid_state+1){ngo--;}

  //debug
#ifdef TEST_FILTER
  for(icount=1;icount<=ngo;icount++){
    aka = (double)(kastore[(icount+koff)]);
    akb = (double)(kbstore[(icount+koff)]);
    akc = (double)(kcstore[(icount+koff)]);
    xk = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
    yk = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
    zk = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
    g2 = xk*xk+yk*yk+zk*zk;
    g4 = g2*g2;
    g  = sqrt(g2);
    ak2[icount] = g2;
  }
  int ngo_temp = ngo;
  ngo = 0;
#endif

  for(icount=1;icount<=ngo;icount++){

/*======================================================================*/
/* I) Get the k vectors                                                 */

    aka = (double)(kastore[(icount+koff)]);
    akb = (double)(kbstore[(icount+koff)]);
    akc = (double)(kcstore[(icount+koff)]);
    xk = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
    yk = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
    zk = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
    g2 = xk*xk+yk*yk+zk*zk;
    g4 = g2*g2;
    g  = sqrt(g2);
    ak2[icount] = g2;

/*======================================================================*/
/* II) If break point number one or you are just starting out calculate */
/*     the helpful vectors                                              */
    if(ibreak1[(icount+koff)]==1||icount==1){
      for(ipart=1;ipart<=natm_use;ipart++){
	atemp = ewd_scr_x[ipart];
	btemp = ewd_scr_y[ipart];
	ctemp = ewd_scr_z[ipart];
	arg = (aka*atemp+akb*btemp+akc*ctemp)*tpi;
	helr[ipart] = cos(arg);
	heli[ipart] = sin(arg);
      }//endfor
    }//endif

/*======================================================================*/
/* III) Get the external potential                                      */
/*               (interaction of electron with particles)               */

    if(ipseud_opt==1){
      for(itype=1;itype<=natm_typ;itype++){
	index_atm[itype] = (itype-1)*nsplin_g*(n_ang_max+1)*n_rad_max
			   +loc_opt[itype]*nsplin_g*n_rad_max;
      }//endfor
      get_vpsnow(index_atm,nsplin_g,gmin_spl,dg_spl,g,
		 vps0,vps1,vps2,vps3,vtemp,iatm_typ,natm_typ,natm_use,1); 
/*-------------------------*/
/* charge correction */
      for(i=1; i<= natm_use; i++){
        vtemp[i] += -fpi*(q_tmp[i]-q_pseud[iatm_typ[i]])/g2;
      }
    }
    else{  //q_temp is q
      get_vpslong(natm_use,vtemp,g2,q_tmp,alpha_conv_dual,pivol);    
    }//endif
 /*----------------------------------------------------------------------*/
 /* Cluster boundary condition correction                                */
    if(iperd!=3&&idens_opt==0){
      for(ipart=1;ipart<=natm_use;ipart++){
        vtemp[ipart] -= q[ipart]*clus_corr_r[icount];
      }// endfor
    }// endif cluster boundary conditions
/*----------------------------------------------------------------------*/
    if(cp_ptens==1&&idens_opt==0){
      if(ipseud_opt == 1){
        get_vpsnow(index_atm,nsplin_g,gmin_spl,dg_spl,g,
	    	   dvps0,dvps1,dvps2,dvps3,dvtemp,iatm_typ,natm_typ,natm_use,1);
        for(i=1;i<=natm_use;i++){
	  dvtemp[i] += 2.0*fpi*(q_tmp[i]-q_pseud[iatm_typ[i]])/(g2*g2);
        }// endfor
     } 
     else{
       get_dvpslong(natm_use,dvtemp,g2,q_tmp,alpha_conv_dual,pivol);
     }//endif ipseud_opt
    }//endif cp_ptens

    vextr[icount]  =  ddot1(natm_use,helr,1,vtemp,1)*rvol;
    vexti[icount]  = -ddot1(natm_use,heli,1,vtemp,1)*rvol;

    if(cp_ptens==1&&idens_opt==0){
     dvextr[icount] =  ddot1(natm_use,helr,1,dvtemp,1)*rvol;
     dvexti[icount] = -ddot1(natm_use,heli,1,dvtemp,1)*rvol;
    }//endif

/*======================================================================*/
/* IV) Get the real and imag parts of the structure factor              */

   if(idens_opt==0){//create charge weighted structure factor
     sumr = ddot1(natm_use,helr,1,q,1); 
     sumi = ddot1(natm_use,heli,1,q,1);
     smag = sumr*sumr+sumi*sumi;
   }//endif

/*======================================================================*/
/* VII) If break point two, increment the helpful vectors                 */

    if(ibreak2[(icount+koff)]==1){
      for(ipart=1;ipart<=natm_use;ipart++){
        temp = helr[ipart];
        helr[ipart] = helr[ipart]*cossc[ipart] - heli[ipart]*sinsc[ipart];
        heli[ipart] = heli[ipart]*cossc[ipart] + temp*sinsc[ipart];
      }//endfor
    }//endif
  }//endfor:icount loop over k vectors
/*======================================================================*/
/* VIII) g=0 term (local pseudopotential) including term for CBCs       */
  
  //debug
#ifdef TEST_FILTER
  ngo = ngo_temp;
#endif

  if((myid_state+1)==np_states){
    if(ipseud_opt==1){
      ak2[(ngo+1)] = 0.0;
      for(ipart=1;ipart<=natm_use;ipart++){
        vtemp[ipart] = gzvps[iatm_typ[ipart]];
      }//endfor
      vextr[(ngo+1)] =  dsum1(natm_use,vtemp,1)*rvol;
      vexti[(ngo+1)] = 0.0;
    }
    else{ //large sparse grid
      vextr[(ngo+1)] = 0.0;
      bgr  = dsum1(natm_use,q_tmp,1); 
      bgr  = bgr*M_PI/(alpha_conv_dual*alpha_conv_dual*vol);
      vextr[(ngo+1)] = bgr; 
    }//endif

    if(iperd!=3&&idens_opt==0) {
      vextr[(ngo+1)] -= dsum1(natm_use,q_tmp,1)*clus_corr_r[(ngo+1)]*rvol;
    }//endif

    *pseud_hess_loc = vextr[(ngo+1)];

    if(iperd>0&&iperd!=3&&idens_opt==0){
      q_sum1 = dsum1(natm_use,q,1);
      if(iperd==2){
         vrecip += 0.5*q_sum1*q_sum1*rvol*(clus_corr_r[(ngo+1)]
                 + M_PI/(alp_clus*alp_clus));
      } 
      else{
        vrecip += 0.5*q_sum1*q_sum1*clus_corr_r[(ngo+1)]*rvol;
      }// endif iperd 
    }//endif
  }//endif

/*======================================================================*/
/* IX) Copy and store vext */
 
  if(myid_state==np_states-1){
    memcpy(&vextr_loc[1],&vextr[1],(ngo+1)*sizeof(double));
    memcpy(&vexti_loc[1],&vexti[1],(ngo+1)*sizeof(double));
  }
  else{
    memcpy(&vextr_loc[1],&vextr[1],ngo*sizeof(double));
    memcpy(&vexti_loc[1],&vexti[1],ngo*sizeof(double));
  }
  
  /*
  for(i=1;i<=ngo+1;i++){ 
    printf("11111111111 vext %lg %lg\n",vextr[i],vexti[i]);
    //printf("11111111111 vexti %lg %lg\n",vexti[i],vexti_loc[i]);
  }
  */
  

/*======================================================================*/
/* IX) Collect the forces */

/*======================================================================*/
/* X) Add in the surface term the dumb way */

  if(np_states > 1) Barrier(comm);

/*======================================================================*/
/* XI) Collect diagonal piece of the pressure tensor and assign vrecip  */

/*======================================================================*/
/* XII) Finally, store the final value of vrecip */

/*======================================================================*/
    }/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void getNlPotPvFatmSCF(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
		       CELL *cell,CPCOEFFS_INFO *cpcoeffs_info,CPSCR *cpscr,
		       EWD_SCR *ewd_scr,CPOPTS *cpopts,
		       PSEUDO *pseudo,ATOMMAPS *atommaps,
		       double *cp_enl_ret, int np_nlmax,double *pvten)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* This is the wrapper to calculate the H|phi> without calculating K-S   */
/* potential. This part comes from control_cp_eext_recip. This part are  */
/* all contributions from non-local pp. Nuclei forces and related terms  */
/* are not calculated here						 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int i,l,m,i_shift,ipart,is,iii,ioff,lp1;
  int iatm;
  int ind_loc,nl_max,irad,jrad;
  int nl_chan_max;
  double rvol_cp,vol_cp,cp_enl;
  double p11,p22,p33,p12,p13,p23;

/* Local pointers */
  int npart                = clatoms_info->natm_tot;
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
	  }//endfor
/*-----------------------------------------------------------------------*/
/* ii) Sum the contributions over the 2l+1 directions and the states    */

	  sumNlPot(npart,nstate_up,np_nlmax,nl_chan_max,
	 	   np_nl[lp1],l,np_nl_rad_str[lp1][jrad],
	 	   irad,jrad,
		   ip_nl,vnorm_now,vnlreal_up,vnlimag_up,&cp_enl);
	  if(cp_lsda==1){
	    sumNlPot(npart,nstate_dn,np_nlmax,nl_chan_max,
		     np_nl[lp1],l,np_nl_rad_str[lp1][jrad],
		     irad,jrad,
		     ip_nl,vnorm_now,vnlreal_dn,vnlimag_dn,&cp_enl);
	  }//endif
        }//endfor radial channels
      }//endfor: radial channels
    }//endif: l channel open
  }//endfor: l channels

/*======================================================================*/
/* III) Assign the non-local energy  and add it to the pvten            */

  *cp_enl_ret = cp_enl;

/*======================================================================*/
  }/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void sumNlPot(int npart,int nstate,int np_nlmax,
              int nl_chan_max,int np_nl,int l,int np_nl_rad_str,
              int irad, int jrad,
              int *ip_nl,double *vnorm_now,
              double *vnlreal,double *vnlimag,double *cp_enl_ret)
/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int i,m,i_shift,ipart,is,iii,ioff,ltemp;
  int ind_loc,ind_lm;
  int ioff_i,ioff_j,ind_loc_i,ind_loc_j;
  int hess_ind;
  double p11,p22,p33,p12,p13,p23,cp_enl_now;
  double cp_enl = *cp_enl_ret;
/*==========================================================================*/
/* I) Loop over the 2*l+1 directions and sum the nl contributions           */
  for(m = 1;m<=(2*l+1);m++){
    ind_lm = m + l*l;
    for(is=1;is<=nstate;is++){
/*----------------------------------------------------------------------*/
/*  i) Get the contrib to the non-local energy                          */
      cp_enl_now = 0.0;
      ioff_i = (is-1)*np_nlmax +
               +(ind_lm-1)*nstate*np_nlmax
               +(irad-1)*nl_chan_max*nstate*np_nlmax;
      ioff_j = (is-1)*np_nlmax +
               +(ind_lm-1)*nstate*np_nlmax
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
    }//endfor is
  }//endfor m
/*==========================================================================*/
/* II) Set the return values                                               */

  *cp_enl_ret = cp_enl;

/*======================================================================*/
  }/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void controlEwdNonlocFilter(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                      CPCOEFFS_INFO *cpcoeffs_info,CPCOEFFS_POS *cpcoeffs_pos,
                      CELL *cell, PTENS *ptens, CPEWALD *cpewald,
                      CPSCR *cpscr, PSEUDO *pseudo, EWD_SCR *ewd_scr,
                      CPOPTS *cpopts, ATOMMAPS *atommaps,
                      COMMUNICATE *communicate,FOR_SCR *for_scr)
/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  double arg;
  double aka,akb,akc,xk,yk,zk,atemp,btemp,ctemp;
  double xtemp,ytemp,ztemp;
  double g2,tpi,pi,g,fpi;
  double dx,dy,dz;
  double sx,sy,sz;
  double asx,asy,asz;

  double temp;
  double xtrans,ytrans,ztrans;
  int nl_chan_max,irad,jrad;
  int ipart,iii,i,iatm;
  int nl_max,np_nlmax;
  int icount;
  YLM_CONS ylm_cons;
  double ylmr[2];

/*         Local Pointer declarations                                   */
  double *creal_up   = cpcoeffs_pos->cre_up;
  double *cimag_up   = cpcoeffs_pos->cim_up;
  double *creal_dn   = cpcoeffs_pos->cre_dn;
  double *cimag_dn   = cpcoeffs_pos->cim_dn;
  int icoef_orth_up  = cpcoeffs_pos->icoef_orth_up;
  int icoef_form_up  = cpcoeffs_pos->icoef_form_up;
  int ifcoef_form_up = cpcoeffs_pos->ifcoef_form_up;
  int icoef_orth_dn  = cpcoeffs_pos->icoef_orth_dn;
  int icoef_form_dn  = cpcoeffs_pos->icoef_form_dn;
  int ifcoef_form_dn = cpcoeffs_pos->ifcoef_form_dn;
  int nstate_up      = cpcoeffs_info->nstate_up_proc;
  int nstate_dn      = cpcoeffs_info->nstate_dn_proc;
  int ncoef          = cpcoeffs_info->ncoef;
  int cp_hess_calc   = cpopts->cp_hess_calc;
  int cp_lsda        = cpopts->cp_lsda;

  int npart          = clatoms_info->natm_tot;
  double *x          = clatoms_pos->x;
  double *y          = clatoms_pos->y;
  double *z          = clatoms_pos->z;
  double *hmat_cp    = cell->hmat_cp;
  double *hmati_cp   = cell->hmati_cp;
  double *hmat_big  = cell->hmat;
  double *hmati_big  = cell->hmati;
  double *cp_box_center     = cell->cp_box_center;
  double *cp_box_center_rel = cell->cp_box_center_rel;
  int natm_typ       = atommaps->natm_typ;
  int *iatm_typ      = atommaps->iatm_atm_typ;
  int *iatm_typ_nl   = atommaps->iatm_atm_typ_nl;
  int *index_atm     = for_scr->iexcl;

  int nktot_sm       = cpewald->nktot_sm;
  int *kastore_sm    = cpewald->kastr_sm;
  int *kbstore_sm    = cpewald->kbstr_sm;
  int *kcstore_sm    = cpewald->kcstr_sm;
  double *ak2_sm     = cpewald->ak2_sm;
  int *ibreak1_sm    = cpewald->ibrk1_sm;
  int *ibreak2_sm    = cpewald->ibrk2_sm;
  double *cossc      = ewd_scr->cossc;
  double *sinsc      = ewd_scr->sinsc;
  double *helr       = ewd_scr->helr;
  double *heli       = ewd_scr->heli;
  double *vtemp      = ewd_scr->fx2;
  double *vtemp_now  = ewd_scr->vtemp_now;
  double *dvtemp     = ewd_scr->q;
  double *ewd_scr_x  = ewd_scr->x;
  double *ewd_scr_y  = ewd_scr->y;
  double *ewd_scr_z  = ewd_scr->z;
  int n_ang_max      = pseudo->n_ang_max;
  int n_ang_max_kb   = pseudo->n_ang_max_kb;
  int n_rad_max      = pseudo->n_rad_max;
  int *ip_nl         = pseudo->ip_nl;
  int *np_nl         = pseudo->np_nl;
  double *gzvps      = pseudo->gzvps;
  double *gzvps0     = pseudo->gzvps0;

  int np_nonloc_cp_box_kb  = pseudo->np_nonloc_cp_box_kb;
  int *nrad_max_l       = pseudo->nrad_max_l;
  int **np_nl_rad_str   = pseudo->np_nl_rad_str;
  double *vnlreal_up = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn = cpscr->cpscr_nonloc.vnlim_dn;

  int myid_state     = communicate->myid_state;
  int np_states      = communicate->np_states;


/*======================================================================*/
/* 0) Check the forms                                                   */

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
  }//endfor
  nl_chan_max = (nl_max + 1)*(nl_max + 1);
/*======================================================================*/
/* III) Determine the maximum number of atoms in any                    */
/*       open angular momentum channel                                  */
  np_nlmax = 1;
  for(i = 1;i<=(nl_max+1);i++)np_nlmax = MAX(np_nlmax,np_nl[i]);
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

   controlNlmatFilter(clatoms_info,cpcoeffs_info,cpcoeffs_pos,
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


/*======================================================================*/
  }/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void controlNlmatFilter(CLATOMS_INFO *clatoms_info,
                   CPCOEFFS_INFO *cpcoeffs_info,
                   CPCOEFFS_POS *cpcoeffs_pos,
                   CPSCR *cpscr,CPOPTS *cpopts,PSEUDO *pseudo,
                   EWD_SCR *ewd_scr,ATOMMAPS *atommaps,
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

  double *vnlreal_up       = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up       = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn       = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn       = cpscr->cpscr_nonloc.vnlim_dn;

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

/*======================================================================*/
/* I) Get the ylm(g)                                                   */

  get_ylm(xk,yk,zk,g,ylmr,ylmi,dylmr_gx,dylmi_gx,dylmr_gy,dylmi_gy,
          dylmr_gz,dylmi_gz,ylm_cons);

  nl_chan_max = (nl_max + 1)*(nl_max + 1);
/*======================================================================*/
/* II) Calculate the nl-pseudoponential matrix elements by looping over */
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
	  index_atm[itype] =  (itype_nl-1)*nsplin_g*(n_ang_max+1)*n_rad_max
			     +l*nsplin_g*n_rad_max+(irad-1)*nsplin_g;
	}//endfor
	get_vpsnow(index_atm,nsplin_g,gmin_spl,dg_spl,g,
		   vps0,vps1,vps2,vps3,vtemp,iatm_typ_nl_rev,
		   natm_typ_nl,np_nonloc_cp_box_kb,
		   np_nl_rad_str[lp1][irad]);
	i_shift  = l*npart;
	for(ipart=np_nl_rad_str[lp1][irad];ipart<=np_nl[lp1];ipart++){
	  ktemp             =  ipart+i_shift;
	  ltemp             =  ip_nl_rev[ktemp];
	  helr_now[ipart]   =  helr[ltemp]*vtemp[ltemp];
	  heli_now[ipart]   =  heli[ltemp]*vtemp[ltemp];
	}//endfor
/*---------------------------------------------------------------------*/
/* ii) Loop over m components of the channel and get the matrix elements */
/*     vtemp,dvtemp, */
	sgn_l = 1.0;
	if((l%2)==1)sgn_l = -1.0;
	for(m=1;m<=(2*l+1);m++){
	  ind_lm = m + l*l;
	  getNlmatFilter(ncoef,ismcount,nstate_up,ind_lm,irad,
	    np_nlmax,np_nl[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
	    helr_now,heli_now,creal_up,cimag_up,
	    vnlreal_up,vnlimag_up,xk,yk,zk,
	    ylmr[ind_lm],ylmi[ind_lm],sgn_l,scr1,scr2);
	  if(cp_lsda==1){
	     getNlmatFilter(ncoef,ismcount,nstate_dn,ind_lm,irad,
	       np_nlmax,np_nl[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
	       helr_now,heli_now,creal_dn,cimag_dn,
	       vnlreal_dn,vnlimag_dn,xk,yk,zk,
	       ylmr[ind_lm],ylmi[ind_lm],sgn_l,scr1,scr2);
	  }//endif cp_lsda
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

void getNlmatFilter(int ncoef,int ismcount,int nstate,int ind_lm,int irad,
               int np_nlmax,int np_nl,int np_nl_rad_str,int nl_chan_max,
               double *helr,double *heli,
               double *creal,double *cimag,
               double *vnlreal,double *vnlimag,
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
    ioff_v    = (is-1)*np_nlmax +
                (ind_lm-1)*nstate*np_nlmax
               +(irad-1)*nl_chan_max*nstate*np_nlmax;
    switch(isgn_l){
       case 1:
         for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
           cbyhelr[ipart]     =  cre_ind*helr[ipart] - cim_ind*heli[ipart];
           cbyheli[ipart]     =  cre_ind*heli[ipart] + cim_ind*helr[ipart];
         }// endfor : loop over particles
         for(i=np_nl_rad_str;i<=np_nl;i++){
           tmp1        = -cbyheli[i]*tylmi;
           tmp2        =  cbyheli[i]*tylmr;
           vnlreal[(i+ioff_v)]    += tmp1;
           vnlimag[(i+ioff_v)]    += tmp2;
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
    }/*end switch: sgn of m in Ylm*/
  }/* endfor : loop over states */
/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void getYlmOnly(double xk,double yk,double zk,double g,
             double *ylmr,double *ylmi,YLM_CONS *ylm_cons)
/*========================================================================*/
  {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
  double y00r,y00i,y10r,y10i,y11r,y11i,y20r,y20i,y21r,y21i;
  double y22r,y22i,y30r,y30i,y31r,y31i,y32r,y32i,y33r,y33i;
  double ctheta,stheta,cphi,sphi,c2phi,s2phi,c3phi,s3phi;
  double xydist;
  double pre_phi_x,pre_phi_y,pre_theta_x,pre_theta_y,pre_theta_z;
  int iii;

 /* Local pointers */
  double rt_fpi     = ylm_cons->rt_fpi;
  double rt_thrfpi  = ylm_cons->rt_thrfpi;
  double rt_threpi  = ylm_cons->rt_threpi;
  double hrt_fivfpi = ylm_cons->hrt_fivfpi;
  double rt_fiftepi = ylm_cons->rt_fiftepi;
  double hrt_sevfpi = ylm_cons->hrt_sevfpi;
  double hrt_toepi  = ylm_cons->hrt_toepi;
  double hrt_ohffpi = ylm_cons->hrt_ohffpi;
  double hrt_tfepi  = ylm_cons->hrt_tfepi;
/*==========================================================================*/
/* I) Calculate polar angles of the vector g                                */

  ctheta = zk/g;
  stheta = sqrt(xk*xk + yk*yk)/g;
  xydist = sqrt(xk*xk + yk*yk);
  if(xydist==0){
    cphi = 1.0;
    sphi = 0.0;
  }else{
    cphi = xk/xydist;
    sphi = yk/xydist;
  }/*endif*/
  c2phi = cphi*cphi - sphi*sphi;
  s2phi = 2.0*cphi*sphi;
  c3phi = cphi*c2phi - sphi*s2phi;
  s3phi = cphi*s2phi + c2phi*sphi;

/*==========================================================================*/
/* I.i) l=0                                                                 */

  y00r    = rt_fpi;
  y00i    = 0.0;
  ylmr[1] = y00r;
  ylmi[1] = y00i;


/*==========================================================================*/
/* I.ii) l=1 (phi derivatives have stheta divided out)                      */

  y10r        =  rt_thrfpi*ctheta;
  y10i        = 0.0;
  y11r        =  rt_threpi*stheta*cphi;
  y11i        =  rt_threpi*stheta*sphi;
  ylmr[2] = y10r;
  ylmi[2] = y10i;
  ylmr[3] = y11r;
  ylmi[3] = y11i;
  ylmr[4] = y11r;
  ylmi[4] = -y11i;

/*==========================================================================*/
/* I.iii) l=2 (phi derivatives have stheta divided out)                     */

  y20r    = hrt_fivfpi*(3.0*ctheta*ctheta - 1.0);
  y20i    = 0.0;
  y21r    = rt_fiftepi*stheta*ctheta*cphi;
  y21i    = rt_fiftepi*stheta*ctheta*sphi;
  y22r    = 0.50*rt_fiftepi*stheta*stheta*c2phi;
  y22i    = 0.50*rt_fiftepi*stheta*stheta*s2phi;
  ylmr[5] = y20r;
  ylmi[5] = y20i;
  ylmr[6] = y21r;
  ylmi[6] = y21i;
  ylmr[7] = y21r;
  ylmi[7] = -y21i;
  ylmr[8] = y22r;
  ylmi[8] = y22i;
  ylmr[9] = y22r;
  ylmi[9] = -y22i;

/*==========================================================================*/
/* I.iv) l=3 (phi derivatives have stheta divided out)                      */

  y30r        = hrt_sevfpi*(5.0*ctheta*ctheta*ctheta - 3.0*ctheta);
  y30i        = 0.0;

  y31r        = hrt_toepi*stheta*(5.0*ctheta*ctheta - 1.0)*cphi;
  y31i        = hrt_toepi*stheta*(5.0*ctheta*ctheta - 1.0)*sphi;
  y32r        = hrt_ohffpi*stheta*stheta*ctheta*c2phi;
  y32i        = hrt_ohffpi*stheta*stheta*ctheta*s2phi;
  y33r        = hrt_tfepi*stheta*stheta*stheta*c3phi;
  y33i        = hrt_tfepi*stheta*stheta*stheta*s3phi;

  ylmr[10] = y30r;
  ylmi[10] = y30i;
  ylmr[11] = y31r;
  ylmi[11] = y31i;
  ylmr[12] = y31r;
  ylmi[12] = -y31i;
  ylmr[13] = y32r;
  ylmi[13] = y32i;
  ylmr[14] = y32r;
  ylmi[14] = -y32i;
  ylmr[15] = y33r;
  ylmi[15] = y33i;
  ylmr[16] = y33r;
  ylmi[16] = -y33i;

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void allocRealNl(CP *cp,CLASS *class)
/*========================================================================*/
  {/*begin routine*/
/*************************************************************************/
/* Prepare the real space nonlocal pp. Allocate necessary memory.        */
/*************************************************************************/
/*========================================================================*/
/*             Local variable declarations                                */
  PSEUDO *pseudo = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  PARA_FFT_PKG3D *cpParaFftPkg3d = &(cp->cp_para_fft_pkg3d_lg);
  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  COMMUNICATE *communicate = &(cp->communicate);
  
  int nfft = cpParaFftPkg3d->nfft;
  int numGrid = nfft/2;
  int numAtom = clatoms_info->natm_tot;
  int iPart;
  int numGridMax;
  int numThreads = communicate->numThreads;
  int *numGridNlppMap = pseudoReal->numGridNlppMap;
  
  pseudoReal->forceRealNlpp = (double*)cmalloc(numGrid*sizeof(double));

  numGridMax = numGridNlppMap[0];
  for(iPart=0;iPart<numAtom;iPart++){
    if(numGridNlppMap[iPart]>numGridMax)numGridMax = numGridNlppMap[iPart];
  }
  pseudoReal->numGridMax = numGridMax;
  //printf("numThreads %i iThread %i\n",numThreads,iThread);

  pseudoReal->forceTemp = (double*)cmalloc(numAtom*numGridMax*sizeof(double));
  pseudoReal->wfNbhd = (double*)cmalloc(numThreads*numGridMax*sizeof(double));

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/

