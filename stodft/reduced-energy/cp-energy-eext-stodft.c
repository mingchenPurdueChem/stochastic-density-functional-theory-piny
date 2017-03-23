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
    kastore   = ewald->kastr;
    kbstore   = ewald->kbstr;
    kcstore   = ewald->kcstr;
    ibreak1   = ewald->ibrk1;
    ibreak2   = ewald->ibrk2;
    vextr     = cpscr->cpscr_loc.vextr;
    vexti     = cpscr->cpscr_loc.vexti;
    vextr_loc = cpscr->cpscr_loc.vextr_loc;
    vexti_loc = cpscr->cpscr_loc.vexti_loc;
    dvextr    = cpscr->cpscr_loc.dvextr;
    dvexti    = cpscr->cpscr_loc.dvexti;
    rhocr     = cpscr->cpscr_rho.rhocr_up;
    rhoci     = cpscr->cpscr_rho.rhoci_up;
    ak2       = cpewald->ak2;
    nktot     = ewald->nktot;
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
    
  memcpy(&vextr_loc[1],&vextr[1],(ngo+1)*sizeof(double));
  memcpy(&vexti_loc[1],&vextr[1],(ngo+1)*sizeof(double));

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
		   ip_nl,vnorm_now,vnlreal_up,vnlimag_up,,&cp_enl);
	  if(cp_lsda==1){
	    sumNlPot(npart,nstate_dn,np_nlmax,nl_chan_max,
		     np_nl[lp1],l,np_nl_rad_str[lp1][jrad],
		     irad,jrad,
		     ip_nl,vnorm_now,vnlreal_dn,vnlimag_dn,,&cp_enl);
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
/*==========================================================================*/
/* II) Set the return values                                               */

  *cp_enl_ret = cp_enl;

/*======================================================================*/
  }/*end routine*/
/*======================================================================*/


