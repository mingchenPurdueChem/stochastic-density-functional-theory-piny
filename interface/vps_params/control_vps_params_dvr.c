/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_vps_params                           */
/*                                                                          */
/* This reads in and sets up the electron-atom interaction pseudopotential  */ 
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_vps_params_entry.h"
#include "../proto_defs/proto_search_entry.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_vps_params_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


void make_vps_splin_dvr(char *vps_file,int loc_opt,int natm_typ,
                        int iatm_typ, int *pnsplin_mall_dvr, PSEUDO *pseudo,
                        int *pnsplin_r_tot, double *gzvps,double *gzvps0,
                        int iopt_nl_filter, double phi0_0, double phi0_1, 
                        double phi0_2,double alpha, double beta, int myid,
                        MPI_Comm comm, int num_proc, double grid_h, int iperd)

/*==========================================================================*/
/*               Begin subprogram:                                          */
{/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int n_ang         =  pseudo->n_ang[iatm_typ];
  int ivps_label    =  pseudo->ivps_label[iatm_typ];
  int n_ang_max     =  pseudo->n_ang_max;
  int *nsplin_r     =  pseudo->nsplin_r;
  double *dr_spl    =  pseudo->dr_spl;
  double *rmin      =  pseudo->rmin;

  /*
  double *rmax      =  pseudo->rmax;
  double *z1        =  pseudo->z_1;
  double *z2        =  pseudo->z_2;
  double *alp1      =  pseudo->alp_1;
  double *alp2      =  pseudo->alp_2;
  */

  double rcut_nl_now = pseudo->rcut_nl[iatm_typ];
  double alpha_conv_dual = pseudo->alpha_conv_dual;

  double *vps0,*vps1,*vps2,*vps3;
  double *dvps0,*dvps1,*dvps2,*dvps3;

  int nsplin_g      =  pseudo->nsplin_g;
  double dg_spl     =  pseudo->dg_spl;
  double gmax_spl   =  pseudo->gmax_spl;
  double gmin_spl   =  pseudo->gmin_spl;
  double gmin_true  =  pseudo->gmin_spl; /* DY */
  double *vpsnorm   =  pseudo->vpsnorm;

  double r_max;                  /* Num: Maximum spline radius        */
  double z_1,z_2;                /* Num: Charges for long range piece */
  double alpha_1,alpha_2;        /* Num: Ewald alpha's for long range */
  double dr;                     /* Num: delta r for spline           */
  double amat;                   /* Num: Useful temporary for spline  */
  double v_now;                  /* Num: current value of pp          */
  double rphi_now;               /* Num: current value of pseudo wf   */
  double zpol;                   /* Num: Polarization charge          */
  double gamma;                  /* Num: Useful temporary for long range */
  double dummy;                  /* Num: Really dum for pressure tensor*/
  int i,iii;                     /* Num: Generic counter              */

  int n_ang_now;                 /* Num: Number of angular momentum
                                               components            */
  int iang;                      /* Num: Angular momentum counter    */
  int ishift_now;                /* Num: Angular momentum shift      */
  int iang_now;                  /* Num: Angular momentum counter    */
  int ioff_nr;                   /* Num: Offset for real space spline */
  int ioff_ng;                   /* Num: Offset for reciprocal space spline */

  double *phi0;                  /* Num: r=0 value for pseudo wave func */
  double *v_phi,*v_loc,*r;       /* Lst: temporary spline arrays     */
  double *v_rphi, *v_g, *rphi;
  int ir,ig;                        /* Num: counter for above arrays    */
  int nr;                        /* Num: Length of above arrays      */
  int nsplin_r_loc;
  int ishift_norm;
  FILE *fp_vps_file;             /* File pointer to pseudo pot file  */

  double tmp;
  int nsplin_r_tot = *pnsplin_r_tot;
  double *g;

  double gmin_dvr = 0.001;
  double dg_dvr;
  int n_rad_max = 1;
  int n_rad_max_sq = n_rad_max*n_rad_max;

 /* FOR ERFC FUNCTION */
  double pi=M_PI;
  double erfc1,erfc2;
  double e1 = 0.2041422096422003, e2 = 0.1997535956961481;
  double e3 = 0.2213176596405576, e4 = 0.03360430734640255;
  double e5 = 0.4732592578721755, e6 =-0.509078520069735;
  double e7 = 0.6772631491947646, e8 =-0.369912979092217;
  double e9 = 0.06965131976970335;
  double p=0.3614;

  double tt,eee,ralp;

 /* FOR SPLINING IN G SPACE */
  int ilong;
  int cp_ptens_calc = 0; /* no pressure tensor for DVR yet */
  int cp_dual_grid_opt = 0;

  /*
  double alpha = pseudo->nl_alp[iatm_typ];
  double beta  = pseudo->nl_beta[iatm_typ];
  */
  double g_cut,big_R;
  double *sf;

/*==========================================================================*/
/*  0) Open the electron-atom pseudopotential file for reading              */

   if(myid==0){
     fp_vps_file = cfopen(vps_file,"r");
   }/*endif*/

/*==========================================================================*/
/*  I) Allocate g-space vector for spline                                   */

   g      = (double *) cmalloc(nsplin_g*sizeof(double))-1;
   v_g    = (double *) cmalloc(nsplin_g*sizeof(double))-1;
   sf     = (double *) cmalloc(nsplin_g*sizeof(double))-1;

   if(iperd==3){
     dg_dvr = (gmax_spl - gmin_spl)/((double)nsplin_g);

     for(i=1;i <= nsplin_g;i++){
       g[i] = dg_dvr*((double) (i-1)) + gmin_spl;
     }
   }else if(iperd==0){

     dg_dvr = (gmax_spl - 0.0)/((double)nsplin_g);

     for(i=1;i <= nsplin_g;i++){
       g[i] = dg_dvr*((double) (i));
     }/*endfor*/

     g_cut=pi/grid_h*alpha;

     if(myid==0){
       printf("Local pseudopotential Fourier filtering for atom_type %d\n",
               iatm_typ);
       printf("pi/h = %.10g,  g_cut= %.10g\n\n",pi/grid_h, g_cut);
     }

     for(i=1;i <= nsplin_g;i++){
       if(g[i] > g_cut){
         sf[i]=exp(-beta*(g[i]/g_cut-1.0)*(g[i]/g_cut-1.0));
       }else{
         sf[i]=1.0;
       }
     }/*endfor*/

   }else{
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("   WRONG PERIODICITY in make_vps_splin_dvr  \n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     exit(1);
   }

/*==========================================================================*/
/*  II) Set up a local pseudo potential                                     */

   if(ivps_label <= 4) {

     phi0 = (double *) cmalloc(3*sizeof(double))-1;
     phi0[1] = phi0_0;
     phi0[2] = phi0_1;
     phi0[3] = phi0_2;

     if(myid==0){
       if(fscanf(fp_vps_file,"%d %lf %d\n",&nr,&r_max,&n_ang_now) != 3)
                      {vps_read_error(vps_file);}
       if(fscanf(fp_vps_file,"%lf %lf %lf %lf\n",
                              &z_1,&alpha_1,&z_2,&alpha_2) != 4)
                      {vps_read_error(vps_file);}

       if(fscanf(fp_vps_file,"%lf %lf\n",&zpol,&gamma) != 2)
                 {vps_read_error(vps_file);}
     }/*endif*/

     if( (n_ang_now < n_ang) && (myid==0) && (ivps_label!=5) ) {
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dude, given the number of angular momentum\n");
       printf("components you've specified in the\n");
       printf("pseudopotential file %s\n",vps_file);
       printf("proceeding with the simulation would be a\n");
       printf("pointless exercise leading to most bogus results\n");
       printf("%d vs %d\n",n_ang_now,n_ang);
       putchar('\n');
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/*endif*/
     /*==============MA===============*/
     if((myid==0)&&(fabs(zpol) > 1.0e-8)){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("The value of zpol is non-zero in file %s \n",vps_file);
       printf("polarization corrections are not implemented for the \n");
       printf("DVR basis set. \n");
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/*endif*/

     if(num_proc>1){
       Bcast(&(nr),1,MPI_INT,0,comm);
       Bcast(&(r_max),1,MPI_DOUBLE,0,comm);
       Bcast(&(n_ang_now),1,MPI_INT,0,comm);
       Bcast(&(z_1),1,MPI_DOUBLE,0,comm);
       Bcast(&(alpha_1),1,MPI_DOUBLE,0,comm);
       Bcast(&(z_2),1,MPI_DOUBLE,0,comm);
       Bcast(&(alpha_2),1,MPI_DOUBLE,0,comm);
       Bcast(&(zpol),1,MPI_DOUBLE,0,comm);
       Bcast(&(gamma),1,MPI_DOUBLE,0,comm);
     }/*endif*/

     /* Assign DVR Pseudo array values */
     nsplin_r_loc = nr;
     nsplin_r[iatm_typ] = (nr > nsplin_g ? nr : nsplin_g);

     dr_spl[iatm_typ] = r_max/nr;
     rmin[iatm_typ]   = 0.00;

     /*
     rmax[iatm_typ]   = r_max;
     z1[iatm_typ]     = z_1;
     alp1[iatm_typ]   = alpha_1;
     z2[iatm_typ]     = z_2;
     alp2[iatm_typ]   = alpha_2;
     */

     nsplin_r_tot += (n_ang_max + 1)*(nsplin_r[iatm_typ])*(natm_typ); 
      
     if(nsplin_r_tot > (*pnsplin_mall_dvr) ){
       if(myid == 0){ 
         printf("Allocating more memory for dvr pseudopotentials\n");
       }
       pseudo->vps0 = (double *) crealloc(&(pseudo->vps0[1]),nsplin_r_tot*sizeof(double))-1;

       pseudo->vps1 = (double *) crealloc(&(pseudo->vps1[1]),nsplin_r_tot*sizeof(double))-1;
       pseudo->vps2 = (double *) crealloc(&(pseudo->vps2[1]),nsplin_r_tot*sizeof(double))-1;
       pseudo->vps3 = (double *) crealloc(&(pseudo->vps3[1]),nsplin_r_tot*sizeof(double))-1;

       pseudo->dvps0 = (double *) crealloc(&(pseudo->dvps0[1]),nsplin_r_tot*sizeof(double))-1;
       pseudo->dvps1 = (double *) crealloc(&(pseudo->dvps1[1]),nsplin_r_tot*sizeof(double))-1;
       pseudo->dvps2 = (double *) crealloc(&(pseudo->dvps2[1]),nsplin_r_tot*sizeof(double))-1;
       pseudo->dvps3 = (double *) crealloc(&(pseudo->dvps3[1]),nsplin_r_tot*sizeof(double))-1;

       *pnsplin_mall_dvr = nsplin_r_tot;
     }
   }/*endif ivps_label */


   vps0      =  pseudo->vps0;
   vps1      =  pseudo->vps1;
   vps2      =  pseudo->vps2;
   vps3      =  pseudo->vps3;

   dvps0     =  pseudo->dvps0;
   dvps1     =  pseudo->dvps1;
   dvps2     =  pseudo->dvps2;
   dvps3     =  pseudo->dvps3;

   /* allocate more memory for NL PP */

   v_phi =  (double *) cmalloc(nr*sizeof(double))-1;
   v_rphi = (double *) cmalloc(nr*sizeof(double))-1;
   v_loc  = (double *) cmalloc(nr*sizeof(double))-1;
   r      = (double *) cmalloc(nr*sizeof(double))-1;
   rphi   = (double *) cmalloc(nr*sizeof(double))-1;

   /*Allocate real space arrays   */

   dr = r_max/((double) nr);

   for(i=1;i <= nr;i++){
     r[i] = ((double ) (i-1))*dr;
   }/*endfor*/

/*==========================================================================*/
/* III) Find Local  pseudo potential */

   if(ivps_label < 5){ /* NOT GOEDECKER TYPE */
     if(myid==0){
       for(iang=1;iang <= (loc_opt + 1);iang++) {
         for(ir=1;ir <= nr;ir++) {
           if(fscanf(fp_vps_file,"%lf %lf\n",&v_now,&rphi_now) != 2)
                           {vps_read_error(vps_file);}
           v_loc[ir] = v_now;
         }/* endfor */
       } /* endfor */
     }/*endif*/
   }/*endif*/

   if(num_proc>1){
     Bcast(&(v_loc[1]),nr,MPI_DOUBLE,0,comm);
   }

   ioff_nr = 0;
   for(i=1; i< iatm_typ; i++){ 
     ioff_nr += (n_ang_max+1)*nsplin_r[i];
   }
   ioff_nr += loc_opt*nsplin_r[iatm_typ];
   ishift_now = ioff_nr;

   for(ir=1;ir <= nr;ir++) { 
     vps0[ir+ishift_now] = v_loc[ir]; 
   }


   if(iperd==0){

     if(myid==0) printf("SPLINING IN REAL SPACE FOR Local PP\n");

     ilong = 0;
     iang_now = 0;

     for(ir=1;ir <= nr;ir++) { v_rphi[ir] = vps0[ir+ishift_now]*r[ir]; }

     slow_bess_vps(v_rphi,nr,dr,r,v_g,&dummy,
                   nsplin_g,g,gmin_true,
                   z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
                   gzvps,&gzvps0[1],cp_ptens_calc,cp_dual_grid_opt,
                   alpha_conv_dual);

     for(ig=1;ig<=nsplin_g;ig++){
       v_g[ig]=v_g[ig]*sf[ig];
     }

     slow_invs_bess_vps(v_g,nsplin_g,dg_dvr,g,&(vps0[ishift_now]),
                         nr,r,iang_now);

      /* add long range contribution */

     for(ir=2;ir <= nr;ir++) {
       vps0[ir+ishift_now] -= (z_1*gerf(alpha_1*r[ir])+z_2*gerf(alpha_2*r[ir]))/r[ir];
     }
     vps0[1+ishift_now] -= (z_1*gerf(alpha_1*r[2])+z_2*gerf(alpha_2*r[2]))/r[2];


     spline_fit(&(vps0)[ishift_now],&(vps1)[ishift_now],
               &(vps2)[ishift_now],&(vps3)[ishift_now],r,nsplin_r_loc);

   }else if (iperd==3){

     if(myid==0) printf("SPLINING IN RECIPROCAL SPACE For Local PP\n");

     for(ir=1;ir <= nr;ir++) {v_rphi[ir] = v_loc[ir]*r[ir];}

     ilong      = 1;  /* 0 = no long range, 1 = long range  */
     iang_now   = 0;

     slow_bess_vps(v_rphi,nr,dr,r,&(vps0)[ishift_now],&dummy,
                   nsplin_g,g,gmin_true,
                   z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
                   gzvps,gzvps0,cp_ptens_calc,cp_dual_grid_opt,
                   alpha_conv_dual);

     spline_fit(&(vps0)[ishift_now],&(vps1)[ishift_now],
                &(vps2)[ishift_now],&(vps3)[ishift_now],g,nsplin_g);


   }else{
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("   WRONG PERIODICITY in make_vps_splin_dvr   \n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     exit(1);
   }

/*==========================================================================*/
/* IV) Set up a null local pseudo potential                                 */

   if(ivps_label == 4) {

     ioff_nr = 0;
     for(i=1; i< iatm_typ; i++){ 
       ioff_nr += nsplin_r[i];
     }
     ioff_nr *= (n_ang_max + 1);
     ishift_now = ioff_nr;

     for(i=1;i <= nsplin_r[iatm_typ]; i++) {
       vps0[i+ishift_now] = 0.0;
       vps1[i+ishift_now] = 0.0;
       vps2[i+ishift_now] = 0.0;
       vps3[i+ishift_now] = 0.0;
     }/*endfor*/

   }/* endif:null potential */


/*======================================================================*/
/* V) set-up Fourier filtering parameter */

   if(iopt_nl_filter == 1){

     if(iperd==3){
       dg_dvr = (gmax_spl - 0.0)/((double)nsplin_g);

       for(i=1;i <= nsplin_g;i++){
         g[i] = dg_dvr*((double) (i));
       }/*endfor*/
     }

     g_cut = pi/grid_h*alpha; 
   
     for(i=1;i <= nsplin_g;i++){
       if(g[i] > g_cut){
         sf[i]=exp(-beta*(g[i]/g_cut-1.0)*(g[i]/g_cut-1.0));
       }else{
         sf[i]=1.0;
       }
     }
     if(myid==0){
       printf("Non-local pseudopotential Fourier filtering for atom_type %d\n",
               iatm_typ);
       printf("pi/h = %.10g,  g_cut= %.10g\n\n",pi/grid_h, g_cut);
     }
   }else{
     if(n_ang !=0 && myid==0){
       printf("$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
       printf("You have chosen DVR basis, but non-local filtering\n");
       printf("option is not turned on for atm_typ = %d.\n",iatm_typ);
       printf("Make it sure that you provied correct pseudopotential\n");
       printf("value at r=0 in pi_md.vps\n");
       printf("$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
     }
   }

/*==========================================================================*/
/* VI) nonlocal potential */

   /*Only KB nonlocal potential for now */

   if(ivps_label == 2 || ivps_label == 3 || ivps_label == 5) {
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf(" Only KB type nonlocal pseudopotential has been\n");
       printf(" implemented for DVR basis\n");
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }else{

     /* A) rewind */

     if(myid==0){

       rewind(fp_vps_file);

       if(fscanf(fp_vps_file,"%d %lf %d\n",&nr,&r_max,&n_ang_now) != 3){
         vps_read_error(vps_file);
       }
       if(fscanf(fp_vps_file,"%lf %lf %lf %lf\n",
                                     &z_1,&alpha_1,&z_2,&alpha_2)!= 4) {
         vps_read_error(vps_file);
       }

       if(fscanf(fp_vps_file,"%lf %lf\n",&zpol,&gamma) != 2){
         vps_read_error(vps_file);
       }
     }/*endif myid*/

     /* B) Spline projection operator*/

     ishift_norm = (iatm_typ - 1)*(pseudo->n_ang_max+1);
     vpsnorm[(loc_opt + 1)+ishift_norm] = 0.0;

     for(iang = 1; iang <= n_ang + 1; iang++) {
       amat = 0.0;

       if(iopt_nl_filter==0){ /*no-filtering*/
         if(myid==0){
           /* r=0 point (first point)*/
           if(fscanf(fp_vps_file,"%lf %lf\n",&v_now,&rphi_now) != 2){
             vps_read_error(vps_file); 
           }
           v_phi[1] = (v_now-v_loc[1])*phi0[iang];
           /* the rest of r points */
           for(ir=2; ir <= nr; ir++) {
             if(fscanf(fp_vps_file,"%lf %lf\n",&v_now,&rphi_now) != 2){
               vps_read_error(vps_file); 
             }
             v_phi[ir] = (v_now-v_loc[ir])*rphi_now/r[ir];
             amat += dr*v_phi[ir]*rphi_now*r[ir];
           } 
         }/*endif myid*/

         if(num_proc>1){
          Bcast(&(v_phi[1]),nr,MPI_DOUBLE,0,comm);
          Bcast(&(amat),1,MPI_DOUBLE,0,comm);
         }

         if(iang != loc_opt+1) {
           vpsnorm[((iang-1)+1)+ishift_norm] = (1.0/amat);
           ioff_nr = 0;
           for(i=1; i< iatm_typ; i++){ 
             ioff_nr += nsplin_r[i];
           }
           ioff_nr *= (n_ang_max + 1);
           ishift_now = ioff_nr + (iang-1)*nsplin_r[iatm_typ];
           for(ir=1; ir <= nr; ir++){ 
             vps0[ishift_now + ir] = v_phi[ir];
           }
           spline_fit(&(vps0[ishift_now]),&(vps1[ishift_now]),
                      &(vps2[ishift_now]),&(vps3[ishift_now]),
                      r,nsplin_r_loc);
         }/*endfor loc_opt*/
       }/*end if no-filtering*/

       if(iopt_nl_filter==1){
         if(myid==0){
           for(ir=1; ir <= nr; ir++) {
             if(fscanf(fp_vps_file,"%lf %lf\n",&v_now,&rphi_now) != 2){
               vps_read_error(vps_file); 
             }
             v_rphi[ir] = (v_now-v_loc[ir])*rphi_now;
             rphi[ir] = rphi_now;
             amat += rphi_now*v_rphi[ir]*dr;
           } /* endfor */
         }/*endif*/

         if(num_proc>1){
           Bcast(&(v_rphi[1]),nr,MPI_DOUBLE,0,comm);
           Bcast(&(rphi[1]),nr,MPI_DOUBLE,0,comm);
         }/*endif*/

         if(iang != loc_opt+1) {
           ioff_nr = 0;
           for(i=1; i< iatm_typ; i++){ 
             ioff_nr += nsplin_r[i];
           }
           ioff_nr *= (n_ang_max + 1);
           ishift_now = ioff_nr + (iang-1)*nsplin_r[iatm_typ];

           ilong = 0;
           iang_now = iang-1;

           slow_bess_vps(v_rphi,nr,dr,r,v_g,&dummy, nsplin_g,g,gmin_true,
                         z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
                       gzvps,&gzvps0[1],cp_ptens_calc,cp_dual_grid_opt,
                       alpha_conv_dual);

           for(ig=1;ig<=nsplin_g;ig++){
             v_g[ig]=v_g[ig]*sf[ig];
           }

           slow_invs_bess_vps(v_g,nsplin_g,dg_dvr,g,&(vps0[ishift_now]),
                              nr,r,iang_now);

           amat=0.0;
           for(ir=1;ir<=nr;ir++){
             if(r[ir] > rcut_nl_now-0.5 && r[ir] <=rcut_nl_now) {
             /* cut-off oscillation */
               big_R=(r[ir]-(rcut_nl_now-0.5))/0.5;
               vps0[ishift_now+ir]*=(1.0+big_R*big_R*(2.0*big_R-3.0));

             }else if(r[ir] > rcut_nl_now){
               vps0[ishift_now+ir]=0.0;
             }
             amat += dr*vps0[ishift_now+ir]*rphi[ir]*r[ir];
           }

           vpsnorm[((iang-1)+1)+ishift_norm] = (1.0/amat);

           spline_fit(&(vps0[ishift_now]),&(vps1[ishift_now]),
                      &(vps2[ishift_now]),&(vps3[ishift_now]),
                      r,nsplin_r_loc);

         }/* endif:loc_opt */
       }/*endif filtering */
     }/* endfor:channels */
   } /*endif KB*/

/*==========================================================================*/
/*   VII) Free memory                                                       */

    cfree(&(phi0[1]));
    cfree(&(v_phi[1]));
    cfree(&(v_rphi[1]));
    cfree(&(v_loc[1]));
    cfree(&(r[1]));
    cfree(&g[1]);
    cfree(&(v_g[1]));
    cfree(&(rphi[1]));
    cfree(&(sf[1]));

/*==========================================================================*/
/*   VIII) Close the file and done                                          */

   if(myid==0){
     fclose(fp_vps_file);
   }/*endif*/

  *pnsplin_r_tot = nsplin_r_tot;

/*-----------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


 void slow_invs_bess_vps(double v_g[],int nsplin_g,double dg,double g[],
                         double fv_rphi[],int nr,double r[],int iang)


/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */

   int ig,ir,i;
   double phi;
   double g2,rj0,rj1,rj2,rj3,rzero,arg,fpi,pi,tpi,tpisqdg;


/*==========================================================================*/
/* II) Get some constants                                                    */


   pi  = M_PI;
   tpi = 2.0*pi;
   fpi = 4.0*pi;
   tpisqdg = dg/(2.0*pi*pi);
   for(ir=1;ir <= nr;ir++) {
     fv_rphi[ir] = 0.0;
   }

/*==========================================================================*/


   switch(iang) {
/*==========================================================================*/
/* I) L=0 Term */

    case 0:

     rzero=0.0;

     for(ig=1;ig<=nsplin_g;ig++){rzero+=tpisqdg*g[ig]*g[ig]*v_g[ig];}
     fv_rphi[1]=rzero;

     for(ir=2;ir <= nr; ir++) {
       for(ig=1;ig <= nsplin_g; ig++) {
        arg = r[ir]*g[ig];
        rj0 = (sin(arg)/arg);
        fv_rphi[ir]  += tpisqdg*rj0*v_g[ig]*g[ig]*g[ig];
       } /* endfor */
      } /* endfor */

    break;

/*==========================================================================*/
/* II) L=1 Term */

    case 1:

/*--------------------------------------------------------------------------*/
/*    i) g ne 0                                                             */


      fv_rphi[1]=0.0;

      for(ir=2;ir <= nr; ir++) {
        for(ig=1;ig <= nsplin_g; ig++) {

          arg = r[ir]*g[ig];
          rj1 = ((sin(arg)/arg - cos(arg))/arg);
          fv_rphi[ir] += tpisqdg*rj1*v_g[ig]*g[ig]*g[ig];
        } /* endfor */
      } /* endfor */

    break;

/*==========================================================================*/
/* III) L=2 Term */
    case 2:
/*--------------------------------------------------------------------------*/
/*    i) g ne 0                                                             */

      fv_rphi[1]=0.0;

      for(ir=2;ir <= nr; ir++) {
        for(ig=1;ig <= nsplin_g; ig++) {
          arg = r[ir]*g[ig];
          rj2 = (((3.0/(arg*arg)-1.0)*sin(arg)-3.0*cos(arg)/arg)/arg);
          fv_rphi[ir] += tpisqdg*rj2*v_g[ig]*g[ig]*g[ig];
        } /* endfor */
      } /* endfor */

      break;

/*==========================================================================*/
/* IV) L=3 Term */
    case 3:
/*--------------------------------------------------------------------------*/
/*    i) g ne 0     */

      fv_rphi[1]=0.0;

      for(ir=2;ir <= nr; ir++) {
        for(ig=1;ig <= nsplin_g; ig++) {

          arg = r[ir]*g[ig];
          rj3 = (((15.0/(arg*arg) - 6.0)*sin(arg)/arg +
                 (1.0 - 15.0/(arg*arg))*cos(arg))/arg);
          fv_rphi[ir] += tpisqdg*rj3*v_g[ig]*g[ig]*g[ig];
        } /* endfor */
      } /* endfor */

    break;

   } /* end switch */

/*-----------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/







