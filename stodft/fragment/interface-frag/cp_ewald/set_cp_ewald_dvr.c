/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_ewald (DVR)                              */
/*                                                                          */
/* This reads in and sets up the electron-atom interaction pseudopotential  */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_cp_ewald_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"

#define CP_EMAGIC 0.00050

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calc_kmaxv(int *kmaxv, int *kmax_cp, int *dvr_gd_nx, int *dvr_gd_ny,
                int *dvr_gd_nz, double *ecut_cp, double deth, double kmax_ewd,
                int myid)

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/

   int kmin;
   int n1,n2,n3;

   double rtwoth,rvol23,ecut;

/*==========================================================================*/
/* I) No truncation for DVR basis   */

   n1 = (*dvr_gd_nx-1)/2;
   n2 = (*dvr_gd_ny-1)/2;
   n3 = (*dvr_gd_nz-1)/2;

   kmin = n1;
   kmin = (kmin < n2 ? kmin:n2);
   kmin = (kmin < n3 ? kmin:n3);

   rtwoth = 2.0/3.0;
   rvol23 = 1.0/pow(deth,rtwoth);
   *ecut_cp = M_PI * .5 * M_PI * (double) (kmin * kmin)*rvol23;

   kmaxv[1] = n1;
   kmaxv[2] = n2;
   kmaxv[3] = n3;

   kmax_cp[1] = kmaxv[1];
   kmax_cp[2] = kmaxv[2];
   kmax_cp[3] = kmaxv[3];

/*==========================================================================*/
/* II) Compare Ewald to CP cutoff  */

   ecut = M_PI * .5 * M_PI * (double) (kmax_ewd * kmax_ewd)*rvol23;

   if (*ecut_cp * .5 < ecut) {
     if(myid==0){
       printf("@@@@@@@@@@@@@@@@@@@@_WARNING_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Ewald cutoff greater than the approximate cp cutoff\n");
       printf(": %g vs %g .\n", ecut*2.0,*ecut_cp);
       printf("This is probably OK with DVR basis, but just make\n"); 
       printf("it sure that ewalds parameters are OK\n");
       printf("@@@@@@@@@@@@@@@@@@@@_WARNING_@@@@@@@@@@@@@@@@@@@@\n");
     }
   }/*endif*/

/*-------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void setkvec3d_all(int nktot,int *kmaxv,double *hmatik,
                   int *kastore, int *kbstore, int *kcstore,
                   int *ibreak1, int *ibreak2, 
                   double *gmin_spl, double *gmax_spl)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

    int iii;
    int i1, i2, i3;

    int i, kbmin, kcmin, kbmax, kcmax, kamax, ka, kb, kc;
    double xk, yk, zk;
    int icount;
    double aka, akb, akc;
    double tpi;
    double try;

/*==========================================================================*/
/* Count the K-vectors */

    tpi = 2.0*M_PI;
    for(i=1;i<=nktot;i++){
      ibreak1[i] = 0;
      ibreak2[i] = 0;
    }/*endfor*/
    icount = 0;
    (*gmin_spl) = 1.0e10;
    (*gmax_spl) = 0.0;

    kamax = kmaxv[3];
    for (ka = -kamax; ka <= kamax; ++ka) {
      aka = (double) ka;
      kbmin = -kmaxv[2];
      kbmax = kmaxv[2];

      for (kb = kbmin; kb <= kbmax; ++kb) {
        akb = (double) kb;
        kcmin = -kmaxv[1];
        kcmax = kmaxv[1];

        for (kc = kcmin; kc <= kcmax; ++kc) {

          if(!((ka==0)&&(kb==0)&&(kc==0))){
            ++icount;
            aka = (double) ka;
            akb = (double) kb;
            akc = (double) kc;
            xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
            yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
            zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
            kastore[icount] = kc;
            kbstore[icount] = kb;
            kcstore[icount] = ka;
            if (kc == kcmin) {ibreak1[icount] = 1; }
            if (kc <  kcmax) {ibreak2[icount] = 1; }
            try = sqrt(xk * xk + yk * yk + zk * zk);
            (*gmin_spl) = MIN(try,(*gmin_spl));
            (*gmax_spl) = MAX(try,(*gmax_spl));
          }
    } } }

    kastore[nktot+1] = 0;
    kbstore[nktot+1] = 0;
    kcstore[nktot+1] = 0;

    if(nktot!= icount){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("      Mismatch number of kvectors\n");
        printf("        %d vs %d\n",icount,nktot);
        printf("      Contact technical support\n");
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
    }

/*--------------------------------------------------------------------------*/
} /* setkvec3d_all */
/*==========================================================================*/
/*=======================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*=======================================================================*/

void set_cpmass_dvr(CPCOEFFS_INFO *cpcoeffs_info, double *cmass_tau_def,
                    int *icmass_unif)

/*=======================================================================*/
/*            Begin subprogram:                                          */
    {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations:                               */

  int ncoef       = cpcoeffs_info->ncoef;
  int grid_nx     = cpcoeffs_info->grid_nx;
  int grid_ny     = cpcoeffs_info->grid_ny;
  int grid_nz     = cpcoeffs_info->grid_nz;

  double *cmass   = cpcoeffs_info->cmass;
  double scale_fact = cpcoeffs_info->scale_fact;

  double cmass0,tpi;        
  int i,j,k;            
  int icount;
  int index_x,index_y,index_z;

/*=======================================================================*/
/*  A) Set up CP masses                                                  */

  *icmass_unif = 1;
  tpi = 2.0*M_PI;

  *cmass_tau_def /= TIME_CONV;

  cmass0 = *cmass_tau_def;     
  cmass0 = 4.0*(*cmass_tau_def)*(*cmass_tau_def)*CP_EMAGIC;

  icount = 0;
  for(k=1;k<= grid_nz;k++) {
    index_z = (k-1)*grid_nz + k;
    for(j=1;j<= grid_ny;j++) {
      index_y = (j-1)*grid_ny + j;
      for(i=1;i<= grid_nx;i++) {
        index_x = (i-1)*grid_nx + i;
        icount++;
        cmass[icount] = cmass0*scale_fact;
  }}}

/*-----------------------------------------------------------------------*/
} /* end routine */
/*=======================================================================*/


