 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/
 /*                                                                          */
 /*                         PI_MD:                                           */
 /*             The future of simulation technology                          */
 /*             ------------------------------------                         */
 /*                   Routine: set_dvr_grid.c                                */
 /*                                                                          */
 /* The routine to initialize DVR related arrays                             */
 /*                                                                          */
 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_parse_local.h"
#include "../proto_defs/proto_math.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calculate_cp_nfree_dvr(CP *cp)

/*========================================================================*/
/*             Begin Routine                                              */
{/*Begin subprogram: */
/*========================================================================*/
/*             Local variable declarations                                */

  int ncons,ncons_dn;

/*------------------------------------------------------------------------*/
/* Number of up degrees of freedom */

  if(cp->cpopts.cp_norb <= 1) {
    ncons = (cp->cpcoeffs_info.nstate_up)*
             ((cp->cpcoeffs_info.nstate_up)+1)/2;
  }/* endif norb off or norb full_ortho */

  /*norm only*/
  if(cp->cpopts.cp_norb == 2) {ncons = cp->cpcoeffs_info.nstate_up; } 
  if(cp->cpopts.cp_norb == 3) {ncons = 0; }    /* No constraint */
  cp->cpcoeffs_info.cp_nfree = (double) (cp->cpcoeffs_info.ncoef*
                                         cp->cpcoeffs_info.nstate_up-ncons);
  cp->cpcoeffs_info.cp_nfree_up = cp->cpcoeffs_info.cp_nfree;

/*------------------------------------------------------------------------*/
/* Number of down degrees of freedom */

  if(cp->cpopts.cp_lsda == 1) {
    if(cp->cpopts.cp_norb <= 1) {
      ncons = (cp->cpcoeffs_info.nstate_up)*
               ((cp->cpcoeffs_info.nstate_up)+1)/2+
               (cp->cpcoeffs_info.nstate_dn)*
               ((cp->cpcoeffs_info.nstate_dn)+1)/2;
      ncons_dn = (cp->cpcoeffs_info.nstate_dn)*
                  ((cp->cpcoeffs_info.nstate_dn)+1)/2;
    }/* endif norb off or norb full_ortho */
    if(cp->cpopts.cp_norb == 2) {
      ncons = cp->cpcoeffs_info.nstate_up + cp->cpcoeffs_info.nstate_dn;
      ncons_dn = cp->cpcoeffs_info.nstate_dn;
    }/* endif norm only */
    if(cp->cpopts.cp_norb == 3) {
      ncons = 0;
      ncons_dn = 0;
    }/* endif no constraint */
    cp->cpcoeffs_info.cp_nfree = (double)(cp->cpcoeffs_info.ncoef*
                                         (cp->cpcoeffs_info.nstate_up+
                                          cp->cpcoeffs_info.nstate_dn)-ncons);
    cp->cpcoeffs_info.cp_nfree_dn = (double)(cp->cpcoeffs_info.ncoef*
                                            cp->cpcoeffs_info.nstate_dn-ncons_dn);
  } /* endif */

/*------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

void move_to_cell_center(CLASS *class)
{
  int ip,ip_now,i,j,k;
  double dx,dy,dz;
  int ip_start        = class->clatoms_info.pi_beads_proc_st;
  int ip_end          = class->clatoms_info.pi_beads_proc_end;
  int natm_tot        = class->clatoms_info.natm_tot;
  double *x,*y,*z;

  for(ip=ip_start; ip <=ip_end;ip++){
    x = class->clatoms_pos[ip].x;
    y = class->clatoms_pos[ip].y;
    z = class->clatoms_pos[ip].z;
    dx=0.0;
    dy=0.0;
    dz=0.0;
    for(i=1;i<=natm_tot;i++){
      dx += x[i];
      dy += y[i];
      dz += z[i];
    }
    dx /= ((double)natm_tot);
    dy /= ((double)natm_tot);
    dz /= ((double)natm_tot);

    for(i=1;i<=natm_tot;i++){
      x[i] -= dx;
      y[i] -= dy;
      z[i] -= dz;
    }
  }/*endfor*/

}

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void gen_Ts_dvr(DVR_MATRIX *dvr_matrix,CPCOEFFS_INFO *cpcoeffs_info,
                CELL *cell,int iperd,int myid)

/*========================================================================*/
/*             Begin Routine                                              */
 { /*Begin subprogram: */
/*========================================================================*/
/*             Local variable declarations                                */

  int grid_nx  = cpcoeffs_info->grid_nx;
  int grid_ny  = cpcoeffs_info->grid_ny;
  int grid_nz  = cpcoeffs_info->grid_nz;

  int cubic_box_flag = cell->cubic_box_flag;
  double *hmat       = cell->hmat;
  double boxl_x,boxl_y,boxl_z;

  double *Tsx   = dvr_matrix->Tsx;
  double *Tsy   = dvr_matrix->Tsy;
  double *Tsz   = dvr_matrix->Tsz;

  int grid_nx2 = (grid_nx - 1)/2;
  int grid_ny2 = (grid_ny - 1)/2;
  int grid_nz2 = (grid_nz - 1)/2;

  double  pi   = M_PI;
  double fpi2  = 4.0*pi*pi;

  double nogx1 = (double)(grid_nx + 1);
  double nogy1 = (double)(grid_ny + 1);
  double nogz1 = (double)(grid_nz + 1);
  double pd2nx = pi/(2.0*nogx1);
  double pd2ny = pi/(2.0*nogy1);
  double pd2nz = pi/(2.0*nogz1);
  double pn_x = pi/(double)grid_nx;
  double pn_y = pi/(double)grid_ny;
  double pn_z = pi/(double)grid_nz;

  double pdl2_x,pdl2_y,pdl2_z;
  double tdd_x,tdd_y,tdd_z;

  int i,j,k,k1;
  int ind;
  double dk;
  double dx,dy,dz;

  if(myid == 0){
    printf("\n");PRINT_LINE_STAR;
    printf("Setting up DVR related matrices\n");
    PRINT_LINE_DASH;printf("\n");

    printf("  SETTING UP KINETIC ENERGY MATRICES\n");

  }

/*===========================================================================*/
/* 0) Define constants  */

  if(cubic_box_flag == 1){
    boxl_x = hmat[1];
    boxl_y = hmat[5];
    boxl_z = hmat[9];
    dx = boxl_x/((double) grid_nx);
    dy = boxl_y/((double) grid_ny);
    dz = boxl_z/((double) grid_nz);
  }else{
    boxl_x = sqrt(hmat[1]*hmat[1] + hmat[4]*hmat[4] + hmat[7]*hmat[7]);
    boxl_y = sqrt(hmat[2]*hmat[2] + hmat[5]*hmat[5] + hmat[8]*hmat[8]);
    boxl_z = sqrt(hmat[3]*hmat[3] + hmat[6]*hmat[6] + hmat[9]*hmat[9]);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("DVR has not been implemented for a NON-CUBIC Box \n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }

  pdl2_x  =  fpi2/(boxl_x*boxl_x);
  pdl2_y  =  fpi2/(boxl_y*boxl_y);
  pdl2_z  =  fpi2/(boxl_z*boxl_z);
  tdd_x   = -pdl2_x*(double)grid_nx2*((double)grid_nx2 + 1.0)/3.0;
  tdd_y   = -pdl2_y*(double)grid_ny2*((double)grid_ny2 + 1.0)/3.0;
  tdd_z   = -pdl2_z*(double)grid_nz2*((double)grid_nz2 + 1.0)/3.0;

/*========================================================================*/
/* I) set up kinetic energy matrices  */

  switch(iperd) {
    case 3: /* periodic T matrix */
      for(i=1;i<= grid_nx;i++) {
        ind = (i-1)*grid_nx + i;  /* Diagonal Elements */
        Tsx[ind] = tdd_x;
        for (j=1;j <= grid_nx;j++) {
          if(i!=j) {
            k = i-j;
            dk = (double)k;
            ind = (i-1)*grid_nx + j;  /* OFF-Diagonal Elements */
            Tsx[ind] = pdl2_x*cos(pn_x*dk)/(2.0*sin(pn_x*dk)*sin(pn_x*dk));
            if(k%2==0) Tsx[ind] *= -1.0;
          }/*endif*/
        }/*endfor*/
      }/*endfor*/

      for(i=1;i<= grid_ny;i++) {
        ind = (i-1)*grid_ny + i;  /* Diagonal Elements */
        Tsy[ind] = tdd_y;
        for (j=1;j<= grid_ny;j++) {
          if(i!=j) {
            k = i-j;
            dk = (double)k;
            ind = (i-1)*grid_ny + j;  /* OFF-Diagonal Elements */
            Tsy[ind] = pdl2_y*cos(pn_y*dk)/(2.0*sin(pn_y*dk)*sin(pn_y*dk));
            if(k%2==0) Tsy[ind] *= -1.0;
          }/*endif*/
        }/*endfor*/
      }/*endfor*/

      for(i=1;i<= grid_nz;i++) {
        ind = (i-1)*grid_nz + i;  /* Diagonal Elements */
        Tsz[ind] = tdd_z;
        for (j=1;j <= grid_nz;j++) {
          if(i!=j) {
            k = i-j;
            dk = (double)k;
            ind = (i-1)*grid_nz + j;  /* OFF-Diagonal Elements */
            Tsz[ind] = pdl2_z*cos(pn_z*dk)/(2.0*sin(pn_z*dk)*sin(pn_z*dk));
            if(k%2==0) Tsz[ind] *= -1.0;
          }/*endif*/
        }/*endfor*/
      }/*endfor*/

    break;

    case 0: /* cluster T matrices */

      for(i=1;i<= grid_nx;i++) {
        ind = (i-1)*grid_nx + i;  /* Diagonal Elements */
        Tsx[ind] = -pi*pi/(3.0*dx*dx);
        for (j=1;j <= grid_nx;j++) {
          if(i!=j) {
            k = i-j;
            dk = (double)k;
            ind = (i-1)*grid_nx + j;  /* OFF-Diagonal Elements */
            Tsx[ind] = 2.0/(dk*dk*dx*dx);
            if(k%2==0) Tsx[ind] *= -1.0;
          }/*endif*/
        }/*endfor*/
      }/*endfor*/

      for(i=1;i<= grid_ny;i++) {
        ind = (i-1)*grid_ny + i;  /* Diagonal Elements */
        Tsy[ind] = -pi*pi/(3.0*dy*dy);;
        for (j=1;j<= grid_ny;j++) {
          if(i!=j) {
            k = i-j;
            dk = (double)k;
            ind = (i-1)*grid_ny + j;  /* OFF-Diagonal Elements */
            Tsy[ind] = 2.0/(dk*dk*dy*dy);
            if(k%2==0) Tsy[ind] *= -1.0;
          }/*endif*/
        }/*endfor*/
      }/*endfor*/

      for(i=1;i<= grid_nz;i++) {
        ind = (i-1)*grid_nz + i;  /* Diagonal Elements */
        Tsz[ind] = -pi*pi/(3.0*dz*dz);
        for (j=1;j <= grid_nz;j++) {
          if(i!=j) {
            k = i-j;
            dk = (double)k;
            ind = (i-1)*grid_nz + j;  /* OFF-Diagonal Elements */
            Tsz[ind] =  2.0/(dk*dk*dz*dz);
            if(k%2==0) Tsz[ind] *= -1.0;
          }/*endif*/
        }/*endfor*/
      }/*endfor*/

    break;

    if(myid==0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("Periodicity other than 0 or 3 is not allowed yet \n");
      printf("for CP calculations with DVR basis\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    }/*endif*/
    fflush(stdout);
    exit(1);

  }/*end switch*/

/*------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void gen_dgrid_dvr(DVR_MATRIX *dvr_matrix, CPCOEFFS_INFO *cpcoeffs_info,
                   CELL *cell,int iperd,int myid)

/*========================================================================*/
/*             Begin Routine                                              */
 { /*Begin subprogram: */
/*========================================================================*/
/*             Local variable declarations                                */

  int grid_nx  = cpcoeffs_info->grid_nx;
  int grid_ny  = cpcoeffs_info->grid_ny;
  int grid_nz  = cpcoeffs_info->grid_nz;

  int cubic_box_flag = cell->cubic_box_flag;
  double boxl_x,boxl_y,boxl_z,dx,dy,dz;

  double *dGx   = dvr_matrix->dGx;
  double *dGy   = dvr_matrix->dGy;
  double *dGz   = dvr_matrix->dGz;

  double pi  = M_PI;
  double pnx = pi/(double)grid_nx;
  double pny = pi/(double)grid_ny;
  double pnz = pi/(double)grid_nz;

  double nogx1 = (double)(grid_nx+1);
  double nogy1 = (double)(grid_ny+1);
  double nogz1 = (double)(grid_nz+1);

  double pd2nx = pi/(2.0*nogx1);
  double pd2ny = pi/(2.0*nogy1);
  double pd2nz = pi/(2.0*nogz1);

  int i,j,k,k1;
  int ind;
  double dk;

  if(myid == 0){
    printf("  SETTING UP MATRICES FOR THE GRADIENTS OF DVR FUNCTIONS\n");
  }

/*========================================================================*/

  if(cubic_box_flag == 1){
    boxl_x = cell->hmat[1];
    boxl_y = cell->hmat[5];
    boxl_z = cell->hmat[9];
    dx = boxl_x/((double) grid_nx);
    dy = boxl_y/((double) grid_ny);
    dz = boxl_z/((double) grid_nz);
  }else{
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("DVR has not been implemented for a NON-CUBIC Box \n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

/*========================================================================*/
/* I) DERIVATIVE MATRIX FOR PBC */

  switch(iperd){
    case 3: 
      for(i=1;i<=grid_nx;i++) {
        ind = (i-1)*grid_nx + i; /*DIAGONAL ELEMENTS*/
        dGx[ind] = 0.0;
        for(j=1;j<=grid_nx;j++) {
          if(i!=j) {
            k = j-i;
            ind = (i-1)*grid_nx + j; /*OFF-DIAGONAL ELEMENTS*/
            dGx[ind] = M_PI/(sin(pnx*k)*boxl_x);
            if(k%2!=0) dGx[ind] *= -1.0;
          }/*endif*/
        }/*endfor*/
      }/*endfor*/

      for(i=1;i<=grid_ny;i++) {
        ind = (i-1)*grid_ny + i; /*DIAGONAL ELEMENTS*/
        dGy[ind] = 0.0;
        for(j=1;j<=grid_ny;j++) {
          if(i!=j) {
            k = j-i;
            ind = (i-1)*grid_ny + j; /*OFF-DIAGONAL ELEMENTS*/
            dGy[ind] = pi/(sin(pny*k)*boxl_y);
            if(k%2!=0) dGy[ind] *= -1.0;
          }/*endif*/
        }/*endfor*/
      }/*endfor*/

      for(i=1;i<=grid_nz;i++) {
        ind = (i-1)*grid_nz + i; /* DIAGONAL ELEMENTS*/
        dGz[ind] = 0.0;
        for(j=1;j<=grid_nz;j++) {
          if(i!=j) {
            k = j-i;
            ind = (i-1)*grid_nz + j; /*OFF-DIAGONAL ELEMENTS*/
            dGz[ind] = pi/(sin(pnz*k)*boxl_z);
            if(k%2!=0) dGz[ind] *= -1.0;
          }/*endif*/
        }/*endfor*/
      }/*endfor*/

    break;

    case 0:
      for(i=1;i<=grid_nx;i++) {
        ind = (i-1)*grid_nx + i; /*DIAGONAL ELEMENTS*/
        dGx[ind] = 0.0;
        for(j=1;j<=grid_nx;j++) {
          if(i!=j) {
            k = j-i;
            dk = ((double) k);
            ind = (i-1)*grid_nx + j; /*OFF-DIAGONAL ELEMENTS*/
            dGx[ind] = cos(pi*dk)/(dk*dx);
          }/*endif*/
        }/*endfor*/
      }/*endfor*/

      for(i=1;i<=grid_ny;i++) {
        ind = (i-1)*grid_ny + i; /*DIAGONAL ELEMENTS*/
        dGy[ind] = 0.0;
        for(j=1;j<=grid_ny;j++) {
          if(i!=j) {
            k = j-i;
            dk = ((double) k);
            ind = (i-1)*grid_ny + j; /*OFF-DIAGONAL ELEMENTS*/
            dGy[ind] = cos(pi*dk)/(dk*dy);
          }/*endif*/
        }/*endfor*/
      }/*endfor*/

      for(i=1;i<=grid_nz;i++) {
        ind = (i-1)*grid_nz + i; /* DIAGONAL ELEMENTS*/
        dGz[ind] = 0.0;
        for(j=1;j<=grid_nz;j++) {
          if(i!=j) {
            k = j-i;
            dk = ((double) k);
            ind = (i-1)*grid_nz + j; /*OFF-DIAGONAL ELEMENTS*/
            dGz[ind] = cos(pi*dk)/(dk*dz);
          }/*endif*/
        }/*endfor*/
      }/*endfor*/

    break;
 }/*end switch*/

/*------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void gen_dvr_fbr(DVR_MATRIX *dvr_matrix, CPCOEFFS_INFO *cpcoeffs_info,
                 CELL *cell, CPSCR_WAVE *cpscr_wave,
                 int iperd,int myid, double cmass_cut_def)

/*========================================================================*/
/*             Begin Routine                                              */
 { /*Begin subprogram: */
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int grid_nx  = cpcoeffs_info->grid_nx;
  int grid_ny  = cpcoeffs_info->grid_ny;
  int grid_nz  = cpcoeffs_info->grid_nz;
  int grid_nx2 = (grid_nx-1)/2;
  int grid_ny2 = (grid_ny-1)/2;
  int grid_nz2 = (grid_nz-1)/2;
  int cubic_box_flag = cell->cubic_box_flag;

  double boxl_x,boxl_y,boxl_z;

  double *Tsx   = dvr_matrix->Tsx;
  double *Tsy   = dvr_matrix->Tsy;
  double *Tsz   = dvr_matrix->Tsz;

  double *TRx   = dvr_matrix->TRx;
  double *TRy   = dvr_matrix->TRy;
  double *TRz   = dvr_matrix->TRz;
  double *TIx   = dvr_matrix->TIx;
  double *TIy   = dvr_matrix->TIy;
  double *TIz   = dvr_matrix->TIz;

  double *Tfbr_x = dvr_matrix->Tfbr_x;
  double *Tfbr_y = dvr_matrix->Tfbr_y;
  double *Tfbr_z = dvr_matrix->Tfbr_z;

  double *kin_mat,*scr1,*scr2;

  double grid_x,grid_y,grid_z,dx,dy,dz;
  double pi;
  double tpdlx,tpdly,tpdlz;
  double wght_x,wght_y,wght_z;

  int i,j,k,n,id,ngrid_max,ierr,ind;
  int job = 1;

  if(myid == 0){
    printf("  SETTING UP DVR-FBR TRANSFORMATION MATRICES\n");
  }

/*=========================================================================*/
/* 0) define constants */

  if(cubic_box_flag == 1){
    boxl_x = cell->hmat[1];
    boxl_y = cell->hmat[5];
    boxl_z = cell->hmat[9];
  }else{
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("DVR has not been implemented for a NON-CUBIC Box \n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

  if(grid_nx > grid_ny) {
    ngrid_max = grid_nx;
  }else{
    ngrid_max = grid_ny;
  }
  if(grid_nz > ngrid_max) ngrid_max = grid_nz;

  pi  = M_PI;
  tpdlx = 2.0*pi/boxl_x;
  tpdly = 2.0*pi/boxl_y;
  tpdlz = 2.0*pi/boxl_z;

  wght_x= 1.0/sqrt((double) grid_nx);
  wght_y= 1.0/sqrt((double) grid_ny);
  wght_z= 1.0/sqrt((double) grid_nz);

/*========================================================================*/
/* I) transformation matrices */

   dx = boxl_x/((double) grid_nx);
   dy = boxl_y/((double) grid_ny);
   dz = boxl_z/((double) grid_nz);

  switch(iperd) {
    case 3 :
      id=0;
      for (n=1;n<=grid_nx;n++){
        k=-grid_nx2+n-1;
        for (j=1;j<=grid_nx;j++){
          id=id+1;
          grid_x=-0.5*boxl_x+(j-1)*dx+dx/2.0;
          TRx[id]=cos(tpdlx*k*grid_x)*wght_x;
          TIx[id]=sin(tpdlx*k*grid_x)*wght_x;
        }

        if (abs(k) > cmass_cut_def) {
          Tfbr_x[n]= ((double) tpdlx*tpdlx*k*k)*0.5;
        }
        else{
          Tfbr_x[n]= ((double)tpdlx*tpdlx*cmass_cut_def*cmass_cut_def)*0.5;
        }
      }

      id=0;
      for (n=1;n<=grid_ny;n++){
        k=-grid_ny2+n-1;
        for (j=1;j<=grid_ny;j++){
          id=id+1;
          grid_y=-0.5*boxl_y+(j-1)*dy+dy/2.0;
          TRy[id]=cos(tpdly*k*grid_y)*wght_y;
          TIy[id]=sin(tpdly*k*grid_y)*wght_y;
        }
        if (abs(k) > cmass_cut_def) {
          Tfbr_y[n]=((double) tpdly*tpdly*k*k)*0.5;
        }
        else{
          Tfbr_y[n]=((double) tpdly*tpdly*cmass_cut_def*cmass_cut_def)*0.5;
        }
      }

      id=0;
      for (n=1;n<=grid_nz;n++){
        k=-grid_nz2+n-1;
        for (j=1;j<=grid_nz;j++){
          id=id+1;
          grid_z=-0.5*boxl_z+(j-1)*dz+dz/2.0;
          TRz[id]=cos(tpdlz*k*grid_z)*wght_z;
          TIz[id]=sin(tpdlz*k*grid_z)*wght_z;
        }
        if (abs(k) > cmass_cut_def) {
          Tfbr_z[n]=((double) tpdlz*tpdlz*k*k)*0.5;
        }
        else{
          Tfbr_z[n]=((double) tpdlz*tpdlz*cmass_cut_def*cmass_cut_def)*0.5;
        }
      }

    break;

    case 0:

      kin_mat = cpscr_wave->cre_up;
      scr1    = cpscr_wave->cim_up;
      scr2    = &(cpscr_wave->cim_up[ngrid_max]);
      for(i=1;i<=grid_nx;i++) {
        for(j=1;j<=grid_nx;j++) {
          ind = (i-1)*grid_nx + j;
          kin_mat[ind] = -Tsx[ind];
        }
      }

      RS(&grid_nx,&grid_nx,&(kin_mat[1]),&(Tfbr_x[1]),&job,&(TRx[1]),
         &(scr1[1]),&(scr2[1]),&ierr);

      for(i=1;i<=grid_ny;i++) {
        for(j=1;j<=grid_ny;j++) {
          ind = (i-1)*grid_ny + j;
          kin_mat[ind] = -Tsy[ind];
        }
      }

      RS(&grid_ny,&grid_ny,&(kin_mat[1]),&(Tfbr_y[1]),&job,&(TRy[1]),
         &(scr1[1]),&(scr2[1]),&ierr);

      for(i=1;i<=grid_nz;i++) {
        for(j=1;j<=grid_nz;j++) {
          ind = (i-1)*grid_nz + j;
          kin_mat[ind] = -Tsz[ind];
        }
      }

      RS(&grid_nz,&grid_nz,&(kin_mat[1]),&(Tfbr_z[1]),&job,&(TRz[1]),
         &(scr1[1]),&(scr2[1]),&ierr);

      ind = (int)(cmass_cut_def);
      if(ind > grid_nx) ind = grid_nx;
      id = ind+1;
      if(id > grid_nx) id = grid_nx;
      for(i=1;i<=ind;i++){
        Tfbr_x[i]=Tfbr_x[id];
      }

      ind = (int)(cmass_cut_def);
      if(ind > grid_ny) ind = grid_ny;
      id = ind+1;
      if(id > grid_ny) id = grid_ny;
      for(i=1;i<=ind;i++){
        Tfbr_y[i]=Tfbr_y[id];
      }

      ind = (int)(cmass_cut_def);
      if(ind > grid_nz) ind = grid_nz;
      id = ind+1;
      if(id > grid_nz) id = grid_nz;
      for(i=1;i<=ind;i++){
        Tfbr_z[i]=Tfbr_z[id];
      }
    break;

  } /*end case */

/*======================================================================*/
}/* end of subprogram gen_dvr_fbr*/
/*======================================================================*/

