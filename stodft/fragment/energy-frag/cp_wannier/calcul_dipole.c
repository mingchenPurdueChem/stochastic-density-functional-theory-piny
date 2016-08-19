/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: calcul_dipole                                */
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
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_wannier_cpcon_local.h"
#include "../proto_defs/proto_wannier_cpcon_entry.h"

#define Debye 4.8033  /* 1e*Angstrom=4.8033 Debyes */

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calcul_dipole(CLASS *class, GENERAL_DATA *general_data,CP *cp,
                   double *W, double ***Z_real, double ***Z_imag)

/*===========================================================================*/
{ /*begin routine */
/*============================================================================*/
/*            Local variable declarations                                    */
/*===========================================================================*/

#include "../typ_defs/typ_mask.h"

  int myid         = cp->communicate.myid_state;
  int num_proc     = cp->communicate.np_states;
  MPI_Comm comm    = cp->communicate.comm_states;

  int nstate_up    = cp->cpcoeffs_info.nstate_up;
  int nstate       = nstate_up;
  int tot_beads    = cp->cpcoeffs_info.pi_beads;
  int natm         = class->clatoms_info.natm_tot;

  double vol_cp    = general_data->cell.vol_cp;
  double *hmati_cp = general_data->cell.hmati_cp;
  double *hmat_cp  = general_data->cell.hmat_cp;
  double *x        = class->clatoms_pos[1].x;
  double *y        = class->clatoms_pos[1].y;
  double *z        = class->clatoms_pos[1].z;
  double *q        = class->clatoms_info.q;

  double **dipole_cp = cp->electronic_properties.cp_dip;
  double *HMat_1     = cp->electronic_properties.HMat_1;
  double *HMat_2     = cp->electronic_properties.HMat_2;
  double *HMat_3     = cp->electronic_properties.HMat_3;
  double *GMat_1     = cp->electronic_properties.GMat_1;
  double *GMat_2     = cp->electronic_properties.GMat_2;
  double *GMat_3     = cp->electronic_properties.GMat_3;

  int    iopt_cp_pw  = cp->cpcoeffs_info.iopt_cp_pw;
  int    iopt_cp_dvr = cp->cpcoeffs_info.iopt_cp_dvr;

  int proc_beads,st_beads,end_beads;
  int ip,is,js,i,j,iflag;
  double a1,a2,a3;
  double P_elec_x,P_elec_y,P_elec_z,P_nucl_x,P_nucl_y,P_nucl_z,P_tot_x,P_tot_y,P_tot_z;

  double Pi = M_PI;

  zomplex DET[2];               /* the determinant = Det[0]*10^Det[1]      */
  
  zomplex **S1, **S2, **S3; 
  zomplex *WORK;                /* needed to communicate between zsifa and zsidi */

  int IPVT[nstate],INFO;        /* INFO != 0 signals possible error              */
  int N=nstate;
  int JOB=10;                   /* 10 signifies calculation of the determinant ONLY */
  int LDA=nstate;

/*--------------------------------------------------------------------------*/
/* 0) Simulation check and initialization */

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
  proc_beads=st_beads=end_beads=1;

  S1 = cmall_zomp2(0,nstate-1,0,nstate-1);
  S2 = cmall_zomp2(0,nstate-1,0,nstate-1);
  S3 = cmall_zomp2(0,nstate-1,0,nstate-1);
  WORK = (zomplex *) malloc((size_t)(nstate*sizeof(zomplex)));

  for (i=0;i<nstate;i++){
    for (j=0;j<nstate;j++){
      S1[i][j].re=0.0;
      S1[i][j].im=0.0;
      S2[i][j].re=0.0;
      S2[i][j].im=0.0;
      S3[i][j].re=0.0;
      S3[i][j].im=0.0;
    }
  }

/*--------------------------------------------------------------------------*/
  for (ip=1; ip<=proc_beads; ip++){  /* for each bead*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* I) compute the z_tensor   */


    if(iopt_cp_pw == 1){
      comp_tensor_z(general_data,cp,W,Z_real,Z_imag);
    }

    if(iopt_cp_dvr == 1){
      comp_tensor_z_dvr(general_data,cp,Z_real,Z_imag);
    }

    if(myid==0){
/*-------------------------------------------------------------------------------*/
/* II) Root processor compute the S-matrix                                        */
/* Note: Z-maxtir use e^{igr}, but S-maxtir use e^{-igr}                         */

      for (is=0;is<nstate;is++) {
        for (js=is;js<nstate;js++){
          S1[js][is].re= Z_real[1][is+1][js+1];
          S1[is][js].re= Z_real[1][is+1][js+1];
          S1[js][is].im= -Z_imag[1][is+1][js+1];
          S1[is][js].im= -Z_imag[1][is+1][js+1];
          S2[js][is].re= Z_real[2][is+1][js+1];
          S2[is][js].re= Z_real[2][is+1][js+1];
          S2[js][is].im= -Z_imag[2][is+1][js+1];
          S2[is][js].im= -Z_imag[2][is+1][js+1];
          S3[js][is].re= Z_real[3][is+1][js+1];
          S3[is][js].re= Z_real[3][is+1][js+1];
          S3[js][is].im= -Z_imag[3][is+1][js+1];
          S3[is][js].im= -Z_imag[3][is+1][js+1];
        }
      }

/*-------------------------------------------------------------------------------*/
/* III) Root processor compute the determinant of S-matrix using Linpack routines  */
/* We need to compute Im Ln det S. If we write det S, which is a complex number, */
/* in a polar form, Im Ln is just the phase angle of det S                       */

      ZSIFA(&S1[0][0],&LDA,&N,&IPVT[0],&INFO);
      if (INFO != 0){
        printf("Divide by Zero Error in zsidi for X-component!\n");
        fflush(stdout);
        Finalize();
        exit(1);
      }
      ZSIDI(&S1[0][0],&LDA,&N,&IPVT[0],&DET[0],&WORK[0],&JOB);

      P_elec_x =atan(-DET[0].im/DET[0].re);
      if (P_elec_x<=0.0) P_elec_x+=Pi;
      P_elec_x*=2.0;

      ZSIFA(&S2[0][0],&LDA,&N,&IPVT[0],&INFO);
      if (INFO!=0){
        printf("Divide by Zero Error in zsidi for Y-component!\n");
        fflush(stdout);
        Finalize(); exit(1);
      }
      ZSIDI(&S2[0][0],&LDA,&N,&IPVT[0],&DET[0],&WORK[0],&JOB);

      P_elec_y=atan(-DET[0].im/DET[0].re);
      if (P_elec_y<=0.0) P_elec_y+=Pi;
      P_elec_y*=2.0;

      ZSIFA(&S3[0][0],&LDA,&N,&IPVT[0],&INFO);
      if (INFO!=0){
        printf("Divide by Zero Error in zsidi for Z-component!\n");
        fflush(stdout);
        Finalize();
        exit(1);
      }
      ZSIDI(&S3[0][0],&LDA,&N,&IPVT[0],&DET[0],&WORK[0],&JOB);

      P_elec_z=atan(-DET[0].im/DET[0].re);
      if (P_elec_z<=0.0) P_elec_z+=Pi;
      P_elec_z*=2.0;

/*------------------------------------------------------------------------------------*/
/* IV)Root processor compute the dipole moment of Nuclei                               */

    /*------------------------------------------------------------------------*/
    /* A. First, compute the h ang G vectors                                  */
    /* Recall: basis vectors for g-space, G = 2pi x [h_1,h_2,h_3]^{-T}, where */
    /* each column vector of matrix G=[G_1,G_2,G_3] is the basis for g-space  */
    /*------------------------------------------------------------------------*/

      GMat_1[1]=(2.0*Pi)/BOHR*hmati_cp[1];
      GMat_1[2]=(2.0*Pi)/BOHR*hmati_cp[2];
      GMat_1[3]=(2.0*Pi)/BOHR*hmati_cp[3];

      GMat_2[1]=(2.0*Pi)/BOHR*hmati_cp[4];
      GMat_2[2]=(2.0*Pi)/BOHR*hmati_cp[5];
      GMat_2[3]=(2.0*Pi)/BOHR*hmati_cp[6];

      GMat_3[1]=(2.0*Pi)/BOHR*hmati_cp[7];
      GMat_3[2]=(2.0*Pi)/BOHR*hmati_cp[8];
      GMat_3[3]=(2.0*Pi)/BOHR*hmati_cp[9];

      HMat_1[1]=BOHR*hmat_cp[1];
      HMat_1[2]=BOHR*hmat_cp[4];
      HMat_1[3]=BOHR*hmat_cp[7];

      HMat_2[1]=BOHR*hmat_cp[2];
      HMat_2[2]=BOHR*hmat_cp[5];
      HMat_2[3]=BOHR*hmat_cp[8];

      HMat_3[1]=BOHR*hmat_cp[3];
      HMat_3[2]=BOHR*hmat_cp[6];
      HMat_3[3]=BOHR*hmat_cp[9];

    /*---------------------------------------------------------------------*/
    /* B. Express the positon of each nuclei in terms of cell matrix       */
    /* i.e r = a x h_1 + b x h_2 + c x h_3 where h_1,h_2,h_3 are           */
    /*     the comlum vector of cell matrix                                */
    /* (a,b,c) are easily computed via  a = G_1*r/2pi                      */
    /* Then compute the dipole moment of nuclei                            */
    /*---------------------------------------------------------------------*/

      P_nucl_x=0.0;
      P_nucl_y=0.0;
      P_nucl_z=0.0;
      for (i=1;i<=natm;i++){
        a1=BOHR*(x[i]*GMat_1[1]+y[i]*GMat_1[2]+z[i]*GMat_1[3]);
        a2=BOHR*(x[i]*GMat_2[1]+y[i]*GMat_2[2]+z[i]*GMat_2[3]);
        a3=BOHR*(x[i]*GMat_3[1]+y[i]*GMat_3[2]+z[i]*GMat_3[3]);
        P_nucl_x+=q[i]*a1;
        P_nucl_y+=q[i]*a2;
        P_nucl_z+=q[i]*a3;
      }

      iflag = 1;
      while (iflag) {
        if (0<=P_nucl_x && P_nucl_x<=2.0*Pi) iflag=0;
        else if (2.0*Pi<P_nucl_x) P_nucl_x -= 2.0*Pi;
        else P_nucl_x += 2.0*Pi;
      }

      iflag = 1;
      while (iflag) {
        if (0<=P_nucl_y && P_nucl_y<=2.0*Pi) iflag=0;
        else if (2.0*Pi<P_nucl_y) P_nucl_y -= 2.0*Pi;
        else P_nucl_y += 2.0*Pi;
      }

      iflag = 1;
      while (iflag) {
        if (0<=P_nucl_z && P_nucl_z<=2.0*Pi) iflag=0;
        else if (2.0*Pi<P_nucl_z) P_nucl_z -= 2.0*Pi;
        else P_nucl_z += 2.0*Pi;
      }

    /*------------------------------------------------------------*/
    /* C. Now, compute the Total dipole with the correct factor   */
    /*------------------------------------------------------------*/

      P_tot_x = P_nucl_x-P_elec_x;
      P_tot_y = P_nucl_y-P_elec_y;
      P_tot_z = P_nucl_z-P_elec_z;


      iflag=1;
      while (iflag) {
        if (0<=P_tot_x && P_tot_x<=2.0*Pi) iflag=0;
        else if (2.0*Pi<P_tot_x) P_tot_x -= 2.0*Pi;
        else P_tot_x += 2.0*Pi;
      }
      if(P_tot_x > Pi) P_tot_x -= 2.0*Pi;

      iflag=1;
      while (iflag) {
        if (0<=P_tot_y && P_tot_y<=2.0*Pi) iflag=0;
        else if (2.0*Pi<P_tot_y) P_tot_y -= 2.0*Pi;
        else P_tot_y += 2.0*Pi;
      }
      if(P_tot_y > Pi) P_tot_y -= 2.0*Pi;

      iflag=1;
      while (iflag) {
        if (0<=P_tot_z && P_tot_z<=2.0*Pi) iflag=0;
        else if (2.0*Pi<P_tot_z) P_tot_z -= 2.0*Pi;
        else P_tot_z += 2.0*Pi;
      }
      if(P_tot_z > Pi) P_tot_z -= 2.0*Pi;

      dipole_cp[ip+st_beads-1][1]=
          (HMat_1[1]*P_tot_x+HMat_2[1]*P_tot_y+HMat_3[1]*P_tot_z)*Debye/(2.0*Pi);
      dipole_cp[ip+st_beads-1][2]=
          (HMat_1[2]*P_tot_x+HMat_2[2]*P_tot_y+HMat_3[2]*P_tot_z)*Debye/(2.0*Pi);
      dipole_cp[ip+st_beads-1][3]=
          (HMat_1[3]*P_tot_x+HMat_2[3]*P_tot_y+HMat_3[3]*P_tot_z)*Debye/(2.0*Pi);
    }/*endif myid=0*/

/*-----------------------------------------------------------------------------*/
  }/* end for ip: bead */
/*-----------------------------------------------------------------------------*/

  if(num_proc > 1){
    Barrier(comm);
  }

/*===========================================================================*/
} /*end routine */
/*============================================================================*/

