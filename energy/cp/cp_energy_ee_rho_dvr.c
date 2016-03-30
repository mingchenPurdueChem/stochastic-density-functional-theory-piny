/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: control_energy_ee_rho.c                        */
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

#define  CHECK_HARTREE_OFF

/*==========================================================================*/
/* calculate the Hartree energy and forces on the coefficients for DVR CP   */
/*==========================================================================*/

void calc_hartree_dvr_pbc(CPOPTS *cpopts,CPEWALD *cpewald,double *rhocr,double *rhoci,
                          double *zfft,double *zfft_tmp, double *fcre,double *fcim,
                          double *dvrfc_up, double *dvrc_up,double *dvrfc_dn,
                          double *dvrc_dn, double *eh_ret,int nstate_up, int nstate_dn,
                          EWALD *ewald,CELL *cell,COMMUNICATE *communicate,
                          int nkf1,int nkf2,int nkf3, 
                          int nfft,int nfft_proc,int *fft_map,int *fft_map_rev,
                          PARA_FFT_PKG3D *cp_para_fft_pkg3d_sm)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"


   int       cp_lsda         =    cpopts->cp_lsda;

   double    *ak2            =    cpewald->ak2;
   int       *kastore        =    ewald->kastr;
   int       *kbstore        =    ewald->kbstr;
   int       *kcstore        =    ewald->kcstr;

   int       nfft2           =    nfft/2;
   int       nfft2_proc      =    nfft_proc/2;

   double    *hmat           =    cell->hmat;

   int       myid_state      =    communicate->myid_state;
   int       np_states       =    communicate->np_states;
   MPI_Comm  comm            =    communicate->comm_states;

   double *fcre_scr,*fcim_scr,*fscr;
   int *irecv,*idispl;

/*     Local Variable declarations         */

   double ghfact,pi,fpi,tpi,dfact,vol,eh;
   double aka,akb,akc,xk,yk,zk;
   int i,j,k,is,ind;
   int index,irem,idiv,sign,ngo;

/*--------------------------------------------------------------------*/
/*  0) Set local variables */

   ngo = ( (myid_state + 1) == np_states ? nfft2_proc - 1: nfft2_proc);

   irecv  = (int *) malloc((np_states)*sizeof(int))-1;
   idispl = (int *) malloc((np_states)*sizeof(int))-1;

   idiv    = (nkf2*nkf3)/np_states;
   irem    = (nkf2*nkf3)%np_states;

   for(i=1; i <= np_states; i++){
     irecv[i] = ( (i-1) < irem ? (idiv+1)*nkf1 : (idiv)*nkf1);
   }

   idispl[1] = 0;
   for(i=2; i <= np_states; i++){
    idispl[i] = idispl[(i-1)] + irecv[(i-1)];
   }

   pi     = M_PI; fpi = 4.0*pi; tpi = 2.0*pi;
   vol    = getdeth(hmat);
   dfact  = vol/((double) nfft2);

/*----------------------------------------------------------------------------*/
/* I) Compute hartree energy and force contribution */

   eh   = 0.0;

   for(i=1;i<= ngo; i++){
     index = i;
     ghfact     = fpi/(ak2[i]*vol);
     eh       += 0.50*ghfact*(rhocr[index]*rhocr[index] + rhoci[index]*rhoci[index]);

     /* the += is used since fc(re/im) also contains the local long range piece */
     fcre[index]   += ghfact*rhocr[index];
     fcim[index]   += ghfact*rhoci[index];
   }/*endfor*/

   *eh_ret += eh;

/*--------------------------------------------------------------------------*/
/* II) Gather forces across processors */

  if(np_states > 1){ 

    Gatherv(&fcre[1],nfft2_proc,MPI_DOUBLE,&zfft[1],
            &irecv[1],&idispl[1],MPI_DOUBLE,0,comm);

    if(myid_state == 0){
      fcre_scr = zfft_tmp;
      for(i=1; i <= nfft2; i++){
        index = fft_map[i];
        fcre_scr[index] = zfft[i];
      }
    }

    Gatherv(&fcim[1],nfft2_proc,MPI_DOUBLE,&zfft[1],
            &irecv[1],&idispl[1],MPI_DOUBLE,0,comm);

    if(myid_state == 0){
      fcim_scr = zfft_tmp + nfft2;
      for(i=1; i <= nfft2; i++){
        index = fft_map[i];
        fcim_scr[index] = zfft[i];
      }

      for(i=1; i <= nfft2; i++){
        zfft[2*i-1] = fcre_scr[i];
        zfft[2*i]   = fcim_scr[i];
      }
    }/*endif myid */
  }else{ /*Serial*/  

    fcre_scr = zfft_tmp;

    for(i=1; i <= nfft2; i++){
      index = fft_map[i];
      fcre_scr[index] = fcre[i];
    }
    for(i=1; i <= nfft2; i++){ 
      zfft[2*i-1] = fcre_scr[i]; 
    }
 
    fcim_scr = zfft_tmp + nfft2;

    for(i=1; i <= nfft2; i++){
      index = fft_map[i];
      fcim_scr[index] = fcim[i];
    }
    for(i=1; i <= nfft2; i++){ 
      zfft[2*i] = fcim_scr[i]; 
    }

  }/* endif SERIAL */

/*-----------------------------------------------------------------------*/
/* III) Fourier Transform the force to real-space */
 
  /* Always serial FFT for current implementation */
  if(myid_state == 0){ 

    sign = 1;
    para_fft_gen3d_dvr_bck(zfft,zfft_tmp,cp_para_fft_pkg3d_sm,-1);

    /* unpacking the force */
    fscr = fcre;

    for(i=1; i <= nfft2/2; i++){
      index = fft_map[i];
      fscr[i] = zfft[2*index-1];
    }

    /* G=0 Term located at end of array in reciprocal space */
    /*   maps to middle of array in real space */

    index = fft_map[nfft2];
    fscr[i] = zfft[2*index-1];

    /* G array shifted down by one to move G=0 term to end of array*/
    /*  want fft_map element i-1 for this half of array */

    for(i=nfft2/2 + 2 ; i<= nfft2; i++){
      index = fft_map[i-1];
      fscr[i] = zfft[2*index-1];
    }
    for(i=1; i <= nfft2; i++){ 
      zfft[fft_map_rev[i]] = fscr[i]; 
    }
    for(i=1; i <= nfft2; i++){ 
      fscr[i] = zfft[i]; 
    }
  }/*endif myid*/

/*------------------------------------------------------------------------------*/
/* IV) Broadcast the force and finish it */

  if(np_states > 1){  
    Scatterv(&fscr[1],&irecv[1],&idispl[1],MPI_DOUBLE,
             &zfft[1],nfft2_proc,MPI_DOUBLE,0,comm);
  }

  for(is=1; is <= nstate_up; is++){
    for(i=1 ; i<= nfft2_proc; i++){
      ind = (is-1)*nfft2_proc + i;
      dvrfc_up[ind] -= (2.0*dvrc_up[ind]*zfft[i]*dfact);
    }
  }
  if(cp_lsda==1){
    for(is=1; is <= nstate_dn; is++){
      for(i=1 ; i<= nfft2_proc; i++){
        ind = (is-1)*nfft2_proc + i;
        dvrfc_dn[ind] -= (2.0*dvrc_dn[ind]*zfft[i]*dfact);
      }
    }
  }

/*---------------------------------------------------*/
/* V) Free locally assigned memory */

  free(&irecv[1]);
  free(&idispl[1]);

/*====================================================================*/
 }/*end routine*/
/*====================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calc_ekin_dvr_par(double *dvrc,double *dvrfc,DVR_MATRIX *dvr_matrix,
                       CPCOEFFS_INFO *cpcoeffs_info,int nstate,
                       double *eke_ret,double *zfft,double *zfft_tmp,double vol,
                       COMMUNICATE *communicate,
                       PARA_FFT_PKG3D *cp_para_fft_pkg)

/*==========================================================================*/
 {/*begin routine*/
/*==========================================================================*/


#include "../typ_defs/typ_mask.h"

  /*   Local Pointer declarations   */

  double *Tsx     = dvr_matrix->Tsx;
  double *Tsy     = dvr_matrix->Tsy;
  double *Tsz     = dvr_matrix->Tsz;

  int  nkf1 = cp_para_fft_pkg->nkf1;
  int  nkf2 = cp_para_fft_pkg->nkf2;
  int  nkf3 = cp_para_fft_pkg->nkf3;

  int myid                  = communicate->myid_state;
  int nproc                 = communicate->np_states;
  MPI_Comm comm             = communicate->comm_states;

  int  nfft2                = cp_para_fft_pkg->nfft/2;
  int  nfft2_proc           = cp_para_fft_pkg->nfft_proc/2;

  int nstate_proc      = cpcoeffs_info->nstate_up_proc;
  int *ioff_up_st      = cpcoeffs_info->ioff_up;
  double dfact         = cpcoeffs_info->scale_fact;

/*   Local Variable declarations   */

  int ia,ib,ic,iap,ibp,icp,is,i;
  int ioff,index1,index2;
  double eke,sum;

  eke=0.0;

/*=========================================================================*/
/* I) Kinetic energy involves a multiplication with a tensor product */

  for(is=1; is <= nstate_proc; is++){

  /* X-contribution */

    ioff=0;
    for(ic=1; ic <= nkf3; ic++){
      for(ib=1; ib <= nkf2; ib++){
        for(ia=1; ia <= nkf1; ia++){
          ioff++;
          sum=0.0;
          for(iap=1; iap <= nkf1; iap++){
            index1 = (is-1)*nfft2 + (ic-1)*nkf2*nkf1 + (ib-1)*nkf1 + iap ;
            index2 = (ia-1)*nkf1 + iap;
            sum+=Tsx[index2]*dvrc[index1];
          }/* iap */
          zfft[ioff]=sum;
        }/* ia */
      }/* ib  */
    }/* ic */

  /*rearrange the vector*/

    ioff=0;
    for (ic=1; ic<= nkf3;ic++){
      for(ia=1; ia <= nkf1; ia++){
        for(ib=1; ib <= nkf2; ib++){
          ioff++;
          index1=(is-1)*nfft2+(ic-1)*nkf1*nkf2+(ib-1)*nkf1+ia;
          zfft_tmp[ioff]=dvrc[index1];
        }
      }
    }


  /* Y-contribution */

    ioff=0;
    for(ic=1; ic <= nkf3; ic++){
      for(ib=1; ib <= nkf2; ib++){
        for(ia=1; ia <= nkf1; ia++){
          ioff++;
          sum=0.0;
          for(ibp=1; ibp <= nkf2; ibp++){
            index1 = (ic-1)*nkf1*nkf2 + (ia-1)*nkf2 + ibp ;
            index2 = (ib-1)*nkf2 + ibp;
            sum+=Tsy[index2]*zfft_tmp[index1];
          }/* iap */
          zfft[ioff]+=sum;
        }/* ia */
      }/* ib  */
    }/* ic */

  /* rearrange the vector*/

    ioff=0;
    for (ib=1;ib<=nkf2;ib++){
      for (ia=1;ia<=nkf1;ia++){
        for (ic=1;ic<=nkf3;ic++){
          ioff++;
          index1=(is-1)*nfft2+(ic-1)*nkf1*nkf2+(ib-1)*nkf1+ia;
          zfft_tmp[ioff]=dvrc[index1];
        }
      }
    }

    /* Z-contribution */

    ioff=0;
    for (ic=1;ic<=nkf3;ic++){
      for (ib=1;ib<=nkf2;ib++){
        for (ia=1;ia<=nkf1;ia++){
          ioff++;
          sum=0.0;
          for (icp=1;icp<=nkf3;icp++){
            index1=(ib-1)*nkf1*nkf3+(ia-1)*nkf3+icp;
            index2=(ic-1)*nkf3+icp;
            sum+=Tsz[index2]*zfft_tmp[index1];
          }
          zfft[ioff]+=sum;
        }
      }
    }

    /*add to the force*/

    for(i=1;i<=nfft2;i++){
      dvrfc[(is-1)*nfft2+i]+=zfft[i]*dfact;
    }

    /*add to the kinetic energy*/

    sum=0.0;
    for(i=1;i<=nfft2;i++){
      sum+=dvrc[(is-1)*nfft2+i]*zfft[i];
    }
    eke+=sum;

  }/*end for is*/


  /* return the kinetic energy contribution*/

  *eke_ret = -0.5*eke*dfact;

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

void calc_hartree_dvr_clus (CPCOEFFS_INFO *cpcoeffs_info, CPSCR_WAVE *cpscr_wave,
                            CP_DVR_CLUS *cp_dvr_clus, double *rho,
                            double *fcre, double *fcre_tmp,
                            double *dvrfc,double *dvrc, double *eh_ret,
                            int nfft2_proc, int nstate, COMMUNICATE *communicate)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"

  int i,j,k,l,m,n,is,ioff;
  double pi = M_PI;
  double fpi = 4.0*pi;
  double eh,sum;

  int   myid_state       =    communicate->myid_state;
  int   np_states        =    communicate->np_states;
  MPI_Comm comm_states   =    communicate->comm_states;

  int grid_nx = cpcoeffs_info->grid_nx;
  int grid_ny = cpcoeffs_info->grid_ny;
  int grid_nz = cpcoeffs_info->grid_nz;
  int nfft2 = grid_nx*grid_ny*grid_nz;

  double dfact = cpcoeffs_info->scale_fact;

  int id,ip;
  double **del_u_x, *u, *t,*coef;
  double dx,box_x,pidx,delta,xi,xj;

  int *irecv,*idispl;
  int idiv,irem;

  irecv  = (int *) malloc((np_states)*sizeof(int))-1;
  idispl = (int *) malloc((np_states)*sizeof(int))-1;

/*===============================================================*/
/* I) In Parallel: AllGather Density to each processor. */

  if( np_states > 1){

    idiv    = (grid_ny*grid_nz)/np_states;
    irem    = (grid_ny*grid_nz)%np_states;
    for(i=0; i <irem; i++){ irecv[i+1] = (idiv+1)*grid_nx;  }
    for(i=irem ; i < np_states; i++){ irecv[i+1] = idiv*grid_nx;  }

    idispl[1] = 0;
    for(i=2; i <= np_states; i++){
      idispl[i] = idispl[i-1] + irecv[i-1];
    }

    Barrier(comm_states);

    Allgatherv(&rho[1],nfft2_proc,MPI_DOUBLE,&fcre_tmp[1],
               &irecv[1],&idispl[1],MPI_DOUBLE,0,comm_states);

  }else{
    for(i=1; i <= nfft2; i++){fcre_tmp[i] = rho[i]; }
  }

/*=========================================================================*/
/* II) Compute Hartree potential */

  for(i=1;i<=nfft2;i++){
    fcre_tmp[i] = fcre_tmp[i]*sqrt(dfact);
  }

  get_hartree_pot_clus(fcre,fcre_tmp,nfft2,cpcoeffs_info,cpscr_wave,cp_dvr_clus);

/*==========================================================================*/
/* III) Check the accuracy of Hartree potential */

#ifdef CHECK_HARTREE_ON

  if(np_states > 1){
    Barrier(comm_states);
    Reduce(&(fcre[1]),&(fcre_tmp[1]),nfft2,MPI_DOUBLE,MPI_SUM,0,comm_states);
  }else{
    for(i=1;i<=nfft2;i++){
      fcre_tmp[i] = fcre[i];
    }
  }

  if(myid_state ==0){

    del_u_x = cmall_mat(1,grid_nx,1,grid_nx);
    u = (double *)cmalloc(grid_nx*grid_ny*grid_nz*sizeof(double))-1;
    t = (double *)cmalloc(grid_nx*grid_ny*grid_nz*sizeof(double))-1;
    coef = (double *)cmalloc(4*sizeof(double));

    dx = pow(dfact, 1.0/3.0);

    coef[0] = -490.0/180.0;
    coef[1] =  270.0/180.0;
    coef[2] =  -27.0/180.0;
    coef[3] =    2.0/180.0;

    for(i=1;i<=nfft2;i++){
      u[i]=0.0;
    }

    for(l=4;l<=grid_nz-4;l++){
      for(m=4;m<=grid_ny-4;m++){
        for(n=4;n<=grid_nx-4;n++){
          id = (l-1)*grid_nx*grid_ny+(m-1)*grid_nx+n;
          ip=(l-1)*grid_nx*grid_ny+(m-1)*grid_nx;
          for(i=-3;i<=3;i++){
            k= abs(i);
            u[id] += fcre_tmp[ip+n+i]*coef[k];
          }
        }
      }
    }

    id=0;
    for(i=1;i<=grid_nx;i++){
      for(k=1;k<=grid_nz;k++){
        for(j=1;j<=grid_ny;j++){
          id=id+1;
          ip=(k-1)*grid_nx*grid_ny+(j-1)*grid_nx+i;
          t[id]=fcre_tmp[ip];
        }
      }
    }

    for(l=4;l<=grid_nz-4;l++){
      for(m=4;m<=grid_ny-4;m++){
        for(n=4;n<=grid_nx-4;n++){
          id = (l-1)*grid_nx*grid_ny+(m-1)*grid_nx+n;
          ip = (n-1)*grid_ny*grid_nz+(l-1)*grid_ny;
          for(j=-3;j<=3;j++){
            k=abs(j);
            u[id] += t[ip+m+j]*coef[k];
          }
        }
      }
    }

    id=0;
    for(j=1;j<=grid_ny;j++){
      for(i=1;i<=grid_nx;i++){
        for(k=1;k<=grid_nz;k++){
          id=id+1;
          ip = (i-1)*grid_ny*grid_nz+(k-1)*grid_ny+j;
          fcre_tmp[id]=t[ip];
        }
      }
    }

    for(l=4;l<=grid_nz-4;l++){
      for(m=4;m<=grid_ny-4;m++){
        for(n=4;n<=grid_nx-4;n++){
          id = (l-1)*grid_nx*grid_ny+(m-1)*grid_nx+n;
          ip = (n-1)*grid_nz+(m-1)*grid_nx*grid_nz;
          for(k=-3;k<=3;k++){
            i=abs(k);
            u[id] += fcre_tmp[ip+l+k]*coef[i];
          }
        }
      }
    }

    for(i=1;i<=nfft2;i++){
      t[i] = u[i]/(dx*dx);
    }
  }/*endif myid==0*/

  if(np_states > 1){
    Barrier(comm_states);

    Allgatherv(&rho[1],nfft2_proc,MPI_DOUBLE,&fcre_tmp[1],
               &irecv[1],&idispl[1],MPI_DOUBLE,0,comm_states);

  }else{
    for(i=1; i <= nfft2; i++){fcre_tmp[i] = rho[i]; }
  }

  if(myid_state==0){

    sum=0.0;
    for(k=4;k<=grid_nz-4;k++){
      for(j=4;j<=grid_ny-4;j++){
        for(i=4;i<=grid_nx-4;i++){
          ip=(k-1)*grid_nx*grid_ny+(j-1)*grid_nx+i;
          sum += (t[ip]+fpi*fcre_tmp[ip])*(t[ip]+fpi*fcre_tmp[ip]);
        }
      }
    }
    sum = sum/((double)((grid_nx-6)*(grid_ny-6)*(grid_nz-6)));
    sum = sqrt(sum);

    printf("RMS error for Hartree potential = %.15g\n",sum);

    cfree_mat(del_u_x,1,grid_nx,1,grid_nx);
    cfree(&(u[1]));
    cfree(&(t[1]));
    cfree(&(coef[1]));
  }/*end if myid==0*/

#endif

/*===========================================================================*/
/* IV) Reduce Scatter Hartree potential to all processors */

  if(np_states > 1){

    for(i=1;i<=nfft2;i++){
      fcre_tmp[i]=fcre[i];
    }

    Barrier(comm_states);

    Reduce_scatter(&fcre_tmp[1],&fcre[1],&irecv[1], MPI_DOUBLE,
                   MPI_SUM,comm_states);
  }


/*===========================================================================*/
/* V) compute Hartree energy and force */

  eh = 0.0;
  for(i=1; i <= nfft2_proc; i++){
    eh += fcre[i]*rho[i];
  }

  *eh_ret += 0.5*eh*dfact;

  /* printf("HARTREE E= %.10g %.10g\n",0.5*eh*dfact, *eh_ret); */

  sum=0.0;
  for(is = 1;is<= nstate; is++){
    ioff = (is-1)*nfft2_proc;
    for(i=1;i<= nfft2_proc;i++){
      dvrfc[ioff + i] -= 2.0*fcre[i]*dvrc[ioff+i]*dfact;
    }
  }

/*--------------------------------------------------------------*/
/* VI) Free locally assigned memory */

  free(&irecv[1]);
  free(&idispl[1]);

/*==============================================================*/
}/*end routine*/
/*==============================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_hartree_pot_clus(double *fcre, double *fcre_tmp, int ncoef,
                          CPCOEFFS_INFO *cpcoeffs_info,
                          CPSCR_WAVE *cpscr_wave, CP_DVR_CLUS *cp_dvr_clus)

/*========================================================================*/
{ /*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int i,j,k,id,is,ip,it,l,m,n,idx,idy,idz;
  double sum,pi,twosqrtpi;

  int grid_nx  = cpcoeffs_info->grid_nx;
  int grid_ny  = cpcoeffs_info->grid_ny;
  int grid_nz  = cpcoeffs_info->grid_nz;
  int ngrid = grid_nx*grid_ny*grid_nz;

  double *Fx   = cp_dvr_clus->Fx;
  double *Fy   = cp_dvr_clus->Fy;
  double *Fz   = cp_dvr_clus->Fz;

  double *vr1= cpscr_wave->zfft;
  double *vr2= &cpscr_wave->zfft[ncoef];
  double *vi1= cpscr_wave->zfft_tmp;
  double *vi2= &cpscr_wave->zfft_tmp[ncoef];

  double *tq_all = cp_dvr_clus->tq_all;
  double *wq_all = cp_dvr_clus->wq_all;
  double *fmat_x = cp_dvr_clus->fmat_x;
  double *fmat_y = cp_dvr_clus->fmat_y;
  double *fmat_z = cp_dvr_clus->fmat_z;
  int    my_num_t = cp_dvr_clus->my_num_t;

/*=======================================================================*/
/* I) Initialization */

  pi=M_PI;
  twosqrtpi = 2.00/sqrt(pi);

  for (i=1;i<=ngrid;i++){
     vi1[i] =0.0;
  }

/*-----------------------------------------------*/
/*        SUM OVER QUADRATURE POINTS IN t        */
/*-----------------------------------------------*/

  idx=0;
  idy=0;
  idz=0;
  for(it=1;it<=my_num_t;it++){

    for(i=1;i<=grid_nx;i++){
      for(j=1;j<=i;j++){
        idx += 1;
        Fx[(i-1)*grid_nx+j] = fmat_x[idx];
        if(j != i) Fx[(j-1)*grid_nx+i] = fmat_x[idx];
      }
    }
    for(i=1;i<=grid_ny;i++){
      for(j=1;j<=i;j++){
        idy += 1;
        Fy[(i-1)*grid_ny+j] = fmat_y[idy];
        if(j != i) Fy[(j-1)*grid_ny+i] = fmat_y[idy];
      }
    }
    for(i=1;i<=grid_nz;i++){
      for(j=1;j<=i;j++){
        idz += 1;
        Fz[(i-1)*grid_nz+j] = fmat_z[idz];
        if(j != i) Fz[(j-1)*grid_nz+i] = fmat_z[idz];
      }
    }

    /* II) transform along the x coordinate (inner most index) */

    id=0;
    for (i=1;i<=grid_nz;i++) {
      for (j=1;j<=grid_ny;j++) {
        ip=(i-1)*grid_nx*grid_ny+(j-1)*grid_nx;
        for (n=1;n<=grid_nx;n++) {
          id=id+1;
          is=(n-1)*grid_nx;
          sum=0.0;
          for (k=1;k<=grid_nx;k++) {
            sum=sum+Fx[is+k]*fcre_tmp[ip+k];
          }
          vr1[id]=sum;
        }
      }
    }

    /* III) transform along the y coordinate */

    id=0;
    for (i=1;i<=grid_nz;i++) {
      for (m=1;m<=grid_ny;m++) {
        is=(m-1)*grid_ny;
        for (n=1;n<=grid_nx;n++) {
          id=id+1;
          sum=0.0;
          ip=(i-1)*grid_nx*grid_ny+n;
          for (j=1;j<=grid_ny;j++) {
            sum=sum+Fy[is+j]*vr1[ip+(j-1)*grid_nx];
          }
          vr2[id]=sum;
        }
      }
    }

    /* IV) shuffling vector elements so that z becomes the inner most index */

    id=0;
    for (m=1;m<=grid_ny;m++) {
      for (n=1;n<=grid_nx;n++) {
        for (i=1;i<=grid_nz;i++) {
          id=id+1;
          ip=(i-1)*grid_nx*grid_ny+(m-1)*grid_nx+n;
          vr1[id]=vr2[ip];
        }
      }
    }

    /* V) transform along the z coordinates*/

    id=0;
    for (m=1;m<=grid_ny;m++) {
      for (n=1;n<=grid_nx;n++) {
        ip=(m-1)*grid_nx*grid_nz+(n-1)*grid_nz;
        for (l=1;l<=grid_nz;l++) {
          id=id+1;
          sum=0.0;
          is=(l-1)*grid_nz;
          for (i=1;i<=grid_nz;i++) {
            sum=sum+Fz[is+i]*vr1[ip+i];
          }
          vr2[id]=sum;
        }
      }
    }

    /* VI) accumulating contribution from each t quadrature point */

    for(i=1;i<=ngrid;i++){
      vi1[i] += vr2[i]*wq_all[it];
    }

  }/*end for it*/

/*=====================================================================*/
/* VII) reorder*/

  id=0;
  for (i=1;i<=grid_nz;i++) {
    for (m=1;m<=grid_ny;m++) {
      for (n=1;n<=grid_nx;n++) {
        id=id+1;
        ip=(m-1)*grid_nz*grid_nx+(n-1)*grid_nz+i;
        fcre[id]=vi1[ip]*twosqrtpi;
      }
    }
  }

/*==============================================================*/
}/*end routine*/
/*==============================================================*/


