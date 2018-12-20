/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: cp_energy_ee_rho_threads_state.c               */
/*                                                                          */
/* This routine calculate the density, coefficient forces with openMP	    */
/* parallel on states.							    */
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
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

//#define REAL_PP_DEBUG
//#define TEST_DM

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  cp_rho_calc:                                                            */
/*      calculates the density (up,down or total)                           */
/*      using the orbitals. Double packing of ffts and other                */
/*      goodies (most studdly coding)                                       */
/*      Even the option to do full g space (cp_rho_calc_full_g)             */
/*      or the studdly hybrid option in (cp_rho_calc_hybrid)                */
/*==========================================================================*/

void cp_rho_calc_hybrid_threads_state(CPEWALD *cpewald,CPSCR *cpscr,
                        CPCOEFFS_INFO *cpcoeffs_info,EWALD *ewald,
                        CELL *cell,double *ccreal, double *ccimag,
                        int icoef_form,int icoef_orth,
                        double *rhocr ,double *rhoci,double *rho,
                        double *rhocr_dens_cp_box,double *rhoci_dens_cp_box,
                        double *del_rho_x, double *del_rho_y,
                        double *del_rho_z,  
                        double *del2_rho,int nstate,int ncoef,
                        int cp_gga,int cp_dual_grid_opt,
                        int n_interp_pme_dual,
                        COMMUNICATE *communicate, 
                        PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg,
                        PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_lg,
                        PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box,
                        PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_dens_cp_box,
                        PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm,CP *cp,
			CLASS *class,GENERAL_DATA *general_data)
      
/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"

/* local variables                                                  */

  int iii,ioff,ioff2;
  int is,i,j,k,iupper;
  double vol_cp,rvol_cp;
  double temp_r,temp_i;
  int realSparseOpt = cpewald->realSparseOpt;

/*  Assign local pointers                                           */
  double *zfft           =    cpscr->cpscr_wave.zfft;
  double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
  double **zfft_threads  =    cpscr->cpscr_wave.zfft_threads;
  double **zfft_tmp_threads = cpscr->cpscr_wave.zfft_tmp_threads;
  double *rho_scr;
  double *rho_scr_threads = NULL;
  int numThreads = cp_para_fft_pkg3d_lg->numThreads;
  int iThread;

  double *hmati_cp       =    cell->hmati_cp;
  double *hmat_cp        =    cell->hmat_cp;

  double dbox_rat        =    cpewald->dbox_rat;
  double *bw_r           =    cpscr->cpscr_dual_pme.bw_r;
  double *bw_i           =    cpscr->cpscr_dual_pme.bw_i;

  int cp_elf_calc_frq    =    cpcoeffs_info->cp_elf_calc_frq; 

  int   myid_state       =    communicate->myid_state;
  int   np_states        =    communicate->np_states;
  int   laplacian_on     =    cpcoeffs_info->cp_laplacian_on;

  //int *recv_counts_coef;
  //int ncoef_l,ncoef_l_use,ncoef_l_proc,nfft_proc,nfft,nfft2,nfft2_proc;
  
  int  *recv_counts_coef =    cp_para_fft_pkg3d_lg->recv_counts_coef;
  int   ncoef_l          =    cp_para_fft_pkg3d_lg->ncoef;
  int   ncoef_l_use      =    cp_para_fft_pkg3d_lg->ncoef_use;
  int   ncoef_l_proc     =    cp_para_fft_pkg3d_lg->ncoef_proc;
  int   nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
  int   nfft             =    cp_para_fft_pkg3d_lg->nfft;
  int   nfft2            =    nfft/2;
  int   nfft2_proc       =    nfft_proc/2;

  int   nfft_dens_cp_box,nfft2_dens_cp_box;
  int   nfft_proc_dens_cp_box,nfft2_proc_dens_cp_box;
  int   ncoef_l_dens_cp_box;
#ifdef DEBUG_LSDA
  double *rhocr_dens_cp_box;
  double *rhoci_dens_cp_box;
#endif

  int *mapFFTWSm = cp_sclr_fft_pkg3d_sm->mapFFTW;
  int *mapConFFTWSm = cp_sclr_fft_pkg3d_sm->mapConFFTW;
  fftw_plan fftwPlan3DForward,fftwPlan3DBackward;
  fftw_complex *fftw3DForwardIn,*fftw3DForwardOut;
  fftw_complex *fftw3DBackwardIn,*fftw3DBackwardOut;

  double integral,int_tmp;
  int    *recv_counts_coef_dens_cp_box;
  int fftw3dFlag = cpcoeffs_info->fftw3dFlag;
  int ind,indt;

 //debug
  int icoef;
  int fftIndTrans;
  int fftInd;
  //int nkf1,nkf2,nkf3;
  int nkf1    = cp_para_fft_pkg3d_lg->nkf1;
  int nkf2    = cp_para_fft_pkg3d_lg->nkf2;
  int nkf3    = cp_para_fft_pkg3d_lg->nkf3;
  int kc,kb,ka;
  int *kastr_sm = cpewald->kastr_sm;
  int *kbstr_sm = cpewald->kbstr_sm;
  int *kcstr_sm = cpewald->kcstr_sm;
  int *kastr = ewald->kastr;
  int *kbstr = ewald->kbstr;
  int *kcstr = ewald->kcstr;
  double sum;
  double time_st,time_end;
#ifdef TEST_DM
  int iigrid,jjgrid,kkgrid,llgrid;
  //double *dm = (double*)cmalloc(numThreads*nkf1*nkf1*sizeof(double));
  //double *dm_reduce = (double*)cmalloc(nkf1*nkf1*sizeof(double));
  double **dm_complete = (double**)cmalloc(numThreads*sizeof(double*));
  for(i=0;i<numThreads;i++){
    dm_complete[i] = (double*)calloc(nkf1*nkf1*nkf2*nkf3);
  }
  double *dm_reduce = (double*)calloc(nkf1*nkf1*nkf2*nkf3);
#endif
  
  /*
  if(realSparseOpt==0){
    nkf1    = cp_para_fft_pkg3d_lg->nkf1;
    nkf2    = cp_para_fft_pkg3d_lg->nkf2;
    nkf3    = cp_para_fft_pkg3d_lg->nkf3;
    recv_counts_coef =    cp_para_fft_pkg3d_lg->recv_counts_coef;
    ncoef_l          =    cp_para_fft_pkg3d_lg->ncoef;
    ncoef_l_use      =    cp_para_fft_pkg3d_lg->ncoef_use;
    ncoef_l_proc     =    cp_para_fft_pkg3d_lg->ncoef_proc;
    nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
    nfft             =    cp_para_fft_pkg3d_lg->nfft;
    nfft2            =    nfft/2;
    nfft2_proc       =    nfft_proc/2;
  }
  else{
    nkf1    = cp_para_fft_pkg3d_sparse->nkf1;
    nkf2    = cp_para_fft_pkg3d_sparse->nkf2;
    nkf3    = cp_para_fft_pkg3d_sparse->nkf3;
    recv_counts_coef =    cp_para_fft_pkg3d_sparse->recv_counts_coef;
    ncoef_l          =    cp_para_fft_pkg3d_sparse->ncoef;
    ncoef_l_use      =    cp_para_fft_pkg3d_sparse->ncoef_use;
    ncoef_l_proc     =    cp_para_fft_pkg3d_sparse->ncoef_proc;
    nfft_proc        =    cp_para_fft_pkg3d_sparse->nfft_proc;
    nfft             =    cp_para_fft_pkg3d_sparse->nfft;
    nfft2            =    nfft/2;
    nfft2_proc       =    nfft_proc/2;
  }
  */


  MPI_Comm comm_states   =    communicate->comm_states;

  if(cp_dual_grid_opt >= 1){
    nfft_dens_cp_box       =  cp_para_fft_pkg3d_dens_cp_box->nfft;
    nfft_proc_dens_cp_box  =  cp_para_fft_pkg3d_dens_cp_box->nfft_proc;
    nfft2_dens_cp_box      =  nfft_dens_cp_box/2;
    nfft2_proc_dens_cp_box =  nfft_proc_dens_cp_box/2;


    ncoef_l_dens_cp_box  =  cp_para_fft_pkg3d_dens_cp_box->ncoef;

#ifdef DEBUG_LSDA
    rhocr_dens_cp_box = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
    rhoci_dens_cp_box = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
#endif

    recv_counts_coef_dens_cp_box = 
              cp_para_fft_pkg3d_dens_cp_box->recv_counts_coef;

  }/*endif cp_dual_grid_opt*/

  if(np_states > 1){
    rho_scr = cpscr->cpscr_rho.v_ks_up; 
  } else {
    rho_scr = rho;
  }

  //rho_scr_threads = (double*)calloc((numThreads*nfft2+1),sizeof(double));
  rho_scr_threads = (double*)cmalloc((numThreads*nfft2+1)*sizeof(double));
  for(i=1;i<=numThreads*nfft2;i++){
    rho_scr_threads[i] = 0.0;
  }
  
/* ================================================================= */
/*0) Check the form of the coefficients                              */

  if(icoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in orthogonal form    \n");
    printf("on state processor %d in cp_rho_calc_hybrid   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1 && icoef_form == 1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in normal (not transposed) \n");
    printf("form on state processor %d in cp_rho_calc_hybrid  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/


/* ================================================================= */
/*1) zero density and gradients if necessary                         */

  if(cp_dual_grid_opt >= 1){
    for(i=1; i<= nfft2_dens_cp_box ; i++){
      rho_scr[i] = 0.0;
    }/*endfor*/
  }else{
    for(i=1; i<= nfft2 ; i++){
      rho_scr[i] = 0.0;
    }/*endfor*/
  }/*endif cp_dual_grid_opt*/

/*==========================================================================*/
/*==========================================================================*/
/*2)sum the density in real space two states at a time                      */
/*  This is done at state level and uses the scalar packages!               */

  iupper = nstate;
  if((nstate % 2) != 0){ 
    iupper = nstate-1; 
  }/* endif */

  time_st = omp_get_wtime();
  omp_set_num_threads(numThreads);
#ifdef TEST_DM
  #pragma omp parallel private(iThread,is,ioff,ioff2,iigrid,jjgrid)
#else
  #pragma omp parallel private(iThread,is,ioff,ioff2)
#endif
  {
    //iThread = 0;
    iThread = omp_get_thread_num();
    #pragma omp for
    for(is = 1; is <= iupper; is = is + 2){
      ioff = (is-1)*ncoef;
      ioff2 = (is)*ncoef;

/*--------------------------------------------------------------------------*/
/*I) pack the complex zfft array with two wavefunctions (real) 
    the wavefunctions are reperesented in spherically cuttof 
    half g space                                                            */

      //if(realSparseOpt==0){
      if(fftw3dFlag==0){
	dble_pack_coef(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
		       zfft_threads[iThread],cp_sclr_fft_pkg3d_sm);
      }
      else{
	dble_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
	    zfft_threads[iThread],cp_sclr_fft_pkg3d_sm);
      }
      //}
      /*
      else{
        if(fftw3dFlag==0){
          dble_pack_coef(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
                         zfft_threads[iThread],cp_sclr_fft_pkg3d_sparse);
        }
        else{
          dble_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
              zfft_threads[iThread],cp_sclr_fft_pkg3d_sparse);
        }
      }
      */
/*--------------------------------------------------------------------------*/
/*II) fourier transform the two wavefunctions to real space
     convention exp(-igr)                                                   */
      
      //if(iThread==0)time_st = omp_get_wtime();
      //if(realSparseOpt==0){
      if(fftw3dFlag==0){
	para_fft_gen3d_fwd_to_r(zfft_threads[iThread],zfft_tmp_threads[iThread],cp_sclr_fft_pkg3d_sm);
      }
      else{
	para_fft_gen3d_fwd_to_r_fftw3d(zfft_threads[iThread],cp_sclr_fft_pkg3d_sm,iThread);
      }
      //}
      /*
      else{
        if(fftw3dFlag==0){
          para_fft_gen3d_fwd_to_r(zfft_threads[iThread],zfft_tmp_threads[iThread],cp_sclr_fft_pkg3d_sparse);
        }
        else{
          para_fft_gen3d_fwd_to_r_fftw3d(zfft_threads[iThread],cp_sclr_fft_pkg3d_sparse,iThread);

        }
      }
      */
      
      //printf("zfft_threads[0][2] %lg\n",zfft_threads[0][3]);
      /*
      if(iThread==0){
        time_end = omp_get_wtime();
        cp_sclr_fft_pkg3d_sm->cputime1 += time_end-time_st;
      }
      */
      /*
      for(kc=0;kc<nkf3;kc++){
	for(kb=0;kb<nkf2;kb++){
	  for(ka=0;ka<nkf1;ka++){
	    ind = kc*nkf2*nkf1+kb*nkf1+ka;
	    indt = ka*nkf2*nkf3+kb*nkf3+kc;
	    printf("zfft_threads %i %i %i %lg\n",kc,kb,ka,zfft_threads[iThread][indt*2+1]);
	  }
	}
      }

      fflush(stdout);
      exit(0);
      */

/*--------------------------------------------------------------------------*/
/* III) add the square of the two wave functions to the density(real space) */

      if(cp_dual_grid_opt >= 1){
	sum_rho(zfft_threads[iThread],&rho_scr_threads[iThread*nfft2],cp_sclr_fft_pkg3d_dens_cp_box); 
      }else{
	//if(realSparseOpt==0){
	sum_rho(zfft_threads[iThread],&rho_scr_threads[iThread*nfft2],cp_sclr_fft_pkg3d_lg);
	//}
	/*
	else{
	  sum_rho(zfft_threads[iThread],&rho_scr_threads[iThread*nfft2],cp_sclr_fft_pkg3d_sparse);
        }
	*/
#ifdef TEST_DM
        /*
        if(fftw3dFlag==1){
          for(iigrid=0;iigrid<nkf1;iigrid++){
            for(jjgrid=0;jjgrid<nkf1;jjgrid++){
              dm[iThread*nkf1*nkf1+iigrid*nkf1+jjgrid] += zfft_threads[iThread][2*iigrid*nkf2*nkf3+1]*zfft_threads[iThread][2*jjgrid*nkf2*nkf3+1];
              dm[iThread*nkf1*nkf1+iigrid*nkf1+jjgrid] += zfft_threads[iThread][2*iigrid*nkf2*nkf3+2]*zfft_threads[iThread][2*jjgrid*nkf2*nkf3+2];
            }
          }
        }
        else{
          for(iigrid=0;iigrid<nkf1;iigrid++){
            for(jjgrid=0;jjgrid<nkf1;jjgrid++){
              dm[iThread*nkf1*nkf1+iigrid*nkf1+jjgrid] += zfft_threads[iThread][2*iigrid+1]*zfft_threads[iThread][2*jjgrid+1];
              dm[iThread*nkf1*nkf1+iigrid*nkf1+jjgrid] += zfft_threads[iThread][2*iigrid+2]*zfft_threads[iThread][2*jjgrid+2];
            }
          }
        }
        */
        if(fftw3dFlag==1){
          for(iigrid=0;iigrid<nkf1;iigrid++){
            for(jjgrid=0;jjgrid<nfft2;jjgrid++){
              dm_complete[iThread][iigrid*nkf1*nkf2*nkf3+jjgrid] += zfft_threads[iThread][2*iigrid*nkf2*nkf3+1]*zfft_threads[iThread][2*jjgrid+1];
              dm_complete[iThread][iigrid*nkf1*nkf2*nkf3+jjgrid] += zfft_threads[iThread][2*iigrid*nkf2*nkf3+2]*zfft_threads[iThread][2*jjgrid+2];
            }
          }
        }
#endif
      }//endif cp_dual_grid_opt
    }//endfor is
  }//end omp

  
  for(iThread=0;iThread<numThreads;iThread++){
    for(i=1;i<=nfft2;i++){
      rho_scr[i] += rho_scr_threads[iThread*nfft2+i];
    }
#ifdef TEST_DM
    for(i=0;i<nkf1;i++){
      for(j=0;j<nkf1;j++){
        dm_reduce[i*nkf1+j] += dm[iThread*nkf1*nkf1+i*nkf1+j];
      }
    }
#endif
  }

#ifdef TEST_DM
  vol_cp  = getdeth(hmat_cp);
  rvol_cp = 1.0/vol_cp;
  for(i=0;i<nkf1;i++){
    for(j=0;j<nkf1;j++){
      printf("11111111 dm %.8lg\n",dm_reduce[i*nkf1+j]*rvol_cp);
    }
  }
  fflush(stdout); 
  exit(0);  
#endif

/*--------------------------------------------------------------------------*/
/*IV) if there is an odd number of states, use single pack
      for the last state (i caen't be chaengin' the laws of physics
      captn, i've got to be usin the single pack!!!!!!)                     */

  if((nstate % 2 ) != 0){
    ioff = (nstate-1)*ncoef;
    //if(realSparseOpt==0){
    if(fftw3dFlag==0){
      //printf("ioff %i creal %lg cimag %lg\n",ioff+1,ccreal[ioff+1],ccimag[ioff+1]);
      sngl_pack_coef(&ccreal[ioff],&ccimag[ioff],zfft_threads[0],cp_sclr_fft_pkg3d_sm);
    }
    else{
      sngl_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],zfft_threads[0],cp_sclr_fft_pkg3d_sm);
    }
    //}
    /*
    else{
      if(fftw3dFlag==0){
        //printf("ioff %i creal %lg cimag %lg\n",ioff+1,ccreal[ioff+1],ccimag[ioff+1]);
        sngl_pack_coef(&ccreal[ioff],&ccimag[ioff],zfft_threads[0],cp_sclr_fft_pkg3d_sparse);
      }
      else{
        sngl_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],zfft_threads[0],cp_sclr_fft_pkg3d_sparse);
      }
    }
    */

/*--------------------------------------------------------------------------*/
/* V) fourier transform the last wavefunction to real space
       convention exp(-igr)                                                 */

    //if(realSparseOpt==0){
    if(fftw3dFlag==0){
      para_fft_gen3d_fwd_to_r(zfft_threads[0],zfft_tmp_threads[0],cp_sclr_fft_pkg3d_sm);
    }
    else{
      para_fft_gen3d_fwd_to_r_fftw3d(zfft_threads[0],cp_sclr_fft_pkg3d_sm,0);
    }
    //}
    /*
    else{
      if(fftw3dFlag==0){
        para_fft_gen3d_fwd_to_r(zfft_threads[0],zfft_tmp_threads[0],cp_sclr_fft_pkg3d_sparse);
      }
      else{
        para_fft_gen3d_fwd_to_r_fftw3d(zfft_threads[0],cp_sclr_fft_pkg3d_sparse,0);
      }
    }
    */

/*--------------------------------------------------------------------------*/
/*VI) add the square of the last wave function to the density(real space)   */

    if(cp_dual_grid_opt >= 1){ 
      sum_rho(zfft_threads[0],rho_scr,cp_sclr_fft_pkg3d_dens_cp_box); 
    }else{
      //if(realSparseOpt==0){
      sum_rho(zfft_threads[0],rho_scr,cp_sclr_fft_pkg3d_lg);
      //}
      /*
      else{
	sum_rho(zfft_threads[0],rho_scr,cp_sclr_fft_pkg3d_sparse);
      }
      */
    }//endif cp_dual_grid_opt
  }//endif

  //if(fftw3dFlag==0)printf("rrrrho_scr %lg\n",rho_scr[2]);
  //else printf("rrrrho_scr %lg %lg\n",rho_scr[2],rho_scr[5185]);

  time_end = omp_get_wtime();
  //cp_sclr_fft_pkg3d_sm->cputime1 += time_end-time_st;

/*=========================================================================*/
/*=========================================================================*/
/*  3) get density in g space                                              */
/*  I) pack it up                                                          */
/*  rho_scr = rho in scalar                                                */

  time_st = omp_get_wtime();

  if(cp_dual_grid_opt >= 1){ 
    sngl_pack_rho(zfft,rho_scr,cp_sclr_fft_pkg3d_dens_cp_box); 
  }else{ // now the leading dimension of zfft is x
    //if(realSparseOpt==0){
    sngl_pack_rho(zfft,rho_scr,cp_sclr_fft_pkg3d_lg);
    //}
    /*
    else{
      sngl_pack_rho(zfft,rho_scr,cp_sclr_fft_pkg3d_sparse);
    }
    */
  //if(fftw3dFlag==0)sngl_pack_rho(zfft,rho_scr,cp_sclr_fft_pkg3d_lg); 
  //else sngl_pack_rho_fftw3d(zfft,rho_scr,cp_sclr_fft_pkg3d_lg);
  }/*endif cp_dual_grid_opt*/
  
  /*
  for(i=0;i<nfft2;i++){
    printf("000000000 zfft %lg %lg\n",zfft[2*i+1],zfft[2*i+2]);
  }
  */
  

/*--------------------------------------------------------------------------*/
/*  II) back transform to g-space  convention exp(igr)                      */
 
 /*
 int fragFlag = 0;
 if(fftw3dFlag==1)fragFlag = 1;
 fftw3dFlag = 0;
 */
 

  if(cp_dual_grid_opt >= 1){ 
    para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_dens_cp_box); 
  }else{
    //if(realSparseOpt==0){
    if(fftw3dFlag==0){
      para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_lg);
    }
    else{ //It's fine that zfft is x leading 
      para_fft_gen3d_bck_to_g_fftw3d(zfft,cp_sclr_fft_pkg3d_lg,0); 
    }
    //}
    /*
    else{
      if(fftw3dFlag==0){
        para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sparse);
      }
      else{ //It's fine that zfft is x leading 
        para_fft_gen3d_bck_to_g_fftw3d(zfft,cp_sclr_fft_pkg3d_sparse,0);
      }      
    }
    */
  }/*endif cp_dual_grid_opt*/

/*==========================================================================*/
/*  III) unpack the density                                                 */ 

  /*  
  if(fragFlag==0){
    for(i=1;i<=nfft2;i++){
      printf("1111 zfft %i %lg %lg\n",2*i+1,zfft[2*i+1],zfft[2*i+2]);
    }
    //fflush(stdout);
    //exit(0);
  }
  */

  if(cp_dual_grid_opt == 0){
    if(np_states == 1){
      //if(realSparseOpt==0){
      if(fftw3dFlag==0)sngl_upack_coef(rhocr,rhoci,zfft,cp_sclr_fft_pkg3d_lg);
      else sngl_upack_coef_fftw3d(rhocr,rhoci,zfft,cp_sclr_fft_pkg3d_lg);
      //}
      /*
      else{
        if(fftw3dFlag==0)sngl_upack_coef(rhocr,rhoci,zfft,cp_sclr_fft_pkg3d_sparse);
        else sngl_upack_coef_fftw3d(rhocr,rhoci,zfft,cp_sclr_fft_pkg3d_sparse);
      }//endif realSparseOpt
      */
    }else{
      //if(realSparseOpt==0){
      if(fftw3dFlag==0){
	sngl_upack_coef(zfft_tmp,&zfft_tmp[ncoef_l],zfft,cp_sclr_fft_pkg3d_lg);
      }
      else{
	sngl_upack_coef_fftw3d(zfft_tmp,&zfft_tmp[ncoef_l],zfft,cp_sclr_fft_pkg3d_lg);
      }
      //}
      /*
      else{
        if(fftw3dFlag==0){
          sngl_upack_coef(zfft_tmp,&zfft_tmp[ncoef_l],zfft,cp_sclr_fft_pkg3d_sparse);
        }
        else{
          sngl_upack_coef_fftw3d(zfft_tmp,&zfft_tmp[ncoef_l],zfft,cp_sclr_fft_pkg3d_sparse);
        }
      }
      */
    }//endif
  }else{
    if(np_states == 1){
      sngl_upack_coef(rhocr_dens_cp_box,rhoci_dens_cp_box,zfft,
                      cp_sclr_fft_pkg3d_dens_cp_box);
    }else{
      sngl_upack_coef(zfft_tmp,&zfft_tmp[ncoef_l_dens_cp_box],zfft,
                      cp_sclr_fft_pkg3d_dens_cp_box);
    }//endif
  }//endif cp_dual_grid_opt

  //printf("rhokkkkkkkk %lg %lg\n",rhocr[1],rhoci[1]);

/*==========================================================================*/
/* VII) Reduce rho in g space and get all of it in real space               */
/*      Now in g-level parallel and need parallel packages                  */

 if(np_states > 1){
   if(cp_dual_grid_opt == 0){
      Reduce_scatter(&zfft_tmp[1],&rhocr[1],&recv_counts_coef[1],
                     MPI_DOUBLE,MPI_SUM,comm_states);
      Reduce_scatter(&zfft_tmp[(ncoef_l+1)],&rhoci[1],&recv_counts_coef[1],
                     MPI_DOUBLE,MPI_SUM,comm_states);  

      //if(realSparseOpt==0){
      sngl_pack_coef(rhocr,rhoci,zfft,cp_para_fft_pkg3d_lg);
      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_lg);
      sngl_upack_rho(zfft,rho,cp_para_fft_pkg3d_lg); 
      //}
      /*
      else{
        sngl_pack_coef(rhocr,rhoci,zfft,cp_para_fft_pkg3d_sparse);
        para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_sparse);
        sngl_upack_rho(zfft,rho,cp_para_fft_pkg3d_sparse);
      }
      */
      // now rho is z leading
      Barrier(comm_states);
   }else{ //we dont need this
     Reduce_scatter(&zfft_tmp[1],&rhocr_dens_cp_box[1],
                    &recv_counts_coef_dens_cp_box[1],
                    MPI_DOUBLE,MPI_SUM,comm_states);
     Reduce_scatter(&zfft_tmp[(ncoef_l_dens_cp_box+1)],
                     &rhoci_dens_cp_box[1],&recv_counts_coef_dens_cp_box[1],
                     MPI_DOUBLE,MPI_SUM,comm_states);  

      sngl_pack_coef(rhocr_dens_cp_box,rhoci_dens_cp_box,zfft,
                     cp_para_fft_pkg3d_dens_cp_box);

      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_dens_cp_box);

      sngl_upack_rho(zfft,rho,cp_para_fft_pkg3d_dens_cp_box); 
      // now rho is z leading
      
   }/*endif cp_dual_grid_opt*/
 }
 else{ //transpose for sequential case
    if(fftw3dFlag==1){
      for(kc=0;kc<nkf3;kc++){
	for(kb=0;kb<nkf2;kb++){
	  for(ka=0;ka<nkf1;ka++){
	    ind = kc*nkf2*nkf1+kb*nkf1+ka+1;
	    indt = ka*nkf2*nkf3+kb*nkf3+kc+1;
	    zfft[ind] = rho_scr[indt];
	  }
	}
      }
      memcpy(&rho_scr[1],&zfft[1],nfft2*sizeof(double)); // now rho is z leading
    }
 }//endif np_states
 //printf("rrrrrho %lg %lg\n",rho[2],rho[3]);

/*===========================================================================*/
/* IF DUALED put rho real space onto the large grid and fft it to g space    */

 if(cp_dual_grid_opt >= 1){
/* sending density*vol_cp on small grid  */
     control_spread_rho(cpscr,rho,cell,dbox_rat,np_states,
                        n_interp_pme_dual,
                        cp_para_fft_pkg3d_dens_cp_box,
                        cp_para_fft_pkg3d_lg,cp_dual_grid_opt);  

  if(np_states == 1){
    //if(realSparseOpt==0){
    para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_lg); 
    sngl_upack_coef(rhocr,rhoci,zfft,cp_sclr_fft_pkg3d_lg);
    //}
    /*
    else{
      para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sparse);
      sngl_upack_coef(rhocr,rhoci,zfft,cp_sclr_fft_pkg3d_sparse);
    }
    */
  }else{
    //if(realSparseOpt==0){
    para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_lg); 
    sngl_upack_coef(rhocr,rhoci,zfft,cp_para_fft_pkg3d_lg);
    //}
  }/*endif np_states*/
 }/*endif cp_dual_grid_opt*/


/*--------------------------------------------------------------------*/
/*  Post-processing for pme grid need to multiply by complex weight factor */
 
  if((n_interp_pme_dual > 1) && (cp_dual_grid_opt == 2)){
    for(i=1; i<= ncoef_l_use; i++){
     temp_r   =  (rhocr[i]*bw_r[i] - rhoci[i]*bw_i[i]); 
     temp_i   =  (rhocr[i]*bw_i[i] + rhoci[i]*bw_r[i]); 
     rhocr[i] =  temp_r;
     rhoci[i] =  temp_i;
    }/*endfor*/

    if((myid_state+1) == np_states){rhocr[ncoef_l_proc]*=bw_r[ncoef_l_proc];}

   }/*endif pme grid */

/*===========================================================================*/
/* IV) finish the density in real space by dividing by the volume           */
/* DUALED SYSTEMS only keep the real space density on the cp_grid            */

     vol_cp  = getdeth(hmat_cp);
     rvol_cp = 1.0/vol_cp;
     //printf("rvol_cp %.16lg\n",rvol_cp);

     if(cp_dual_grid_opt >= 1){
      for(i=1 ; i<= nfft2_proc_dens_cp_box;i++){
         rho[i] *= rvol_cp;
      }/*endfor*/
     }else{
      for(i=1 ; i<= nfft2_proc;i++){
         rho[i] *= rvol_cp;
      }/*endfor*/
     }/*endif cp_dual_grid_opt*/


  if(fftw3dFlag==100){
    //double sum = 0.0;
    sum = 0.0;
    FILE *fp_rho = fopen("rho_bm","w");
    if(np_states == 1){
      for(kc=1;kc<=nkf3;kc++){
        for(kb=1;kb<=nkf2;kb++){
          for(ka=1;ka<=nkf1;ka++){
            i = (ka-1) + (kb-1)*nkf1 + (kc-1)*nkf1*nkf2 + 1;
            fprintf(fp_rho,"%i %i %i %.16g\n",kc,kb,ka,rho[i]);
            sum += rho[i];
          }//endfor
        }//endfor
      }//endfor
    }
    //printf("rho sum test %lg\n",sum);
    fclose(fp_rho);
    //fflush(stdout);
    //exit(0);
  }

/*==============================================================*/
/* VII) if doing gradient corrections, get gradient of density  */

  if((cp_gga == 1 || cp_elf_calc_frq > 0)) {
   if(cp_dual_grid_opt >= 1){
    control_grad_rho(cpewald,cpscr,ewald,rhocr_dens_cp_box,rhoci_dens_cp_box,
                     del_rho_x,del_rho_y,del_rho_z,del2_rho,
                     hmati_cp,vol_cp,laplacian_on,cp_dual_grid_opt,
                     cp_para_fft_pkg3d_dens_cp_box);
   }else{
    //if(realSparseOpt==0){
    control_grad_rho(cpewald,cpscr,ewald,rhocr,rhoci,
	             del_rho_x,del_rho_y,del_rho_z,del2_rho,
	             hmati_cp,vol_cp,laplacian_on,cp_dual_grid_opt,
	             cp_para_fft_pkg3d_lg);
    //}
    /*
    else{
      control_grad_rho(cpewald,cpscr,ewald,rhocr,rhoci,
                       del_rho_x,del_rho_y,del_rho_z,del2_rho,
                       hmati_cp,vol_cp,laplacian_on,cp_dual_grid_opt,
                       cp_para_fft_pkg3d_sparse);
    }
    */
    // now rho and gradiant are z leading
   }/*endif cp_dual_grid_opt*/
  }/*endif*/

  cfree(rho_scr_threads);

  time_end = omp_get_wtime();

  //cp_sclr_fft_pkg3d_sm->cputime2 += time_end-time_st;


/*==============================================================*/
}/*end routine*/
/*==============================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* calculate the force on the coeff's (up or down, both) */
/* given an appropriate v_ks (up or down, both).         */
/*==========================================================================*/

void coef_force_calc_hybrid_threads_state(CPEWALD *cpewald,int nstate,
                             double *ccreal,double *ccimag, 
                             double *fccreal,double  *fccimag,
                             double *cre_scr,double *cim_scr,
                             double *cp_hess_re,double *cp_hess_im,
                             double **zfft_threads,double **zfft_tmp_threads,
                             double *v_ks,double *v_ks_tau,double *ak2_sm,
                             double *eke_ret,double *pvten_cp,
                             int cp_ptens_calc,double *hmati,
                             COMMUNICATE *communicate,
                             int icoef_form,int icoef_orth,int ifcoef_form,
                             int cp_tau_functional,int cp_min_on,
                             PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm,
			     CP *cp, CLASS *class,GENERAL_DATA *general_data)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                               */
#include "../typ_defs/typ_mask.h"

  STAT_AVG *stat_avg = &(general_data->stat_avg);
  CLATOMS_POS *clatoms_pos = &(class->clatoms_pos[1]);
  PSEUDO *pseudo = &(cp->pseudo);
  PSEUDO_REAL *pseudoReal = &(pseudo->pseudoReal);
  //PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sparse = &(cp->cp_sclr_fft_pkg3d_sparse);

  int is,i,j,k,iupper;
  int igrid,jgrid;
  double tpi;
  double aka,akb,akc,xk,yk,zk,cfact;
  double eke;
  int ioff,ncoef1,ioff2;
  int iii,iis,nis;
  //int realSparseOpt = cpewald->realSparseOpt;

  int nfft       = cp_sclr_fft_pkg3d_sm->nfft;
  int nfft2   = nfft/2;
  int ncoef      = cp_sclr_fft_pkg3d_sm->ncoef;
  int nkf1	 = cp_sclr_fft_pkg3d_sm->nkf1;
  int nkf2       = cp_sclr_fft_pkg3d_sm->nkf2;
  int nkf3       = cp_sclr_fft_pkg3d_sm->nkf3;
  int myid_state = communicate->myid_state;
  int np_states  = communicate->np_states;
  int pseudoRealFlag = cp->pseudo.pseudoReal.pseudoRealFlag;
  int nlppForceOnly = cp->pseudo.pseudoReal.nlppForceOnly;
  int forceCalcFlag = pseudoReal->forceCalcFlag;
  int numThreads = cp_sclr_fft_pkg3d_sm->numThreads;
  int iThread;


/*            Local pointers                                       */

  int  *kastore_sm    =  cpewald->kastr_sm;
  int  *kbstore_sm    =  cpewald->kbstr_sm;
  int  *kcstore_sm    =  cpewald->kcstr_sm; 

#define DEBUG_OFF
#ifdef DEBUG
  int icount;
  double c_g,g2,anorm,sum,vol,cre_now,cim_now;
  double dx,x_pos,y_pos,z_pos,phase_r,phase_i,arg;
  FILE *fp;
#endif

  double sum_check,sum_check_tmp;
  MPI_Comm comm_states = communicate->comm_states;

  int fftw3dFlag = cpewald->fftw3dFlag;
  int numAtomTot = class->clatoms_info.natm_tot;
  double *energyNl = (double *)cmalloc(numThreads*sizeof(double));
  double *fxThreads = (double *)cmalloc(numThreads*numAtomTot*sizeof(double));
  double *fyThreads = (double *)cmalloc(numThreads*numAtomTot*sizeof(double));
  double *fzThreads = (double *)cmalloc(numThreads*numAtomTot*sizeof(double));
  double *fx = clatoms_pos->fx;
  double *fy = clatoms_pos->fy;
  double *fz = clatoms_pos->fz;
  double *ak2Kinetic = cpewald->ak2Kinetic;

  double time_st,time_end;

/* ================================================================= */
/*0) Check the form of the coefficients                              */

  if(icoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in orthogonal form    \n");
    printf("on state processor %d in coef_force_calc_hybrid\n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1) 
   if( ifcoef_form==1 || icoef_form == 1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs/coef forces must be in normal (not transposed) \n");
    printf("form on state processor %d in coef_force_calc_hybrid  \n",
           myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/*=================================================================*/
/*  Find the upper state limit                                     */

  //reset variables
  /*
  if(realSparseOpt==1){
    nfft = cp_sclr_fft_pkg3d_sparse->nfft;
    nfft2 = nfft/2;
    nkf1 = cp_sclr_fft_pkg3d_sparse->nkf1;
    nkf2 = cp_sclr_fft_pkg3d_sparse->nkf2;
    nkf3 = cp_sclr_fft_pkg3d_sparse->nkf3;
  }
  */

  ncoef1 = ncoef - 1;
  iupper = nstate;
  if(nstate % 2 == 1){
     iupper = nstate - 1;
  }

  //debug
#ifdef REAL_PP_DEBUG  
  for(i=1;i<=ncoef;i++){
    fccreal[i] = 0.0;
    fccimag[i] = 0.0;
  }
#endif
  /*
  for(i=1;i<=ncoef;i++){
    printf("wfkkkkkkkkk %.16lg %.16lg\n",ccreal[i],ccimag[i]);
  }
  fflush(stdout);
  exit(0);
  */
  for(i=0;i<numThreads;i++)energyNl[i] = 0.0;
  for(i=0;i<numThreads*numAtomTot;i++){
    fxThreads[i] = 0.0;
    fyThreads[i] = 0.0;
    fzThreads[i] = 0.0;    
  }

/*=================================================================*/
/*  get the forces on the coefs of each state                      */

  time_st = omp_get_wtime();

  omp_set_num_threads(numThreads);
  #pragma omp parallel private(iThread,is,ioff,ioff2)
  {
    iThread = omp_get_thread_num();
    #pragma omp for
    for(is=1 ; is<= iupper; is+=2 ){
      ioff = (is-1)*ncoef;
      ioff2 = (is)*ncoef;

/*==========================================================================*/
/* 1) get the wave functions in real space two at a time                    */
/*   I) double pack the complex zfft array with two real wavefunctions      */
     
      //if(realSparseOpt==0){
      if(fftw3dFlag==0){
	dble_pack_coef(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
		    zfft_threads[iThread],cp_sclr_fft_pkg3d_sm);
      }
      else{
	dble_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
			   zfft_threads[iThread],cp_sclr_fft_pkg3d_sm);
      }
      //}else{
      /*
        if(fftw3dFlag==0){
          dble_pack_coef(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
                      zfft_threads[iThread],cp_sclr_fft_pkg3d_sparse);
        }
        else{
          dble_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],&ccreal[ioff2],&ccimag[ioff2],
                             zfft_threads[iThread],cp_sclr_fft_pkg3d_sparse);
        }
      }//endif realSparseOpt
      */

/*--------------------------------------------------------------------------*/
/* II) fourier transform the wavefunctions to real space                    */
/*      convention exp(-igr)                                                */

      //if(realSparseOpt==0){
      if(fftw3dFlag==0){
	para_fft_gen3d_fwd_to_r(zfft_threads[iThread],zfft_tmp_threads[iThread],cp_sclr_fft_pkg3d_sm);
      }
      else{
	para_fft_gen3d_fwd_to_r_fftw3d(zfft_threads[iThread],cp_sclr_fft_pkg3d_sm,iThread);
	// x leading
      }
      //}
      /*
      else{
	if(fftw3dFlag==0){
	  para_fft_gen3d_fwd_to_r(zfft_threads[iThread],zfft_tmp_threads[iThread],cp_sclr_fft_pkg3d_sparse);
	}
	else{
	  para_fft_gen3d_fwd_to_r_fftw3d(zfft_threads[iThread],cp_sclr_fft_pkg3d_sparse,iThread);
	  // x leading
	}
      }
      */

/*==========================================================================*/
/* 2) get v|psi> in g space and store it in zfft                            */
/*   I) get  v|psi> in real space                                           */
   
      if(pseudoRealFlag==1){
	memcpy(&zfft_tmp_threads[iThread][1],&zfft_threads[iThread][1],nfft*sizeof(double));
	if(nlppForceOnly==0)cp_vpsi(zfft_threads[iThread],v_ks,nfft);  
	cp->pseudo.pseudoReal.energyCalcFlag = 1;
	controlEnergyNlppReal(cp,class,general_data,zfft_tmp_threads[iThread],
			      zfft_threads[iThread],1,&energyNl[iThread],
			      &fxThreads[iThread*numAtomTot],
			      &fyThreads[iThread*numAtomTot],
			      &fzThreads[iThread*numAtomTot],
			      cp_sclr_fft_pkg3d_sm);
      }
      else{
	cp_vpsi(zfft_threads[iThread],v_ks,nfft);
      }
      /*
      for(i=1;i<=nfft2;i++){
        printf("vvvvvvvvks %lg\n",v_ks[i]);
      }
      fflush(stdout);
      exit(0);
      */
      
      
/*--------------------------------------------------------------------------*/
/*  II) fourier transform  to g-space                                       */
/*     convention exp(igr)                                                  */

      //if(realSparseOpt==0){
      if(fftw3dFlag==0){
	para_fft_gen3d_bck_to_g(zfft_threads[iThread],zfft_tmp_threads[iThread],cp_sclr_fft_pkg3d_sm);
      }
      else{
	para_fft_gen3d_bck_to_g_fftw3d(zfft_threads[iThread],cp_sclr_fft_pkg3d_sm,iThread);
      }
      //}
      /*
      else{
        if(fftw3dFlag==0){
          para_fft_gen3d_bck_to_g(zfft_threads[iThread],zfft_tmp_threads[iThread],cp_sclr_fft_pkg3d_sparse);
        }
        else{
          para_fft_gen3d_bck_to_g_fftw3d(zfft_threads[iThread],cp_sclr_fft_pkg3d_sparse,iThread);
        }
      }
      */

/*==========================================================================*/
/* 3) get forces on coefficients by double unpacking the array zfft         */
    
      //if(realSparseOpt==0){
	if(fftw3dFlag==0){
	  dble_upack_coef_sum(&fccreal[ioff],&fccimag[ioff],
			  &fccreal[ioff2],&fccimag[ioff2],
			  zfft_threads[iThread],cp_sclr_fft_pkg3d_sm);
	}
	else{
	  dble_upack_coef_sum_fftw3d(&fccreal[ioff],&fccimag[ioff],
			  &fccreal[ioff2],&fccimag[ioff2],
			  zfft_threads[iThread],cp_sclr_fft_pkg3d_sm);

	}
      //}
      /*
      else{
        if(fftw3dFlag==0){
          dble_upack_coef_sum(&fccreal[ioff],&fccimag[ioff],
                          &fccreal[ioff2],&fccimag[ioff2],
                          zfft_threads[iThread],cp_sclr_fft_pkg3d_sparse);
        }
        else{
          dble_upack_coef_sum_fftw3d(&fccreal[ioff],&fccimag[ioff],
                          &fccreal[ioff2],&fccimag[ioff2],
                          zfft_threads[iThread],cp_sclr_fft_pkg3d_sparse);

        }
      }//endif realSparseOpt
      */
    }//endfor is
  }//end omp

  
  //if(fftw3dFlag==0){
#ifdef REAL_PP_DEBUG  
    for(i=1;i<=ncoef;i++){
      printf("forceeeee %i %.16lg %.16lg\n",i,fccreal[i],fccimag[i]);
    }
    fflush(stdout);
    exit(0);
#endif
  //}
 
  if(pseudoRealFlag==1&&forceCalcFlag==1){
    for(iThread=0;iThread<numThreads;iThread++){
      for(i=0;i<numAtomTot;i++){
	fx[i+1] += fxThreads[iThread*numAtomTot+i];
        fy[i+1] += fyThreads[iThread*numAtomTot+i];
        fz[i+1] += fzThreads[iThread*numAtomTot+i];
      }
    }
  } 
  if(pseudoRealFlag==1){
    for(iThread=0;iThread<numThreads;iThread++){
      stat_avg->cp_enl += energyNl[iThread];
    }
  }


/*==========================================================================*/
/*==========================================================================*/
/* 4) if there is an odd number of states, go through                       */
/*      the same procedure using sng_packs                                  */

  energyNl[0] = 0.0;
  for(i=0;i<numAtomTot;i++){
    fxThreads[i] = 0.0;
    fyThreads[i] = 0.0;
    fzThreads[i] = 0.0;
  }

  if(nstate % 2 != 0){
    is = nstate;
    ioff = (is -1)*ncoef;

/*--------------------------------------------------------------------------*/
/*   I) sngl pack                                                           */
    
    //if(realSparseOpt==0){  
      if(fftw3dFlag==0){
	sngl_pack_coef(&ccreal[ioff],&ccimag[ioff],zfft_threads[0],cp_sclr_fft_pkg3d_sm);
      }
      else{
	sngl_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],zfft_threads[0],
			       cp_sclr_fft_pkg3d_sm);
      }
    //}
    /*
    else{
      if(fftw3dFlag==0){
        sngl_pack_coef(&ccreal[ioff],&ccimag[ioff],zfft_threads[0],cp_sclr_fft_pkg3d_sparse);
      }
      else{
        sngl_pack_coef_fftw3d(&ccreal[ioff],&ccimag[ioff],zfft_threads[0],
                               cp_sclr_fft_pkg3d_sparse);
      }
    }
    */

/*--------------------------------------------------------------------------*/
/* II) fourier transform the wavefunctions to real space                    */
/*      convention exp(-igr)                                                */

    //if(realSparseOpt==0){ 
    if(fftw3dFlag==0){
      para_fft_gen3d_fwd_to_r(zfft_threads[0],zfft_tmp_threads[0],cp_sclr_fft_pkg3d_sm);
    }
    else{ 
      para_fft_gen3d_fwd_to_r_fftw3d(zfft_threads[0],cp_sclr_fft_pkg3d_sm,0);
    }
    //}
    /*
    else{
      if(fftw3dFlag==0){
        para_fft_gen3d_fwd_to_r(zfft_threads[0],zfft_tmp_threads[0],cp_sclr_fft_pkg3d_sparse);
      }
      else{
        para_fft_gen3d_fwd_to_r_fftw3d(zfft_threads[0],cp_sclr_fft_pkg3d_sparse,0);
      }
    }
    */

/*==========================================================================*/
/* 5) get v|psi> in g space and store it in zfft                            */
/*   I) get  v|psi> in real space                                           */

    if(pseudoRealFlag==1){
      memcpy(&zfft_tmp_threads[0][1],&zfft_threads[0][1],nfft*sizeof(double));
      if(nlppForceOnly==0)cp_vpsi(zfft_threads[0],v_ks,nfft);
      cp->pseudo.pseudoReal.energyCalcFlag = 1;
      controlEnergyNlppReal(cp,class,general_data,zfft_tmp_threads[0],
			    zfft_threads[0],0,&energyNl[0],&fxThreads[0],
                            &fyThreads[0],&fzThreads[0],cp_sclr_fft_pkg3d_sm);
    }
    else{
      cp_vpsi(zfft_threads[0],v_ks,nfft);
    }

/*--------------------------------------------------------------------------*/
/*   II) fourier transform the result back to g-space */
/*     convention exp(igr)  */
    //if(realSparseOpt==0){
      if(fftw3dFlag==0){
	para_fft_gen3d_bck_to_g(zfft_threads[0],zfft_tmp_threads[0],cp_sclr_fft_pkg3d_sm);
      }
      else{
	para_fft_gen3d_bck_to_g_fftw3d(zfft_threads[0],cp_sclr_fft_pkg3d_sm,0);
      }
    //}
    /*
    else{
      if(fftw3dFlag==0){
        para_fft_gen3d_bck_to_g(zfft_threads[0],zfft_tmp_threads[0],cp_sclr_fft_pkg3d_sparse);
      }
      else{
        para_fft_gen3d_bck_to_g_fftw3d(zfft_threads[0],cp_sclr_fft_pkg3d_sparse,0);
      }
    }
    */

/*==========================================================================*/
/* 6) get forces on coefficients by single unpacking the array zfft         */

    //if(realSparseOpt==0){
    if(fftw3dFlag==0){
      sngl_upack_coef_sum(&fccreal[ioff],&fccimag[ioff],zfft_threads[0],
			  cp_sclr_fft_pkg3d_sm);
    }
    else{
      sngl_upack_coef_sum_fftw3d(&fccreal[ioff],&fccimag[ioff],zfft_threads[0],
			    cp_sclr_fft_pkg3d_sm);
    }
    //}
    /*
    else{
      if(fftw3dFlag==0){
        sngl_upack_coef_sum(&fccreal[ioff],&fccimag[ioff],zfft_threads[0],
                            cp_sclr_fft_pkg3d_sparse);
      }
      else{
        sngl_upack_coef_sum_fftw3d(&fccreal[ioff],&fccimag[ioff],zfft_threads[0],
                              cp_sclr_fft_pkg3d_sparse);
      }
    }
    */

  }/* endif: odd number of states*/

  if(pseudoRealFlag==1&&forceCalcFlag==1){
    for(i=0;i<numAtomTot;i++){
      fx[i+1] += fxThreads[i];
      fy[i+1] += fyThreads[i];
      fz[i+1] += fzThreads[i];
    }
  }
  stat_avg->cp_enl += energyNl[0];

  time_end = omp_get_wtime();

  //cp_sclr_fft_pkg3d_sm->cputime3 += time_end-time_st;

/*==========================================================================*/
/* 7) If there is an electron KE density dependent functional, calculate    */
/*    this contribution to the force                                        */

  time_st = omp_get_wtime();
  
  if(cp_tau_functional==1)
    coef_force_tau_fun_hybrid(cpewald,nstate,ccreal,ccimag,fccreal,fccimag,
                              cre_scr,cim_scr,zfft_threads[0],zfft_tmp_threads[0],
			      v_ks_tau,ak2_sm,
                              pvten_cp,cp_ptens_calc,hmati,communicate,
                              icoef_form,icoef_orth,ifcoef_form,
                              cp_sclr_fft_pkg3d_sm);

/*==========================================================================*/
/* 8) If doing minimization, Fourier transform the KS potential to g-space  */
/*    and unpack it into the diagonal Hessian                               */

  if(cp_min_on > 0){
    for(i=1;i<=ncoef;i++){ cp_hess_re[i] = cp_hess_im[i] = 0.0;}
    //if(realSparseOpt==0){
      if(fftw3dFlag==0){
	sngl_pack_rho(zfft_threads[0],v_ks,cp_sclr_fft_pkg3d_sm);
	para_fft_gen3d_bck_to_g(zfft_threads[0],zfft_tmp_threads[0],cp_sclr_fft_pkg3d_sm);
	sngl_upack_coef(cp_hess_re,cp_hess_im,zfft_threads[0],cp_sclr_fft_pkg3d_sm);
	//printf("hhhess %lg %lg\n",cp_hess_re[1],cp_hess_im[1]);
      }
      else{
	// if using fftw3d, v_ks is x leading, we don't need transpose
	sngl_pack_rho(zfft_threads[0],v_ks,cp_sclr_fft_pkg3d_sm);
	/*
	for(i=0;i<nkf3;i++){
	  for(j=0;j<nkf2;j++){
	    for(k=0;k<nkf1;k++){
	      igrid = i*nkf2*nkf1+j*nkf1+k; //z leading
	      jgrid = k*nkf2*nkf3+j*nkf3+i; //x leading
	      zfft_tmp_threads[0][2*igrid+1] = zfft_threads[0][2*jgrid+1];
	      zfft_tmp_threads[0][2*igrid+2] = zfft_threads[0][2*jgrid+2];
	    }//endfor k
	  }//endfor j
	}//endfor i
	memcpy(&(zfft_threads[0][1]),&(zfft_tmp_threads[0][1]),nfft*sizeof(double));
	*/
	para_fft_gen3d_bck_to_g_fftw3d(zfft_threads[0],cp_sclr_fft_pkg3d_sm,0);
	
	sngl_upack_coef_fftw3d(cp_hess_re,cp_hess_im,zfft_threads[0],cp_sclr_fft_pkg3d_sm);
	//printf("hhhess %lg %lg\n",cp_hess_re[1],cp_hess_im[1]);
      }
    //}
    /*
    else{
      if(fftw3dFlag==0){
        sngl_pack_rho(zfft_threads[0],v_ks,cp_sclr_fft_pkg3d_sm);
        para_fft_gen3d_bck_to_g(zfft_threads[0],zfft_tmp_threads[0],cp_sclr_fft_pkg3d_sparse);
        sngl_upack_coef(cp_hess_re,cp_hess_im,zfft_threads[0],cp_sclr_fft_pkg3d_sparse);
        //printf("hhhess %lg %lg\n",cp_hess_re[1],cp_hess_im[1]);
      }
      else{
        // if using fftw3d, v_ks is x leading, we don't need transpose
        sngl_pack_rho(zfft_threads[0],v_ks,cp_sclr_fft_pkg3d_sm);
        para_fft_gen3d_bck_to_g_fftw3d(zfft_threads[0],cp_sclr_fft_pkg3d_sparse,0);
       
        sngl_upack_coef_fftw3d(cp_hess_re,cp_hess_im,zfft_threads[0],cp_sclr_fft_pkg3d_sparse);
        //printf("hhhess %lg %lg\n",cp_hess_re[1],cp_hess_im[1]);
      }
    }
    */
  }// endif cp_min_on


/*==========================================================================*/
/* 9) calculate the kinetic energy term and add its contribution to the force*/

  //printf("I'm here kinetic energy!\n");
  tpi = 2.0*M_PI;
  eke = 0.0;
  //printf("ak2_sm[1] %lg\n",ak2_sm[1]);
  /*
  for(is=1;is<=nstate;is++){
    ioff = (is-1)*ncoef;
    for(i=1; i<= ncoef1 ; i++){
      iis = ioff + i;
      if(isnan(ccreal[iis])==1)printf("is %i i %i ccreal[iis] %lg\n",is,i,ccreal[iis]);
      if(isnan(ccimag[iis])==1)printf("is %i i %i ccimag[iis] %lg\n",is,i,ccreal[iis]);
      if(isnan(ccimag[iis])==1)printf("is %i i %i ccimag[iis] %lg\n",is,i,ccreal[iis]);
    }
  }
  */

  #pragma omp parallel for private(is,ioff,i,iis,nis) reduction(+:eke)
  for(is=1 ; is<= nstate ; is++){
    ioff = (is-1)*ncoef;
    for(i=1; i<= ncoef1 ; i++){
      iis = ioff + i;
      fccreal[iis] -= 2.0*ak2Kinetic[i]*ccreal[iis];
      fccimag[iis] -= 2.0*ak2Kinetic[i]*ccimag[iis];
      eke += (2.0*ak2Kinetic[i]*(ccreal[iis]*ccreal[iis] + ccimag[iis]*ccimag[iis]));
    }/*endfor i*/
   nis = is*ncoef;
   fccimag[nis] = 0.0;
  }/*endfor*/

  eke *= .50;
  *eke_ret = eke;

/*================================================================================*/
/* 10) If doing minimization, calculat kinetic contribution to diagonal Hessian   */

  if(cp_min_on){    
    for(i=1; i<= ncoef ; i++){
      cp_hess_re[i] *= -1;
      cp_hess_im[i] *= -1;
    }/* endfor */
    for(i=1; i<= ncoef1 ; i++){
      //cp_hess_re[i] += 2.0*ak2_sm[i];
      //cp_hess_im[i] += 2.0*ak2_sm[i];
      cp_hess_re[i] += 2.0*ak2Kinetic[i];
      cp_hess_im[i] += 2.0*ak2Kinetic[i];      
    }/* endfor */
  }/* endif cp_min_on */

/*===================================================================*/
/* Get a debug tool for this system size                             */

/*==========================================================================*/
/* 9) calculate kinetic contribution to pressure tensor                     */

  if(cp_ptens_calc == 1) {
   tpi = 2.0*M_PI;
   for(is=1; is<= nstate; is++){
     ioff = (is-1)*ncoef;
     for(i=1; i<= ncoef1; i++){
       iis = ioff + i;
       aka = (double)kastore_sm[i];
       akb = (double)kbstore_sm[i];
       akc = (double)kcstore_sm[i];

       xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
       yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
       zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;

       cfact = 2.0*(ccreal[iis]*ccreal[iis] + ccimag[iis]*ccimag[iis]);

       pvten_cp[1] += xk*xk*cfact;
       pvten_cp[2] += xk*yk*cfact;
       pvten_cp[3] += xk*zk*cfact;

       pvten_cp[4] += xk*yk*cfact;
       pvten_cp[5] += yk*yk*cfact;
       pvten_cp[6] += yk*zk*cfact;

       pvten_cp[7] += xk*zk*cfact;
       pvten_cp[8] += yk*zk*cfact;
       pvten_cp[9] += zk*zk*cfact;

     }/*endfor*/
   }/*endfor*/
  }/*endif */

/*==========================================================================*/
/* 11) Free local memories						    */

  cfree(energyNl);
  cfree(fxThreads);
  cfree(fyThreads);
  cfree(fzThreads);

  time_end = omp_get_wtime();

  //cp_sclr_fft_pkg3d_sm->cputime4 += time_end-time_st;

/* fine fertig terminado finito */
/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/


