/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                 Module: read_coef_dvr.c                                  */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_coords_cp_entry.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coef_alloc_init_dvr(CP *cp,int cp_min_on, double *tot_memory)

/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*    Local Variables   */

#include "../typ_defs/typ_mask.h"

  int myid = cp->communicate.myid;
  int np_states = cp->communicate.np_states;
  int i,ip,nread,is;

  int par_size_up,par_size_dn;
  int ncoef_up_tot,ncoef_dn_tot;
  int nstate,nstate2;
  double *dvrc,*dvrvc,*dvrfc;
  int *ioff_up,*ioff_upt,*ioff_dn,*ioff_dnt;
  double mem_test;

  int cp_lda         = cp->cpopts.cp_lda;
  int cp_lsda        = cp->cpopts.cp_lsda;
  int cp_norb        = cp->cpopts.cp_norb;


  int pi_beads_proc  = cp->cpcoeffs_info.pi_beads_proc;

  int nstate_up      = cp->cpcoeffs_info.nstate_up;
  int nstate_dn      = cp->cpcoeffs_info.nstate_dn;
  int nstate_up_proc = cp->cpcoeffs_info.nstate_up_proc;
  int nstate_dn_proc = cp->cpcoeffs_info.nstate_dn_proc;
  int nstate_max_up            = cp->cp_comm_state_pkg_dvr_up.nstate_max;
  int nstate_max_dn            = cp->cp_comm_state_pkg_dvr_dn.nstate_max;

  int ncoef_up       = cp->cpcoeffs_info.ncoef;
  int ncoef_dn       = cp->cpcoeffs_info.ncoef;
  int ncoef_up_proc  = cp->cp_comm_state_pkg_dvr_up.nstate_ncoef_proc;
  int ncoef_dn_proc  = cp->cp_comm_state_pkg_dvr_dn.nstate_ncoef_proc;

  int nstate_ncoef_proc_max_up = cp->cp_comm_state_pkg_dvr_up.nstate_ncoef_proc_max;
  int nstate_ncoef_proc_max_dn = cp->cp_comm_state_pkg_dvr_dn.nstate_ncoef_proc_max;
  int nstate_up_proc_max = cp->cp_comm_state_pkg_dvr_up.nstate_proc_max;
  int nstate_dn_proc_max = cp->cp_comm_state_pkg_dvr_dn.nstate_proc_max;
 


/*==========================================================================*/
/* 0) Calculate the sizes */

  par_size_up  = nstate_max_up*(nstate_ncoef_proc_max_up);
  ncoef_up_tot = nstate_up_proc*ncoef_up;
  ncoef_up_tot = MAX(ncoef_up_tot,par_size_up);

  par_size_dn = nstate_max_dn*(nstate_ncoef_proc_max_dn);
  ncoef_dn_tot = MAX(nstate_dn_proc*ncoef_dn,1);
  ncoef_dn_tot = MAX(ncoef_dn_tot,par_size_dn);

  nstate  = MAX(nstate_up,nstate_dn);
  nstate2 = nstate*nstate;

  if(cp_lda==1){ncoef_dn=0;ncoef_dn_tot=1;}

/*==========================================================================*/
/* I) Malloc the variables */

  for(i=1;i<=pi_beads_proc;i++){
    cp->cpcoeffs_pos_dvr[i].dvrc_up =(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
    cp->cpcoeffs_pos_dvr[i].dvrc_dn =(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
    cp->cpcoeffs_pos_dvr[i].dvrvc_up=(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
    cp->cpcoeffs_pos_dvr[i].dvrvc_dn=(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
    cp->cpcoeffs_pos_dvr[i].dvrfc_up=(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
    cp->cpcoeffs_pos_dvr[i].dvrfc_dn=(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;

    cp->cpcoeffs_pos_dvr[i].ksmat_up    =(double *)cmalloc(nstate2*sizeof(double))-1;
    cp->cpcoeffs_pos_dvr[i].ksmat_dn    =(double *)cmalloc(nstate2*sizeof(double))-1;
    cp->cpcoeffs_pos_dvr[i].ksmat_eig_up=(double *)cmalloc(nstate*sizeof(double))-1;
    cp->cpcoeffs_pos_dvr[i].ksmat_eig_dn=(double *)cmalloc(nstate*sizeof(double))-1;
    cp->cpcoeffs_pos_dvr[i].norbmat_up  =(double *)cmalloc(nstate2*sizeof(double))-1;
    cp->cpcoeffs_pos_dvr[i].norbmat_dn  =(double *)cmalloc(nstate2*sizeof(double))-1;
    cp->cpcoeffs_pos_dvr[i].norbmati_up =(double *)cmalloc(nstate2*sizeof(double))-1;
    cp->cpcoeffs_pos_dvr[i].norbmati_dn =(double *)cmalloc(nstate2*sizeof(double))-1;
    cp->cpcoeffs_pos_dvr[i].ovmat_eigv_up=(double *)
                                  cmalloc(nstate2*sizeof(double))-1;
    cp->cpcoeffs_pos_dvr[i].ovmat_eigv_dn=(double *)
                                  cmalloc(nstate2*sizeof(double))-1;
  }/*endfor*/

  cp->cpcoeffs_info.ioff_up  = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_dn  = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_upt = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_dnt = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpopts.occ_up          = (double *) cmalloc(nstate*sizeof(double))-1;
  cp->cpopts.occ_dn          = (double *) cmalloc(nstate*sizeof(double))-1;
  cp->cpopts.rocc_sum_up = (double *) cmalloc(nstate*nstate*sizeof(double))-1;
  cp->cpopts.rocc_sum_dn = (double *) cmalloc(nstate*nstate*sizeof(double))-1;

  if(myid==0){
    mem_test = (pi_beads_proc*(ncoef_up_tot + ncoef_dn_tot)*3*sizeof(double)
             + nstate*4*sizeof(int) + nstate2*8*sizeof(double)
             + nstate*4*sizeof(double)) *1.0e-6;
  }/*endif for myid==0*/

/*==========================================================================*/
/* II) Assign the offsets */

  ioff_up  = cp->cpcoeffs_info.ioff_up;
  ioff_upt = cp->cpcoeffs_info.ioff_upt;
  ioff_dn  = cp->cpcoeffs_info.ioff_dn;
  ioff_dnt = cp->cpcoeffs_info.ioff_dnt;

  for(is=1;is<=nstate;is++){ioff_up[is]=(is-1)*ncoef_up;}
  for(is=1;is<=nstate;is++){ioff_dn[is]=(is-1)*ncoef_dn;}

  if(np_states==1){
    for(is=1;is<=nstate;is++){ioff_upt[is]=(is-1)*ncoef_up;}
    for(is=1;is<=nstate;is++){ioff_dnt[is]=(is-1)*ncoef_dn;}
  }else{
    for(is=1;is<=nstate;is++){ioff_upt[is]=(is-1)*ncoef_up_proc;}
    for(is=1;is<=nstate;is++){ioff_dnt[is]=(is-1)*ncoef_dn_proc;}
  }/*endif*/

/*========================================================================*/
/* III) Initialize coeficient arrays and form flags                       */


  for(ip=1;ip<=pi_beads_proc;ip++){
    cp->cpcoeffs_pos_dvr[ip].icoef_form_up  = 0;
    cp->cpcoeffs_pos_dvr[ip].ivcoef_form_up = 0;
    cp->cpcoeffs_pos_dvr[ip].ifcoef_form_up = 1;
    cp->cpcoeffs_pos_dvr[ip].icoef_orth_up  = 1;
    if(cp_norb>0){cp->cpcoeffs_pos_dvr[ip].icoef_orth_up = 0;}
    cp->cpcoeffs_pos_dvr[ip].ivcoef_orth_up  = 1;
    if(cp_norb>0){cp->cpcoeffs_pos_dvr[ip].ivcoef_orth_up = 0;}

    dvrc  = cp->cpcoeffs_pos_dvr[ip].dvrc_up;
    dvrvc = cp->cpcoeffs_pos_dvr[ip].dvrvc_up;
    dvrfc = cp->cpcoeffs_pos_dvr[ip].dvrfc_up;
    nread = ncoef_up_tot;
    for(i=1;i<=nread;i++){
      dvrc[i]  = 0.0;
      dvrvc[i] = 0.0;
      dvrfc[i] = 0.0;
    }/*endfor*/
    if(cp_lsda==1){
      cp->cpcoeffs_pos_dvr[ip].icoef_form_dn  = 0;
      cp->cpcoeffs_pos_dvr[ip].ivcoef_form_dn = 0;
      cp->cpcoeffs_pos_dvr[ip].ifcoef_form_dn = 1;
      cp->cpcoeffs_pos_dvr[ip].icoef_orth_dn  = 1;
      if(cp_norb>0){cp->cpcoeffs_pos_dvr[ip].icoef_orth_dn = 0;}
      cp->cpcoeffs_pos_dvr[ip].ivcoef_orth_dn  = 1;
      if(cp_norb>0){cp->cpcoeffs_pos_dvr[ip].ivcoef_orth_dn = 0;}

      dvrc  = cp->cpcoeffs_pos_dvr[ip].dvrc_dn;
      dvrvc = cp->cpcoeffs_pos_dvr[ip].dvrvc_dn;
      dvrfc = cp->cpcoeffs_pos_dvr[ip].dvrfc_dn;

      nread = ncoef_dn_tot;
      for(i=1;i<=nread;i++){
        dvrc[i]  = 0.0;
        dvrvc[i] = 0.0;
        dvrfc[i] = 0.0;
      }/*endfor*/
    }/*endif:lsda*/

  }/*endfor:pi_bead*/

/*------------------------------------------------------------------------*/
/* For cp minimization */

  if(cp_min_on==1){
    par_size_up  = nstate_max_up*(nstate_ncoef_proc_max_up);
    ncoef_up_tot = nstate_up_proc_max*ncoef_up;
    ncoef_up_tot = MAX(ncoef_up_tot,par_size_up);

    par_size_dn = nstate_max_dn*(nstate_ncoef_proc_max_dn);
    ncoef_dn_tot = MAX(nstate_dn_proc_max*ncoef_dn,1);
    ncoef_dn_tot = MAX(ncoef_dn_tot,par_size_dn);

    if(cp_lda==1){ncoef_dn=0;ncoef_dn_tot=1;}

    for(i=1;i<=pi_beads_proc;i++){
      cp->cpcoeffs_pos_dvr[i].fpc_up =(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
      cp->cpcoeffs_pos_dvr[i].fpc_dn =(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
    }

    if(myid==0){
      mem_test += (pi_beads_proc*(ncoef_up_tot + ncoef_dn_tot)*sizeof(double))*1.0e-6;
    }
  }

  if(myid==0){
    *tot_memory += mem_test;
    printf("CP wave func allocation: %g Mbytes; Total memory: %g Mbytes\n",
           mem_test,*tot_memory);
  }

/*-----------------------------------------------------------------------*/
}/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coef_fetch_coefs_dvr(CP *cp, FILE *fp_dnameci, char *dnameci,
                               int initial_spread_opt)

/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*  Local Variables */
#include "../typ_defs/typ_mask.h"
  int upper,pflag,ip,is,igo,i;
  double *cre_up,*cre_dn;
  int  isoff,ipoff,iii;
  int ihaveit,iwantit;
  int myid_send,myid_recv;

/* Local Pointers */
  int myid               = cp->communicate.myid;
  int nproc              = cp->communicate.np;
  MPI_Comm world         = cp->communicate.world;
  int pi_beads           = cp->cpcoeffs_info.pi_beads;
  int nstate_up          = cp->cpcoeffs_info.nstate_up;
  int nstate_dn          = cp->cpcoeffs_info.nstate_dn;
  int istate_up_st       = cp->cpcoeffs_info.istate_up_st;
  int istate_up_end      = cp->cpcoeffs_info.istate_up_end;
  int istate_dn_st       = cp->cpcoeffs_info.istate_dn_st;
  int istate_dn_end      = cp->cpcoeffs_info.istate_dn_end;
  int pi_beads_st        = cp->cpcoeffs_info.pi_beads_proc_st;
  int pi_beads_end       = cp->cpcoeffs_info.pi_beads_proc_end;
  double *cre_up_tmp     = cp->cpscr.cpscr_wave.zfft;
  double *cre_dn_tmp     = cp->cpscr.cpscr_wave.zfft_tmp;

  int ncoef              = cp->cpcoeffs_info.ncoef;
  int cp_lsda            = cp->cpopts.cp_lsda;
  int cp_lda             = cp->cpopts.cp_lda;

  int ibinary,n;
  char *c_array1,*c_array2;
  float cre_dum;

  ibinary = cp->cpopts.iread_coef_binary;

  if( (myid == 0) && (ibinary == 1)){
    c_array1 = (char *) cmalloc(MAXWORD*sizeof(char ));
    c_array2 = (char *) cmalloc(MAXWORD*sizeof(char ));
  }/*endif*/

/*==========================================================================*/
/* 0) Start reading */

  if(myid==0){
    printf("Reading in DVR coefficients\n");
   if(ibinary == 0){
     readtoendofline(fp_dnameci);
   }else{
     fread(c_array1,sizeof(char),MAXWORD,fp_dnameci);
   }/*endif ibinary */
  }/*endif for myid==0*/
  if(nproc>1){Barrier(world);}

/*==========================================================================*/
/*  I) Read Up states                                                    */

  upper = pi_beads;
  pflag = 1;
  if(initial_spread_opt == 1){upper = 1;pflag=2;}

  fflush(stdout);  if(nproc>1){Barrier(world);}
  for(ip=1;ip<=upper;ip++){
    for(is=1;is<=nstate_up;is++){
      igo=0;
      if((is>=istate_up_st) && (is<=istate_up_end) &&
        (ip>=pi_beads_st)  && (ip<=pi_beads_end)){
        if(myid != 0){
          Ssend(&myid,1,MPI_INT,0,0,world);
        }/*endif*/
        igo=1;
      }/* endif */
      if(igo==0 && myid==0){
        Recv(&myid_recv,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,world);
      }

      if(myid==0){
        for(i=1;i<=ncoef;i++){
          if(ibinary == 0){
            fscanf(fp_dnameci,"%lf",&(cre_up_tmp[i]));
            readtoendofline(fp_dnameci);
          }else{
            n = 1;
            fread(&(cre_dum),sizeof(float),n,fp_dnameci);
            cre_up_tmp[i] = (double) cre_dum;
          }/*endif ibinary*/
        }/*endfor: read*/
        if(igo==0){
          Ssend(&(cre_up_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
        }/* endif igo*/
      } /* endif */

      if(igo==1){
        if(myid != 0){
          Recv(&(cre_up_tmp[1]),ncoef,MPI_DOUBLE,0,MPI_ANY_TAG,world);
        }/* endif */
        isoff=(is-istate_up_st)*ncoef;
        ipoff=(ip-pi_beads_st+1);
        cre_up = cp->cpcoeffs_pos_dvr[ipoff].dvrc_up;
        for(i=1;i<=ncoef;i++){
          cre_up[(i+isoff)] = cre_up_tmp[i];
        }/* endfor */
      }/* endif igo */
      if(nproc>1){Barrier(world);}
    }/*endfor:states*/
  }/* endfor ip*/


/*==========================================================================*/
/*  II) Read Dn states                                                    */

  if(cp_lsda == 1){

    if(myid==0){
      if(ibinary == 0){
        readtoendofline(fp_dnameci);
      }else{
        fread(c_array1,sizeof(char),MAXWORD,fp_dnameci);
      }/*endif ibinary*/
    }/*endif myid*/
    for(ip=1;ip<=upper;ip++){
      for(is=1;is<=nstate_dn;is++){
        igo=0;
        if((is>=istate_dn_st) && (is<=istate_dn_end) &&
           (ip>=pi_beads_st) && (ip<=pi_beads_end)){
          if(myid != 0){Ssend(&myid,1,MPI_INT,0,0,world);}
          igo=1;
        }/* endif */
        if(igo==0 && myid==0){
          Recv(&myid_recv,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,world);
        }
        if(myid==0){
          for(i=1;i<=ncoef;i++){
            if(ibinary == 0){
              fscanf(fp_dnameci,"%lf",&(cre_dn_tmp[i]));
              readtoendofline(fp_dnameci);
            }else{
              n = 1;
              fread(&(cre_dum),sizeof(float),n,fp_dnameci);
              cre_dn_tmp[i] = (double) cre_dum;
            }/*endif ibinary*/
          }/*endfor: ncoef*/
          if(igo==0){
            Ssend(&(cre_dn_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
          }/* endif igo*/
        } /* endif */
        if(igo==1){
          if(myid != 0){
            Recv(&(cre_dn_tmp[1]),ncoef,MPI_DOUBLE,0,MPI_ANY_TAG,world);
          }/* endif */
          isoff=(is-istate_dn_st)*ncoef;
          ipoff=(ip-pi_beads_st+1);
          cre_dn = cp->cpcoeffs_pos_dvr[ipoff].dvrc_dn;
          for(i=1;i<=ncoef;i++){
            cre_dn[(i+isoff)] = cre_dn_tmp[i];
          }/* endfor */
        }/* endif igo */
        if(nproc>1){Barrier(world);}
      }/*endfor:states*/
    }/* endfor ip*/
  }/* endif cp_lsda */


/*==========================================================================*/
/* III) Spread Up states */

   if(upper < pi_beads){

    for(ip=2;ip<=pi_beads;ip++){
     for(is=1;is<=nstate_up;is++){
      ihaveit = 0; iwantit = 0;
      if(pi_beads_st==1 && ((is >= istate_up_st) && (is <= istate_up_end))){
        ihaveit = myid+1;
        isoff=(is-istate_up_st)*ncoef;
        cre_up = cp->cpcoeffs_pos_dvr[1].dvrc_up;
        for(i=1;i<=ncoef;i++){
          cre_up_tmp[i] = cre_up[(i+isoff)];
        }/* endfor */
      }/* endif */
      if((is>=istate_up_st) && (is<=istate_up_end) &&
          (ip>=pi_beads_st) && (ip<=pi_beads_end)){
        iwantit = myid+1;
      }/* endif */
      myid_send=0; myid_recv=0;
      if(nproc>1){
       Allreduce(&ihaveit,&myid_send,1,MPI_INT,MPI_SUM,0,world);
       Allreduce(&iwantit,&myid_recv,1,MPI_INT,MPI_SUM,0,world);
       myid_send--; myid_recv--;
      } /* endif */
      if(myid_send != myid_recv){
       if(myid == myid_send){
        Ssend(&(cre_up_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
       }
       if(myid == myid_recv){
         Recv(&(cre_up_tmp[1]),ncoef,MPI_DOUBLE,myid_send,MPI_ANY_TAG,world);
       }
      } /* endif */
      if(iwantit==myid+1){
        isoff=(is-istate_up_st)*ncoef;
        ipoff=(ip-pi_beads_st+1);
        cre_up = cp->cpcoeffs_pos_dvr[ipoff].dvrc_up;
        for(i=1;i<=ncoef;i++){
          cre_up[(i+isoff)] = cre_up_tmp[i];
        }/* endfor */
      }/* endif iwantit*/
      if(nproc>1){Barrier(world);}
     }/* endfor is */
    }/* endfor ip */

   }/* endif upper<pi_beads*/

/*==========================================================================*/
/* III) Spread Dn states    */

   if(upper < pi_beads && cp_lsda == 1 && nstate_dn != 0){

     for(ip=2;ip<=pi_beads;ip++){
      for(is=1;is<=nstate_dn;is++){
       ihaveit = 0; iwantit = 0;
       if(pi_beads_st==1 && ((is >= istate_dn_st) && (is <= istate_dn_end))){
         ihaveit = myid+1;
         isoff=(is-istate_dn_st)*ncoef;
         cre_dn = cp->cpcoeffs_pos_dvr[1].dvrc_dn;
         for(i=1;i<=ncoef;i++){
           cre_dn_tmp[i] = cre_dn[(i+isoff)];
         }/* endfor */
       }/* endif */
       if((is>=istate_dn_st) && (is<=istate_dn_end) &&
           (ip>=pi_beads_st) && (ip<=pi_beads_end)){
         iwantit = myid+1;
       }/* endif */
       myid_send=0; myid_recv=0;
       if(nproc>1){
        Allreduce(&myid_send,&ihaveit,1,MPI_INT,MPI_SUM,0,world);
        Allreduce(&myid_recv,&iwantit,1,MPI_INT,MPI_SUM,0,world);
        myid_send--; myid_recv--;
       } /* endif */
       if(myid_send != myid_recv){
        if(myid == myid_send){
          Ssend(&(cre_dn_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
        }
        if(myid == myid_recv){
          Recv(&(cre_dn_tmp[1]),ncoef,MPI_DOUBLE,myid_send,MPI_ANY_TAG,world);
        }
       } /* endif */
       if(iwantit==myid+1){
         isoff=(is-istate_dn_st)*ncoef;
         ipoff=(ip-pi_beads_st+1);
         cre_dn = cp->cpcoeffs_pos_dvr[ipoff].dvrc_dn;
         for(i=1;i<=ncoef;i++){
           cre_dn[(i+isoff)] = cre_dn_tmp[i];
         }/* endfor */
       }/* endif iwantit*/
       if(nproc>1){Barrier(world);}
      }/* endfor is */
     }/* endfor ip */

   }/* endif upper<pi_beads and lsda*/

/* free locally assigned memory */
  if( (myid == 0) && (ibinary == 1)){
    cfree(c_array1);
    cfree(c_array2);
  }/*endif*/

/*-----------------------------------------------------------------------*/
  } /* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coef_fetch_vcoefs_dvr(CP *cp, FILE *fp_dnameci, char *dnameci)

/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int upper,pflag,ip,is,igo,i;
  double *vcre_up,*vcre_dn;
  int  isoff,ipoff;
  int ihaveit,iwantit;
  int myid_send,myid_recv;

/* Local Pointers */
  int myid               = cp->communicate.myid;
  int nproc              = cp->communicate.np;
  MPI_Comm world         = cp->communicate.world;
  int pi_beads           = cp->cpcoeffs_info.pi_beads;
  int nstate_up          = cp->cpcoeffs_info.nstate_up;
  int nstate_dn          = cp->cpcoeffs_info.nstate_dn;
  int istate_up_st       = cp->cpcoeffs_info.istate_up_st;
  int istate_up_end      = cp->cpcoeffs_info.istate_up_end;
  int istate_dn_st       = cp->cpcoeffs_info.istate_dn_st;
  int istate_dn_end      = cp->cpcoeffs_info.istate_dn_end;
  int pi_beads_st        = cp->cpcoeffs_info.pi_beads_proc_st;
  int pi_beads_end       = cp->cpcoeffs_info.pi_beads_proc_end;
  double *cre_up_tmp     = cp->cpscr.cpscr_wave.zfft;
  double *cre_dn_tmp     = cp->cpscr.cpscr_wave.zfft_tmp;

  int ncoef              = cp->cpcoeffs_info.ncoef;
  int cp_lsda            = cp->cpopts.cp_lsda;
  int cp_lda             = cp->cpopts.cp_lda;

  char *c_array1,*c_array2;
  int n;  /* for binary read */
  int ibinary;

  float cre_dum;

  ibinary = cp->cpopts.iread_coef_binary;

  if( (myid == 0) && (ibinary == 1)){
    c_array1 = (char *) cmalloc(MAXWORD*sizeof(char ));
    c_array2 = (char *) cmalloc(MAXWORD*sizeof(char ));
  }/*endif*/

/*==========================================================================*/
/* 0) Announce to the screen */

   if(myid==0){
     printf("Reading DVR coefficient velocities\n");
     if(ibinary == 0){
      readtoendofline(fp_dnameci);
     }else{
       fread(c_array1,sizeof(char),MAXWORD,fp_dnameci);
     }/*endif ibinary */
   }/*endif for myid==0*/

/*==========================================================================*/
/*  I) Read up states velocities                                            */

  for(ip=1;ip<=pi_beads;ip++){
    if(nproc>1){Barrier(world);}
      for(is=1;is<=nstate_up;is++){
        igo=0;
        if((is>=istate_up_st) && (is<=istate_up_end) &&
          (ip>=pi_beads_st) && (ip<=pi_beads_end)){
          if(myid != 0){
            Ssend(&myid,1,MPI_INT,0,0,world);
          }
          igo=1;
        }/* endif */
      if(igo==0 && myid==0){
        Recv(&myid_recv,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,world);
      }
      if(myid==0){
        for(i=1;i<=ncoef;i++){
          if(ibinary == 0){
            fscanf(fp_dnameci,"%lf",&(cre_up_tmp[i]));
            readtoendofline(fp_dnameci);
          }else{
            n = 1;
            fread(&(cre_dum),sizeof(float),n,fp_dnameci);
            cre_up_tmp[i] = (double) cre_dum;
          }/*endif ibinary */
        }/*endfor: write*/
        if(igo==0){
          Ssend(&(cre_up_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
        }/* endif igo*/
      } /* endif */
      if(igo==1){
        if(myid != 0){
          Recv(&(cre_up_tmp[1]),ncoef,MPI_DOUBLE,0,MPI_ANY_TAG,world);
        }/* endif */
        isoff=(is-istate_up_st)*ncoef;
        ipoff=(ip-pi_beads_st+1);
        vcre_up =   cp->cpcoeffs_pos_dvr[ipoff].dvrvc_up;
        for(i=1;i<=ncoef;i++){
          vcre_up[(i+isoff)] = cre_up_tmp[i];
        }/* endfor */
      }/* endif igo */
      if(nproc>1){Barrier(world);}

    }/*endfor:states*/
  }/* endfor ip*/

/*-----------------------------------------------------------------------*/
/*  A) Read dn states                                                    */

  if(cp_lsda == 1){

    if(myid==0){
      if(ibinary == 0){
        readtoendofline(fp_dnameci);
      }else{
        fread(c_array1,sizeof(char),MAXWORD,fp_dnameci);
      }/*endif ibinary*/
    }/*endif myid*/

    for(ip=1;ip<=pi_beads;ip++){
      for(is=1;is<=nstate_dn;is++){
        igo=0;
        if((is>=istate_dn_st) && (is<=istate_dn_end) &&
          (ip>=pi_beads_st) && (ip<=pi_beads_end)){
          if(myid != 0){
            Ssend(&myid,1,MPI_INT,0,0,world);
          }
          igo=1;
        }/* endif */
        if(igo==0 && myid==0){
          Recv(&myid_recv,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,world);
        }
        if(myid==0){
          for(i=1;i<=ncoef;i++){
            if(ibinary == 0){
              fscanf(fp_dnameci,"%lf",&(cre_dn_tmp[i]));
              readtoendofline(fp_dnameci);
            }else{
              n = 1;
              fread(&(cre_dum),sizeof(float),n,fp_dnameci);
              cre_dn_tmp[i] = (double) cre_dum;
            }/*endif ibinary*/
          }/*endfor: write*/
          if(igo==0){
            Ssend(&(cre_dn_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
          }/* endif igo*/
        } /* endif */
        if(igo==1){
          if(myid != 0){
            Recv(&(cre_dn_tmp[1]),ncoef,MPI_DOUBLE,0,MPI_ANY_TAG,world);
          }/* endif */
          isoff=(is-istate_dn_st)*ncoef;
          ipoff=(ip-pi_beads_st+1);
          vcre_dn =   cp->cpcoeffs_pos_dvr[ipoff].dvrvc_dn;
          for(i=1;i<=ncoef;i++){
            vcre_dn[(i+isoff)] = cre_dn_tmp[i];
          }/* endfor */
        }/* endif igo */
        if(nproc>1){Barrier(world);}
      }/*endfor:states*/
    }/* endfor ip*/
  }/* endif lsda */

/*-----------------------------------------------------------------------*/
}/* end routine */
/*==========================================================================*/




