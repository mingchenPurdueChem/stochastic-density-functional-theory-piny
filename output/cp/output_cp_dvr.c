/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: output_cp                                    */
/*                                                                          */
/* This subprogram provides output for a CP on the                          */ 
/* ground state Born-Oppenheimer GGA-LDA surface                            */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_output_cp_entry.h"
#include "../proto_defs/proto_output_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_output_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_dump_coef_cp_dvr(FILE *fp_dnamec,CP *cp,CLASS *class,int ibinary)

/*==========================================================================*/

{/* begin routine */
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

 int  ncoef         = cp->cpcoeffs_info.ncoef;
 int  istate_up_st  = cp->cpcoeffs_info.istate_up_st;
 int  istate_up_end = cp->cpcoeffs_info.istate_up_end;
 int  nstate_up     = cp->cpcoeffs_info.nstate_up;
 int  nstate_up_proc= cp->cpcoeffs_info.nstate_up_proc;

 int  istate_dn_st  = cp->cpcoeffs_info.istate_dn_st;
 int  istate_dn_end = cp->cpcoeffs_info.istate_dn_end;
 int  nstate_dn     = cp->cpcoeffs_info.nstate_dn;
 int  nstate_dn_proc= cp->cpcoeffs_info.nstate_dn_proc;
 int cp_lsda        = cp->cpopts.cp_lsda;
 int nrecv,isoff;
 int is,igo;
 int i,n;

 float cre_dum,cim_dum;

 char *c_array1,*c_array2;  /*for binary write */
 int csize = MAXWORD;

 MPI_Status stat;

/* Local pointers */
 double *cre_up_tmp     = cp->cpscr.cpscr_wave.zfft;
 double *cre_dn_tmp     = cp->cpscr.cpscr_wave.zfft_tmp;

 int myid       = class->communicate.myid;
 int np_states  = class->communicate.np_states;
 int nproc      = class->communicate.np;
 MPI_Comm world = class->communicate.world;

/*=======================================================================*/

  if( (myid == 0) && (ibinary == 1)){
    c_array1 = (char *) cmalloc(csize*sizeof(char ));
    c_array2 = (char *) cmalloc(csize*sizeof(char ));

    /* initialize array */
    for(i=0; i< csize-1; i++){
     c_array1[i] = ' ';
     c_array2[i] = ' ';
    }
     c_array1[csize -1] = '\0';
     c_array2[csize -1] = '\0';

  }/*endif */

/*=========================================================================*/
/* Write the coefficients                                                  */
/*-------------------------------------------------------------------------*/
/* Up states                                                               */

  if(myid==0){
    if(ibinary == 0){
      fprintf(fp_dnamec,"Up_State_Real       \n");
    }else{
      strcpy(c_array1,"Up_State_Real ");
      fwrite(c_array1,sizeof(char),csize,fp_dnamec);
    }/*endif ibinary */
  }/*endif myid*/

  if(nproc>1){Barrier(world);}

  for(is=1;is<=nstate_up;is++){
    igo=0;
    if((is>=istate_up_st) && (is<=istate_up_end)){
      isoff=(is-istate_up_st)*ncoef;
      for(i=1;i<=ncoef;i++){
        cre_up_tmp[i] = cp->cpcoeffs_pos_dvr[1].dvrc_up[(i+isoff)];
      }/* endfor */
      igo = 1;
      if(myid!=0){
        Ssend(&(cre_up_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
      }/*endif*/
    }/* endif */

    if(myid==0&&igo==0){
      Recv(&(cre_up_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
           world);
    }/* endif */
    if(myid==0){
      n = 1;
      for(i=1;i<=ncoef;i++){
        if(ibinary == 0){
          fprintf(fp_dnamec,"%.10g \n",cre_up_tmp[i]);
        }else{
          cre_dum = (float)cre_up_tmp[i];
          fwrite(&(cre_dum),sizeof(float),n,fp_dnamec);
        }/*endif ibinary */
      }/*endfor: write*/
    }/* endif myid */

    if(nproc>1){Barrier(world);}

  }/*endfor:states*/

/*-------------------------------------------------------------------------*/
/* Dn states                                                               */

  if(cp->cpopts.cp_lsda==1){
    if(myid==0){
      if(ibinary == 0){
        fprintf(fp_dnamec,"Dn_State_Real  \n");
      }else{
        strcpy(c_array1,"Dn_State_Real ");
        fwrite(c_array1,sizeof(char),csize,fp_dnamec);
      }/*endif ibinary */
    }/*endif myid*/

    if(nproc>1){Barrier(world);}

    for(is=1;is<=nstate_dn;is++){
      igo=0;
      if((is>=istate_dn_st) && (is<=istate_dn_end)){
        isoff=(is-istate_dn_st)*ncoef;
        for(i=1;i<=ncoef;i++){
          cre_dn_tmp[i] = cp->cpcoeffs_pos_dvr[1].dvrc_dn[(i+isoff)];
        }/* endfor */
        igo = 1;
        if(myid!=0){
          Ssend(&(cre_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
        }/*endif*/
      }/* endif */
      if(myid==0&&igo==0){
        Recv(&(cre_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
               world);
      }/* endif */
      if(myid==0){
        n = 1;
        for(i=1;i<=ncoef;i++){
          if(ibinary == 0){
            fprintf(fp_dnamec,"%.10g \n",(float)cre_dn_tmp[i]);
          }else{
            cre_dum = (float) cre_dn_tmp[i];
            fwrite(&(cre_dum),sizeof(double),n,fp_dnamec);
          }/*endif ibinary */
        }/*endfor: write*/
      }/* endif myid */
      if(nproc>1){Barrier(world);}
    }/*endfor:states*/
  }/*endif:lsda*/

/* free locally assigned memory */
/*
  if((myid == 0) && (ibinary == 1)){
    free(&c_array1[1]);
    free(&c_array2[1]);
  }
*/

/*==========================================================================*/
 }/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_dump_vcoef_cp_dvr(FILE *fp_dnamec,CP *cp,CLASS *class,
                             GENERAL_DATA *general_data,int ibinary)

/*==========================================================================*/

{/* begin routine */
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

 int  ncoef         = cp->cpcoeffs_info.ncoef;
 int  istate_up_st  = cp->cpcoeffs_info.istate_up_st;
 int  istate_up_end = cp->cpcoeffs_info.istate_up_end;
 int  nstate_up     = cp->cpcoeffs_info.nstate_up;
 int  nstate_up_proc= cp->cpcoeffs_info.nstate_up_proc;

 int  istate_dn_st  = cp->cpcoeffs_info.istate_dn_st;
 int  istate_dn_end = cp->cpcoeffs_info.istate_dn_end;
 int  nstate_dn     = cp->cpcoeffs_info.nstate_dn;
 int  nstate_dn_proc= cp->cpcoeffs_info.nstate_dn_proc;
 int cp_lsda        = cp->cpopts.cp_lsda;
 int nrecv,isoff;
 int is,igo;
 int i,n;
 float cre_dum,cim_dum;

 char *c_array1,*c_array2;  /*for binary write */
 int csize = MAXWORD;

 MPI_Status stat;

/* Local pointers */
 double *cre_up_tmp     = cp->cpscr.cpscr_wave.zfft;
 double *cre_dn_tmp     = cp->cpscr.cpscr_wave.zfft_tmp;

 int myid       = class->communicate.myid;
 int np_states  = class->communicate.np_states;
 int nproc      = class->communicate.np;
 MPI_Comm world = class->communicate.world;


  if( (myid == 0) && (ibinary == 1)){
    c_array1 = (char *) cmalloc(csize*sizeof(char ));
    c_array2 = (char *) cmalloc(csize*sizeof(char ));

    /* initialize array */
    for(i=0 ; i< csize-1; i++){
      c_array1[i] = ' ';
      c_array2[i] = ' ';
    }
    c_array1[csize-1] = '\0';
    c_array2[csize-1] = '\0';
  }/*endif*/

/*==========================================================================*/
/* Write the coefficient velocities                                         */

/*-------------------------------------------------------------------------*/
/* Up states                                                               */

  if((general_data->simopts.cp+general_data->simopts.cp_wave)==1){

    if(myid==0){
      if(ibinary == 0){
        fprintf(fp_dnamec,"V_Up_State_Real   \n");
      }else{
        strcpy(c_array1,"V_UP_State_Real");
        fwrite(c_array1,sizeof(char),csize,fp_dnamec);
      }/*endif ibinary */
    }/*endif myid*/

    if(nproc>1){Barrier(world);}

    for(is=1;is<=nstate_up;is++){
      igo=0;
      if((is>=istate_up_st) && (is<=istate_up_end)){
        isoff=(is-istate_up_st)*ncoef;
        for(i=1;i<=ncoef;i++){
          cre_up_tmp[i] = cp->cpcoeffs_pos_dvr[1].dvrvc_up[(i+isoff)];
        }/* endfor */
        igo = 1;
        if(myid!=0){
          Ssend(&(cre_up_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
        }/*endif*/
      }/* endif */

      if(myid==0&&igo==0){
        Recv(&(cre_up_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
               world);
      }/* endif */
      if(myid==0){
        for(i=1;i<=ncoef;i++){
          if(ibinary == 0){
            fprintf(fp_dnamec,"%.10g \n",cre_up_tmp[i]);
          }else{
            n = 1;
            cre_dum = (float) cre_up_tmp[i];
            fwrite(&(cre_dum),sizeof(float),n,fp_dnamec);
          }/*endif ibinary*/
        }/*endfor: write*/
      }/* endif myid */
      if(nproc>1){Barrier(world);}
    }/*endfor:states*/

/*-------------------------------------------------------------------------*/
/* Dn states                                                               */

    if(cp->cpopts.cp_lsda==1){
      if(myid==0){
        if(ibinary == 0){
          fprintf(fp_dnamec,"V_Dn_State_Real   \n");
        }else{
          strcpy(c_array1,"V_Dn_State_Real");
          fwrite(c_array1,sizeof(char),csize,fp_dnamec);
        }/*endif ibinary*/
      }/*endif myid*/

      if(nproc>1){Barrier(world);}

      for(is=1;is<=nstate_dn;is++){
        igo=0;
        if((is>=istate_dn_st) && (is<=istate_dn_end)){
          isoff=(is-istate_dn_st)*ncoef;
          for(i=1;i<=ncoef;i++){
            cre_dn_tmp[i] = cp->cpcoeffs_pos_dvr[1].dvrvc_dn[(i+isoff)];
          }/* endfor */
          igo = 1;
          if(myid!=0){
            Ssend(&(cre_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
          }/*endif*/
        }/* endif */

        if(myid==0&&igo==0){
          Recv(&(cre_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
               world);
        }/* endif */
        if(myid==0){
          for(i=1;i<=ncoef;i++){
            if(ibinary == 0){
              fprintf(fp_dnamec,"%.10g \n",cre_dn_tmp[i]);
            }else{
              n = 1;
              cre_dum = (float) cre_dn_tmp[i];
              fwrite(&(cre_dum),sizeof(float),n,fp_dnamec);
            }/*endif ibinary*/
          }/*endfor: write*/
        }/* endif myid */
        if(nproc>1){Barrier(world);}
      }/*endfor:states*/
    }/*endif:lsda*/

  }/*endif:write the velocities*/

/* free locally assigned memory */
/*
  if( (myid == 0) && (ibinary == 1)){
    free(&c_array1[1]);
    free(&c_array2[1]);
  }
*/

/*==========================================================================*/
}/* end routine */
/*==========================================================================*/

