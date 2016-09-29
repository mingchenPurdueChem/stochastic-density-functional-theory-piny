/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: read_coord                                   */
/*                                                                          */
/* This subprogram reads atm-atm_NHC vol-vol_NHC input for a MD on a        */ 
/* LD-classical potential energy surface (LD-PES)                           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_coords_entry.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_pimd_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void readHmatFrag(CLASS *classMini,GENERAL_DATA *generalDataMini,
		CP *cpMini,int cp_dual_grid_opt_on, double *dbox_rat, 
		int *box_rat)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */

  SIMOPTS *simopts		= &(generalDataMini->simopts);
  CELL *cell			= &(generalDataMini->cell);
  CLATOMS_INFO *clatoms_info	= &(classMini->clatoms_info);
  ATOMMAPS *atommaps		= &(classMini->atommaps);

  int ii,iii,upper;
  int i,j,ip;              /* Num: For loop counter                 */
  
  double sx,sy,sz;
  double vol,vol_cp,area;

/*  Local pointers */
  double *hmat           = cell->hmat;
  double *hmati          = cell->hmati;
  double *hmat_ewd       = cell->hmat_ewd;
  double *hmat_ewd_cp    = cell->hmat_ewd_cp;
  double *hmat_cp        = cell->hmat_cp;
  double *hmati_cp       = cell->hmati_cp;
  double *cp_box_center  = cell->cp_box_center;
  double *cp_box_center_rel  = cell->cp_box_center_rel;

  int hmat_int_typ       = cell->hmat_int_typ;
  int hmat_cons_typ      = cell->hmat_cons_typ;
  int iperd              = cell->iperd;

  int ensemble_flag;
  ensemble_flag          = generalDataMini->ensopts.nve 
                         + generalDataMini->ensopts.nvt
                         + generalDataMini->ensopts.npt_i;


/*========================================================================*/
/*  I) Assign Volumes */
  for(i=1;i<=9;i++)hmat_cp[i] = hmat[i];

  gethinv(hmat,hmati,&(vol),iperd);
  gethinv(hmat_cp,hmati_cp,&(vol_cp),iperd);

  general_data->cell.vol        = vol;
  general_data->cell.vol_cp     = vol_cp;
  general_data->cell.vol0       = vol;
  general_data->baro.vol        = vol;
  general_data->par_rahman.vol  = vol;
  general_data->stat_avg.vol    = vol;
  general_data->baro.x_lnv      = log(vol)/3.0;

   for(i=1;i<=9;i++){hmat_ewd_cp[i] = hmat_cp[i];}
   for(i=1;i<=9;i++){hmat_ewd[i] = hmat[i];}

  area = hmat[1]*hmat[5] - hmat[2]*hmat[4];
  general_data->cell.area       = area;
  general_data->baro.area       = area;
  general_data->par_rahman.area = area;

  *dbox_rat = hmat[1]/hmat_cp[1];
  *box_rat  = (int)(*dbox_rat);

/*========================================================================*/
/*  II) Determine if Box is Cubic for Fast Imaging in Period.c */

  if( (hmat[4] == 0.0) && (hmat[7] == 0.0) &&
      (hmat[2] == 0.0) && (hmat[8] == 0.0) &&
      (hmat[3] == 0.0) && (hmat[6] == 0.0) &&
      (ensemble_flag != 0)){
    cell->cubic_box_flag = 1;
  }else{
    cell->cubic_box_flag = 0;
  }/*endif*/

/*========================================================================*/
/*  III) Check cell                                                      */

  check_cell(&(general_data->cell),cp_dual_grid_opt_on,*dbox_rat,dnamei);

/*========================================================================*/
/* convert cp_box_center from xtal coordinates to cartesian coordinates   */

  sx = cp_box_center[1];
  sy = cp_box_center[2];
  sz = cp_box_center[3];

  cp_box_center[1] =sx*hmat[1]+sy*hmat[4]+sz*hmat[7];
  cp_box_center[2] =sx*hmat[2]+sy*hmat[5]+sz*hmat[8];
  cp_box_center[3] =sx*hmat[3]+sy*hmat[6]+sz*hmat[9];

  sx = 0.5;
  sy = 0.5;
  sz = 0.5;

  cp_box_center_rel[1] =sx*hmat_cp[1]+sy*hmat_cp[4]+sz*hmat_cp[7];
  cp_box_center_rel[2] =sx*hmat_cp[2]+sy*hmat_cp[5]+sz*hmat_cp[8];
  cp_box_center_rel[3] =sx*hmat_cp[3]+sy*hmat_cp[6]+sz*hmat_cp[9];

/*========================================================================*/
   }/* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void check_box_center(CELL *cell,PARA_FFT_PKG3D *para_fft_pkg3d_lg,int myid)
/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */

 int nkf1 = para_fft_pkg3d_lg->nkf1;
 int nkf2 = para_fft_pkg3d_lg->nkf2;
 int nkf3 = para_fft_pkg3d_lg->nkf3;

 double *hmati = cell->hmati;
 double *cp_box_center = cell->cp_box_center;
 double sx,sy,sz;
 double x,y,z;
 double ax,ay,az;
 double eps = 1.0e-7;
 double M_EPS = 1.0e-10;

/*========================================================================*/

  x =  cp_box_center[1];
  y =  cp_box_center[2];
  z =  cp_box_center[3];

  sx = x*(hmati)[1]+y*(hmati)[4]+z*(hmati)[7];
  sy = x*(hmati)[2]+y*(hmati)[5]+z*(hmati)[8];
  sz = x*(hmati)[3]+y*(hmati)[6]+z*(hmati)[9];

  ax = sx*(double)nkf1 + M_EPS;
  ay = sy*(double)nkf2 + M_EPS;
  az = sz*(double)nkf3 + M_EPS;

  ax -= (int)ax;
  ay -= (int)ay;
  az -= (int)az;

 if( (fabs(ax) > eps) || (fabs(ay) > eps) || (fabs(az) > eps)){
   if(myid == 0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The cp box center must be located on a grid point \n");
      printf("Acceptable error %lg  error is %lg %lg %lg \n",eps,
              fabs(ax),fabs(ay),fabs(az));
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
   }/*endif myid*/
 }

/*========================================================================*/
   }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void check_cell(CELL *cell,int cp_dual_grid_opt_on,double dbox_rat, char *dnamei)
/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */
/*----------------------------------------------------------------------*/
  int ii;

  double eps = 1.0e-8;
  int icheck_dual_flag1 = 0;    
  int icheck_dual_flag2 = 0;    
  int icheck_dual_flag4 = 0;    
  int icheck_dual_flag8 = 0;    

  double rat_a,rat_b,rat_c;
  double eps_a,eps_b,eps_c;

/*----------------------------------------------------------------------*/
/*  Local pointers */

  int hmat_int_typ           = cell->hmat_int_typ;
  int hmat_cons_typ          = cell->hmat_cons_typ;
  int iperd                  = cell->iperd;

  double *hmat               = cell->hmat;
  double *hmat_cp            = cell->hmat_cp;
  double *cp_box_center      = cell->cp_box_center;
  double *cp_box_center_rel  = cell->cp_box_center_rel;


/*----------------------------------------------------------------------*/
/*I) Check cell symmetry                                                */

  if(hmat_int_typ==1){
     if(  (hmat[2] != 0.0) || (hmat[3] != 0.0) || (hmat[6] != 0.0) ){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("The Cell matrix must be upper triangular\n");
       printf("in file \"%s\"\n",dnamei);
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/*end if*/
  }/*endif*/


  if(iperd==0 || iperd==4){
    if( (hmat[4] != 0.0) || (hmat[7] != 0.0) ||
        (hmat[2] != 0.0) || (hmat[8] != 0.0) ||
        (hmat[3] != 0.0) || (hmat[6] != 0.0) ||
        (hmat[1] == 0.0) || (hmat[5] == 0.0) ||
        (hmat[9] == 0.0) ){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The Cell matrix must be orthorhombic in clusters\n");
      printf("in file \"%s\"\n",dnamei);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*end if*/
  }/*endif*/

  if(iperd ==2){
    if( (hmat[7] != 0.0) || (hmat[8] != 0.0) ||
        (hmat[3] != 0.0) || (hmat[6] != 0.0) ||
        (hmat[9] == 0.0) ){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The Cell matrix must contain no c-coupling in 2d\n");
      printf("in file \"%s\"\n",dnamei);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*end if*/
  }/*endif*/

  if(hmat_cons_typ==1){
    if( (hmat[2] != 0.0) || (hmat[3] != 0.0) ||
        (hmat[4] != 0.0) || (hmat[6] != 0.0) || 
        (hmat[7] != 0.0) || (hmat[8] != 0.0) ){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The Cell matrix should contain no off-diagonal \n"); 
      printf("coupling with the orthorhombic constraint \n");
      printf("in file \"%s\"\n",dnamei);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);      
    }/* endif */ 
  }/* endif */ 

  if(hmat_cons_typ==2){
    if( (hmat[3] != 0.0) || (hmat[6] != 0.0) ||
        (hmat[7] != 0.0) || (hmat[8] != 0.0) ){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The Cell matrix should contain no c-coupling with \n");
      printf("the monoclinic constraint in file \"%s\"\n",dnamei);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);      
    }/* endif */ 
  }/* endif */ 

/*-------------------------------------------------------------------*/
/* Check dual box                                                    */

/* Check to see that the cp box center is input in crystal coordinates  */
 if(cp_dual_grid_opt_on >= 1){
      if(  ((cp_box_center[1] <  0.0 ) || (cp_box_center[1] >= 1.0 ))
         ||((cp_box_center[2] <  0.0 ) || (cp_box_center[2] >= 1.0 ))         
         ||((cp_box_center[3] <  0.0 ) || (cp_box_center[3] >= 1.0 )) ){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The cp_box_center must be input in crystal     \n");
      printf("coordinates for dualed systems in file \"%s\"\n",dnamei);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
      }/*endif*/
 }


 if(cp_dual_grid_opt_on >= 1){
 /* Check to see that both hmat_cp and hmat are ORTHORHOMIBIC */
    if( (hmat_cp[4] != 0.0) || (hmat_cp[7] != 0.0) ||
        (hmat_cp[2] != 0.0) || (hmat_cp[8] != 0.0) ||
        (hmat_cp[3] != 0.0) || (hmat_cp[6] != 0.0) ||
        (hmat_cp[1] == 0.0) || (hmat_cp[5] == 0.0) ||
        (hmat_cp[9] == 0.0) ){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The CP Cell matrix must be orthorhombic for dualed systems\n");
      printf("in file \"%s\"\n",dnamei);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*end if*/

    if( (hmat[4] != 0.0) || (hmat[7] != 0.0) ||
        (hmat[2] != 0.0) || (hmat[8] != 0.0) ||
        (hmat[3] != 0.0) || (hmat[6] != 0.0) ||
        (hmat[1] == 0.0) || (hmat[5] == 0.0) ||
        (hmat[9] == 0.0) ){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The Cell matrix must be orthorhombic for dualed systems\n");
      printf("in file \"%s\"\n",dnamei);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*end if*/
 }

/* Integer multiple check for proportional grid option */
 if(cp_dual_grid_opt_on == 1){
    eps = 1.0e-4;
    for(ii=1; ii<=9; ii+= 4){
      if(fabs(hmat_cp[ii]*1.0 - hmat[ii]) > eps){
        icheck_dual_flag1 = 1;
      }
      if(fabs(hmat_cp[ii]*2.0 - hmat[ii]) > eps){
        icheck_dual_flag2 = 1;
      }
      if(fabs(hmat_cp[ii]*4.0 - hmat[ii]) > eps){
        icheck_dual_flag4 = 1;
      }
      if(fabs(hmat_cp[ii]*8.0 - hmat[ii]) > eps){
        icheck_dual_flag8 = 1;
      }
    }/*endfor ii*/
    if( (icheck_dual_flag1 + icheck_dual_flag2
       + icheck_dual_flag4 + icheck_dual_flag8) == 4){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The box matrix must be 2 4 or 8 times the size of the cp box\n");
      printf("in file \"%s\"\n",dnamei);
      printf("when using the proportional dual grid option \n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }
  }/*cp_dual_grid_opt_on*/

  if(cp_dual_grid_opt_on == 2){
   if( (hmat_cp[1] > hmat[1]) || (hmat_cp[5] > hmat[5]) ||
       (hmat_cp[9] > hmat[9])){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The cp box matrix must be smaller than the \n");
      printf("classical box in file \"%s\"\n",dnamei);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
   }/*endif*/
  }/*endif cp_dual_grid_opt*/

/* check that all sides for PME dual opt have the same proportionality */
  if(cp_dual_grid_opt_on == 2){

    rat_a = hmat[1]/hmat_cp[1];  
    rat_b = hmat[5]/hmat_cp[5];  
    rat_c = hmat[9]/hmat_cp[9];  

    eps_a  = fabs(rat_a - dbox_rat);
    eps_b  = fabs(rat_b - dbox_rat);
    eps_c  = fabs(rat_c - dbox_rat);

    if(eps_a > eps || eps_b > eps || eps_c > eps){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("You are performing a dual_grid mixed simulation using \n");
       printf("the incommensuerate option for the grids  \n");
       printf("the ratio of the a,b,c edges must be the same \n");
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
      exit(1);
     }
  }

  if(cp_dual_grid_opt_on == 2){
    if(dbox_rat <= 2.0 ){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You are performing a dual_grid mixed simulation using \n");
       printf("the incommensuerate option for the grids with a  \n");
       printf("box ratio of %lg \n",dbox_rat);
       printf("Are you certain this is what you would like to do?   \n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    }
  }

/*========================================================================*/
}/* end routine */
/*==========================================================================*/
