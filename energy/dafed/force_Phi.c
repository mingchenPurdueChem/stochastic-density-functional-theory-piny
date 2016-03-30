#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_dafed_energy.h"
#include "../proto_defs/proto_math.h"

/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void force_Phi(DAFED_INFO *dinfo,DAFED *dafed, CLATOMS_POS *clatoms_pos)
{

  double array_ab[3], array_bc[3], array_cd[3], crs_1[3], crs_2[3], crs_12[3];
  double *clatoms_x  = clatoms_pos->x;
  double *clatoms_y  = clatoms_pos->y;
  double *clatoms_z  = clatoms_pos->z;
  int index_a, index_b, index_c, index_d;
  int i,j;
  int *atm_list = dafed->atm_list;
  double *Fx    = dafed->Fx;
  double *Fy    = dafed->Fy;
  double *Fz    = dafed->Fz;
  double s      = dafed->s;
  double ks     = dafed->ks;
  double m_crs_1, m_crs_2, m_v2, phi, phi_cos, phi_sin,pre,diff;
  double prod12,prod11,prod22,prod123,prod33;
  double a_ab_d[3], a_bc_d[3], a_cd_d[3], crs_1_d[3], crs_2_d[3], m_1[3], m_2[3];
  //double cos_x,cos_y,cos_z,sin_x,sin_y,sin_z;
  //double delta       = 0.00000001;
  //double test_ab[3],test_bc[3],test_cd[3],phi_test_1, phi_test_2, pe_test_1, pe_test_2, f_test, f_ana;
  double dev_tri_local[4][3][2]; //index of atom *x/y/z*cos/sin
  
  index_a = atm_list[0];
  index_b = atm_list[1];
  index_c = atm_list[2];
  index_d = atm_list[3];
  //printf("index_a %i index_b %i index_c %i index_d %i\n",index_a,index_b,index_c,index_d);
  array_ab[0] = clatoms_x[index_b] - clatoms_x[index_a];
  array_ab[1] = clatoms_y[index_b] - clatoms_y[index_a];
  array_ab[2] = clatoms_z[index_b] - clatoms_z[index_a];

  array_bc[0] = clatoms_x[index_c] - clatoms_x[index_b];
  array_bc[1] = clatoms_y[index_c] - clatoms_y[index_b];
  array_bc[2] = clatoms_z[index_c] - clatoms_z[index_b];

  array_cd[0] = clatoms_x[index_d] - clatoms_x[index_c];
  array_cd[1] = clatoms_y[index_d] - clatoms_y[index_c];
  array_cd[2] = clatoms_z[index_d] - clatoms_z[index_c];
  cross_product(array_ab,array_bc,crs_1);
  cross_product(array_bc,array_cd,crs_2);

  //m_crs_1 = 1.0/sqrt(dot(crs_1,crs_1));
  //m_crs_2 = 1.0/sqrt(dot(crs_2,crs_2));
  //m_v2    = 1.0/sqrt(dot(array_bc,array_bc));
  prod11  = 1.0/dot(crs_1,crs_1);
  prod22  = 1.0/dot(crs_2,crs_2);
  prod33  = 1.0/dot(array_bc,array_bc);
  prod12  = sqrt(prod11*prod22);
  prod123 = prod12*sqrt(prod33);


  phi_cos = dot(crs_1,crs_2)*prod12;
  cross_product(crs_1,crs_2,crs_12);
  phi_sin = dot(crs_12,array_bc)*prod123;

  phi = atan2(phi_sin,phi_cos);
  if(s<-3.5){dafed->s = phi - 0.001;s = dafed->s;}
  //printf("s %lg phi %lg\n",s,phi);
  dafed->q = phi;
  //if(Phi[ipt][i].s<=-180.0) Phi[ipt][i].s += 360.0;
  //if(Phi[ipt][i].s>180.0) Phi[ipt][i].s -= 360.0;
  diff = phi - s;
  if(diff<=-M_PI) diff += M_PI*2.0;
  if(diff>=M_PI) diff -= M_PI*2.0;
  dafed->pot_harm = 0.5*ks*diff*diff;
  dafed->Fs   = ks*diff;

  pre = -ks*diff;

  force_trifun(dev_tri_local,array_ab,array_bc,array_cd,crs_1,crs_2,crs_12,
               prod11,prod22,prod33,prod123,prod12,phi_cos,phi_sin);

     Fx[0] = pre*(-phi_sin*dev_tri_local[0][0][0]+phi_cos*dev_tri_local[0][0][1]);
     Fy[0] = pre*(-phi_sin*dev_tri_local[0][1][0]+phi_cos*dev_tri_local[0][1][1]);
     Fz[0] = pre*(-phi_sin*dev_tri_local[0][2][0]+phi_cos*dev_tri_local[0][2][1]);
     Fx[1] = pre*(-phi_sin*dev_tri_local[1][0][0]+phi_cos*dev_tri_local[1][0][1]);
     Fy[1] = pre*(-phi_sin*dev_tri_local[1][1][0]+phi_cos*dev_tri_local[1][1][1]);
     Fz[1] = pre*(-phi_sin*dev_tri_local[1][2][0]+phi_cos*dev_tri_local[1][2][1]);
     Fx[2] = pre*(-phi_sin*dev_tri_local[2][0][0]+phi_cos*dev_tri_local[2][0][1]);
     Fy[2] = pre*(-phi_sin*dev_tri_local[2][1][0]+phi_cos*dev_tri_local[2][1][1]);
     Fz[2] = pre*(-phi_sin*dev_tri_local[2][2][0]+phi_cos*dev_tri_local[2][2][1]);
     Fx[3] = pre*(-phi_sin*dev_tri_local[3][0][0]+phi_cos*dev_tri_local[3][0][1]);
     Fy[3] = pre*(-phi_sin*dev_tri_local[3][1][0]+phi_cos*dev_tri_local[3][1][1]);
     Fz[3] = pre*(-phi_sin*dev_tri_local[3][2][0]+phi_cos*dev_tri_local[3][2][1]);

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
inline void cross_product(double *a, double *b, double *result)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
    result[0] = a[1]*b[2]-a[2]*b[1];
    result[1] = a[2]*b[0]-a[0]*b[2];
    result[2] = a[0]*b[1]-a[1]*b[0];
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
inline double dot(double *dot_1, double *dot_2)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/

    return dot_1[0]*dot_2[0]+dot_1[1]*dot_2[1]+dot_1[2]*dot_2[2];
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
double get_force(double *crs_1, double *crs_2, double *crs_1_d, 
                 double *crs_2_d, double m_crs_1, double m_crs_2)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/

     double dot_1[3], dot_2[3];
     double force = 0.0;

/*==========================================================================*/

     dot_1[0] = crs_2[0]/m_crs_2;
     dot_1[1] = crs_2[1]/m_crs_2;
     dot_1[2] = crs_2[2]/m_crs_2;

     dot_2[0] = crs_1_d[0]/m_crs_1 + crs_1[0]/pow(m_crs_1,3)*dot(crs_1,crs_1_d);
     dot_2[1] = crs_1_d[1]/m_crs_1 + crs_1[1]/pow(m_crs_1,3)*dot(crs_1,crs_1_d);
     dot_2[2] = crs_1_d[2]/m_crs_1 + crs_1[2]/pow(m_crs_1,3)*dot(crs_1,crs_1_d);


     force += dot(dot_1,dot_2);
         
     dot_1[0] = crs_1[0]/m_crs_1;
     dot_1[1] = crs_1[1]/m_crs_1;
     dot_1[2] = crs_1[2]/m_crs_1;

     dot_2[0] = crs_2_d[0]/m_crs_2 + crs_2[0]/pow(m_crs_2,3)*dot(crs_2,crs_2_d);
     dot_2[1] = crs_2_d[1]/m_crs_2 + crs_2[1]/pow(m_crs_2,3)*dot(crs_2,crs_2_d);
     dot_2[2] = crs_2_d[2]/m_crs_2 + crs_2[2]/pow(m_crs_2,3)*dot(crs_2,crs_2_d);

     force += dot(dot_1,dot_2);

     //printf("force %lg\n",force);
     return force;
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
double get_phi(double *array_ab,double *array_bc, double *array_cd)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
     double phi, m_crs_1, m_crs_2, m_v2, phi_cos, phi_sin;
     double crs_1[3], crs_2[3], crs_12[3];     

     cross_product(array_ab,array_bc,crs_1);
     cross_product(array_bc,array_cd,crs_2);

     m_crs_1 = sqrt(dot(crs_1,crs_1));
     m_crs_2 = sqrt(dot(crs_2,crs_2));
     m_v2    = sqrt(dot(array_bc,array_bc));


     phi_cos = dot(crs_1,crs_2)/(m_crs_1*m_crs_2);
     cross_product(crs_1,crs_2,crs_12);
     phi_sin = dot(crs_12,array_bc)/(m_crs_1*m_crs_2*m_v2);

     phi = atan2(phi_sin,phi_cos);

     return phi;
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void force_trifun(double dev_tri_local[][3][2],double *array_ab, double *array_bc,
                    double *array_cd, double *crs_1, double *crs_2,double *crs_12,
                    double prod11, double prod22, double prod33,double prod123,double prod12,
		    double phi_cos, double phi_sin)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/

    //double a_ab_d[3],a_bc_d[3],a_cd_d[3],crs_1_d[3],crs_2_d[3],m_1[3],m_2[3];
    double temp1[3],temp2[3],temp3[3],m_1[3];
    double cos11 = phi_cos*prod11;
    double sin11 = phi_sin*prod11;
    double cos22 = phi_cos*prod22;
    double sin22 = phi_sin*prod22;
    double sin33 = phi_sin*prod33;
    double force;


/*-------------------------------------------------------------------------*/
/* force of x_a                                                            */

     temp1[0] = 0.0;//temp1 = crs_1_d
     temp1[1] = array_bc[2];
     temp1[2] = -array_bc[1];

     dev_tri_local[0][0][0] = dot(crs_2,temp1)*prod12-dot(crs_1,temp1)*cos11;
     cross_product(temp1,crs_2,m_1);
     dev_tri_local[0][0][1] = dot(m_1,array_bc)*prod123-dot(crs_1,temp1)*sin11;

/*-------------------------------------------------------------------------*/
/* force of y_a                                                            */

     temp1[0] = -array_bc[2];//temp1 = crs_1_d
     temp1[1] = 0.0;
     temp1[2] = array_bc[0];

     dev_tri_local[0][1][0] = dot(crs_2,temp1)*prod12-dot(crs_1,temp1)*cos11;
     cross_product(temp1,crs_2,m_1);
     dev_tri_local[0][1][1] = dot(m_1,array_bc)*prod123-dot(crs_1,temp1)*sin11;

/*-------------------------------------------------------------------------*/
/* force of z_a                                                            */

     temp1[0] = array_bc[1];//temp1 = crs_1_d
     temp1[1] = -array_bc[0];
     temp1[2] = 0.0;

     dev_tri_local[0][2][0] = dot(crs_2,temp1)*prod12-dot(crs_1,temp1)*cos11;
     cross_product(temp1,crs_2,m_1);
     dev_tri_local[0][2][1] = dot(m_1,array_bc)*prod123-dot(crs_1,temp1)*sin11;

/*-------------------------------------------------------------------------*/
/* force of x_b                                                            */
     temp1[0] = -1.0;//temp1 = a_bc_d
     temp1[1] = 0.0;
     temp1[2] = 0.0;

     temp2[0] = 0.0;//temp2 = crs_1_d
     temp2[1] = -array_bc[2]-array_ab[2];
     temp2[2] = array_bc[1]+array_ab[1];

     temp3[0] = 0.0;//temp3 = crs_2_d
     temp3[1] = array_cd[2];
     temp3[2] = -array_cd[1];

     dev_tri_local[1][0][0] = dot(crs_2,temp2)*prod12-dot(crs_1,temp2)*cos11+dot(crs_1,temp3)*prod12-dot(crs_2,temp3)*cos22;
     force = 0.0;
     cross_product(temp2,crs_2,m_1);
     force += dot(m_1,array_bc)*prod123-dot(crs_1,temp2)*sin11;
     cross_product(crs_1,temp3,m_1);
     force += dot(m_1,array_bc)*prod123-dot(crs_2,temp3)*sin22+dot(crs_12,temp1)*prod123-dot(array_bc,temp1)*sin33;
     dev_tri_local[1][0][1] = force;

/*-------------------------------------------------------------------------*/
/* force of y_b                                                            */

     temp1[0] = 0.0;//temp1 = a_bc_d
     temp1[1] = -1.0;
     temp1[2] = 0.0;

     temp2[0] = array_bc[2]+array_ab[2];//temp2 = crs_1_d
     temp2[1] = 0.0;
     temp2[2] = -array_bc[0]-array_ab[0];

     temp3[0] = -array_cd[2];//temp3 = crs_2_d
     temp3[1] = 0.0;
     temp3[2] = array_cd[0];

     dev_tri_local[1][1][0] = dot(crs_2,temp2)*prod12-dot(crs_1,temp2)*cos11+dot(crs_1,temp3)*prod12-dot(crs_2,temp3)*cos22;
     force = 0.0;
     cross_product(temp2,crs_2,m_1);
     force += dot(m_1,array_bc)*prod123-dot(crs_1,temp2)*sin11;
     cross_product(crs_1,temp3,m_1);
     force += dot(m_1,array_bc)*prod123-dot(crs_2,temp3)*sin22+dot(crs_12,temp1)*prod123-dot(array_bc,temp1)*sin33;
     dev_tri_local[1][1][1] = force;

/*-------------------------------------------------------------------------*/
/* force of z_b                                                            */

     temp1[0] = 0.0;//temp1 = a_bc_d
     temp1[1] = 0.0;
     temp1[2] = -1.0;

     temp2[0] = -array_bc[1]-array_ab[1];//temp2 = crs_1_d
     temp2[1] = array_bc[0]+array_ab[0];
     temp2[2] = 0.0;

     temp3[0] = array_cd[1];//temp3 = crs_2_d
     temp3[1] = -array_cd[0];
     temp3[2] = 0.0;

     dev_tri_local[1][2][0] = dot(crs_2,temp2)*prod12-dot(crs_1,temp2)*cos11+dot(crs_1,temp3)*prod12-dot(crs_2,temp3)*cos22;
     force = 0.0;
     cross_product(temp2,crs_2,m_1);
     force += dot(m_1,array_bc)*prod123-dot(crs_1,temp2)*sin11;
     cross_product(crs_1,temp3,m_1);
     force += dot(m_1,array_bc)*prod123-dot(crs_2,temp3)*sin22+dot(crs_12,temp1)*prod123-dot(array_bc,temp1)*sin33;
     dev_tri_local[1][2][1] = force;

/*-------------------------------------------------------------------------*/
/* force of x_c                                                            */

     temp1[0] = 1.0;//temp1 = a_bc_d
     temp1[1] = 0.0;
     temp1[2] = 0.0;

     temp2[0] = 0.0;//temp2 = crs_1_d
     temp2[1] = array_ab[2];
     temp2[2] = -array_ab[1];

     temp3[0] = 0.0;////temp3 = crs_2_d
     temp3[1] = -array_bc[2]-array_cd[2];
     temp3[2] = array_bc[1]+array_cd[1];

     dev_tri_local[2][0][0] = dot(crs_2,temp2)*prod12-dot(crs_1,temp2)*cos11+dot(crs_1,temp3)*prod12-dot(crs_2,temp3)*cos22;
     force = 0.0;
     cross_product(temp2,crs_2,m_1);
     force += dot(m_1,array_bc)*prod123-dot(crs_1,temp2)*sin11;
     cross_product(crs_1,temp3,m_1);
     force += dot(m_1,array_bc)*prod123-dot(crs_2,temp3)*sin22+dot(crs_12,temp1)*prod123-dot(array_bc,temp1)*sin33;
     dev_tri_local[2][0][1] = force;

/*-------------------------------------------------------------------------*/
/* force of y_c                                                            */

     temp1[0] = 0.0;//temp1 = a_bc_d
     temp1[1] = 1.0;
     temp1[2] = 0.0;

     temp2[0] = -array_ab[2];//temp2 = crs_1_d
     temp2[1] = 0.0;
     temp2[2] = array_ab[0];

     temp3[0] = array_bc[2]+array_cd[2];//temp3 = crs_2_d
     temp3[1] = 0.0;
     temp3[2] = -array_bc[0]-array_cd[0];

     dev_tri_local[2][1][0] = dot(crs_2,temp2)*prod12-dot(crs_1,temp2)*cos11+dot(crs_1,temp3)*prod12-dot(crs_2,temp3)*cos22;
     force = 0.0;
     cross_product(temp2,crs_2,m_1);
     force += dot(m_1,array_bc)*prod123-dot(crs_1,temp2)*sin11;
     cross_product(crs_1,temp3,m_1);
     force += dot(m_1,array_bc)*prod123-dot(crs_2,temp3)*sin22+dot(crs_12,temp1)*prod123-dot(array_bc,temp1)*sin33;
     dev_tri_local[2][1][1] = force;

/*-------------------------------------------------------------------------*/
/* force of z_c                                                            */

     temp1[0] = 0.0;//temp1 = a_bc_d
     temp1[1] = 0.0;
     temp1[2] = 1.0;

     temp2[0] = array_ab[1];//temp2 = crs_1_d
     temp2[1] = -array_ab[0];
     temp2[2] = 0.0;

     temp3[0] = -array_bc[1]-array_cd[1];//temp3 = crs_2_d
     temp3[1] = array_bc[0]+array_cd[0];
     temp3[2] = 0.0;

     dev_tri_local[2][2][0] = dot(crs_2,temp2)*prod12-dot(crs_1,temp2)*cos11+dot(crs_1,temp3)*prod12-dot(crs_2,temp3)*cos22;
     force = 0.0;
     cross_product(temp2,crs_2,m_1);
     force += dot(m_1,array_bc)*prod123-dot(crs_1,temp2)*sin11;
     cross_product(crs_1,temp3,m_1);
     force += dot(m_1,array_bc)*prod123-dot(crs_2,temp3)*sin22+dot(crs_12,temp1)*prod123-dot(array_bc,temp1)*sin33;
     dev_tri_local[2][2][1] = force;

/*-------------------------------------------------------------------------*/
/* force of x_d                                                            */

     temp1[0] = 0.0;//temp1 = crs_2_d
     temp1[1] = array_bc[2];
     temp1[2] = -array_bc[1];

     dev_tri_local[3][0][0] = dot(crs_1,temp1)*prod12-dot(crs_2,temp1)*cos22;
     cross_product(crs_1,temp1,m_1);
     dev_tri_local[3][0][1] = dot(m_1,array_bc)*prod123-dot(crs_2,temp1)*sin22;

/*-------------------------------------------------------------------------*/
/* force of y_d                                                            */

     temp1[0] = -array_bc[2];//temp1 = crs_2_d
     temp1[1] = 0.0;
     temp1[2] = array_bc[0];

     dev_tri_local[3][1][0] = dot(crs_1,temp1)*prod12-dot(crs_2,temp1)*cos22;
     cross_product(crs_1,temp1,m_1);
     dev_tri_local[3][1][1] = dot(m_1,array_bc)*prod123-dot(crs_2,temp1)*sin22;

/*-------------------------------------------------------------------------*/
/* force of z_d                                                            */
     temp1[0] = array_bc[1];//temp1 = crs_2_d
     temp1[1] = -array_bc[0];
     temp1[2] = 0.0;

     dev_tri_local[3][2][0] = dot(crs_1,temp1)*prod12-dot(crs_2,temp1)*cos22;
     cross_product(crs_1,temp1,m_1);
     dev_tri_local[3][2][1] = dot(m_1,array_bc)*prod123-dot(crs_2,temp1)*sin22;

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/




    
    


 
