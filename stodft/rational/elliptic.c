////////////////////////////////////////////////////////////////////////////////
// File: complete_elliptic_integral_first_kind.c                              //
// Routine(s):                                                                //
//    Complete_Elliptic_Integral_First_Kind                                   //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>       // required for fabs(), fabsl(), sqrtl(), and M_PI_2
#include <float.h>      // required for LDBL_EPSILON, DBL_MAX

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"

#include "../proto_defs/proto_rational_elliptic.h"

////////////////////////////////////////////////////////////////////////////////
// double Complete_Elliptic_Integral_First_Kind(char arg, double x)           //
//                                                                            //
//  Description:                                                              //
//     The complete elliptic integral of the first kind is the integral from  //
//     0 to pi / 2 of the integrand                                           //
//                   dtheta / sqrt( 1 - k^2 sin^2(theta) ).                   //
//     The parameter k is called the modulus.  This integral is even in k.    //
//     The modulus, k, must satisfy |k| <= 1.  If k = 0 then the integral     //
//     can be readily evaluated.  If |k| = 1, then the integral is infinite.  //
//     Otherwise it must be approximated.                                     //
//                                                                            //
//     In practise the arguments the complete elliptic function of the first  //
//     kind are also given as F(pi/2 \ alpha) or F(pi/2 | m) where the angle  //
//     alpha, called the modular angle, satisfies k = sin(alpha) and the      //
//     argument m = k^2 is simply called the parameter.                       //
//     In terms of these arguments K = F(pi/2 \ alpha) = F(pi/2, sin(alpha))  //
//     and K = F(pi/2 | m) = F(pi/2, sqrt(m)), where                          //
//             K = Complete_Elliptic_Integral_First_Kind( k ).                //
//                                                                            //
//     Let K(k) be the complete elliptic integral of the second kind where    //
//     k is the modulus and  k' = sqrt(1-k^2) is the complementary modulus.   //
//                                                                            //
//     The common mean method, sometimes called the Gauss transform, is a     //
//     variant of the descending Landen transformation in which two sequences //
//     are formed: Setting a[0] = 1 and g[0] = k', the complementary modulus, //
//     a[i] is the arithmetic average and g[i] is the geometric mean of a[i-1]//
//     and g[i-1], i.e. a[i+1] = (a[i] + g[i])/2 and g[i+1] = sqrt(a[i]*g[i]).//
//     The sequences satisfy the inequalities g[0] < g[1] < ... < a[1] < a[0].//
//     Further, lim g[n] = lim a[n].                                          //
//     The value of the complete elliptic integral of the first kind is       //
//     (pi/2) lim (1/G[n]) as n -> inf.                                       //
//                                                                            //
//  Arguments:                                                                //
//     char    arg                                                            //
//                The type of argument of the second argument of F():         //
//                  If arg = 'k', then x = k, the modulus of F(pi/2,k).       //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(pi/2 \ alpha), alpha in radians.          //
//                  If arg = 'm', then x = m, the parameter of F(pi/2 | m).   //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the elliptic function F(pi/2,k),     //
//                F(pi/2 \ alpha) or F(pi/2 | m) corresponding to the value   //
//                of 'arg'.  Note that if arg = 'k', then | x | <= 1 and if   //
//                arg = 'm', then 0 <= x <= 1.                                //
//                                                                            //
//  Return Value:                                                             //
//     The value of the complete elliptic integral of the first kind for the  //
//     given modulus, modular angle, or parameter.  Note that if |k| = 1,     //
//     or m = 1 or a = (+/-) pi/2 then the integral is infinite and DBL_MAX   //
//     is returned.                                                           //
//                                                                            //
//  Example:                                                                  //
//     double K;                                                              //
//     double m, k, a;                                                        //
//                                                                            //
//     ( code to initialize a )                                               //
//                                                                            //
//     k = sin(a);                                                            //
//     m = k * k;                                                             //
//     K = Complete_Elliptic_Integral_First_Kind( 'a', a );                   //
//     printf("K(alpha) = %12.6f where angle(radians) = %12.6f\n",K, a);      //
//     K = Complete_Elliptic_Integral_First_Kind( 'k', k );                   //
//     printf("K(k) = %12.6f where k = %12.6f\n",K, k);                       //
//     K = Complete_Elliptic_Integral_First_Kind( 'm', m );                   //
//     printf("K(m) = %12.6f where m = %12.6f\n",K, m);                       //
////////////////////////////////////////////////////////////////////////////////

static const long double M_PI_2 =  1.5707963267948966192313216916397514L; // pi/2
static const long double PI_2 =  1.5707963267948966192313216916397514L; // pi/2

double Complete_Elliptic_Integral_First_Kind(char arg, double x)
{
   long double k;          // modulus 
   long double m;          // parameter 
   long double a;          // average
   long double g;          // geometric mean
   long double a_old;      // previous average
   long double g_old;      // previous geometric mean

   if ( x == 0.0 ) return M_PI_2;

   switch (arg) {
      case 'k': k = fabsl((long double) x);
                m = k * k;
                break;
      case 'm': m = (long double) x;
                k = sqrtl(fabsl(m));
                break;
      case 'a': k = sinl((long double)x);
                m = k * k;
                break;
      default:  k = fabsl((long double) x);
                m = k * k;
   }

   if ( m == 1.0 ) return DBL_MAX;

   a = 1.0L;
   g = sqrtl(1.0L - m);
   while (1) {
      g_old = g;
      a_old = a;
      a = 0.5L * (g_old + a_old);
      g = sqrtl(g_old * a_old);
      if ( fabsl(a_old - g_old) <= (a_old * LDBL_EPSILON) ) break;
   }
   return (double) (PI_2 / g); 
}




////////////////////////////////////////////////////////////////////////////////
// File: jacobi_elliptic_functions.c                                          //
// Routine(s):                                                                //
//    Jacobi_sn_cn_dn                                                         //
//    Jacobi_cs_ds_ns                                                         //
//    Jacobi_sc_dc_nc                                                         //
//    Jacobi_sd_cd_nd                                                         //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void Jacobi_sn_cn_dn(double u, char arg, double x, double* sn, double* cn, //
//                                                                double* dn) //
//                                                                            //
//  Description:                                                              //
//     This function calculates the Jacobi elliptic functions: sn, cn, and dn //
//     which are defined in terms of the inverse to Legendre's elliptic       //
//     integral of the first kind.  If one sets u to the integral from        //
//     0 to phi, |phi| < pi / 2, of the integrand                             //
//                   dtheta / sqrt( 1 - k^2 sin^2(theta) )                    //
//     then         sn(u,k) = sin(phi),                                       //
//                  cn(u,k) = cos(phi),                                       //
//     and          dn(u,phi) = sqrt( 1 - k^2 sin^2(phi) )                    //
//                                                                            //
//      If the modulus k > 1, then the Jacobi real modulus transforms are     //
//      applied:                                                              //
//                          sn(x,k) = sn(kx,1/k)/k                            //
//                          cn(x,k) = dn(kx,1/k)                              //
//                          dn(x,k) = cn(kx,1/k).                             //
//                                                                            //
//      If the modulus k is purely imaginary, m < 0, then the transforms are  //
//      applied:                                                              //
//                          sn(x,ik) = sn(k'x,k/k')/k'                        //
//                          cn(x,ik) = dn(k'x,k/k')                           //
//                          dn(x,ik) = cn(k'x,k/k').                          //
//  Arguments:                                                                //
//     double  u                                                              //
//                The argument of the Jacobi elliptic functions.              //
//     char    arg                                                            //
//                The type of argument of the second argument of am():        //
//                  If arg = 'k', then x = k, the modulus of F(phi,k).        //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(phi \ alpha), alpha in radians.           //
//                  If arg = 'm', then x = m, the parameter of F(phi | m).    //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the jacobi functions sn(u,x),        //
//                cn(u,x), or dn(u,x) corresponding to the second argument    //
//                of the elliptic integral of the first kind F(phi,x).        //
//                'x' may the the modulus, modular angle, or parameter        //
//                depending on the value of 'arg'.                            //
//     double* sn                                                             //
//                The address of the value the Jacobi elliptic function sn    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* cn                                                             //
//                The address of the value the Jacobi elliptic function cn    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* dn                                                             //
//                The address of the value the Jacobi elliptic function dn    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//                                                                            //
//  Return Value:                                                             //
//     The values of sn, cn, and dn are returned via the argument list.       //
//                                                                            //
//  Example:                                                                  //
//     double u, m, k, a;                                                     //
//     double sn, cn, dn;                                                     //
//                                                                            //
//     ( code to initialize u and a )                                         //
//                                                                            //
//     k = sin(a);                                                            //
//     m = k * k;                                                             //
//     Jacobi_sn_cn_dn( u, 'a', a, &sn, &cn, &dn);                            //
//     printf("sn cn dn %12,6f, %12.6f %12.6f, u %12.6f, angle(radians)       //
//                                                 %12.6f\n",sn,cn,dn,u,a);   //
//     Jacobi_sn_cn_dn( u, 'k', k, &sn, &cn, &dn);                            //
//     printf("sn cn dn %12,6f, %12.6f %12.6f, u %12.6f, modulus              //
//                                                 %12.6f\n",sn,cn,dn,u,k);   //
//     Jacobi_sn_cn_dn( u, 'm', m, &sn, &cn, &dn);                            //
//     printf("sn cn dn %12,6f, %12.6f %12.6f, u %12.6f, parameter            //
//                                                 %12.6f\n",sn,cn,dn,u,m);   //
////////////////////////////////////////////////////////////////////////////////
                
void Jacobi_sn_cn_dn(double u, char arg, double x, double* sn, double* cn,
                                                                    double* dn)
{
   double k;
   double m;
   double phi;
   int n;

   switch (arg) {
      case 'a': k = sin(x);
                m = k * k;
                break;
      case 'm': m = x;
                k = sqrt(fabs(m));
                break;
      default:  k = fabs(x);
                m = k * k;
   }

                   // Check special cases m = 1 or m = 0 //

   if ( m == 1.0 ) {
      *sn = tanh(u);
      *cn = 1.0 / cosh(u);
      *dn = *cn;
      return;
   }
   if ( m == 0.0 ) {
      *sn = sin(u);
      *cn = cos(u);
      *dn = 1.0;
      return;
   }

                  // If m > 1, perform the transformation: //
                  //        sn(u,k) = sn(k*u, 1/k) /k,     //
                  //        cn(u,k) = dn(k*u, 1/k),        //
                  //        dn(u,k) = cn(k*u, 1/k).        //

   if (m > 1.0) {
      Jacobi_sn_cn_dn(k * u, 'k', 1.0 / k, sn, dn, cn);
      *sn /= k;
      return;
   }

                  // If m < 0, perform the transformation: //
//sn(u,m) = [sn(sqrt(1-m)*u,-m/(1-m)) / dn(sqrt(1-m)*u,-m/(1-m))] / sqrt(1-m) //
//cn(u,m) = [cn(sqrt(1-m)*u,-m/(1-m)) / dn(sqrt(1-m)*u,-m/(1-m))]             //
//dn(u,m) = [1.0 / dn(sqrt(1-m)*u,-m/(1-m))]                                  //

   if (m < 0.0) {
      Jacobi_sd_cd_nd(sqrt(1.0 - m) * u, 'm', -m / (1.0 - m), sn, cn, dn);
      *sn /= sqrt(1.0 - m);
      return;
   }
 
   phi = Jacobi_am(u, arg, x);
   *sn = sin(phi);
   *cn = cos(phi);
   *dn = sqrt(1.0 - m * (*sn * *sn));

   return; 
}

////////////////////////////////////////////////////////////////////////////////
// void Jacobi_cs_ds_ns(double u, char arg, double x, double* cs, double* ds,
//                                                                double* ns) //
//                                                                            //
//  Description:                                                              //
//     This function calculates the Jacobi elliptic functions: cs, ds, and    //
//     ns which are defined in terms of the Jacobi elliptic functions sn, cn  //
//     and dn (see above)                                                     //
//                        cs(u,k) = cn(u,k) / sn(u,k)                         //
//                        ds(u,k) = dn(u,k) / sn(u,k)                         //
//     and                ns(u,k) = 1 / sn(u,k).                              //
//                                                                            //
//  Arguments:                                                                //
//     double  u                                                              //
//                The argument of the Jacobi elliptic functions.              //
//     char    arg                                                            //
//                The type of argument of the second argument of a Jacobi     //
//                elliptic function.                                          //
//                  If arg = 'k', then x = k, the modulus of F(phi,k).        //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(phi \ alpha), alpha in radians.           //
//                  If arg = 'm', then x = m, the parameter of F(phi | m).    //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the jacobi functions cs(u,x),        //
//                ds(u,x), or ns(u,x) corresponding to the second argument    //
//                of the elliptic integral of the first kind F(phi,x).        //
//                'x' may the the modulus, modular angle, or parameter        //
//                depending on the value of 'arg'.                            //
//     double* cs                                                             //
//                The address of the value the Jacobi elliptic function cs    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* ds                                                             //
//                The address of the value the Jacobi elliptic function ds    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* ns                                                             //
//                The address of the value the Jacobi elliptic function ns    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//                                                                            //
//  Return Value:                                                             //
//     The values of cs, ds, and ns are returned via the argument list.       //
//                                                                            //
//  Example:                                                                  //
//     double u, m, k, a;                                                     //
//     double cs, ds, ns;                                                     //
//                                                                            //
//     ( code to initialize u, and a )                                        //
//                                                                            //
//     k = sin(a);                                                            //
//     m = k * k;                                                             //
//     Jacobi_cs_ds_ns( u, 'a', a, &cs, &ds, &ns);                            //
//     printf("cs ds ns %12,6f, %12.6f %12.6f, u %12.6f, angle(radians)       //
//                                                 %12.6f\n",cs,ds,ns,u,a);   //
//     Jacobi_cs_ds_ns( u, 'k', k, &cs, &ds, &ns);                            //
//     printf("cs ds ns %12,6f, %12.6f %12.6f, u %12.6f, modulus              //
//                                                 %12.6f\n",cs,ds,ns,u,k);   //
//     Jacobi_cs_ds_ns( u, 'm', m, &cs, &ds, &ns);                            //
//     printf("cs ds ns %12,6f, %12.6f %12.6f, u %12.6f, parameter            //
//                                                 %12.6f\n",cs,ds,ns,u,m);   //
////////////////////////////////////////////////////////////////////////////////
                
void Jacobi_cs_ds_ns(double u, char arg, double x, double* cs, double* ds,
                                                                    double* ns)
{
   double sn, cn, dn;

   Jacobi_sn_cn_dn( u, arg, x, &sn, &cn, &dn);

   if (sn == 0.0) {
      *cs = DBL_MAX;
      *ds = DBL_MAX;
      *ns = DBL_MAX;
   } else {
      *cs = cn / sn;
      *ds = dn / sn;
      *ns = 1.0 / sn;
   }

   return; 
}

////////////////////////////////////////////////////////////////////////////////
// void Jacobi_sc_dc_nc(double u, char arg, double x, double* sc, double* dc, //
//                                                                double* nc) //
//                                                                            //
//  Description:                                                              //
//     This function calculates the Jacobi elliptic functions: sc, dc, and    //
//     nc which are defined in terms of the Jacobi elliptic functions sn, cn  //
//     and dn (see above)                                                     //
//                        sc(u,k) = sn(u,k) / cn(u,k)                         //
//                        dc(u,k) = dn(u,k) / cn(u,k)                         //
//     and                nc(u,k) = 1 / cn(u,k).                              //
//                                                                            //
//  Arguments:                                                                //
//     double  u                                                              //
//                The argument of the Jacobi elliptic functions.              //
//     char    arg                                                            //
//                The type of argument of the second argument of a Jacobi     //
//                elliptic function.                                          //
//                  If arg = 'k', then x = k, the modulus of F(phi,k).        //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(phi \ alpha), alpha in radians.           //
//                  If arg = 'm', then x = m, the parameter of F(phi | m).    //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the amplitude function am(u,x)       //
//                corresponding to the second argument of the elliptic        //
//                integral of the first kind F(phi,x).  'x' may the the       //
//                modulus, modular angle, or parameter depending on the value //
//                of 'arg'.                                                   //
//     double* sc                                                             //
//                The address of the value the Jacobi elliptic function sc    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* dc                                                             //
//                The address of the value the Jacobi elliptic function dc    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* nc                                                             //
//                The address of the value the Jacobi elliptic function nc    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//                                                                            //
//  Return Value:                                                             //
//     The values of sc, dc, and nc are returned via the argument list.       //
//                                                                            //
//  Example:                                                                  //
//     double u, m, k, a;                                                     //
//     double sc, dc, nc;                                                     //
//                                                                            //
//     ( code to initialize u, and a )                                        //
//                                                                            //
//     k = sin(a);                                                            //
//     m = k * k;                                                             //
//     Jacobi_sc_dc_nc( u, 'a', a, &sc, &dc, &nc);                            //
//     printf("sc dc nc %12,6f, %12.6f %12.6f, u %12.6f, angle(radians)       //
//                                                 %12.6f\n",sc,dc,nc,u,a);   //
//     Jacobi_sc_dc_nc( u, 'k', k, &sc, &dc, &nc);                            //
//     printf("sc dc nc %12,6f, %12.6f %12.6f, u %12.6f, modulus              //
//                                                 %12.6f\n",sc,dc,nc,u,k);   //
//     Jacobi_sc_dc_nc( u, 'm', m, &sc, &dc, &nc);                            //
//     printf("sc dc nc %12,6f, %12.6f %12.6f, u %12.6f, parameter            //
//                                                 %12.6f\n",sc,dc,nc,u,m);   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
                
void Jacobi_sc_dc_nc(double u, char arg, double x, double* sc,
                                                      double* dc, double* nc)
{
   double sn, cn, dn;

   Jacobi_sn_cn_dn( u, arg, x, &sn, &cn, &dn);

   if (cn == 0.0) {
      *sc = DBL_MAX;
      *dc = DBL_MAX;
      *nc = DBL_MAX;
   } else {
      *sc = sn / cn;
      *dc = dn / cn;
      *nc = 1.0 / cn;
   }

   return; 
}

////////////////////////////////////////////////////////////////////////////////
// void Jacobi_sd_cd_nd(double u, char arg, double x, double* sd, double* cd, //
//                                                                double* nd) //
//                                                                            //
//  Description:                                                              //
//     This function calculates the Jacobi elliptic functions: sd, cd, and    //
//     nd which are defined in terms of the Jacobi elliptic functions sn, cn  //
//     and dn (see above)                                                     //
//                        sd(u,k) = sn(u,k) / dn(u,k)                         //
//                        cd(u,k) = cn(u,k) / dn(u,k)                         //
//     and                nd(u,k) = 1 / dn(u,k).                              //
//                                                                            //
//  Arguments:                                                                //
//     double  u                                                              //
//                The argument of the Jacobi's elliptic functions.            //
//     char    arg                                                            //
//                The type of argument of the second argument of a Jacobi     //
//                elliptic function.                                          //
//                  If arg = 'k', then x = k, the modulus of F(phi,k).        //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(phi \ alpha), alpha in radians.           //
//                  If arg = 'm', then x = m, the parameter of F(phi | m).    //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the amplitude function am(u,x)       //
//                corresponding to the second argument of the elliptic        //
//                integral of the first kind F(phi,x).  'x' may the the       //
//                modulus, modular angle, or parameter depending on the value //
//                of 'arg'.                                                   //
//     double* sd                                                             //
//                The address of the value the Jacobi's elliptic function sd  //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* sc                                                             //
//                The address of the value the Jacobi's elliptic function cd  //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* nd                                                             //
//                The address of the value the Jacobi's elliptic function nd  //
//                evaluated at u with modulus | modular angle | parameter x.  //
//                                                                            //
//  Return Value:                                                             //
//     The values of sd, cd, and nd are returned via the argument list.       //
//                                                                            //
//  Example:                                                                  //
//     double u, m, k, a;                                                     //
//     double sd, cd, nd;                                                     //
//                                                                            //
//     ( code to initialize u, and a )                                        //
//                                                                            //
//     k = sin(a);                                                            //
//     m = k * k;                                                             //
//     Jacobi_sd_cd_nd( u, 'a', a, &sd, &cd, &nd);                            //
//     printf("sd cd nd %12,6f, %12.6f %12.6f, u %12.6f, angle(radians)       //
//                                                 %12.6f\n",sd,cd,nd,u,a);   //
//     Jacobi_sd_cd_nd( u, 'k', k, &sd, &cd, &nd);                            //
//     printf("sd cd nd %12,6f, %12.6f %12.6f, u %12.6f, modulus              //
//                                                 %12.6f\n",sd,cd,nd,u,k);   //
//     Jacobi_sd_cd_nd( u, 'm', m, &sd, &cd, &nd);                            //
//     printf("sd cd nd %12,6f, %12.6f %12.6f, u %12.6f, parameter            //
//                                                 %12.6f\n",sd,cd,nd,u,m);   //
////////////////////////////////////////////////////////////////////////////////
                
void Jacobi_sd_cd_nd(double u, char arg, double x, double* sd,
                                                      double* cd, double* nd)
{
   double sn, cn, dn;

   Jacobi_sn_cn_dn( u, arg, x, &sn, &cn, &dn);

   if (dn == 0.0) {
      *sd = DBL_MAX;
      *cd = DBL_MAX;
      *nd = DBL_MAX;
   } else {
      *sd = sn / dn;
      *cd = cn / dn;
      *nd = 1.0 / dn;
   }

   return; 
}




////////////////////////////////////////////////////////////////////////////////
// File: jacobi_am.c                                                          //
// Routine(s):                                                                //
//    Jacobi_am                                                               //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// double Jacobi_am(double u, char arg, double x)                             //
//                                                                            //
//  Description:                                                              //
//     Let F(phi,k) = F(phi \ alpha) = F(phi | m) be Legendre's elliptic      //
//     function of the first kind with modulus k, modular angle alpha where   //
//     k = sin(alpha) or parameter m where m = k^2, i.e.                      //
//        F(phi,k) = Integral(0,phi) dtheta / sqrt(1 - k^2 sin^2(theta))      //
//        F(phi \ alpha) = Integral(0,phi) dtheta /                           //
//                                        sqrt(1 - sin^2(alpha) sin^2(theta)) //
//        F(phi | m) = Integral(0,phi) dtheta / sqrt(1 - m sin^2(theta))      //
//                                                                            //
//     This Jacobi elliptic amplitude function, am, is defined as             //
//               am(u,k) = am(u \ alpha) = am(u | m)  = phi                   //
//     where u = F(phi,k) = F(phi \ alpha) = F(phi | m).                      //
//                                                                            //
//     The common mean method, sometimes called the Gauss transform method,   //
//     is a variant of the descending Landen transformation in which two      //
//     sequences are formed: Setting a[0] = 1 and g[0] = 1-m, a[i] is the     //
//     arithmetic average and g[i] is the geometric mean of a[i-1] and g[i-1],//
//     i.e. a[i+1] = (a[i] + g[i])/2 and g[i+1] = sqrt(a[i]*g[i]).  The       //
//     sequences, a[i] and g[i], satisfy the inequalities                     //
//     g[0] < g[1] < ... < a[1] < a[0].  Further, lim g[n] = lim a[n].        //
//                                                                            //
//     Set phi[n] = 2^n a[n] u, the recursively compute phi[n-1] by           //
//        phi[n-1] = [ phi[n] + arcsin( c[n] sin(phi[n]) / a[n] ] / 2         //
//     for until n = 1.  Then am(u,k) = am(u \ alpha) = am(u | m) = phi[0].   //
//                                                                            //
//  Arguments:                                                                //
//     double  u                                                              //
//                The first argument of am(u,x) corresponding to the value of //
//                the elliptic integral of the first kind u = F(am(u,x),x).   //
//     char    arg                                                            //
//                The type of argument of the second argument of am():        //
//                  If arg = 'k', then x = k, the modulus of F(phi,k).        //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(phi \ alpha), alpha in radians.           //
//                  If arg = 'm', then x = m, the parameter of F(phi | m).    //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the amplitude function am(u,x)       //
//                corresponding to the second argument of the elliptic        //
//                integral of the first kind F(phi,x).  'x' may the the       //
//                modulus, modular angle, or parameter depending on the value //
//                of 'arg'.  If 'arg' = 'm', then x must be between 0 and 1   //
//                inclusively and if 'arg' = 'k', then x must be between -1   //
//                and 1 inclusively.                                          //
//                                                                            //
//  Return Value:                                                             //
//     The amplitude am(u,m) in radians.                                      //
//                                                                            //
//  Example:                                                                  //
//     double u, x;                                                           //
//     double am;                                                             //
//     char   arg;                                                            //
//                                                                            //
//     ( code to initialize u, arg, and x )                                   //
//                                                                            //
//     phi = Jacobi_am( u, arg, x );                                          //
////////////////////////////////////////////////////////////////////////////////

#define N 30                // More than sufficient for extended precision
                            // Near m = 1, usually an N of 10 would do.
                
double Jacobi_am(double u, char arg,  double x)
{
   long double a[N+1];
   long double g[N+1];
   long double c[N+1];
   long double two_n;
   long double phi;
   long double k;
   int n;

                        // Check special case x = 0 //
                        // i.e. k = m = alpha = 0.  //

   if ( x == 0.0 ) return u;

   switch (arg) {
      case 'a': k = sinl( fabsl((long double) x) ); break;
      case 'm': k = sqrtl( fabsl((long double) x) ); break;
      default:  k = (long double) fabs(x);
   }

                   // Check special case k = 1 //

   if ( k == 1.0 ) return 2.0 * atan( exp(u) ) - M_PI_2;

         // If k > 1, then perform a Jacobi modulus transformation. //
         // Initialize the sequence of arithmetic and geometric     //
         // means, a = 1, g = k'.                                   //

   a[0] = 1.0L;
   g[0] = sqrtl(1.0L - k * k);
   c[0] = k;
   
   // Perform the sequence of Gaussian transformations of arithmetic and //
   // geometric means of successive arithmetic and geometric means until //
   // the two means converge to a common mean (upto machine accuracy)    //
   // starting with a = 1 and g = k', which were set above.              //
   
   two_n = 1.0L; 
   for (n = 0; n < N; n++) {
      if ( fabsl(a[n] - g[n]) < (a[n] * LDBL_EPSILON) ) break;
      two_n += two_n;
      a[n+1] = 0.5L * (a[n] + g[n]);
      g[n+1] = sqrtl(a[n] * g[n]);
      c[n+1] = 0.5L * (a[n] - g[n]);
   }

         // Prepare for the inverse transformation of phi = x * cm. //

   phi = two_n * a[n] * u;

                      // Perform backward substitution //

   for (; n > 0; n--) phi = 0.5L * ( phi + asinl( c[n] * sinl(phi) / a[n]) );

   return (double) phi; 
}
