double Complete_Elliptic_Integral_First_Kind(char, double);

void Jacobi_sn_cn_dn(double, char, double, double* , double* , double* );
void Jacobi_cs_ds_ns(double, char, double, double* , double* , double* );
void Jacobi_sc_dc_nc(double, char, double, double* , double* , double* );
void Jacobi_sd_cd_nd(double, char, double, double* , double* , double* );

double Jacobi_am(double, char, double);

void test(double );
void filterRational(CP *, CLASS *, GENERAL_DATA *, int );
void calcChemPotRational(CP *,CLASS *,GENERAL_DATA *, int );
void init_zseed( CP *, double complex *, int, double *, double complex *, double complex *);
