double Complete_Elliptic_Integral_First_Kind(char, double);

void Jacobi_sn_cn_dn(double, char, double, double* , double* , double* );
void Jacobi_cs_ds_ns(double, char, double, double* , double* , double* );
void Jacobi_sc_dc_nc(double, char, double, double* , double* , double* );
void Jacobi_sd_cd_nd(double, char, double, double* , double* , double* );

double Jacobi_am(double, char, double);

void filterRational(CP *, CLASS *, GENERAL_DATA *, KOMEGAINFO *, int );
void calcChemPotRational(CP *,CLASS *,GENERAL_DATA *, KOMEGAINFO *, int );

double complex fermi_fun(double complex , double , double , double );
void solve_shifted_eqn_cocg( CP *, CLASS *, GENERAL_DATA *, KOMEGAINFO *, int , int );
void init_zseed( CP *, int );
double calcNumberElecRational(CP *, double);
void applyFilterRational(CP *, CLASS *, GENERAL_DATA *, double, int );

void komega_COCG_init(KOMEGAINFO *, int, int, int, double complex*, double complex*, int, double);
void komega_COCG_finalize(KOMEGAINFO *);
void komega_COCG_update(KOMEGAINFO *, double complex*, double complex*, double complex*, double complex*, int*, int);
void komega_COCG_shiftedeqn(KOMEGAINFO *, double complex *, double complex *, int);
void komega_COCG_seed_switch(KOMEGAINFO *, double complex *, int *);


double complex fermi_fun_g(double complex , double , double , double );
void solve_shifted_eqn_cocg_g( CP *, CLASS *, GENERAL_DATA *, KOMEGAINFO *, int , int );
void filterRational_g(CP *,CLASS *,GENERAL_DATA *, KOMEGAINFO *, int );
void applyFilterRational_g(CP *, CLASS *, GENERAL_DATA *, double , int );
double calcNumberElecRational_g(CP *, double );
void calcChemPotRational_g(CP *,CLASS *,GENERAL_DATA *, KOMEGAINFO *, int );
void init_zseed_g( CP *, int );

void komega_COCG_init_g(KOMEGAINFO *, int, int, int,
                       double *, double *, double *, double *,
                       double complex *, int, double );
void komega_COCG_finalize_g(KOMEGAINFO *);
void komega_COCG_update_g(KOMEGAINFO *,
                       double *, double *, double *, double *,
                       double *, double *, double *, double *,
                       double *, double *, double *, double *,
                       double *, double *, double *, double *,
                       int *, int );
void komega_COCG_shiftedeqn_g(KOMEGAINFO *,
                              double *, double *, double *, double *,
                              double *, double *, double *, double *,
                              int );
void komega_COCG_seed_switch_g(KOMEGAINFO *,
                               double *, double *, double *, double *,
                               int *);

