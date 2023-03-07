double Complete_Elliptic_Integral_First_Kind(char, double);

void Jacobi_sn_cn_dn(double, char, double, double* , double* , double* );
void Jacobi_cs_ds_ns(double, char, double, double* , double* , double* );
void Jacobi_sc_dc_nc(double, char, double, double* , double* , double* );
void Jacobi_sd_cd_nd(double, char, double, double* , double* , double* );

double Jacobi_am(double, char, double);

void filterRational(CP *, CLASS *, GENERAL_DATA *, int );
void calcChemPotRational(CP *,CLASS *,GENERAL_DATA *, int );

double complex fermi_fun(double complex , double , double , double );
void solve_shifted_eqn_cocg( CP *, CLASS *, GENERAL_DATA *, int , int );
void init_zseed( CP *, int );

void komega_COCG_init(KOMEGAINFO *, int, int, int, double complex*, double complex*, int, double);
void komega_COCG_finalize(KOMEGAINFO *);
void komega_COCG_update(KOMEGAINFO *, double complex*, double complex*, double complex*, double complex*, int*);
void komega_COCG_shiftedeqn(KOMEGAINFO *, double complex *, double complex *);
void komega_COCG_seed_switch(KOMEGAINFO *, double complex *, int *);

