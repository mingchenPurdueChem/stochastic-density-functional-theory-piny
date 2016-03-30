
void comp_miller_weight(GENERAL_DATA *, double *);

void comp_tensor_z(GENERAL_DATA *, CP *, double *, double ***, double ***);

void comp_tensor_z_prim(double ***,double ***,int ,int , int ,int ,
                           int , int , double *, double *,double **,
                           double **, double **, double **);

void comp_tensor_z_dvr(GENERAL_DATA *, CP *, double ***, double ***);

void lbfg_min(GENERAL_DATA *, CP *, double *,int , double ***, double ***, double *);

void newton(GENERAL_DATA *, CP *, double *, int , double ***, double ***, double *);

void comp_wannier_center(GENERAL_DATA *,CP *,double *, double ***, double ***,
                          double *, double **, double *);

void comp_omega(GENERAL_DATA *,CP *, double *,double ***, double ***, double *,
                          double *, int );

void comp_domega(GENERAL_DATA *,CP *, double *, double *, double ***, double ***, 
                          double *,int );

void comp_domega_zero(int , int , int , double *,double ***, double ***, double *);

void comp_ddomega_zero(int , int , int , double *,double ***, double ***, double *);

void diagonalize_asymm_mat(double **,int ,double *, double **, double **, double **, 
                           double **, int );

void component_wise_matrix_multiply(int ,int ,int ,double **, double **,
                                     double **, double **,double **, double **);


void order_coef(double *,double *,double **,double **, double **, double **,
                CP *,int ***,int ,int ,int );

void write_wannier_orb(GENERAL_DATA *,CP *);

void write_wannier_orb_dvr(GENERAL_DATA *,CP *);
