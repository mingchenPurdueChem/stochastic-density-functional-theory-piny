void mall_coef(CP *,SIMOPTS *,int);

void assign_coef(CPCOEFFS_INFO *,CPCOEFFS_POS *,double *,double *,double *,
                 double *,int, MPI_Comm, int, int, int,COMMUNICATE *,int,int);

void assign_nhc_coef_vel(CPTHERM_INFO *,CPTHERM_POS *,double *,int,
                         MPI_Comm,int,int );

void bess_trans(double *,int ,double ,double *,
             double *,int ,double *,int ,double *);

void get_gpsi(double ,int ,
              double **,double **,double **,double **,
              double *,double ,double ,int );

void  splin_btrans(int *,double **,double **,
                    double **,double **,
                    double ,double ,double *,
                    double *,int *,char *,int ,
                    int ,int ,MPI_Comm  );

void  fit_spline(double *,double *,double *,
                 double *,double *,int );

void read_wan_cent(CP *,GENERAL_DATA *);


void splin_pseudo_wave(int ,double **,double **, double **, double **,int *,
                       double *, char *, int , int , MPI_Comm);

void get_rpsi(double ,double ,double , double *, double *, double *, int , int ,
              double **, double **,double **, double **, double *, double ,
              double , int , CELL *);
void read_coef_sys_info(CP *, FILE *,int , char *);
void read_coef_alloc_init(CP *, int ,double *);
void read_coef_alloc_init_dvr(CP *, int, double *);
void read_coef_fetch_occs(CP *, FILE *, char *);
void read_coef_fetch_coefs(CP *, FILE *, char *,int);
void read_coef_fetch_coefs_dvr(CP *, FILE *, char *,int);
void read_coef_fetch_vcoefs(CP *, FILE *, char *);
void read_coef_fetch_vcoefs_dvr(CP *, FILE *, char *);
void read_coef_fetch_vnhc(CP *, FILE *, char *);
void read_coef_init_nhc(CP *);
void read_coef_transpose(CP *,int ,int, int);

