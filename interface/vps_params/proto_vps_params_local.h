/*---------------------------------------------------------------*/
/* Set_inter_dict.c */

void set_potvps_dict(DICT_WORD *[],int *, int );

/*---------------------------------------------------------------*/

/* control_vps_params.c  */

void set_vps_params(DICT_WORD [],char *,char *, int *,char *,
                    int *,int *,double *,int *, 
                    int *, int *, int *, int *,
                    double *, double *, int *,
                    double *, double *, double *);


void make_vps_splin(char *,int ,int ,
		    int ,int ,double ,
		    double ,double ,double ,
                    double *,double *,double *,double *,
                    double *,double *,double *,double *,double *,
                    double *,double *,double *,int, int, int, int,
                    int ,int, MPI_Comm,int,
                    int ,double,int ,int ,
                    int *,double *,double *);

void make_vps_splin_dvr(char *,int ,int ,int , int *,PSEUDO *,
                        int *, double *,double *, int, double , double , double ,
                        double, double, int , MPI_Comm , int , double , int );


void slow_bess_vps(double [],int ,double ,double [],
               double [],double *,int ,double [],
               double ,
               double ,double ,double ,double ,
               double ,double ,int ,int ,
               double *,double *,int,int ,double);

void slow_invs_bess_vps(double [],int ,double ,double [],
                         double [],int ,double [],int );


void vps_read_error(char *);

void weight_node_gauss_hermite(int ,double *,double * );

void limit_gauss_hermite(int *,double *,double *);

void newt_eval(double *,int ,double *,double *); 

void get_dvl(int ,double *,double *,double *,
             double *, double *, double *,double ,int );

void make_weight_gen(double *,double *,double *, double *,
                     int ,int ,int ,int ,int );

/*---------------------------------------------------------------*/







































