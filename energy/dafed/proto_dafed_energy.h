void force_dafed_final(CLATOMS_INFO *, CLATOMS_POS *);
void force_bias_num(DAFED_INFO *dinfo,DAFED *dafed);
void force_Phi(DAFED_INFO *,DAFED *,CLATOMS_POS *);
inline void cross_product(double *,double *,double *);
inline double dot(double *,double *);
double get_force(double *,double *,double *,double *,double ,double);
double get_phi(double *,double *,double *);
void force_bias(CLATOMS_INFO *,STAT_AVG *,int);
void force_trifun(double [][3][2],double *,double *,double *,double *,
		  double *,double *,double,double,double,double,
		  double,double,double);
double diff_perodic(double,double);

