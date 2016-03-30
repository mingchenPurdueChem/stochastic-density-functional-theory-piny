
void read_coef(CP *,GENERAL_DATA *,CLASS *,FILENAME_PARSE *,CP_PARSE *,
               double *);

void set_coef_NHC(CPOPTS *,CP_COMM_STATE_PKG *,CP_COMM_STATE_PKG *,
		  int ,CPTHERM_INFO *,CPTHERM_POS *,
                  CP_PARSE *,double *,COMMUNICATE *);


void gen_wave(CLASS *,GENERAL_DATA *,CP *,CP_PARSE *,NAME *);

void gen_wave_dvr(CLASS *,GENERAL_DATA *,CP *cp,CP_PARSE *,NAME *);

void mall_properties(CP *);

void mall_wannier(CP *);
