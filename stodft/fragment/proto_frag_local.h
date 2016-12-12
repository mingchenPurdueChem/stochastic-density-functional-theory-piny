void controlCpMinFrag(CLASS *,BONDED *,GENERAL_DATA *,CP *,ANALYSIS *);


void initFragMol(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void initFragUnitCell(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);

void projRhoMini(CP *,GENERAL_DATA *,CLASS *,
                 CP *,GENERAL_DATA *,CLASS *,int);

void rhoRealCalcDriverFrag(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *,CP *);
void rhoRealCalcDriverNoise(GENERAL_DATA *,CP *,CLASS *,int);

void rhoRealCalcWrapper(GENERAL_DATA *,CP *,CLASS *,double *,double *,int*,int*,double *,int);
void rhoRealCalcFragWrapper(GENERAL_DATA *,CP *,CLASS *,
                        CP *,double *,double *,int*,int*,double *,double *,int);


