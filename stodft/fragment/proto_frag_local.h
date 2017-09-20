void controlCpMinFrag(CLASS *,BONDED *,GENERAL_DATA *,CP *,ANALYSIS *);


void initFragMol(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void initFragUnitCell(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void shiftSystem(int ,int ,int ,int *,int *,int *,double *,int *,double *,int *,int *,int *,
                 int *,double *,double *,double *,CELL *);
void partMolUC(double *,int,int *,COMMUNICATE *,FRAGINFO *);
void mapFragMol(FRAGINFO *,COMMUNICATE *,int,int *,int *);
void reorderMol(FRAGINFO *,int,int,int *,int *,int *,int);
int checkInList(int,int *,int);

void initFragEnergy(CP *,CLASS *,CLASS *,CP *);
void initFragUnitCell(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);

void projRhoMini(CP *,GENERAL_DATA *,CLASS *,
                 CP *,GENERAL_DATA *,CLASS *,int);

void rhoRealCalcDriverFragMol(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *,CP *);
void rhoRealCalcDriverNoise(GENERAL_DATA *,CP *,CLASS *,int);

void rhoRealCalcWrapper(GENERAL_DATA *,CP *,CLASS *,double *,double *,int*,int*,double *,int);
void rhoRealCalcFragWrapper(GENERAL_DATA *,CP *,CLASS *,
                        CP *,double *,double *,int*,int*,double *,double *,int);

void energyCorrect(CP *,GENERAL_DATA *,CLASS *,CP *,CLASS *,int);
void calcKECor(CP *,GENERAL_DATA *,CP *,double *);
void calcKEMatrix(CP *,CP *);
void calcVnlCor(CLASS *, CP *,GENERAL_DATA *,CP *,CLASS *,double *,double *,
                double *,double *);
void calcNonLocalMatrix(CP *,CP *,CLASS *, GENERAL_DATA *);
void getnlPotPvFatmFrag(CLATOMS_INFO *clatoms_info,CLATOMS_POS *,CELL *,CPCOEFFS_INFO *,CPSCR *,
                        EWD_SCR *,CPOPTS *,PSEUDO *,ATOMMAPS *,FRAGINFO *,double *,int,double *);
void sumnlPotPvFatmHessFrag(int,int,int,int,int,int,int,int,int,int *,double *,
                            double *,double *,double *,double *,double *,double *,
                            double *,double *,double *,double *,double *,double *,
                            double *,double *,double *,double *,double *,double *,
                            double *,double *,double *,double *,double *,
                            double *,double *,double *,double *,double *,double *,
                            double *,double *,double *,int,int,double *, double *,
                            double *,double *,double *,double *);
void rhoRealCalcDriverFragUnitCell(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *,CP *);
void embedWfReal(GENERAL_DATA *,CP *,CLASS *,CP *,double *,double *);
void qrWrapper(double *,int,int);
void noiseRealReGen(GENERAL_DATA *,CP *,CLASS *,int);



