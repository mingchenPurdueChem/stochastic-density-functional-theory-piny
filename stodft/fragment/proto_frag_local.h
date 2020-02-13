/*-----------------------------------------------------------------*/
/* control-cp-min-frag.c */

void controlCpMinFrag(CLASS *,BONDED *,GENERAL_DATA *,CP *,ANALYSIS *);
/*-----------------------------------------------------------------*/
/* init-frag.c */
void initFragMol(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void initFragUnitCell(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void shiftSystem(int ,int ,int ,int *,int *,int *,double *,int *,
		 double *,int *,int *,int *,
                 int *,double *,double *,double *,CELL *);
void partMolUC(double *,int,int *,COMMUNICATE *,FRAGINFO *,double *);
void mapFragMol(FRAGINFO *,COMMUNICATE *,int,int *,int *);
void mapFragMolHalf(FRAGINFO *,COMMUNICATE *,int,int *,int *,double *);
void reorderMol(FRAGINFO *,int,int,int *,int *,int *,int*,int);
int checkInList(int,int *,int);
int checkIndexList(int,int *,int);

void initFragEnergy(CP *,CLASS *,CLASS *,CP *);
void initFragUnitCell(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);

/*-----------------------------------------------------------------*/
/* proj-wf.c */
void projRhoMiniMol(CP *,GENERAL_DATA *,CLASS *,
                    CP *,GENERAL_DATA *,CLASS *,int);
void projRhoMiniUnitCell(CP *,GENERAL_DATA *,CLASS *,
                    CP *,GENERAL_DATA *,CLASS *,int);
void combineRhoUC(CP *,GENERAL_DATA *,CLASS *,
                 CP *,GENERAL_DATA *,CLASS *,int);
void combineStoUC(CP *,GENERAL_DATA *,CLASS *,
                  CP *,GENERAL_DATA *,CLASS *,int);
void combineStoUCEnergyWindow(CP *,GENERAL_DATA *,CLASS *,
                              CP *,GENERAL_DATA *,CLASS *,int);

/*-----------------------------------------------------------------*/
/* wf-real-frag.c */
void rhoRealCalcDriverFragMol(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *,CP *);
void rhoRealCalcDriverNoise(GENERAL_DATA *,CP *,CLASS *,int);

void rhoRealCalcWrapper(GENERAL_DATA *,CP *,CLASS *,double *,double *,int*,int*,double *,int);
void rhoRealCalcFragWrapper(GENERAL_DATA *,CP *,CLASS *,
                        CP *,double *,double *,int*,int*,double *,double *,int);
void noiseFilterGen(GENERAL_DATA *,CP *,CLASS *,int);
void noiseFilterRealReGen(CLASS *,GENERAL_DATA *,CP *,double *, double *,double *,int);

/*-----------------------------------------------------------------*/
/* energy-cor.c */
void energyCorrect(CP *,GENERAL_DATA *,CLASS *,CP *,CLASS *,int);
void calcKECorMol(CP *,GENERAL_DATA *,CP *,double *);
void calcKEMatrix(CP *,CP *);
void calcVnlCor(CLASS *, CP *,GENERAL_DATA *,CP *,CLASS *,double *,double *,
                double *,double *);
void calcVnlCorEnergyWindow(CLASS *, CP *,GENERAL_DATA *,
                CP *,CLASS *,double *,double *,double *,double *);
void outputFragForce(CP *,CLASS *,CLASS *);
/*-----------------------------------------------------------------*/
/* frag-nlmat.c */
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
                            double *,double *,double *,double *,int *,int *);
/*-----------------------------------------------------------------*/
/* ke-cor-uc.c */
void calcKECorUC(CP *,GENERAL_DATA *,CLASS *,CP *,double *);
void calcKECorUCEnergyWindow(CP *,GENERAL_DATA *,CLASS *,CP *,double *);
void calcKEMatrixUC(GENERAL_DATA *,CP *,CLASS *,CP *, double *);

void rhoRealCalcDriverFragUnitCell(GENERAL_DATA *generalDataMini,CP *cpMini,CLASS *,CP *);
void embedWfReal(GENERAL_DATA *,CP *,CLASS *,CP *,double *,double *);
void qrWrapper(double *,int,int);
void noiseRealReGen(GENERAL_DATA *,CP *,CLASS *,int);
/*-----------------------------------------------------------------*/
/* frag-nlmat-real.c */
void calcRealNonLocalMatrix(CP *, CP *, CLASS *,GENERAL_DATA *);
void calcVnlRealDot(CP *, CLASS *,GENERAL_DATA *,
                    double *,double *,double *,double *,double *,
                    double *,double *,double *,double *,double *,int);
void calcVnlRealDotState(CP *, CLASS *,GENERAL_DATA *,
                    double *,double *,double *,double *,double *,
                    double *,double *,double *,double *);
void calcMatrixFromDot(CP *, CP *,CLASS *,GENERAL_DATA *,double *,double *,double *,
		       double *,double *,double *,double *,
		       double *,double *,double *,double *,double *,int);
