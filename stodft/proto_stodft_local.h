/*-----------------------------------------------------------------*/
/* Local functions for stochastic dft                              */
/*-----------------------------------------------------------------*/
/* coeff.c			                                   */
void genNewtonHermit(STODFTINFO *,STODFTCOEFPOS *);
void genCoeffNewtonHermit(STODFTINFO *,STODFTCOEFPOS *);
void genSampNewtonHermit(STODFTINFO *,STODFTCOEFPOS *);
void genEigenOrb(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
double calcFitError(STODFTINFO *,STODFTCOEFPOS *);

/*-----------------------------------------------------------------*/
/* filter.c                                                        */
void filterNewtonPolyHerm(CP *,CLASS *,GENERAL_DATA *,int);
void genEnergyMax(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void genEnergyMin(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void calcEnergyChemPot(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void calcTotEnergy(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);


/*-----------------------------------------------------------------*/
/* init.c                                                          */

void initStodft(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void reInitWaveFunMin(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void reInitComm(CP *,CPCOEFFS_POS *);
void stoRealloc(CP *,CPCOEFFS_POS *);
void reallocScratch(CP *,int);

/*-----------------------------------------------------------------*/
/* min-CP-stodft.c                                                 */

void scfStodft(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void genStoOrbital(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void genNoiseOrbital(CP *,CPCOEFFS_POS *);

/*-----------------------------------------------------------------*/
/* normh.c                                                         */

void normHNewtonHerm(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,
		     CLATOMS_POS *,double);
void calcCoefForceWrap(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS *,CLATOMS_POS *);
void calcKSPotExtRecipWrap(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS  *,CLATOMS_POS *);
void calcCoefForceExtRecipWrap(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS  *,CLATOMS_POS *);
void calcKSForceControlWrap(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS  *,CLATOMS_POS *);
void calcCoefForceForceControlWrap(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS  *,CLATOMS_POS *);
void calcCoefForceWrapReduce(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS *,CLATOMS_POS *);


/*-----------------------------------------------------------------*/
/* density-init.c                                                  */

void calcRhoInit(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void calcRhoDetInit(CLASS *,BONDED *,GENERAL_DATA *,CP *,CPCOEFFS_POS *);
void calcRhoStoInit(CLASS *,BONDED *,GENERAL_DATA *,CP *,CPCOEFFS_POS *);

/*-----------------------------------------------------------------*/
/* density.c                                                       */

void rhoCalcRealStoHybrid(CPSCR *,CPCOEFFS_INFO *,CELL *,STODFTINFO *,
        double *, double *,double *,int ,int ,int ,int ,int ,COMMUNICATE *,
        PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,PARA_FFT_PKG3D *);

void calcRhoStoRecipHybrid(CPEWALD *,CPSCR *,CPCOEFFS_INFO *,EWALD *,CELL *,
	STODFTINFO *,double *, double *,int ,int ,double *,double *,double *,double *,
        double *,double *,double *, double *,double *,double *,int ,int ,int ,
        int ,int ,int ,COMMUNICATE *,PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,
        PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,PARA_FFT_PKG3D *);

void rhoCalcRealStoFullg(CPSCR *,CPCOEFFS_INFO *,CELL *,STODFTINFO *,double *, 
	double *,double *,int ,int ,int ,int ,int ,COMMUNICATE *,
        PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,PARA_FFT_PKG3D *);

void calcRhoStoRecipFullg(CPEWALD *,CPSCR *,CPCOEFFS_INFO *,EWALD *,CELL *,
                        double *,double *,double *,double *,double *,
                        double *, double *,double *,double *,int ,int ,
                        int ,COMMUNICATE *,PARA_FFT_PKG3D *,PARA_FFT_PKG3D *);


void calcRhoStoHybrid(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);

/*-----------------------------------------------------------------*/
/* calc-chempot.c                                                  */

void calcChemPotInterp(CP *);
void genChemPotInterpPoints(STODFTINFO *,STODFTCOEFPOS *);
double solveLagrangePolyInterp(int ,double *, double *, double, double *,double *);
double calcLagrangeInterpFun(int ,double ,double * ,double * ,double *);
double calcLagrangeInterpDrv(int ,double ,double * ,double * ,double *);

/*-----------------------------------------------------------------*/
/* diis.c	                                                   */

void genDensityMix(CP *,int);
void updateBank(STODFTINFO *,STODFTCOEFPOS *,double *,double **);
void updateErr(STODFTINFO *,STODFTCOEFPOS *,double *,double **,double **);
void calcDensityDiis(CP *,double **,double **);
void matrixInvSVD(double *,double *,double *,int);


