/*-----------------------------------------------------------------*/
/* Local functions for stochastic dft                              */
/*-----------------------------------------------------------------*/
/* coeff.c			                                   */
void genNewtonHermit(STODFTINFO *,STODFTCOEFPOS *);
void genNewtonHermitTrueChemPot(STODFTINFO *,STODFTCOEFPOS *);
void genCoeffNewtonHermit(STODFTINFO *,STODFTCOEFPOS *);
void genSampNewtonHermit(STODFTINFO *,STODFTCOEFPOS *);
void genEigenOrb(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
double calcFitErrorNewton(STODFTINFO *,STODFTCOEFPOS *);
void calcChebyCoeff(STODFTINFO *,STODFTCOEFPOS *,double,double *);
void testChebyCoeff(STODFTINFO *,STODFTCOEFPOS *,double,double *);
/*-----------------------------------------------------------------*/
/* filter.c                                                        */
void filterNewtonPolyHerm(CP *,CLASS *,GENERAL_DATA *,int);
/*-----------------------------------------------------------------*/
/* calc-spectral-range.c                                           */
void genEnergyMax(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void genEnergyMin(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
/*-----------------------------------------------------------------*/
/* calc-energy.c						   */
void calcEnergyChemPot(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void calcTotEnergy(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void calcKNEEnergyFilterDiag(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void calcTotEnergyFilterDiag(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
/*-----------------------------------------------------------------*/
/* init.c                                                          */
void initStodft(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void reInitWaveFunMin(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void reInitComm(CP *,CPCOEFFS_POS *);
void stoRealloc(CP *,CPCOEFFS_POS *);
void reallocScratch(CP *,int);
void initFilterDiag(CP *);
/*-----------------------------------------------------------------*/
/* min-CP-stodft.c                                                 */
void scfStodftInterp(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void scfStodftCheby(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void scfStodftFilterDiag(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
/*-----------------------------------------------------------------*/
/* gen-stodft-wf.c						   */
void genStoOrbitalInterp(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void genStoOrbitalCheby(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
/*-----------------------------------------------------------------*/
/* gen-noise.c							   */
void genNoiseOrbital(CP *,CPCOEFFS_POS *);
/*-----------------------------------------------------------------*/
/* normh.c                                                         */
void normHNewtonHerm(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,
		     CLATOMS_POS *,double);
void normHCheby(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *,int);
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
void readRho(CLASS *,BONDED *,GENERAL_DATA *,CP *,CPCOEFFS_POS *);
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

void calcRhoStoHybridInterp(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void calcRhoStoHybridCheby(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void outputDensity(CP *,CELL *);
void calcRhoFilterDiagHybrid(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
/*-----------------------------------------------------------------*/
/* calc-chempot.c                                                  */
void calcChemPotInterp(CP *);
void genChemPotInterpPoints(STODFTINFO *,STODFTCOEFPOS *);
double solveLagrangePolyInterp(int ,double *, double *, double, double *,double *);
double calcLagrangeInterpFun(int ,double ,double * ,double * ,double *);
double calcLagrangeInterpDrv(int ,double ,double * ,double * ,double *);
void updateChemPot(STODFTINFO *,STODFTCOEFPOS *);
void adjChemPot(STODFTINFO *,STODFTCOEFPOS *);
/*-----------------------------------------------------------------*/
/* calc-chempot-chebyshev.c                                        */
void calcChemPotCheby(CP *,CLASS *,GENERAL_DATA *,int);
void calcChebyMoments(CP *,CLASS *,GENERAL_DATA *,int);
double calcNumElecCheby(CP *,double,double *);
/*-----------------------------------------------------------------*/
/* diis.c	                                                   */
void genDensityMix(CP *,int);
void updateBank(STODFTINFO *,STODFTCOEFPOS *,double *,double **);
void updateErr(STODFTINFO *,STODFTCOEFPOS *,double *,double *,double **);
void calcDensityDiis(CP *,double **,double **);
void matrixInvSVD(double *,double *,double *,int);
/*-----------------------------------------------------------------*/
/* filter-diag.c                                                   */
void orthDiagDriver(CP *,CLASS *,GENERAL_DATA *,int);
void orthNormStoWf(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
int orthogSVD(int ,int ,double *);
void buildKSMatrix(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void diagKSMatrix(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void diagSymMatWrapper(int ,double *,double *);
void genMatrixMulWrapper(int ,int ,double *,double *,double *);
void calcForceWrapper(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *,
		      double *,double *);

