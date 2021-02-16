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
void genChebyHermit(STODFTINFO *,STODFTCOEFPOS *,int);
void genChebyHermitTrueChemPot(STODFTINFO *,STODFTCOEFPOS *,int);
void calcChebyCoeffWrapper(STODFTINFO *,STODFTCOEFPOS *,int);
void calcChebyCoeff(STODFTINFO *,STODFTCOEFPOS *,double *,double *, int, double);
double calcFitErrorCheby(STODFTINFO *,STODFTCOEFPOS *,int);
double testChebyCoeff(STODFTINFO *,STODFTCOEFPOS *,double *,double *, int, double);
void genCoeffNewtonHermitEntropy(STODFTINFO *,STODFTCOEFPOS *);
int roundFFT(int);
/*-----------------------------------------------------------------*/
/* filter.c                                                        */
void filterNewtonPolyHerm(CP *,CLASS *,GENERAL_DATA *,int);
void filterChebyPolyHerm(CP *,CLASS *,GENERAL_DATA *,int);
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
void calcStoEntropy(CP *);
/*-----------------------------------------------------------------*/
/* init.c                                                          */
void commStodft(CLASS *,BONDED *,GENERAL_DATA *,CP *);
void initStodft(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void reInitWaveFunMin(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void reInitComm(CP *,CPCOEFFS_POS *);
void stoRealloc(CP *,CPCOEFFS_POS *);
void reallocScratch(CP *,int);
void initFilterDiag(CP *);
/*-----------------------------------------------------------------*/
/* min-CP-stodft.c                                                 */
#ifndef FAST_FILTER   
void scfStodftInterp(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void scfStodftCheby(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void scfStodftEnergyWindow(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void scfStodftEnergyWindowFrag(CLASS *,BONDED *,GENERAL_DATA *,
                    CP *, CP *, GENERAL_DATA *, CLASS *,
                    int);
#endif
void scfStodftFilterDiag(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
/*-----------------------------------------------------------------*/
/* gen-stodft-wf.c						   */
void genStoOrbitalInterp(CLASS *,GENERAL_DATA *,CP *,int);
void genStoOrbitalCheby(CLASS *,GENERAL_DATA *,CP *,int);
void genStoOrbitalFake(CLASS *,GENERAL_DATA *,CP *,int );
void genStoOrbitalEnergyWindow(CLASS *,GENERAL_DATA *,CP *,int );
void genStoOrbitalEnergyWindowFake(CLASS *,GENERAL_DATA *,CP *,int);
void genStoOrbitalEnergyWindowFragFake(CLASS *,GENERAL_DATA *,
            CP *,GENERAL_DATA *,CP *,CLASS *,int);
void genStoOrbitalEnergyWindowFrag(CLASS *,GENERAL_DATA *,
            CP *,GENERAL_DATA *,CP *,CLASS *,int);
/*-----------------------------------------------------------------*/
/* gen-noise.c							   */
void genNoiseOrbital(CP *,CPCOEFFS_POS *);
void genNoiseOrbitalReal(CP *,CPCOEFFS_POS *);
/*-----------------------------------------------------------------*/
/* normh.c                                                         */
void normHNewtonHerm(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,
		     CLATOMS_POS *,double);
void normHCheby(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *,int,int);
/*-----------------------------------------------------------------*/
/* density-init.c                                                  */
void calcRhoInit(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void calcRhoDetInit(CLASS *,BONDED *,GENERAL_DATA *,CP *,CPCOEFFS_POS *);
void outputInitDensity(CP *,CELL *);
void calcRhoStoInit(CLASS *,BONDED *,GENERAL_DATA *,CP *,CPCOEFFS_POS *);
void readRho(CLASS *,BONDED *,GENERAL_DATA *,CP *,CPCOEFFS_POS *);
void calcRhoFragInit(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void calcRhoOffInit(CLASS *,BONDED *,GENERAL_DATA *,CP *,CPCOEFFS_POS  *);
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
void outputDensity(CP *,CELL *,int);
void calcRhoFilterDiagHybrid(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void calcRhoStoHybridEnergyWindow(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);

/*-----------------------------------------------------------------*/
/* calc-chempot.c                                                  */
void calcChemPotInterp(CP *);
void genChemPotInterpPoints(STODFTINFO *,STODFTCOEFPOS *);
double solveLagrangePolyInterp(int ,double *, double *, double, double *,double *);
double calcLagrangeInterpFun(int ,double ,double * ,double * ,double *);
double calcLagrangeInterpDrv(int ,double ,double * ,double * ,double *);
void updateChemPot(STODFTINFO *,STODFTCOEFPOS *);
void adjChemPot(STODFTINFO *,STODFTCOEFPOS *);
void calcChemPotMetal(CP *);
void genChemPotEnergyWindows(STODFTINFO *,STODFTCOEFPOS *);
double calcNumElecSmear(int, double, double,double *,int);
/*-----------------------------------------------------------------*/
/* calc-chempot-chebyshev.c                                        */
void calcChemPotCheby(CP *,CLASS *,GENERAL_DATA *,int);
void calcChemPotChebyEWFrag(CP *,CLASS *,GENERAL_DATA *,int);
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
int orthogSVD(int ,int ,double *,int);
void buildKSMatrix(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void diagKSMatrix(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void diagSymMatWrapper(int ,double *,double *);
void genMatrixMulWrapper(int ,int ,double *,double *,double *,int);
void calcForceWrapper(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *,
		      double *,double *);
// Reduced Energy Calculation
/*-----------------------------------------------------------------*/
/* cp-energy-ee-rho-stodft.c                                       */
void coefForceCalcHybridSCF(CPEWALD *,int,double *,double *,
                            double *,double *,double *,double *,
                            double *,double *,double *,double *,double *,
                            int ,double *,COMMUNICATE *,int ,int ,int ,int ,int ,
                            PARA_FFT_PKG3D *,CP *,CLASS *,GENERAL_DATA *);
void coefForceCalcHybridEnergy(CPEWALD *,int,double *,double *,
                            double *,double *,double *,double *,
                            double *,double *,double *,double *,double *,
                            int ,double *,COMMUNICATE *,int ,int ,int ,int ,int ,
                            PARA_FFT_PKG3D *,CP *,CLASS *,GENERAL_DATA *);
void cpGetVksStodft(CPOPTS *,CPSCR *,CPEWALD *,EWALD *,
                COMMUNICATE *,CP_COMM_STATE_PKG *,CP_COMM_STATE_PKG *,
                STAT_AVG *,double *,CELL *, char *,double *,double ,double,
                int,int,int,int,int,int,PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,
                PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,int);
/*-----------------------------------------------------------------*/
/* cp-energy-eext-stodft.c					   */
void controlEwdLocPreScf(CLATOMS_INFO *,CLATOMS_POS *,CELL *, PTENS *, EWALD *,
			 CPEWALD *,CPSCR *,PSEUDO *,EWD_SCR *,CPOPTS *,
			 ATOMMAPS *, double *,double *,COMMUNICATE *,
                         FOR_SCR *,int ,int);
void getNlPotPvFatmSCF(CLATOMS_INFO *,CLATOMS_POS *,CELL *,CPCOEFFS_INFO *,CPSCR *,
                       EWD_SCR *,CPOPTS *,PSEUDO *,ATOMMAPS *,double *,int,double *);
void sumNlPot(int,int,int,int,int,int,int,int,int,
              int *,double *,double *,double *,double *);
void controlEwdNonlocFilter(CLATOMS_INFO *,CLATOMS_POS *,CPCOEFFS_INFO *,CPCOEFFS_POS *,
                      CELL *, PTENS *, CPEWALD *,CPSCR *, PSEUDO *, EWD_SCR *,
                      CPOPTS *, ATOMMAPS *,COMMUNICATE *,FOR_SCR *);
void controlNlmatFilter(CLATOMS_INFO *,CPCOEFFS_INFO *,CPCOEFFS_POS *,
                   CPSCR *,CPOPTS *,PSEUDO *,EWD_SCR *,ATOMMAPS *,
                   int ,int ,int *,double ,double ,double ,double ,
                   int ,YLM_CONS *);
void getNlmatFilter(int,int,int,int,int,int,int,int,int,
               double *,double *,double *,double *,double *,double *,
               double,double,double,double,double,double ,double *,double *);
void getYlmOnly(double,double,double,double,double *,double *,YLM_CONS *);
void allocRealNl(CP *,CLASS *);
/*-----------------------------------------------------------------*/
/* energy-wrapper-scf.c                                            */
void calcLocalPseudoScf(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS  *,CLATOMS_POS *);
void calcNonLocalPseudoScf(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS *,CLATOMS_POS *);
void calcKSPot(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS *,CLATOMS_POS *);
void calcCoefForceScf(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS *,CLATOMS_POS *);
void calcCoefForceEnergy(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS *,CLATOMS_POS *);
void calcCoefForceWrapSCF(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS *,CLATOMS_POS *);
/*-----------------------------------------------------------------*/
/* energy-wrapper-post-scf.c                                       */
void calcCoefForceWrap(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS *,CLATOMS_POS *);
void calcLocExtPostScf(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS  *,CLATOMS_POS *);
void calcNlPseudoPostScf(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS *,CLATOMS_POS *);
void calcCoefForcePosScf(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS *,CLATOMS_POS *);
void calcCoefForceWrapReduce(CLASS *,GENERAL_DATA *,CP *,CPCOEFFS_POS *,CLATOMS_POS *);

// Nuclei force calculation
/*-----------------------------------------------------------------*/
/* calc-nuclei-force.c	                                           */
void calcEnergyForce(CLASS *,GENERAL_DATA *,CP *,BONDED *,CPCOEFFS_POS *,CLATOMS_POS *);
void calcEnergyForceFilterDiag(CLASS *,GENERAL_DATA *,CP *,BONDED *,CPCOEFFS_POS *,CLATOMS_POS *);

/*-----------------------------------------------------------------*/
/* checkpointIO.c */
void checkpointOutput(CP *,GENERAL_DATA *);
void checkpointInput(CP *,GENERAL_DATA *,CLASS *);
void checkpointOutputDist(CP *,GENERAL_DATA *);
void checkpointInputDist(CP *,GENERAL_DATA *,CLASS *);

#ifdef FAST_FILTER   
void scfStodftCheby(CLASS *,BONDED *,GENERAL_DATA *,
                    CP *,CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void scfStodftEnergyWindow(CLASS *,BONDED *,GENERAL_DATA *,
                    CP *,CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void scfStodftInterp(CLASS *,BONDED *,GENERAL_DATA *,
                    CP *,CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void scfStodftEnergyWindowFrag(CLASS *,BONDED *,GENERAL_DATA *,CP *,
                    CLASS *,BONDED *,GENERAL_DATA *,CP *,
                    CP *, GENERAL_DATA *, CLASS *,int);
void genStoOrbitalInterpTest(CLASS *,GENERAL_DATA *,CP *,
            CLASS *,GENERAL_DATA *,CP *,int);
void genStoOrbitalChebyTest(CLASS *,GENERAL_DATA *,
                CP *,CLASS *,GENERAL_DATA *,CP *,int);
void genStoOrbitalEnergyWindowTest(CLASS *,GENERAL_DATA *,
                CP *,CLASS *,GENERAL_DATA *,CP *,int);
void genStoOrbitalEnergyWindowFragTest(CLASS *class,GENERAL_DATA *general_data,
            CP *,GENERAL_DATA *,CP *,CLASS *,CLASS *,GENERAL_DATA *,CP *,int);
void scfStodftEnergyWindowFragTest(CLASS *,BONDED *,GENERAL_DATA *,
                    CP *,CLASS *,BONDED *,GENERAL_DATA *,CP *,
                    CP *, GENERAL_DATA *, CLASS *,int);
void filterNewtonPolyHermFake(CP *,CLASS *,GENERAL_DATA *,int);
void filterChebyPolyHermFake(CP *,CLASS *,GENERAL_DATA *,int,int);
void broadcastWfDet(CP *,CLASS *,GENERAL_DATA *,CP *);
void calcChebyMomentsFake(CP *,CLASS *,GENERAL_DATA *,int);
#endif
