/*-----------------------------------------------------------------*/
/* Local functions for stochastic dft                              */
/*-----------------------------------------------------------------*/
/* coeff.c			                                   */

void genCoeffNewtonHermit(STODFTINFO *,STODFTCOEFPOS *);
void genSampNewtonHermit(STODFTINFO *,STODFTCOEFPOS *);
void genEigenOrb(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);

/*-----------------------------------------------------------------*/
/* filter.c                                                        */
void filterNewtonPolyHerm(CP *,CLASS *,GENERAL_DATA *,int);
void genEnergyMax(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void genEnergyMin(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);

/*-----------------------------------------------------------------*/
/* init.c                                                          */

void initStodft(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void calcRhoInit(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void reInitWaveFunMin(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void reInitComm(CP *,CPCOEFFS_POS *);
void stoRealloc(CP *,CPCOEFFS_POS *);
void reallocScratch(CP *,int);

/*-----------------------------------------------------------------*/
/* min-CP-stodft.c                                                 */

void scfStodft(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void genStoOrbital(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);

/*-----------------------------------------------------------------*/
/* normh.c                                                         */

void normHNewtonHerm(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,
		     CLATOMS_POS *,double);
void calcRhoDet(CLASS *,BONDED *,GENERAL_DATA *,CP *,CPCOEFFS_POS *);
void calcRhoSto(CLASS *,BONDED *,GENERAL_DATA *,CP *,CPCOEFFS_POS *);







