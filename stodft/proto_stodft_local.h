/*-----------------------------------------------------------------*/
/* Local functions for stochastic dft                              */
/*-----------------------------------------------------------------*/
/* coeff.c			                                   */

void genCoeffNewtonHermit(STODFTINFO *,STODFTCOEFPOS *);
void genCoeffNewtonNoHermit(STODFTINFO *,STODFTCOEFPOS *);
void genSampNewtonHermit(STODFTINFO *,STODFTCOEFPOS *);
void genEnergyMax(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void genEnergyMin(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void genEigenOrb(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);

/*-----------------------------------------------------------------*/
/* filter.c                                                        */
void filterNewtonPolyHerm(CP *,CLASS *,GENERAL_DATA *,int);
void filterNewtonPolyNoHerm(CP *,CLASS *,GENERAL_DATA *,int);

/*-----------------------------------------------------------------*/
/* init.c                                                          */

void initStodft(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);

/*-----------------------------------------------------------------*/
/* min-CP-stodft.c                                                 */

void scfStodft(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void genStoOrbital(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);

/*-----------------------------------------------------------------*/
/* normh.c                                                         */

void normHNewtonHerm(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,
		     CLATOMS_POS *,double);
void normHNewtonNoHerm(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,
		       CLATOMS_POS *,double complex);
void calcRhoDeterm(CLASS *,BONDED *,GENERAL_DATA *,CP *,CPCOEFFS_POS *);
void calcRhoSto(CLASS *,BONDED *,GENERAL_DATA *,CP *,CPCOEFFS_POS *);







