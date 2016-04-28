/*-----------------------------------------------------------------*/
/* Local functions for stochastic dft                              */
/*-----------------------------------------------------------------*/
/* coeff.c			                                   */

void genCoeffNewtonHermit(STODFTINFO *,STODFTCOEFPOS *);
void genCoeffNewtonNoHermit(STODFTINFO *,STODFTCOEFPOS *);
void genSampNewtonHermit(STODFTINFO *,STODFTCOEFPOS *);
void genEnergyMax(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);
void genEnergyMin(CP *,CLASS *,GENERAL_DATA *,CPCOEFFS_POS *,CLATOMS_POS *);

/*-----------------------------------------------------------------*/
/* filter.c                                                        */

void filterNewtonPolyHerm(CP *,int,EWALD *,EWD_SCR *,CELL *,CLATOMS_INFO *,
                            CLATOMS_POS *,ATOMMAPS *,STAT_AVG *,PTENS *,SIMOPTS *,
                            FOR_SCR *);
void filterNewtonPolyNoHerm(CP *,int,EWALD *,EWD_SCR *,CELL *,CLATOMS_INFO *,
                            CLATOMS_POS *,ATOMMAPS *,STAT_AVG *,PTENS *,SIMOPTS *,
                            FOR_SCR *);

/*-----------------------------------------------------------------*/
/* init.c                                                          */

void initStodft(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);

/*-----------------------------------------------------------------*/
/* min-CP-stodft.c                                                 */

void scfStodft(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);
void genStoOrbital(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);

/*-----------------------------------------------------------------*/
/* normh.c                                                         */

void normHNewtonHerm(CP *,CPCOEFFS_POS *,CPCOEFFS_INFO *,CELL *,CLATOMS_INFO *,
		 CLATOMS_POS *,EWALD *,EWD_SCR *,ATOMMAPS *,FOR_SCR *,
                 STAT_AVG *,PTENS *,double);

void normHNewtonNoHerm(CP *,CPCOEFFS_POS *,CPCOEFFS_INFO *,CELL *,CLATOMS_INFO *,
                 CLATOMS_POS *,EWALD *,EWD_SCR *,ATOMMAPS *,FOR_SCR *,
                 STAT_AVG *,PTENS *,double complex);
void calcRhoDeterm(CLASS *,BONDED *,GENERAL_DATA *,CP *,CPCOEFFS_POS *);
void calcRhoSto(CLASS *,BONDED *,GENERAL_DATA *,CP *,CPCOEFFS_POS *);







