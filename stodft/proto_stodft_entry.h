/*-----------------------------------------------------------------*/
/* Main control of stochastic dft                                  */
#ifdef FAST_FILTER   
void controlStodftMin(CLASS *,BONDED *,GENERAL_DATA *,CP *,ANALYSIS *,
                      CLASS *,BONDED *,GENERAL_DATA *,CP *,ANALYSIS *);
#else
void controlStodftMin(CLASS *,BONDED *,GENERAL_DATA *,CP *,ANALYSIS *);
#endif
