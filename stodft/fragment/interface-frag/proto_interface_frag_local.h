/*-----------------------------------------------------------------*/
/* Local functions for fragmentation interface                     */
/*-----------------------------------------------------------------*/
/* copy-input.c                                                    */
void copySimParam(GENERAL_DATA *,BONDED *,CLASS *,CP *,GENERAL_DATA *,BONDED *,
                  CLASS *,CP *,CLASS_PARSE *,CP_PARSE *,FILENAME_PARSE *);

/*-----------------------------------------------------------------*/
/* all-control-frag.c	                                           */

void controlInterParamsFrag(GENERAL_DATA *,CLASS *,CP *,BONDED *,CP *,
                SPLINE_PARSE *,FILENAME_PARSE *,CLASS_PARSE *);



void controlMolParamsFrag(CLASS *,GENERAL_DATA *,BONDED *,CP *,CLASS *,
                          GENERAL_DATA *,BONDED *,CP *,CLASS_PARSE *,
                          CP_PARSE *,FREE_PARSE *,FILENAME_PARSE *);

void controlSetMolParamsFrag(CP_PARSE *,CLASS_PARSE *,
                FILENAME_PARSE *,FREE_PARSE *,BONDED *,CLASS *,CP *,
                GENERAL_DATA *,CLASS *class,CP *,GENERAL_DATA *,DICT_MOL *,DICT_WORD *,
                char *,int *,int ,double ,int ,int ,DICT_MOL *);

void controlSetCpEwaldFrag(GENERAL_DATA *,CLASS *,CP *,BONDED *,CP *,CLASS *,
                           GENERAL_DATA *,BONDED *,CP_PARSE *);

void controlFFTPkgFrag(GENERAL_DATA *,CLASS *,CP *,CP *);

void controlVpsParamsFrag(GENERAL_DATA *,CLASS *,CP *,FILENAME_PARSE *,
			    SPLINE_PARSE *,CP_PARSE *);


/*-----------------------------------------------------------------*/
/* all-mall-frag.c                                                 */

void mallMakeListsFrag(CLASS *,GENERAL_DATA *,BONDED *,int;

void controlMallScratchFrag(CLASS *,BONDED *,CP *,GENERAL_DATA *);

void mall_integrator_scr_frag(int ,int ,CLATOMS_INFO *,THERM_INFO *,
			      THERM_INFO *,INT_SCR *,double *,int,MPI_Comm);

void mall_intra_scr_frag(INTRA_SCR *,double *,int ,MPI_Comm);

void mall_ewald_scr_frag(int ,int ,int ,int, int, PART_MESH *,
			 EWD_SCR *,double *,int ,MPI_Comm);

void mall_cp_scr_frag(CPTHERM_INFO *,CPOPTS *,CPEWALD *,CPSCR *,CPCOEFFS_INFO *,
		      PSEUDO *,PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,
		      CP_COMM_STATE_PKG *,CP_COMM_STATE_PKG *,int ,
		      double *,int ,int , MPI_Comm);

void mall_atm_forc_scr_frag(int natm_tot,FOR_SCR *for_scr,int pme_on,
                       int ilnk_lst,int iver_lst, int ,
                       double *,int ,int ,MPI_Comm);

/*-----------------------------------------------------------------*/
/* all-read-frag.c                                                 */

void readCoefFrag(CP *,GENERAL_DATA *,CLASS *,FILENAME_PARSE *,CP_PARSE *,double *);

void readHmatFrag(CLASS *,GENERAL_DATA *,CP *,int , double *,int *);

/*-----------------------------------------------------------------*/
/* gen-wave-frag.c                                                 */

void gen_wave_frag(CLASS *,GENERAL_DATA *,CP *,CP_PARSE *,NAME *);

/*-----------------------------------------------------------------*/
/* init-coord-hmat-fft.c                                           */

void initCoordHmatFFT(GENERAL_DATA *,CLASS *class,CP *,GENERAL_DATA *,
		      CLASS *,CP *);

void passAtomCoord(GENERAL_DATA *,CLASS *,CP *,GENERAL_DATA *,CLASS *,
		    CP *,int ,double *);

void initFFTMap(GENERAL_DATA *,CLASS *,CP *,GENERAL_DATA *,CLASS *,CP *,
		int ,double *);

/*-----------------------------------------------------------------*/
/* init-coord-hmat-fft.c                                           */

void reInitFFT(GENERAL_DATA *,CLASS *,CP *,GENERAL_DATA *,CLASS *,CP *,int);
