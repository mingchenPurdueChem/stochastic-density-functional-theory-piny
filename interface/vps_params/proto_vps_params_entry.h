/*---------------------------------------------------------------*/
/* Control_vps_params.c */

void control_vps_params(PSEUDO *,CELL *,FILENAME_PARSE *,
   		        SPLINE_PARSE *,int ,NAME [],
			double *,int, int,int ,int , 
                        COMMUNICATE *, double, CPCOEFFS_INFO *);
void controlNlppReal(CP *,CLASS *,GENERAL_DATA *,FILENAME_PARSE *);
void mapRealSpaceGrid(CP *, CLASS *, GENERAL_DATA *);
void initRealNlppWf(CP *,CLASS *,GENERAL_DATA *);
