
void cp_wannier_control(CLASS *,GENERAL_DATA *,CP *);

void calcul_dipole(CLASS *, GENERAL_DATA *, CP *, double *, double ***, double ***);

void calcul_wannier(GENERAL_DATA *,CP *, double ***, double ***,double *,
                        double **, int );

void calcul_wannier_dvr(GENERAL_DATA *,CP *, double ***, double ***,double *, 
                        double **, int );

void calcul_initial_wannier(GENERAL_DATA *,CP *);

void calcul_initial_wannier_dvr(GENERAL_DATA *,CP *);

