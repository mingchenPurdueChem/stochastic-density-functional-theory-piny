#=================================================================
#              MAIN_FILES 
#=================================================================


#=================================================================
coeff.o :                $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
			 $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(STODFT_LOC) \
                         $(DCODE)/stodft/coeff.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/coeff.c

control-stodft.o :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(STODFT_LOC) \
                         $(DCODE)/stodft/control-stodft.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/control-stodft.c


filters.o :              $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(STODFT_LOC) \
                         $(DCODE)/stodft/filter.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/filter.c

init.o :                 $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
			 $(ENR_CPCON_ENT) $(ENR_CPCON_LOC) \
                         $(ENR_CP_LOC) $(STODFT_LOC) \
                         $(DCODE)/stodft/init.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/init.c

min-CP-stodft.o :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
			 $(TYP_STAT) $(TYP_PAR) \
			 $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_ENT) $(ENR_CPCON_LOC)\
			 $(VPS_ENT) \
			 $(STODFT_LOC) \
                         $(DCODE)/stodft/min-CP-stodft.c
	$(ECHO)	$@
	$(COBJ)	$(DCODE)/stodft/min-CP-stodft.c

normh.o :                $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(STODFT_LOC) \
                         $(DCODE)/stodft/normh.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/normh.c

density.o :		 $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(STODFT_LOC) \
                         $(DCODE)/stodft/density.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/density.c

density-init.o :         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(STODFT_LOC) \
                         $(DCODE)/stodft/density-init.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/density-init.c	

calc-chempot.o :	 $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(STODFT_LOC) \
                         $(DCODE)/stodft/calc-chempot.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/calc-chempot.c

diis.o :                 $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(STODFT_LOC) \
                         $(DCODE)/stodft/diis.c

	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/diis.c

filter-diag.o :          $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(STODFT_LOC) \
                         $(DCODE)/stodft/filter-diag.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/filter-diag.c

calc-chempot-chebyshev.o:$(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(STODFT_LOC) \
                         $(DCODE)/stodft/calc-chempot-chebyshev.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/calc-chempot-chebyshev.c

calc-spectral-range.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(STODFT_LOC) \
                         $(DCODE)/stodft/calc-spectral-range.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/calc-spectral-range.c

gen-noise.o :            $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_ENT) $(ENR_CPCON_LOC)\
                         $(STODFT_LOC) \
                         $(DCODE)/stodft/gen-noise.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/gen-noise.c

gen-stodft-wf.o :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_ENT) $(ENR_CPCON_LOC)\
                         $(STODFT_LOC) \
                         $(DCODE)/stodft/gen-stodft-wf.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/gen-stodft-wf.c

calc-energy.o :          $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_ENT) $(ENR_CPCON_LOC)\
                         $(STODFT_LOC) \
                         $(DCODE)/stodft/calc-energy.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/calc-energy.c
calc-nuclei-force.o :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_ENT) $(ENR_CPCON_LOC)\
                         $(ENR_CTRL_ENT) $(STODFT_LOC) \
                         $(DCODE)/stodft/calc-nuclei-force.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/calc-nuclei-force.c

cp-energy-ee-rho-stodft.o :  $(STANDARD) $(DEFINES) \
			 $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(STODFT_LOC) \
                         $(DCODE)/stodft/reduced-energy/cp-energy-ee-rho-stodft.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/reduced-energy/cp-energy-ee-rho-stodft.c

cp-energy-eext-stodft.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(STODFT_LOC) \
                         $(DCODE)/stodft/reduced-energy/cp-energy-eext-stodft.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/reduced-energy/cp-energy-eext-stodft.c

energy-wrapper-post-scf.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(STODFT_LOC) \
                         $(DCODE)/stodft/reduced-energy/energy-wrapper-post-scf.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/reduced-energy/energy-wrapper-post-scf.c

energy-wrapper-scf.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(STODFT_LOC) \
                         $(DCODE)/stodft/reduced-energy/energy-wrapper-scf.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/reduced-energy/energy-wrapper-scf.c                     

checkpointIO.o :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) $(TYP_PAR) \
                         $(TYP_CP) $(TYP_PAR) $(ENR_CP_LOC) \
                         $(FRND_ENT) $(MATH) $(COMM_WRAP)\
                         $(STODFT_LOC) \
                         $(DCODE)/stodft/checkpointIO.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/checkpointIO.c

calc-friction.o :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) $(TYP_PAR) \
                         $(TYP_CP) $(TYP_PAR) $(ENR_CP_LOC) \
                         $(FRND_ENT) $(MATH) $(COMM_WRAP)\
                         $(STODFT_LOC) \
                         $(DCODE)/stodft/calc-friction.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/calc-friction.c


cp-energy-eext-fric.o:  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) $(ENR_CP_LOC) $(ENR_CPCON_LOC)\
                         $(FRND_ENT) $(MATH) $(COMM_WRAP)\
                         $(STODFT_LOC) \
                         $(DCODE)/stodft/reduced-energy/cp-energy-eext-fric.c

	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/reduced-energy/cp-energy-eext-fric.c

rational-apprx.o:        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) $(ENR_CP_LOC) $(ENR_CPCON_LOC)\
                         $(FRND_ENT) $(MATH) $(COMM_WRAP)\
                         $(STODFT_LOC) \
                         $(DCODE)/stodft/rational/rational-apprx.c

	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/rational/rational-apprx.c

#------------------------------------------------------------------
