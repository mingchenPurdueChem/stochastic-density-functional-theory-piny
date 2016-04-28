#=================================================================
#              MAIN_FILES 
#=================================================================


#=================================================================
coeff.o :                $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
			 $(STODFT_LOC) \
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
                         $(STODFT_LOC) \
                         $(DCODE)/stodft/filter.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/filter.c

init.o :                 $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(STODFT_LOC) \
                         $(DCODE)/stodft/init.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/init.c

min-CP-stodft.o :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
			 $(TYP_STAT) \
			 $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(ENR_CP_LOC) $(ENR_CPCON_ENT) $(ENR_CPCON_LOC)\
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
                        
#------------------------------------------------------------------
