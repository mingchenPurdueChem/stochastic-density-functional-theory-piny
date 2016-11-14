#=================================================================
#              MAIN_FILES 
#=================================================================


#=================================================================
frag-scf.o :             $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) $(TYP_CP) \
                         $(FRND_ENT) $(MATH) \
			 $(ENR_CPCON_ENT) $(ENR_CPCON_LOC) \
			 $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(COMM_WRAP) \
			 $(INT_CPMIN_ENT) $(FRAG_INT_ENT) \
                         $(DCODE)/stodft/fragment/frag-scf.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/fragment/frag-scf.c

init-frag.o :	         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(COMM_WRAP) $(FRAG_LOC)\
                         $(DCODE)/stodft/fragment/init-frag.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/fragment/init-frag.c

control-cp-min-frag.o :	 $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(COMM_WRAP) $(FRAG_LOC) \
			 $(FRAG_INT_LOC) \
	                 $(DCODE)/stodft/fragment/control-cp-min-frag.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/fragment/control-cp-min-frag.c


parse-frag.o :           $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_PAR) $(TYP_STAT) $(PARSE_ENT) $(PARSE_LOC) \
                         $(SIM_ENT) $(MOL_ENT) $(INTRA_ENT) $(COORD_ENT) \
                         $(COORD_LOC) $(CPEWALD_ENT) $(INTER_ENT) \
                         $(INTER_LOC) $(VPS_ENT) $(LISTS_ENT) $(SCRATCH_ENT) \
                         $(ENR_CPCON_ENT) $(REAL_LOC) $(SAMPL_CLASS_ENT) \
                         $(SAMPL_CP_ENT) $(SAMPL_CP_LOC) $(COORD_CP_ENT) \
                         $(COORD_CP_LOC) $(MATH) $(PIMD_ENT) $(PIMD_LOC) \
                         $(FRND_ENT) $(COMM_ENT) $(COMM_LOC) $(COMM_WRAP) \
                         $(INT_CPMIN_ENT) $(FRAG_INT_LOC) \
                         $(DCODE)/stodft/fragment/interface-frag/parse-frag.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/fragment/interface-frag/parse-frag.c

copy-input.o :           $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(FRAG_INT_LOC) \
                         $(DCODE)/stodft/fragment/interface-frag/copy-input.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/fragment/interface-frag/copy-input.c


all-control-frag.o :	 $(STANDARD) $(DEFINES) \
			 $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
			 $(TYP_STAT) $(TYP_PAR) \
			 $(MATH) $(FRND_ENT) \
			 $(MOL_ENT) $(MOL_LOC) $(INTER_ENT) $(INTER_LOC) \
			 $(INTRA_LOC) $(SCRATCH_ENT) $(HANDLE_ENT) $(COMM_WRAP) \
			 $(CPEWALD_ENT) $(CPEWALD_LOC) $(ENR_CP_LOC) $(VPS_ENT) \
			 $(VPS_LOC) $(FRAG_INT_LOC) \
			 $(DCODE)/stodft/fragment/interface-frag/all-control-frag.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/fragment/interface-frag/all-control-frag.c

all-mall-frag.o	:	 $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
			 $(TYP_STAT) $(TYP_PAR) \
			 $(SCRATCH_ENT) $(SCRATCH_LOC) $(LISTS_ENT) $(LISTS_LOC) \
			 $(REAL_LOC) $(FRND_ENT) $(COMM_WRAP) $(FRAG_INT_LOC) \
			 $(DCODE)/stodft/fragment/interface-frag/all-mall-frag.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/fragment/interface-frag/all-mall-frag.c

all-read-frag.o :	 $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
			 $(TYP_STAT) $(TYP_PAR) \
			 $(COORD_ENT) $(MATH) $(ENR_PIMD_LOC) $(ENR_CPCON_LOC) \
			 $(COORD_CP_ENT) $(COORD_CP_LOC) $(HANDLE_ENT) $(FRND_ENT) \
			 $(COMM_WRAP) $(FRAG_INT_LOC) \
			 $(DCODE)/stodft/fragment/interface-frag/all-read-frag.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/fragment/interface-frag/all-read-frag.c

gen-wave-frag.o :	 $(STANDARD) $(DEFINES) \
			 $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
			 $(TYP_STAT) $(TYP_PAR) \
			 $(FRND_ENT) $(MATH) $(COMM_WRAP) $(ENR_CPCON_LOC) \
			 $(ENR_CP_LOC) $(CPEWALD_LOC) $(COORD_CP_LOC) \
			 $(HANDLE_ENT) $(FRAG_INT_LOC) \
			 $(DCODE)/stodft/fragment/interface-frag/gen-wave-frag.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/fragment/interface-frag/gen-wave-frag.c

init-coord-hmat-fft.o :  $(STANDARD) $(DEFINES) \
			 $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
			 $(TYP_STAT) \
			 $(MATH) $(FRND_ENT) $(COMM_WRAP) $(FRAG_INT_LOC) \
			 $(DCODE)/stodft/fragment/interface-frag/init-coord-hmat-fft.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/fragment/interface-frag/gen-wave-frag.o

reinitFFT.o :		 $(STANDARD) $(DEFINES) \
			 $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
			 $(MATH) $(FRND_ENT) $(COMM_WRAP) $(FRAG_INT_LOC) \
			 $(DCODE)/stodft/fragment/interface-frag/reinitFFT.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/stodft/fragment/interface-frag/reinitFFT.c

#------------------------------------------------------------------
