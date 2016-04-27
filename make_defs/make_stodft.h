#=================================================================
#              MAIN_FILES 
#=================================================================


#=================================================================
:        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(DAFED_LOC) $(DAFED_ENT) \
                         $(DCODE)/dafed/control_dafed.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/dafed/control_dafed.c

integrate_dafed.o :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(DAFED_LOC) $(DAFED_ENT) \
                         $(DCODE)/dafed/integrate_dafed.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/dafed/integrate_dafed.c


dafed_io.o :             $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(DAFED_LOC) $(DAFED_ENT) \
                         $(DCODE)/dafed/dafed_io.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/dafed/dafed_io.c

bias_update.o :          $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(FRND_ENT) $(MATH) \
                         $(DAFED_LOC) $(DAFED_ENT) \
                         $(DCODE)/dafed/bias_update.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/dafed/bias_update.c

energy_control_dafed.o : $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(MATH) \
                         $(DAFED_ENERGY) \
                         $(DCODE)/energy/control/energy_control_dafed.c
	$(ECHO)	$@
	$(COBJ)	$(DCODE)/energy/control/energy_control_dafed.c

force_Phi.o :            $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(MATH) \
                         $(DAFED_ENERGY) \
                         $(DCODE)/energy/dafed/force_Phi.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/energy/dafed/force_Phi.c

force_bias.o :           $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(MATH) \
                         $(DAFED_ENERGY) \
                         $(DCODE)/energy/dafed/force_bias.c
	$(ECHO) $@
	$(COBJ) $(DCODE)/energy/dafed/force_bias.c
                         
force_dafed_final.o :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(MATH) \
                         $(DAFED_ENERGY) \
                         $(DCODE)/energy/dafed/force_dafed_final.c
	$(ECHO)	$@
	$(COBJ)	$(DCODE)/energy/dafed/force_dafed_final.c
                         
#------------------------------------------------------------------
