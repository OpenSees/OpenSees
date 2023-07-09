include ../../../Makefile.def

ifdef HAVESNOPT

OBJS       = SNOPTAnalysis.o \
	SNOPTCLASS.o \
	toyOptfunction.o

else

OBJS = 

endif

################### TARGETS ########################
all: $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core test

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core test
	@$(RM) $(RMFLAGS) test.\$\$\$ 
	@$(RM) $(RMFLAGS) *.ti 
	@$(RM) $(RMFLAGS) tca.map 
	@$(RM) $(RMFLAGS) .inslog*

spotless: clean

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.
