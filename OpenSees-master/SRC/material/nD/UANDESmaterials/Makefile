include ../../../../Makefile.def

OBJS       = SAniSandMS.o\
	     SAniSandMS3D.o\
	     SAniSandMSPlaneStrain.o

all:	$(OBJS)

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



