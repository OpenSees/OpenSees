include ../../../Makefile.def

OBJS       = SimpleContact2D.o  SimpleContact3D.o \
 		BeamContact2D.o BeamContact3D.o BeamEndContact3D.o \
		TclSimpleContactCommand.o 

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak 

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) 
	@$(RM) $(RMFLAGS) *.o
	@$(RM) $(RMFLAGS) .inslog*
	@$(RM) $(RMFLAGS) tca.map

spotless: clean					   

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
