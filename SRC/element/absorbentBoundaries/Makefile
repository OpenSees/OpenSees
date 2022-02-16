include ../../../Makefile.def

OBJS       = LysmerTriangle.o \
	ASDAbsorbingBoundary2D.o \
	ASDAbsorbingBoundary3D.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak 

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) 
	@$(RM) $(RMFLAGS) *.o
	@$(RM) $(RMFLAGS) .inslog* *~
	@$(RM) $(RMFLAGS) tca.map

spotless: clean					   

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
