include ../../../Makefile.def

OBJS       = 	FourNodeQuadUP.o TclFourNodeQuadUPCommand.o \
	BrickUP.o \
	Nine_Four_Node_QuadUP.o \
	Twenty_Eight_Node_BrickUP.o \
	shp3dv.o

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
