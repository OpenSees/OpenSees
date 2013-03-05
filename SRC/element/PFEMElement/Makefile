include ../../../Makefile.def

OBJS = PFEMElement2D.o TclModelBuilder_addPFEMElement.o PFEMMesher2D.o PFEMMesher3D.o TclPFEMCommands.o 

#OBJS = 
all: $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o *.so

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
