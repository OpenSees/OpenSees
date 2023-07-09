include ../../../Makefile.def

# OBJS       = IGAPatch.o 
OBJS       = NurbsDers.o gaussQuadrature.o IGASurfacePatch.o IGAFollowerLoad.o

# Compilation control

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) Makefile.bak $(OBJS) *.o *~ #*# core

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
