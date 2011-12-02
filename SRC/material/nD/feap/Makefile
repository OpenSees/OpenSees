include ../../../../Makefile.def

OBJS       = feapCommon.o \
	TclFeapMaterialCommand.o \
	FeapMaterial01.o matl01.o \
	FeapMaterial02.o matl02.o \
	FeapMaterial03.o matl03.o

all:	$(OBJS) 

# Miscellaneous

tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o
	@$(RM) $(RMFLAGS) .inslog*
	@$(RM) $(RMFLAGS) tca.map

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



