 
include ../../Makefile.def

DOBJ  = dsrc2c.o jcg.o jsi.o rscg.o rssi.o sor.o ssorcg.o ssorsi.o

all: double 
 
double: $(DOBJ)
	@$(RM) $(RMFLAGS) $(ITPACK_LIBRARY)
	@$(AR) $(ARFLAGS) $(ITPACK_LIBRARY) $(DOBJ)
	@$(RANLIB) $(ITPACK_LIBRARY)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core go

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o test go

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
