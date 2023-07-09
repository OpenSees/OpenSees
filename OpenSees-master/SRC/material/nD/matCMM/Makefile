include ../../../../Makefile.def


OBJS       = MaterialCMM.o


all: $(OBJS) 

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



