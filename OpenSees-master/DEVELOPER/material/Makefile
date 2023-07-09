include ../Makefile.def

all:
	echo @pwd;
	@$(CD) cpp; $(MAKE);
	@$(CD) c; $(MAKE);
	@$(CD) fortran; $(MAKE);

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(CD) ./cpp; $(MAKE) wipe;
	@$(CD) ./c; $(MAKE) wipe;
	@$(CD) ./fortran; $(MAKE) wipe;
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so *.dylib

spotless: clean

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.
