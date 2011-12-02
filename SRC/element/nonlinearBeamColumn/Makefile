include ../../../Makefile.def

all:
	@$(CD) element; $(MAKE);
	@$(CD) matrixutil; $(MAKE);
	@$(CD) quadrule; $(MAKE);

# Miscellaneous

tidy:
	@rm -f Makefile.bak *~ #*# core

clean:  tidy
	@rm -f $(OBJS) *.o

spotless: clean
	@rm -f $(PROGRAM) fake

wipe: spotless
	@$(CD) element; $(MAKE) wipe;
	@$(CD) matrixutil; $(MAKE) wipe;
	@$(CD) quadrule; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.





