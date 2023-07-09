include ../../../Makefile.def

all:
	@$(CD) matrixutil; $(MAKE);

# Miscellaneous

tidy:
	@rm -f Makefile.bak *~ #*# core

clean:  tidy
	@rm -f $(OBJS) *.o

spotless: clean
	@rm -f $(PROGRAM) fake

wipe: spotless
	@$(CD) matrixutil; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.





