
include ./Makefile.def

all:
	@cd ./element; $(MAKE);
	@cd ./material; $(MAKE);
	@cd ./recorder; $(MAKE);
	@cd ./system; $(MAKE);
	@cd ./integrator; $(MAKE);

tidy:
	@rm -f Makefile.bak *~ #*# core

clean:  tidy
	@rm -f *.o

spotless: clean

wipe: spotless
	@cd ./element; $(MAKE) wipe;
	@cd ./material; $(MAKE) wipe;
	@cd ./recorder; $(MAKE) wipe;
	@cd ./system; $(MAKE) wipe;
	@cd ./integrator; $(MAKE) wipe;

