include ../Makefile.def

# the example dircetories that 

EXAMPLE1 = $(FE)/../EXAMPLES/Example1
PLANE_FRAME = $(FE)/../EXAMPLES/PlaneFrame
G3 = $(FE)/../EXAMPLES/g3

all:
	@$(CD) $(EXAMPLE1); $(MAKE);
	@$(CD) $(PLANE_FRAME); $(MAKE);

clean:
	@$(CD) $(EXAMPLE1); $(MAKE) clean;
	@$(CD) $(PLANE_FRAME); $(MAKE) clean;

wipe: clean
	@$(CD) $(EXAMPLE1); $(MAKE) wipe;
	@$(CD) $(PLANE_FRAME); $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
