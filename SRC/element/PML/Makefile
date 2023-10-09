include ../../../Makefile.def

OBJS = PML2D.o PML3D.o \
       pml_2d.o pml_3d.o \
	   PML2D_3.o PML2D_5.o PML2D_12.o \
	   PML2DVISCOUS.o

all: $(OBJS) 

# fortran functions for 3D material model
pml_3d.o: pml_3d.f
	$(FC) -fPIC -c pml_3d.f


# fortran functions for 2D material model
pml_2d.o: pml_2d.f
	$(FC) -fPIC -c pml_2d.f

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
