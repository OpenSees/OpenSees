include ../../../../Makefile.def

OBJS = stressDensity.o SDM-UC.o \
	sdm3d.o sdm2d.o MODULE_ALL_INTERFACES_3D.o MODULE_ALL_INTERFACES_2D.o \
	StressDensityModel.o StressDensityModel2D.o StressDensityModel3D.o

all: $(OBJS) 

# fortran functions for 3D material model
sdm3d.o: SDM3D.f90 module_all_interfaces_3d.mod
	$(FC) -fPIC -c SDM3D.f90

module_all_interfaces_3d.mod: MODULE_ALL_INTERFACES_3D.f90 
	$(FC) -fPIC -c MODULE_ALL_INTERFACES_3D.f90

# fortran functions for 2D material model
sdm2d.o: SDM2D.f90 module_all_interfaces_2d.mod
	$(FC) -fPIC -c SDM2D.f90

module_all_interfaces_2d.mod: MODULE_ALL_INTERFACES_2D.f90 
	$(FC) -fPIC -c MODULE_ALL_INTERFACES_2D.f90

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o *.mod

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
