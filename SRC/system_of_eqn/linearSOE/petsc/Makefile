# Makefile for fe objects

include ../../../Makefile.def


PROGRAM       = test


OBJS       = PetscSOE.o PetscSolver.o ShadowPetscSOE.o ActorPetscSOE.o

all:    $(OBJS)

test:    $(PROGRAM) main.o
	$(LINKER) $(LINKFLAGS) $(OBJS) main.o  $(MATRIX_LIBS)  \
        $(SysOfEqn_LIBS) -liostream -lm -ldxml -o $(PROGRAM)


MATRIX_LIBS   = $(HOME)/OOP/matrix/Matrix.o \
        $(HOME)/OOP/matrix/Vector.o \
        $(HOME)/OOP/matrix/ID.o               

SysOfEqn_LIBS = $(FE)/system_of_eqn/SystemOfEqn.o \
		$(FE)/system_of_eqn/linearSOE/LinearSOE.o \
		$(FE)/system_of_eqn/linearSOE/DomainSolver.o \
		$(FE)/system_of_eqn/linearSOE/LinearSOESolver.o \
		$(FE)/system_of_eqn/Solver.o \
		$(FE)/actor/actor/MovableObject.o 

.C.o:
	$(CC++) $(FLAGS) $(INCLUDES) $(PETSC_FLAGS) -c $< -o $@
.c.o:
	$(CC) $(FLAGS) $(INCLUDES) -c $< -o $@
.f.o:      
	/usr/bin/f77 -O -c $< -o $@

# Miscellaneous

tidy:	
	@rm -f Makefile.bak *~ #*# core

clean: tidy
	@rm -f $(OBJS) *.o

spotless: clean
	@rm -f $(PROGRAM) 

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.


