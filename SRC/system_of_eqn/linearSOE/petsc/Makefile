# Makefile for fe objects

include ../../../../Makefile.def

PROGRAM       = test

OBJS       = PetscSOE.o PetscSolver.o PetscSparseSeqSolver.o 

all:    $(OBJS)

LIB_OBJS       = PetscSOE.o PetscSolver.o PetscSparseSeqSolver.o TclPetsc.o SparseGenRowLinSOE.o SparseGenRowLinSolver.o
PETSC_LIB = Petsc.so
PETSC_INC = -I/usr/local/petsc-2.2.1 \
	-I/usr/local/petsc-2.2.1/bmake/linux-gnu \
	-I/usr/local/petsc-2.2.1/include \
        -I/usr/local/petsc-2.2.1/include/mpiuni \
	-DPETSC_USE_BOPT_g -DPETSC_USE_EXTERN_CXX -D_PETSC

lib:
	$(CC++) -fPIC $(INCLUDES) $(PETSC_INC)  -c -Wall PetscSOE.cpp	
	$(CC++) -fPIC $(INCLUDES) $(PETSC_INC)  -c -Wall PetscSolver.cpp	
	$(CC++) -fPIC $(INCLUDES) $(PETSC_INC)  -c -Wall PetscSparseSeqSolver.cpp	
	$(CC++) -fPIC $(INCLUDES) $(PETSC_INC)  -c -Wall TclPetsc.cpp	
	$(CC++) -fPIC $(INCLUDES) $(PETSC_INC)  -c -Wall $(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenRowLinSOE.cpp
	$(CC++) -fPIC $(INCLUDES) $(PETSC_INC)  -c -Wall $(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenRowLinSolver.cpp
	$(CC++) -shared -Wl,-soname,$(PETSC_LIB) -o $(PETSC_LIB) $(LIB_OBJS) \
	  -L/usr/local/petsc-2.2.1/lib/libg/linux-gnu \
	 -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc \
	-rdynamic -Wl,-rpath,/usr/local/petsc-2.2.1/lib/libg/linux-gnu 

test:    main.o $(OBJS)
	$(LINKER) $(LINKFLAGS) $(OBJS) main.o  $(MATRIX_LIBS) $(ERROR_LIBS) $(GRAPH_LIBS)\
        $(SysOfEqn_LIBS) $(MACHINE_NUMERICAL_LIBS) $(PETSC_LIB) -o $(PROGRAM)

ERROR_LIBS = $(FE)/handler/StandardStream.o \
	$(FE)/handler/FileStream.o \
	$(FE)/handler/OPS_Stream.o

GRAPH_LIBS   = $(FE)/graph/graph/Graph.o \
	$(FE)/graph/graph/Vertex.o \
	$(FE)/graph/graph/VertexIter.o \
	$(FE)/tagged/TaggedObject.o \
	$(FE)/tagged/storage/ArrayOfTaggedObjects.o \
	$(FE)/tagged/storage/ArrayOfTaggedObjectsIter.o 	

MATRIX_LIBS   = $(FE)/matrix/Matrix.o \
        $(FE)/matrix/Vector.o \
        $(FE)/matrix/ID.o \
	$(FE)/nDarray/basics.o \
	$(FE)/nDarray/nDarray.o \
	$(FE)/nDarray/BJmatrix.o \
	$(FE)/nDarray/BJvector.o \
	$(FE)/nDarray/BJtensor.o \
	$(FE)/nDarray/straint.o  \
	$(FE)/nDarray/stresst.o

SysOfEqn_LIBS = $(FE)/system_of_eqn/SystemOfEqn.o \
		$(FE)/system_of_eqn/linearSOE/LinearSOE.o \
		$(FE)/system_of_eqn/linearSOE/DomainSolver.o \
		$(FE)/system_of_eqn/linearSOE/LinearSOESolver.o \
		$(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenRowLinSOE.o \
		$(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenRowLinSolver.o \
		$(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenColLinSOE.o \
		$(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenColLinSolver.o \
		$(FE)/system_of_eqn/linearSOE/sparseGEN/SuperLU.o \
		$(FE)/system_of_eqn/Solver.o \
		$(FE)/actor/actor/MovableObject.o 

# Miscellaneous

tidy:	
	@rm -f Makefile.bak *~ #*# core

clean: tidy
	@rm -f $(OBJS) *.o

spotless: clean
	@rm -f $(PROGRAM) 

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.


