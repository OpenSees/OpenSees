include ../../../Makefile.def

OBJS = PFEMElement2D.o PFEMElement3D.o TclModelBuilder_addPFEMElement.o PFEMMesher2D.o PFEMMesher3D.o TclPFEMCommands.o PFEMElement2DCompressible.o PFEMElement2DBubble.o tetgen.o predicates.o triangle.o

all: $(OBJS)

tetgen.o: $(FE)/../OTHER/Tetgen/tetgen.cxx
	$(CC++) -O2 -c $(FE)/../OTHER/Tetgen/tetgen.cxx -o tetgen.o

predicates.o: $(FE)/../OTHER/Tetgen/predicates.cxx
	$(CC++) -O -c $(FE)/../OTHER/Tetgen/predicates.cxx -o predicates.o

triangle.o: $(FE)/../OTHER/Triangle/triangle.c
	$(CC) -O -c $(FE)/../OTHER/Triangle/triangle.c -o triangle.o

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o *.so

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
