include ../../../Makefile.def

OBJS = PFEMElement2D.o PFEMElement2DFIC.o PFEMElement3D.o TclModelBuilder_addPFEMElement.o \
       PFEMElement2DCompressible.o PFEMElement2DQuasi.o PFEMElement2DBubble.o PFEMElement3DBubble.o \
       tetgen.o predicates.o triangle.o PFEMElement2Dmini.o TriangleMeshGenerator.o TetMeshGenerator.o \
       BackgroundMesh.o  Particle.o ParticleGroup.o BackgroundDef.o TaylorHood2D.o TriGaussPoints.o \
       HigherOrder.o MINI.o Mesh.o LineMeshGenerator.o LineMesh.o TriMesh.o TetMesh.o \
       QuadMeshGenerator.o QuadMesh.o PFEMContact2D.o

all: $(OBJS)

tetgen.o: $(FE)/../OTHER/Tetgen/tetgen.cxx
	$(CC++) -O3 -fPIC -c $(FE)/../OTHER/Tetgen/tetgen.cxx -o tetgen.o

predicates.o: $(FE)/../OTHER/Tetgen/predicates.cxx
	$(CC++) -O3 -fPIC -c $(FE)/../OTHER/Tetgen/predicates.cxx -o predicates.o

triangle.o: $(FE)/../OTHER/Triangle/triangle.c
	$(CC) -O3 -c -fPIC $(FE)/../OTHER/Triangle/triangle.c -o triangle.o

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o *.so

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
