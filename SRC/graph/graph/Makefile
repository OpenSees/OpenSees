include ../../../Makefile.def

OBJS       = DOF_Graph.o Vertex.o Graph.o \
	DOF_GroupGraph.o  VertexIter.o


all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
