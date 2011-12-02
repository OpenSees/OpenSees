include ../../../Makefile.def

OBJS       = DOF_Graph.o Vertex.o Graph.o \
	DOF_GroupGraph.o  VertexIter.o


all:         $(OBJS)

test: $(OBJS) main.o Graph.o Vertex.o
	$(LINKER) main.o $(OBJS) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(FE_LIBRARY) \
	-o graph_tst

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
