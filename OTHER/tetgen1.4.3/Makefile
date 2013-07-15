include ../../Makefile.def

all: $(TETGEN_LIBRARY)

obj = tetgen.o predicates.o

tetgen.o: tetgen.cxx
	$(CC++) -O2 -c tetgen.cxx -o tetgen.o

predicates.o: predicates.cxx
	$(CC++) -O -c predicates.cxx -o predicates.o

$(TETGEN_LIBRARY): $(obj)
	@$(ECHO) "creating $(TETGEN_LIBRARY)" 
	@$(RM) $(RMFLAGS) $(TETGEN_LIBRARY)
	@$(AR) $(ARFLAGS) $(TETGEN_LIBRARY) $(obj)
	@$(RANLIB) $(TETGEN_LIBRARY)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*#  core

clean: tidy
	@$(RM) $(RMFLAGS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
