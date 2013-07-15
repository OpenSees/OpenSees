include ../../Makefile.def

all: $(TRIANGLE_LIBRARY)

tri = triangle.o

triangle.o: triangle.c
	$(CC) -O -c triangle.c -o triangle.o

$(TRIANGLE_LIBRARY): $(tri)
	   @$(ECHO) "creating $(TRIANGLE_LIBRARY)" 
	@$(RM) $(RMFLAGS) $(TRIANGLE_LIBRARY)
	@$(AR) $(ARFLAGS) $(TRIANGLE_LIBRARY) $(tri)
	@$(RANLIB) $(TRIANGLE_LIBRARY)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*#  core

clean: tidy
	@$(RM) $(RMFLAGS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
