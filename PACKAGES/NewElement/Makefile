include ../../Makefile.def

LIBRARY         = MyTruss.so

all:    $(LIBRARY)

OBJS = MyTruss.o MyTrussPackage.o 

$(LIBRARY): $(OBJS)
	@$(LINKER) $(LINK_SHARED_FLAGS)  $(OBJS) -o $(LIBRARY)


.cpp.o:
	@$(ECHO) Making $@ from $<
	$(CC++) $(C++SHARED_FLAGS) $(INCLUDES) -c $< -o $@

.C.o:
	@$(ECHO) Making $@ from $<
	$(CC++) $(C++SHARED_FLAGS) $(INCLUDES) -c $< -o $@

.c.o:
	@$(ECHO) Making $@ from $<
	$(CC) $(C_SHARED_FLAGS) -c $<

.f.o:      
	@$(ECHO) Making $@ from $<
	$(FC) $(F_SHARED_FLAGS) -c $< -o $@

.f77.o:      
	@$(ECHO) Making $@ from $<
	$(FC) $(F_SHARED_FLAGS) -c $< -o $@

.f90.o:      
	@$(ECHO) Making $@ from $<
	$(FC90) $(F_SHARED_FLAGS) -c $< -o $@

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core $(LIBRARY)

# DO NOT DELETE THIS LINE -- make depend depends on it.
