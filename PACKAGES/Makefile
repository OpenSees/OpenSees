include ../Makefile.def

OS_FLAG ?= -D_LINUX

ifeq ($(OS_FLAG), -D_LINUX)
all: 
	@$(CD) $(FE)/../PACKAGES/NewElement; $(MAKE);
	@$(CD) $(FE)/../PACKAGES/NewMaterial; $(MAKE);
else
all: dynamicLoadFunctions.o
	@$(CD) $(FE)/../PACKAGES/NewElement; $(MAKE);
	@$(CD) $(FE)/../PACKAGES/NewMaterial; $(MAKE);

dynamicLoadFunctons.o: dynamicLoadFunctions.cpp
	$(CC++) -dynamic $(INCLUDES)  -c -Wall dynamicLoadFunctions.cpp

endif


# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(CD) $(FE)/../PACKAGES/NewElement; $(MAKE) wipe;
	@$(CD) $(FE)/../PACKAGES/NewMaterial; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
