include ../../../../Makefile.def

OBJS = YieldSurfaceSection2d.o \
	   YS_Section2D01.o \
	   YS_Section2D02.o \
           SoilFootingSection2d.o \
	   TclModelBuilderYS_SectionCommand.o

all:         $(OBJS)

test:	$(OBJS) main.o
	$(LINKER) $(LINKFLAGS) main.o $(OBJS) $(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	-o example

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core example main.o

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
