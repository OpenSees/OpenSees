include ../../Makefile.def

ifeq ($(GRAPHICS), UsingOpenGL)

OBJS       = Renderer.o \
	DofColorMap.o PlainMap.o \
	OpenGlDevice.o OpenGlRenderer.o
else

OBJS       = Renderer.o WindowRenderer.o View.o Projection.o Clipping.o \
	Scan.o Viewport.o X11Renderer.o gMatrix.o db.o Device.o \
	VrmlViewer.o DofColorMap.o PlainMap.o X11Device.o WindowDevice.o 

endif

# Compilation control
all:         $(OBJS)


test:  main.o $(OBJS)
	$(LINKER) $(LINKFLAGS) main.o $(OBJS) $(FE_LIBRARY) \
	$(MACHINE_LINKLIBS) $(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY) \
	$(MACHINE_SPECIFIC_LIBS) -o go

triangle:  mainTriangle.o $(OBJS)
	$(LINKER) $(LINKFLAGS) mainTriangle.o $(OBJS) $(FE_LIBRARY) \
	$(MACHINE_LINKLIBS) $(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY) \
	$(MACHINE_SPECIFIC_LIBS) -o go

index:  index.o 
	$(LINKER) $(LINKFLAGS) index.o $(OBJS) $(FE_LIBRARY) \
	$(MACHINE_LINKLIBS) $(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY) \
	$(MACHINE_SPECIFIC_LIBS) $(FE_INCLUDES) -lglut -lpng -o go


# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core go

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o test go

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



