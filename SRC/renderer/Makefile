include ../../Makefile.def


ifeq ($(GRAPHIC_FLAG), -D_AGL)

AGL_OBJS = AGL_Device.o

else

AGL_OBJS = 

endif


ifeq ($(GRAPHICS), UsingOpenGL)

OBJS       = Renderer.o \
	DofColorMap.o PlainMap.o \
	$(AGL_OBJS) OpenGlDevice.o OpenGlRenderer.o
else

OBJS       = Renderer.o WindowRenderer.o View.o Projection.o Clipping.o \
	Scan.o Viewport.o X11Renderer.o gMatrix.o db.o Device.o \
	VrmlViewer.o DofColorMap.o PlainMap.o X11Device.o WindowDevice.o 

endif

ifeq ($(GRAPHICS), NONE)

OBJS = PlainMap.o Renderer.o

endif

# Compilation control
all:         $(OBJS)


test:  main.o $(OBJS)
	$(LINKER) $(LINKFLAGS) main.o $(FE_LIBRARY) \
	$(MACHINE_LINKLIBS) $(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY) \
	$(MACHINE_SPECIFIC_LIBS) -o go

test1:  main1.o 
	$(LINKER) $(LINKFLAGS) main1.o \
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


agl:    aglMain.o
	$(LINKER) $(LINKFLAGS) aglMain.o -L/usr/local/lib -framework AGL -framework OpenGL -framework Carbon

agl1:    aglClipBuffer.o
	$(LINKER) $(LINKFLAGS) aglClipBuffer.o -L/usr/local/lib -framework AGL -framework OpenGL -framework Carbon

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core go

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o test go

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



