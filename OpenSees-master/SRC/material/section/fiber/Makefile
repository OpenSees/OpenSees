include ../../../../Makefile.def

OBJS       = Fiber.o \
	UniaxialFiber2d.o \
	UniaxialFiber3d.o \
	NDFiber2d.o \
	NDFiber3d.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.








