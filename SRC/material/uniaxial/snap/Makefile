include ../../../../Makefile.def

OBJS       = Pinching.o \
	Bilinear.o \
	Bilin.o \
	Bilin02.o \
	Clough.o \
	CloughHenry.o \
	CloughDamage.o \
	PinchingDamage.o \
	Bilin.o \
	TclSnapMaterialCommand.o

all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
