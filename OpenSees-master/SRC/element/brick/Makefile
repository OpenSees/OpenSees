include ../../../Makefile.def

OBJS       = Brick.o \
	Twenty_Node_Brick.o \
	shp3d.o \
	BbarBrick.o \
	BbarBrickWithSensitivity.o \
	TclBrickCommand.o \
	TclTwenty_Node_BrickCommand.o


all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
