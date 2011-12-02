include ../../../Makefile.def

OBJS       = TwentyNodeBrick_u_p_U.o  TclTwentyNodeBrick_u_p_U.o \
	EightNodeBrick_u_p_U.o   TclEightNodeBrick_u_p_U.o \
	FourNodeQuadUP.o TclFourNodeQuadUPCommand.o

            #Tclu_p_UCommand.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak 

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) 
	@$(RM) $(RMFLAGS) *.o
	@$(RM) $(RMFLAGS) .inslog*
	@$(RM) $(RMFLAGS) tca.map

spotless: clean					   

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
