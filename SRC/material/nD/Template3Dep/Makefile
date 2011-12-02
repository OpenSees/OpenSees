#  Written: Boris Jeremic jeremic@ucdavis.edu
#  Created: Aug2000


include ../../../../Makefile.def


OBJS  = Template3Dep.o \
        CDriver.o \
        DP_PS.o	  \
        DP_YS.o	  \
        EL_LEeq.o \
        EL_LEij.o \
        EL_NLEeq.o\
        EL_NLEij.o\
        EL_S.o	  \
        EL_T.o	  \
        EPState.o \
        VM_PS.o	  \
        VM_YS.o	  \
        YS.o	  \
        bmodel00.o\
        bmodel01.o\
        bmodel02.o\
        dpmodel0.o\
        dr_pr.o	  \
        elastic.o \
        epstate.o \
        ma_da.o	  \
        mmodel.o  \
        mrslade.o \
        parabol.o \
        plastbas.o\
        vmises.o  

#        MD_PS.o	  \
#        MD_YS.o	  \
#        MD_EL.o	  \
#        bmodel.o  \
#        ldemodel.o\
#        logarit.o \
#        mmodelld.o\
#        neohook.o \
#        ogden.o	  \
#        simopist.o\


################### TARGETS ########################
all: $(OBJS) 

libmatmod.a: $(OBJS)
	$(AR) -rv $(MATMOD_LIBRARY) $(OBJS)


# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core test

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core 
	@$(RM) $(RMFLAGS) math_tst elast_tst 

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.

