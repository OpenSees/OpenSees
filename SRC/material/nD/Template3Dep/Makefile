#  Written: Boris Jeremic jeremic@ucdavis.edu
#  Created: Aug2000



# #LD = /usr/local/gnu/bin/ld
# MACHINE = LINUX
# 
# #MPI_INCLUDE_DIR = -I/home/zhyang/mpich/include
# MPI_INCLUDE_DIR = -I/usr/local/mpi/mpich/include
# #zhaohui
# MPI_LIB_DIR     = -L/usr/local/mpi/mpich/lib
# #MPI_LIB_DIR     = -L/home/zhyang/mpich/build/LINUX/ch_p4/lib  
# 


include ../../../../Makefile.def



OBJS  = Template3Dep.o \
        TclTemplate3DepCommand.o \
        MatPoint3D.o \
        CAM_PS.o  \
        CAM_YS.o \
        DP_PS.o  \
        DP_YS.o \
        EL_LEeq.o \
        EL_LEij.o \
        EL_NLEp.o\
        EL_NLEeq.o\
        EL_NLEij.o\
        EL_NLEijMD.o\
        EL_S.o	  \
        EL_T.o	  \
        EPState.o \
        VM_PS.o	  \
        VM_YS.o	  \
        YS.o	  \
        MD_PS01.o	  \
        MD_PS.o	  \
        DP_YS01.o	  \
        MD_YS.o	  \
        Tri_a_fail_crit_YS.o \
                RMC01.o      \
                RMC01_PS.o   \
                RMC01_YS.o   

################### TARGETS ########################
all: $(OBJS) 
	$(AR) -rv $(MATMOD_LIBRARY) $(OBJS)

libmatmod.a: $(OBJS)
	$(AR) -rv $(MATMOD_LIBRARY) $(OBJS)




# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core test

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core 
	@$(RM) $(RMFLAGS) math_tst 
	@$(RM) $(RMFLAGS) elast_tst 
	@$(RM) $(RMFLAGS) .inslog*
	@$(RM) $(RMFLAGS) tca.map



spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.

