#  Written: Boris Jeremic jeremic@ucdavis.edu 
#           Zhaohui Yang  zhyang@ucdavis.edu  

#  Created: Aug 2000



# #LD = /usr/local/gnu/bin/ld
# MACHINE = LINUX
# 
# #MPI_INCLUDE_DIR = -I/home/zhyang/mpich/include
# MPI_INCLUDE_DIR = -I/usr/local/mpi/mpich/include
# #zhaohui
# MPI_LIB_DIR     = -L/usr/local/mpi/mpich/lib
# #MPI_LIB_DIR     = -L/home/zhyang/mpich/build/LINUX/ch_p4/lib  
# 

include ../../Makefile.def

3D_El           = 3D_EL
3D_El_Pl        = 3D_EL_PL

all:              $(3D_El)
3D_E:            $(3D_El)
3D_EP:         $(3D_El_Pl)
             

$(3D_El):  main_ElasticIsotropic3D.o
	$(LINKER) $(LINKFLAGS) main_ElasticIsotropic3D.o \
	$(HOME)/FEI/SRC/element/8nbrick/EightNodeBrick.o\
	$(HOME)/FEI/SRC/material/nD/ElasticIsotropic3D.o\
	$(HOME)/FEI/SRC/material/nD/ElasticIsotropicMaterial.o\
	$(HOME)/FEI/SRC/analysis/analysis/Analysis.o\
	$(HOME)/FEI/SRC/analysis/analysis/StaticAnalysis.o\
	$(HOME)/FEI/SRC/analysis/algorithm/equiSolnAlgo/Linear.o\
	$(HOME)/FEI/SRC/analysis/integrator/IncrementalIntegrator.o\
	$(HOME)/FEI/SRC/analysis/integrator/Integrator.o\
	$(HOME)/FEI/SRC/analysis/integrator/StaticIntegrator.o\
	$(HOME)/FEI/SRC/analysis/fe_ele/FE_Element.o\
	$(HOME)/FEI/SRC/material/nD/NDMaterial.o\
	$(HOME)/FEI/SRC/domain/node/Node.o\
	$(MATMOD_LIBRARY) $(BJMISC_LIBRARY) $(NDARRAY_LIBRARY) $(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	-o $(3D_El)
        
$(3D_El_Pl):  main_Template3Dep.o
	$(LINKER) $(LINKFLAGS) main_Template3Dep.o \
	$(HOME)/FEI/SRC/element/8nbrick/EightNodeBrick.o\
	$(HOME)/FEI/SRC/analysis/analysis/Analysis.o\
	$(HOME)/FEI/SRC/analysis/analysis/StaticAnalysis.o\
	$(HOME)/FEI/SRC/analysis/algorithm/equiSolnAlgo/Linear.o\
	$(HOME)/FEI/SRC/analysis/integrator/IncrementalIntegrator.o\
	$(HOME)/FEI/SRC/analysis/integrator/Integrator.o\
	$(HOME)/FEI/SRC/analysis/integrator/StaticIntegrator.o\
	$(HOME)/FEI/SRC/analysis/fe_ele/FE_Element.o\
	$(HOME)/FEI/SRC/material/nD/NDMaterial.o\
	$(HOME)/FEI/SRC/domain/node/Node.o\
	$(MATMOD_LIBRARY) $(BJMISC_LIBRARY) $(NDARRAY_LIBRARY) $(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	-o $(3D_El_Pl)

################### TARGETS ########################
#all: $(OBJS) 

#libmatmod.a: $(OBJS)
#	$(AR) -rv $(MATMOD_LIBRARY) $(OBJS)




# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core test

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core 
	@$(RM) $(RMFLAGS) 3D_EL 3D_EL_PL 

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.

