include ../../../../Makefile.def

OBJS =  MaterialParameter.o \
	ElasticState.o \
	Isotropic_Elastic.o \
	elnp_Elastic.o \
	PressureDependent_Elastic.o \
	DM04_Elastic.o \
	YieldFunction.o \
	VM_YF.o \
	DP_YF.o \
	CC_YF.o \
	DM04_YF.o \
	PlasticFlow.o \
	VM_PF.o \
	DP_PF.o \
	CC_PF.o \
	DM04_PF.o \
	ScalarEvolution.o \
	Linear_Eeq.o \
	CC_Ev.o \
	TensorEvolution.o \
	Linear_Eij.o \
	AF_Eij.o \
	DM04_alpha_Eij.o \
	DM04_z_Eij.o \
	NewTemplate3Dep.o \
	TclNewTemplate3Dep.o  
   

all: $(OBJS)


FT3Dep_test: $(OBJS) FT3Dep_test.o
	$(LINKER) $(LINKFLAGS) FT3Dep_test.o $(OBJS) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(TCL_LIBRARY)\
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	-o FT3Dep_test

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS)  *.o *.bak
	
spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.




