include ../../../../Makefile.def

OBJS =  MaterialParameter.o \
	ElasticState.o \
	Isotropic_Elastic.o \
	elnp_Elastic.o \
	PressureDependent_Elastic.o \
	DM04_Elastic.o \
	YieldFunction.o \
	VM_YF.o \
	RMC_YF.o \
	DP_YF.o \
	CC_YF.o \
	DM04_YF.o \
	PlasticFlow.o \
	VM_PF.o \
	RMC_PF.o \
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

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS)  *.o *.bak
	
spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.




