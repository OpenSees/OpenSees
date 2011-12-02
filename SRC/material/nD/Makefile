include ../../../Makefile.def

OBJS       = NDMaterial.o \
             ElasticIsotropic3D.o \
             ElasticIsotropicMaterial.o \
             ElasticIsotropicPlaneStrain2D.o \
             ElasticIsotropicPlaneStress2D.o \
             ElasticIsotropicPlateFiber.o \
             ElasticIsotropicAxiSymm.o\
             TclModelBuilderNDMaterialCommand.o \
             J2PlaneStrain.o \
             PlateFiberMaterial.o \
             PlaneStressMaterial.o \
             J2PlaneStress.o \
             J2Plasticity.o \
             J2AxiSymm.o \
             J2ThreeDimensional.o \
             J2PlateFiber.o \
             PressureDependentElastic3D.o \
             ElasticIsotropicBeamFiber.o \
             BeamFiberMaterial.o \
             FeapMaterial.o \
             ElasticCrossAnisotropic.o 


TEMPLATE_OBJS =              \
                MatPoint3D.o \
                Template3Dep.o \
                DP_PS.o \
                DP_YS.o \
                DP_YS01.o \
                EL_LEeq.o \
                EL_LEij.o \
                EL_NLEeq.o \
                EL_NLEij.o \
                EL_NLEijMD.o \
                EL_S.o \
                EL_T.o \
                EPState.o \
                MD_PS.o \
                MD_PS01.o \
                MD_YS.o \
                VM_PS.o \
                VM_YS.o \
                YS.o \
                PS.o \
                RMC01.o \
                RMC01_YS.o \
                RMC01_PS.o


all: $(OBJS) 
	@$(CD) $(FE)/material/nD/soil; $(MAKE);
	@$(CD) $(FE)/material/nD/feap; $(MAKE);
	@$(CD) $(FE)/material/nD/cyclicSoil; $(MAKE);
	@$(CD) $(FE)/material/nD/Template3Dep; $(MAKE);  $(AR) -rv $(MATMOD_LIBRARY) $(TEMPLATE_OBJS)
	@$(CD) $(FE)/material/nD/FiniteDeformation; $(MAKE); 

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o
	@$(RM) $(RMFLAGS) .inslog*
	@$(RM) $(RMFLAGS) tca.map

spotless: clean
	@$(CD) $(FE)/material/nD/soil; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/feap; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/cyclicSoil; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/Template3Dep; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/FiniteDeformation; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



