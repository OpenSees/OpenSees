OPS_BASE=../../../../..
EIGENDIR=$(OPS_BASE)/OTHER/eigenAPI
OPS_SRC=$(OPS_BASE)/SRC

CXX=g++
CXXFLAGS=\
   -I$(EIGENDIR)\
   -I$(EIGENDIR)/eigen\
   -I$(OPS_SRC)\
   -I$(OPS_SRC)/matrix\
   -I$(OPS_SRC)/handler\
   -I$(OPS_SRC)/actor\
   -I$(OPS_SRC)/actor/actor\
   -I$(OPS_SRC)/actor/channel\
   -I$(OPS_SRC)/domain \
   -I$(OPS_SRC)/domain/component \
   -I$(OPS_SRC)/tagged \
   -I$(OPS_SRC)/material \
   -I$(OPS_SRC)/material/nD

DEPS = 
OBJ = 
OPS_LIB=$(OPS_BASE)/build
LINKLIBS = $(OPS_LIB)/libOpenSeesLIB.a \
   $(OPS_LIB)/libOpenSeesLIB.a \
   $(OPS_LIB)/libG3.a \
   $(OPS_LIB)/libOPS_InterpTcl.a \
   $(OPS_LIB)/SRC/coordTransformation/libcoordTransformation.a \
   $(OPS_LIB)/SRC/material/uniaxial/libOPS_Material_f.a \
   $(OPS_LIB)/SRC/material/uniaxial/drain/libOPS_Material_Uniaxial_Drain_f.a \
   $(OPS_LIB)/SRC/material/nD/feap/libOPS_Material_nD_Feap_f.a \
   $(OPS_LIB)/SRC/damping/libdamping.a \
   $(OPS_LIB)/SRC/system_of_eqn/linearSOE/sparseSYM/libOPS_SysOfEqn_f.a \
   $(OPS_LIB)/OTHER/ARPACK/libARPACK.a \
   $(OPS_LIB)/OTHER/CSPARSE/libCSPARSE.a \
   $(OPS_LIB)/OTHER/SuperLU_5.1.1/libSUPERLU.a \
   $(OPS_LIB)/OTHER/UMFPACK/libUMFPACK.a \
   $(OPS_LIB)/OTHER/tetgen1.4.3/libtet.a \
   $(OPS_LIB)/OTHER/Triangle/libtriangle.a \
   $(OPS_LIB)/OTHER/AMD/libAMD.a \
   -lopenblas -lm -ldl -lgfortran -lquadmath 


all: test_00_VonMises


test_00_VonMises: test_00_VonMises.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(CPP_OTHER_OBJ) -o $@ $^ $(CFLAGS) $(LINKLIBS) 
	./test_00_VonMises


clean:
	rm *.o test_00_VonMises_YF