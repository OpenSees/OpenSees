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
   -I$(OPS_SRC)/actor/channel

DEPS = 
OBJ = 
LINKLIBS = /home/jaabell/Repositories/OpenSees/build/libOpenSeesLIB.a -lopenblas


all: test_00_VonMises_YF


test_00_VonMises_YF: test_00_VonMises_YF.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(CPP_OTHER_OBJ) -o $@ $^ $(CFLAGS)  $(LINKLIBS)
	./test_00_VonMises_YF


clean:
	rm *.o test_00_VonMises_YF