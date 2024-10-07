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
OBJ = test_elast.o
LINKLIBS = /home/jaabell/Repositories/OpenSees/build/libOpenSeesLIB.a -lopenblas

all: test_elast

# %.o: %.cpp 
# 	$(CXX) $(CXXFLAGS) -c -o $@ $^

test_elast: test_elast.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(CFLAGS) $(LINKLIBS)
	./test_elast


clean:
	rm *.o test_elast