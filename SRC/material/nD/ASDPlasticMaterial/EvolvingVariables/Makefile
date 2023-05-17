OPS_BASE=../../../../..
EIGENDIR=$(OPS_BASE)/OTHER/eigenAPI
OPS_SRC=$(OPS_BASE)/SRC

CXX=g++
CXXFLAGS=  \
   -I$(EIGENDIR)\
   -I$(EIGENDIR)/eigen\
   -I$(OPS_SRC)\
   -I$(OPS_SRC)/matrix\
   -I$(OPS_SRC)/handler\
   -I$(OPS_SRC)/actor\
   -I$(OPS_SRC)/actor/actor\
   -I$(OPS_SRC)/actor/channel

DEPS = 
OBJ = test_variables.o
LINKLIBS = /home/jaabell/Repositories/OpenSees/build/libOpenSeesLIB.a -lopenblas

all: test_variables

# %.o: %.cpp 
# 	$(CXX) $(CXXFLAGS) -c -o $@ $^

test_variables: test_variables.cpp 
	$(CXX) $(CXXFLAGS) $(DEPS) -o $@ $^ $(CFLAGS) $(LINKLIBS)
	./test_variables


clean:
	rm *.o test_variables
