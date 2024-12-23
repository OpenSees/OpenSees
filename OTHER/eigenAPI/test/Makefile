EIGENDIR=../eigen
OPS_SRC=../../../SRC

CXX=g++
CXXFLAGS=\
   -I$(EIGENDIR)\
   -I$(OPS_SRC)\
   -I$(OPS_SRC)/matrix\
   -I$(OPS_SRC)/handler\
   -I$(OPS_SRC)/actor\
   -I$(OPS_SRC)/actor/actor
DEPS = 
OBJ = test01.o test02.o
LINKLIBS = ../../../lib/libOpenSees.a -lopenblas

all: 

%.o: %.c $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ 

test01: test01.o
	$(CXX) -o $@ $^ $(CFLAGS) $(LINKLIBS)
	./test01

test02: test02.o
	$(CXX)  -o $@ $^ $(CFLAGS) $(LINKLIBS)
	./test02

test03: test03.o
	$(CXX)  -o $@ $^ $(CFLAGS) $(LINKLIBS)
	./test03

test04: test04.o
	$(CXX)  -o $@ $^ $(CFLAGS) $(LINKLIBS)
	./test04


test05: test05.o
	$(CXX)  -o $@ $^ $(CFLAGS) $(LINKLIBS)
	./test05

test06: test06.o
	$(CXX)  -o $@ $^ $(CFLAGS) $(LINKLIBS)
	./test06

clean:
	rm *.o test??