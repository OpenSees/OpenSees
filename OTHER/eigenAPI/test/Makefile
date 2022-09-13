EIGENDIR=../eigen

CXX=g++
CXXFLAGS=-I$(EIGENDIR)
DEPS = 
OBJ = test01.o 

%.o: %.c $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ 
test01: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS)

clean:
	rm *.o test??