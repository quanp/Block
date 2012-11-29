CXX=mpic++
CXXFLAGS=-O3 -fPIC

BOOSTDIR=/home/naokin/boost/1.52.0
BOOSTINCLUDE=$(BOOSTDIR)/include
BOOSTLIB=-L$(BOOSTDIR)/lib

INCLUDES= -I. -I$(BOOSTINCLUDE)
LIBS=$(BOOSTLIB) -lboost_serialization -lboost_mpi

FLAGS=$(CXXFLAGS)

OBJ_gaopt = main.o GAOptimize.o

.C.o:
	$(CXX) $(FLAGS) $(INCLUDES) $(LIBS) -c $*.C

all     : gaopt

gaopt: $(OBJ_gaopt)
	$(CXX) $(FLAGS) -o $@ $(OBJ_gaopt) $(LIBS)

clean:
	rm gaopt; rm *.o; rm *.x; rm *.a
#
