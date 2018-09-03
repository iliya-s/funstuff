USE_MPI = no
USE_INTEL = yes
BOOST = /home/iliya/apps/boost/boost_1_68_0/
EIGEN = /home/iliya/apps/eigen/eigen-eigen-b3f3d4950030

FLAGS = -g -Wall -I${BOOST} -I${EIGEN}

ifeq ($(USE_INTEL), yes)
	ifeq ($(USE_MPI), yes)
		CXX = mpiicpc
		LFLAGS = -L${BOOST}/stage/lib -lboost_serialization -lboost_mpi
	else
		CXX = icpc
		LFLAGS = -L${BOOST}/stage/lib -lboost_serialization
		FLAGS += -DSERIAL
	endif
else
	ifeq ($(USE_MPI), yes)
		CXX = mpicxx
		LFLAGS = -L${BOOST}/stage/lib -lboost_serialization -lboost_mpi
	else
		CXX = g++
		LFLAGS = -L${BOOST}/stage/lib -lboost_serialization
		FLAGS += -DSERIAL
	endif
endif

HOSTNAME := $(shell hostname)

SRC_main = main.cpp ao_integrals.cpp

OBJ_main = obj/main.o obj/ao_integrals.o

obj/%.o: %.cpp
	$(CXX) $(FLAGS) $(OPT) -c $< -o $@

all: main

main : $(OBJ_main)
	$(CXX) $(FLAGS) $(OPT) -o main $(OBJ_main) $(LFLAGS)

clean :
	find . -name "*.o"|xargs rm 2>/dev/null; rm -f main >/dev/null 2>&1
