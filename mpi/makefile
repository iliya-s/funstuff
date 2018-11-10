USE_MPI = yes
USE_INTEL = yes
BOOST = /home/iliya/apps/boost/boost_1_68_0

FLAGS = -g -Wall -I${BOOST}

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

SRC_eulers_num = eulers_num.cpp stats.cpp

SRC_pi = pi.cpp

OBJ_eulers_num = obj/eulers_num.o obj/stats.o

OBJ_pi = obj/pi.o

obj/%.o: %.cpp
	$(CXX) $(FLAGS) $(OPT) -c $< -o $@

all: eulers_num pi

eulers_num  : $(OBJ_eulers_num)
	$(CXX) $(FLAGS) $(OPT) -o eulers_num $(OBJ_eulers_num) $(LFLAGS)

pi : $(OBJ_pi)
	$(CXX) $(FLAGS) $(OPT) -o pi $(OBJ_pi) $(LFLAGS)

clean :
	find . -name "*.o"|xargs rm 2>/dev/null; rm -f eulers_num pi >/dev/null 2>&1
