USE_INTEL = yes
EIGEN = /home/iliya/apps/eigen/eigen-eigen-b3f3d4950030

FLAGS = -g -Wall -I${EIGEN}

ifeq ($(USE_INTEL), yes)
		CXX = icpc
else
		CXX = g++
endif

HOSTNAME := $(shell hostname)

SRC_main = main.cpp

OBJ_main = obj/main.o

obj/%.o: %.cpp
	$(CXX) $(FLAGS) $(OPT) -c $< -o $@

all: main

main : $(OBJ_main)
	$(CXX) $(FLAGS) $(OPT) -o main $(OBJ_main) $(LFLAGS)

clean :
	find . -name "*.o"|xargs rm 2>/dev/null; rm -f main >/dev/null 2>&1
