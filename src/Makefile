CC=gcc -std=c99
CXX=g++ -std=c++11
#CC=gcc-4.9 -std=c99
#CXX=g++-4.9 -std=c++11

LDFLAGS=-lpng -ltiff -ljpeg -lm

CFLAGS=-O3 
#CFLAGS =-O3 -DNDEBUG
#CFLAGS=-Wall -Wextra -Werror -g
CXXFLAGS=$(CFLAGS)
OMPFLAGS=-fopenmp

OBJECTS_C=iio.o mask.o xmalloc.o bicubic_interpolation.o elap_recsep.o
OBJECTS_CXX=tvl2_model.o nltv_model.o tvcsad_model.o nltvcsad_model.o tvl2w_model.o nltvcsadw_model.o nltvw_model.o tvcsadw_model.o aux_energy_model.o energy_model.o heuristic_interpolation.o
PROGRAMS=sparse_flow local_faldoi global_faldoi

all: $(PROGRAMS)

%.o : %.c
	$(CC) $(CFLAGS) -c $^

%.o : %.cpp
	$(CXX) $(CFLAGS) $(OMPFLAGS) -c $^

% : %.c $(OBJECTS_C)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) $(OMPFLAGS) -o $@

% : %.cpp $(OBJECTS_C) $(OBJECTS_CXX)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(OMPFLAGS) -o $@

.PHONY:
clean:
	rm -f $(PROGRAMS) $(OBJECTS_C) $(OBJECTS_CXX)

.SECONDARY:
.SUFFIXES:

