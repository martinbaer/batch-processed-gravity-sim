MPICXX?=mpic++
CXX?=g++
CXXFLAGS?=-std=c++11 -O2
# -g is for debugging
NVCC?=nvcc
NVFLAGS?=-O2 --gpu-architecture=sm_35 -Wno-deprecated-gpu-targets

TARGETS = pp_serial bh_serial # bh_mpi test_mpi benchmarking_bh_serial benchmarking_bh_mpi

default: 
	make all

# dependancies
$(TARGETS): helpers.o helpers.h bh_tree.o bh_tree.h

helpers.o: helpers.h

bh_tree.o: bh_tree.h helpers.h

# wildcard rules
%_mpi.o : %_mpi.cpp
	$(MPICXX) $(CXXFLAGS) $(CFLAGS_$(basename $<)) -c $< -o $@

%_mpi : %_mpi.cpp
	$(MPICXX) $(CXXFLAGS) $(filter %.o %.cpp, $^) -o $@

%.o: %.cu
	$(NVCC) $(NVFLAGS) $(NVFLAGS_$(basename $<)) -c $< -o $@

%: %.cu
	$(NVCC) $(NVFLAGS) $(filter %.o %.cu, $^) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CFLAGS_$(basename $<)) -c $< -o $@

%: %.cpp
	$(CXX) $(CXXFLAGS) $(filter %.o %.cpp, $^) -o $@

# rules
all: $(TARGETS)

clean:
	rm -f $(TARGETS) *.o