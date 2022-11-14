MPICXX?=mpic++
CXX?=g++
CXXFLAGS?=-std=c++11 -O2
# -g is for debugging
NVCC?=nvcc
NVFLAGS?=-O2 --gpu-architecture=sm_35 -Wno-deprecated-gpu-targets

TARGETS = pp_serial bh_serial bh_mpi test_mpi benchmarking_bh_serial benchmarking_bh_mpi

default: 
	./shell_scripts/unhide_objects.sh
	make all
	./shell_scripts/hide_objects.sh

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

# # all targets
# TARGETS = Assignment_serial Assignment_openmp Assignment_mpi Assignment_cuda Assignment_avx

# # The first rule in the Makefile is the default target that is made if 'make' is invoked with
# # no parameters.  'all' is a dummy target that will make everything
# default : all

# ## Dependencies

# # all targets depend on the helper programs
# $(TARGETS) : randutil.h randutil.ipp randutil.o eigensolver.h eigensolver.o

# LIBS_Assignment_serial = -larpack
# LIBS_Assignment_avx    = -larpack
# LIBS_Assignment_openmp = -larpack
# LIBS_Assignment_mpi    = -larpack
# LIBS_Assignment_cuda   = -larpack

# CXXFLAGS_Assignment_openmp = -fopenmp
# CXXFLAGS_Assignment_avx = -mavx

# randutil.o : randutil.h randutil.ipp
# eigensolver.o : eigensolver.h

# # wildcard rules
# %_mpi.o : %_mpi.cpp
# 	$(MPICXX) $(CXXFLAGS) $(CFLAGS_$(basename $<)) -c $< -o $@

# %_mpi : %_mpi.cpp
# 	$(MPICXX) $(CXXFLAGS) $(CXXFLAGS_$@) $(filter %.o %.cpp, $^) $(LDFLAGS) $(LIBS_$@) $(LIB) -o $@

# %.o : %.cu
# 	$(NVCC) $(NVFLAGS) $(NVFLAGS_$(basename $<)) -c $< -o $@

# % : %.cu
# 	$(NVCC) $(NVFLAGS) $(NVFLAGS_$@) $(filter %.o %.cu, $^) $(LDFLAGS) $(LIBS_$@) $(LIB) -o $@

# %.o : %.cpp
# 	$(CXX) $(CXXFLAGS) $(CFLAGS_$(basename $<)) -c $< -o $@

# % : %.cpp
# 	$(CXX) $(CXXFLAGS) $(CXXFLAGS_$@) $(filter %.o %.cpp, $^) $(LDFLAGS) $(LIBS_$@) $(LIB) -o $@

# all : $(TARGETS)

# clean:
# 	rm -f $(TARGETS) *.o

# .PHONY: clean default all
