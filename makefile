MPICXX?=mpic++
CXX?=g++
CXXFLAGS?=-std=c++11 -O2
NVCC?=nvcc
NVFLAGS?=-O2 --gpu-architecture=sm_35 -Wno-deprecated-gpu-targets

TARGETS = pp_serial

default : all

# dependancies
$(TARGETS): parse_constants.o parse_constants.h

parse_constants.o: parse_constants.h

# wildcard rules
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
