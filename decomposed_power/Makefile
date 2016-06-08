EXEC = decomposed_power

all: $(EXEC)

# Compiler
CXX       ?= c++
CPPFLAGS  := -O3 -Wall
LIBS      := -lm -lboost_program_options -lhdf5 -lfftw3f -lgsl -lgslcblas

BOOST_DIR ?= #e.g. /opt/local
HDF5_DIR  ?=
GSL_DIR   ?= #e.g. $(HOME)/Research/opt/gcc/gsl
FFTW3_DIR ?=

DIR_PATH   = $(BOOST_DIR) $(GSL_DIR) $(FFTW3_DIR) $(HDF5_DIR)

CXXFLAGS  += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS      += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

OBJS1 := decomposed_power.o
OBJS1 += binary_tree.o nbr_finder.o fft_mesh.o
OBJS1 += halo_file.o
OBJS1 += nearest_nbr_velocity.o power_spectrum3b_v.o histogram.o
OBJS1 += transformation.o density_mesh.o
OBJS1 += hdf5_write.o hdf5_read.o

decomposed_power: $(OBJS1)
	$(CXX) $(OBJS1) $(LIBS) -o $@


.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS1)

