EXEC = smoothed_velocity

all: $(EXEC)

# Compiler
CXX       ?= c++
CPPFLAGS  := -O3 -Wall -I../decomposed_power
LIBS      := -lm -lboost_program_options -lhdf5 -lfftw3f -lgsl -lgslcblas
vpath %.cpp ../decomposed_power
vpath %.h   ../decomposed_power

BOOST_DIR ?= #e.g. /opt/local
HDF5_DIR  ?=
GSL_DIR   ?= #e.g. $(HOME)/Research/opt/gcc/gsl
FFTW3_DIR ?=

DIR_PATH   = $(BOOST_DIR) $(GSL_DIR) $(FFTW3_DIR) $(HDF5_DIR)

CXXFLAGS  += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS      += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

OBJS1 := smoothed_velocity.o kdtree.o
OBJS1 += halo_file.o
#OBJS1 += nearest_nbr_velocity.o
OBJS1 += hdf5_write.o hdf5_read.o

smoothed_velocity: $(OBJS1)
	$(CXX) $(OBJS1) $(LIBS) -o $@


.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS1)

