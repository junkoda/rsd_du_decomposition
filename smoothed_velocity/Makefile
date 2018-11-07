EXEC = smoothed_velocity

all: $(EXEC)

# Compiler
CXX       ?= c++
CPPFLAGS  := -O3 -Wall #-I../decomposed_power
LIBS      := -lm -lboost_program_options -lhdf5 -lfftw3f -lgsl -lgslcblas

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

smoothed_velocity.o: smoothed_velocity.cpp gadget_file2.h halo_file.h hdf5_read.h

gadget_file2.h: ../decomposed_power/gadget_file2.h
	ln -s $< .

halo_file.h: ../decomposed_power/halo_file.h
	ln -s $< .

hdf5_read.h: ../decomposed_power/hdf5_read.h
	ln -s $< .

halo_file.cpp: ../decomposed_power/halo_file.h
	ln -s $< .

hdf5_read.cpp: ../decomposed_power/hdf5_read.cpp
	ln -s $< .


.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS1)

