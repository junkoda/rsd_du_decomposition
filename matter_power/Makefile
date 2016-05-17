#
# Makefile
#

EXEC    = matter_power

all: $(EXEC)

BOOST_DIR ?= #e.g. /opt/local
HDF5_DIR  ?=
FFTW3_DIR ?= 

DIR_PATH   = $(BOOST_DIR) $(HDF5_DIR) $(FFTW3_DIR)

CXXFLAGS  := $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS      := $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

OBJS = matter_power.o power_spectrum3b.o fft_mesh.o histogram.o hdf5_read.o
LIBS      += -lm $(LIBS) -lboost_program_options -lfftw3f -lhdf5

matter_power: $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $@

fft_mesh.o: fft_mesh.cpp fft_mesh.h
hdf5_read.o: hdf5_read.cpp particle.h hdf5_read.h
histogram.o: histogram.cpp histogram.h
matter_power.o: matter_power.cpp fft_mesh.h power_spectrum3b.h particle.h \
  hdf5_read.h
power_spectrum3b.o: power_spectrum3b.cpp histogram.h power_spectrum3b.h \
  particle.h

.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS)
