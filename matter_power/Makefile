#
# Makefile
#

EXEC    = matter_power

all: $(EXEC)

BOOST_DIR ?= #e.g. /opt/local
GSL_DIR   ?= #e.g. $(HOME)/Research/opt/gcc/gsl
FFTW3_DIR ?= 

DIR_PATH   = $(BOOST_DIR) $(GSL_DIR) $(FFTW3_DIR)

CXXFLAGS  := $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS      := $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

OBJS = matter_power.o power_spectrum3b.o fft_mesh.o histogram.o hdf5_read.cpp
LIBS      += -lm $(LIBS) -lboost_program_options -lfftw3f -lhdf5

matter_power: $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $@


.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS)
