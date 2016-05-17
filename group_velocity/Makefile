EXEC = group_velocity
all: $(EXEC)

# Compiler
CC        ?= cc
CXX       ?= c++
CFLAGS    := -std=c99 -O3 -Wall
CPPFLAGS  := -O3 -Wall
LIBS      := -lm -lboost_program_options -lhdf5

# Define paths for BOOST and HDF5 library if necessary.
BOOST_DIR ?= #e.g. /opt/local
HDF5_DIR  ?= #e.g. /Users/jkoda/Research/opt/gcc/hdf5

DIR_PATH   = $(BOOST_DIR) $(HDF5_DIR)

CFLAGS    += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
CPPFLAGS  += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS      += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)



OBJS := group_velocity.o fof.o msg.o read_gadget.o hdf5_write.o


group_velocity: $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $@

.PHONY: clean run dependence
clean :
	rm -f $(EXEC) $(OBJS)

dependence:
	gcc -MM -MG *.c
