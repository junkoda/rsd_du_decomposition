#
# fof
#   

# Compiler
CC        ?= cc -std=c99
CXX       ?= c++
CFLAGS    := -O3 -Wall
CPPFLAGS  := -O3 -Wall
LIBS      := -lm

# Define paths for BOOST and HDF5 library if necessary.
BOOST_DIR ?= #e.g. /opt/local
HDF5_DIR  ?= #e.g. /Users/jkoda/Research/opt/gcc/hdf5

DIR_PATH = $(BOOST_DIR) $(HDF5_DIR)

CFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS   += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

LIBS   += -lboost_program_options -lhdf5

EXEC = group_velocity
all: $(EXEC)

OBJS := group_velocity.o fof.o msg.o read_gadget.o hdf5_write.o

LIBS += -lhdf5

group_velocity: $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $@

.PHONY: clean run dependence
clean :
	rm -f $(EXEC) $(OBJS) $(OBJS2)

dependence:
	gcc -MM -MG *.c
