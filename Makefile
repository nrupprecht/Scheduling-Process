# Makefile snatched from <https://ubuntuforums.org/showthread.php?t=1204739>

COMPILER ?= 1
USE_MPI  ?= 1
DEBUG    ?= 0
NO_SIMD  ?= 0
DO_OPTIMIZE ?= 1

# Add flags to this variable as neccessary
FLAGS = 

MPICC = mpic++

# Compiler options
ifeq ($(USE_MPI), 1)
CC = $(MPICC)
FLAGS += -D USE_MPI=1
COMPILER = 2
else ifeq ($(COMPILER), 0)
CC = icpc
FLAGS += -D _INTEL_=1
else ifeq ($(COMPILER), 1)
CC = g++
FLAGS += -D _CLANG_=1
else
CC = clang++
FLAGS += -D _CLANG_=1
endif 

ifeq ($(DEBUG), 1)
FLAGS += -D DEBUG=1
endif

ifeq ($(NO_SIMD), 1)
FLAGS += -D SIMD_TYPE=0
endif

# Name of the application to compile
APP   = driver
# Directories for source, object, bin files
SRC   = src
OBJ   = obj
BIN   = bin
APPLICATIONS = $(SRC)/applications

EXCLUDE := $(shell find $(APPLICATIONS) -name '*.cpp')

SRCS1   := $(shell find $(SRC) -name '*.cpp')
SRCS    := $(filter-out $(EXCLUDE), $(SRCS1))
SRCDIRS := $(shell find . -name '*.cpp' -exec dirname {} \; | uniq)
OBJS    := $(patsubst %.cpp,$(OBJ)/%.o,$(SRCS))

CSTD     = -std=c++17
DEBUG    = -g
ifeq ($(COMPILER), 0) 
ifeq ($(DO_OPTIMIZE), 1)
OPTIMIZE := -O3 -funroll-loops -restrict -no-inline-max-size #-ffast-math -no-prec-div
endif
else
ifeq ($(DO_OPTIMIZE), 1)
OPTIMIZE := -O3
endif
endif
#OPENCV   := $(shell pkg-config --cflags --libs opencv)
CFLAGS   := -c $(CSTD) $(DEBUG) $(OPTIMIZE) $(FLAGS)

all: $(BIN)/$(APP)

# Executables

$(BIN)/$(APP): obj/driver.o $(OBJS)
	@mkdir -p `dirname $@`
	@echo "Linking $@..."
	@$(CC) -o $@ obj/driver.o $(OBJS)
# General object files

$(OBJ)/%.o: %.cpp
	@mkdir -p `dirname $@`
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@

# Excluded object files

$(OBJ)/driver.o: src/applications/driver.cpp
	@mkdir -p `dirname $@`
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@


.PHONY : clean
clean:
	rm -r -f $(OBJ) $(BIN)/*

.PHONY : forces
forces:
	rm obj/src/interactions/*
	make

.PHONY : creators
creators:
	rm obj/src/creators/*
	make
