# Variables
ARCH ?= sm_35       # Default GPU architecture (modifiable via command line)
HOST_COMP ?= mpicxx   # Host compiler
NVCC = nvcc

# Source and target
SRC = task_gpu.cpp        # Source file
OBJ = $(SRC:.cpp=.o)  # Object file (derived from SRC)
TARGET = program     # Output executable

NVCC_FLAGS = -arch=$(ARCH) -ccbin=$(HOST_COMP) -Xcompiler="-fopenmp"
# Build rules
all: $(TARGET)

$(TARGET): $(OBJ)
        $(NVCC) $(NVCC_FLAGS)  -o $@ $^


clean:
        rm -f $(OBJ) $(TARGET)

# Example of how to set variables at runtime:
# make ARCH=sm_60 HOST_COMP=mpicxx