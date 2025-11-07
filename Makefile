CXX := g++
CXXFLAGS := -O1 -fopenmp -fdiagnostics-color=always -g
INCLUDES := -Iinclude -Iinclude/HDF5
LDFLAGS := -Llibs/HDF5 -lhdf5 -lhdf5_cpp

SRC_DIR := source
OBJ_DIR := build
TARGET := alpha.exe

SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))

all: $(TARGET)

$(TARGET): $(OBJS)
	@echo Linking $@
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@echo Compiling $<
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@


$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: all clean
