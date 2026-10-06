PROJECT := mesh_lod232
CXX := g++
CPPFLAGS += -Ithird_party/eigen3 -Iinclude
CXXFLAGS += -std=c++20 -O2 -Wall -Wextra

BUILD := build
BIN := bin
LIB := $(BUILD)/lib$(PROJECT).a
LIB_SRCS := src/validation.cpp src/mesh.cpp src/quadric.cpp src/simplify.cpp
LIB_OBJS := $(LIB_SRCS:src/%.cpp=$(BUILD)/%.o)

.PHONY: all test run-example clean

all: $(LIB) $(BIN)/example_open_surface $(BIN)/selftest

$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(LIB): $(LIB_OBJS)
	$(AR) rcs $@ $^

$(BIN)/example_open_surface: examples/open_surface.cpp $(LIB) | $(BIN)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -L$(BUILD) -l$(PROJECT) -o $@

$(BIN)/selftest: tests/selftest.cpp $(LIB) | $(BIN)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -L$(BUILD) -l$(PROJECT) -o $@

$(BUILD) $(BIN):
	mkdir -p $@

test: $(BIN)/selftest
	./$(BIN)/selftest

run-example: $(BIN)/example_open_surface
	./$(BIN)/example_open_surface

clean:
	rm -rf $(BUILD) $(BIN)
