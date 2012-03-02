CXX = g++
CXXFLAGS  = -Wall -g -gdwarf-2 -O0
LDFLAGS = -lm -fopenmp

MODULES := Alignment.o Branch.o Matrix.o Node.o Optimizer.o Tree.o helper.o gmm.o
SRC := $(addprefix src/,$(MODULES))
OBJ := $(patsubst src/%.cpp,src/%.o,$(SRC))

BIN = gmm

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(LDFLAGS) $(LIBS) -o $(BIN) $(OBJ)

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) $(BIN) $(OBJ)
