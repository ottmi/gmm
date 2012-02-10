CXX = g++
CXXFLAGS  = -Wall -g -gdwarf-2 -O3 -fopenmp
LFLAGS = -lm -fopenmp

OBJ = Alignment.o Branch.o Matrix.o Node.o Tree.o helper.o gmm.o
BIN = gmm

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(BIN) $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) $(BIN) $(OBJ)
