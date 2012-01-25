CXX = g++
CXXFLAGS  = -Wall -g -O0
LFLAGS = -lm 

OBJ = Alignment.o Branch.o Matrix.o Node.o Tree.o helper.o gmm.o
BIN = gmm

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(BIN) $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) $(BIN) $(OBJ)
