CXX = g++
CXXFLAGS  = -Wall -g -gdwarf-2 -O0
LDFLAGS = -lm -fopenmp

VERSION := $(shell awk '/VERSION/ {print $$3}' src/globals.h| sed 's/\"\(.*\)\"/\1/')
MODULES := Alignment.o Branch.o Matrix.o Node.o Optimizer.o Tree.o helper.o gmm.o
SRC := $(addprefix src/,$(MODULES))
OBJ := $(patsubst src/%.cpp,src/%.o,$(SRC))

BIN = gmm

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(LDFLAGS) $(LIBS) -o $(BIN) $(OBJ)

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

tarball:
	tar --transform "s,^,gmm-$(VERSION)/," -cjf gmm-$(VERSION).tar.bz2 Makefile src/*.cpp src/*.h

clean:
	$(RM) $(BIN) $(OBJ)
