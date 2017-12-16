
CXXFLAGS = -std=c++11 -O2 -Wall -Wextra -Wfatal-errors
CXX = g++

objects = VEC.o MAT.o line_search.o
test = test

.PHONY: all clean

all: $(test)
	./test

$(test): % : $(objects) %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $@.cpp $(objects)

$(objects): %.o : %.h %.cc
	$(CXX) $(CXXFLAGS) -DDEBUG -c $(subst .o,,$@).cc

clean:
	rm *.o
