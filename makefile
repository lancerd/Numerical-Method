
CXXFLAGS = -std=c++11 -O2 -Wall -Wextra -Wfatal-errors
CXX = g++

objects = VEC.o MAT.o line_search.o
test = test_LBFGS test_Steepest_Descent

.PHONY: all clean

all: $(test)
	./test_Steepest_Descent
	./test_LBFGS

$(test): % : $(objects) %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $@.cpp $(objects)

$(objects): %.o : %.h %.cpp
	$(CXX) $(CXXFLAGS) -DDEBUG -c $(subst .o,,$@).cpp

clean:
	rm *.o
