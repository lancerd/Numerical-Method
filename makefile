
all: VEC.o MAT.o line_search.o test_LBFGS.cpp
	g++ -std=c++11 -O2 test_LBFGS.cpp VEC.o MAT.o line_search.o -Wfatal-errors

VEC.o: VEC.h VEC.cpp
	g++ -std=c++11 -O2 VEC.cpp -c -DDEBUG

MAT.o: MAT.h MAT.cpp
	g++ -std=c++11 -O2 MAT.cpp -c -DDEBUG

line_search.o: line_search.h line_search.cpp
	g++ -std=c++11 -O2 line_search.cpp -c
