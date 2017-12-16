#ifndef LINESEARCH
#define LINESEARCH
#include "MAT.h"
#include "VEC.h"

static const double C1 = 1e-4;
static const double C2 = 0.99;
static const double ALPHA_MAX = 10;

struct Problem {
    bool display;
    int MaxIter;
    double tolerance;
    double (*f)(const VEC &);
    VEC (*g)(const VEC &);
    MAT (*H)(const VEC &);
};

struct Result {
    int iterations;
    int fCount;
    double firstorderopt;
};

bool check_ptag(const Problem *);
template <class T> static void display_result_help(const char *, const T &);
void display_result(const VEC &, const Problem * = NULL, const Result * = NULL);

double step_length(const VEC &, const VEC &, const Problem * = NULL,
                   const double = -1, const double = C1, const double = C2,
                   const double = ALPHA_MAX);
VEC Steepest_Descent(const VEC &, const Problem * = NULL, Result * = NULL);
VEC Newton(const VEC &, const Problem * = NULL, Result * = NULL);
VEC LBFGS(const VEC &, const int, const Problem * = NULL, Result * = NULL);
#endif
