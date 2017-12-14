#ifndef LINESEARCH
#define LINESEARCH
#include "MAT.h"
#include "VEC.h"

static const double C1 = 1e-4;
static const double C2 = 0.99;
static const double ALPHA_MAX = 10;

struct Problem {
    bool display;
    bool showIterationNumber;
    int MaxIter;
    double tolerance;
    double (*f)(const VEC &);
    VEC (*g)(const VEC &);
};

bool check_ptag(const Problem *);

double step_length(const VEC &, const VEC &, const Problem * = NULL,
                   const double = -1, const double = C1, const double = C2,
                   const double = ALPHA_MAX);
VEC Steepest_Descent(const VEC &, const Problem * = NULL);
VEC LBFGS(const VEC &, const int, const Problem * = NULL);
#endif
