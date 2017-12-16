#ifndef LINESEARCH
#define LINESEARCH
#include "MAT.h"
#include "VEC.h"

static const double C1 = 1e-4;      // default constant for Armijo condition
static const double C2 = 0.99;      // default constant for curvature condition
static const double ALPHA_MAX = 10; // maximum step length allowed

struct Problem {
    bool display;             // display step by step
    int MaxIter;              // maximum number of iterations
    double tolerance;         // tolerance of 1st order optimality
    double (*f)(const VEC &); // objective function
    VEC (*g)(const VEC &);    // gradient function
    MAT (*H)(const VEC &);    // Hessian function
};

struct Result {
    int iterations;       // number of iterations
    int fCount;           // number of objective function evaluation
    double firstorderopt; // measure of 1st order optimality
};

inline static bool check_ptag(const Problem *);
template <class T> static void display_result_help(const char *, const T &);
void display_result(const VEC &, const Problem * = NULL, const Result * = NULL);

double step_length(const VEC &, const VEC &, const Problem * = NULL,
                   const double = -1, const double = C1, const double = C2,
                   const double = ALPHA_MAX);
VEC Steepest_Descent(const VEC &, const Problem * = NULL, Result * = NULL);
VEC Newton(const VEC &, const Problem * = NULL, Result * = NULL);
VEC LBFGS(const VEC &, const int, const Problem * = NULL, Result * = NULL);
#endif
