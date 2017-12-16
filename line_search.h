#ifndef LINESEARCH
#define LINESEARCH
#include "MAT.h"
#include "VEC.h"

static const double C1 = 1e-4;      // default constant for Armijo condition
static const double C2 = 0.99;      // default constant for curvature condition
static const double ALPHA_MAX = 10; // maximum step length allowed

struct Problem {
    bool display;              // display step by step
    int max_iteration;         // maximum number of iterations
    double tolerance;          // tolerance of 1st order optimality
    double (*f)(const VEC &x); // objective function
    VEC (*g)(const VEC &x);    // gradient function
    MAT (*H)(const VEC &x);    // Hessian function
};

struct Result {
    int iterations;                // number of iterations
    int fCount;                    // number of objective function evaluation
    double first_order_optimality; // measure of 1st order optimality
};

inline static bool check_ptag(const Problem *ptag);
template <class T>
static void display_result_help(const char *name, const T &value);
void display_result(const VEC &x, const Problem *ptag = NULL,
                    const Result *restag = NULL);

double step_length(const VEC &x0, const VEC &p, const Problem *ptag = NULL,
                   const double initial_guess = -1, const double c1 = C1,
                   const double c2 = C2, const double alpha_max = ALPHA_MAX);
VEC Steepest_Descent(const VEC &x0, const Problem *ptag = NULL,
                     Result *restag = NULL);
VEC Newton(const VEC &x0, const Problem *ptag = NULL, Result *restag = NULL);
VEC LBFGS(const VEC &x0, const int m, const Problem *ptag = NULL,
          Result *restag = NULL);
#endif
