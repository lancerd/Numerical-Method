#include "VEC.h"
#include "line_search.h"
#include <cmath>
#include <iostream>

using namespace std;

static const int N = 4;
static const int M = 3;
static const double TOL = 1e-6;
static const double ALPHA = 1;

double f(const VEC &x) {
    double ans = 0;
    for (int i = 0; i < N; i += 2) {
        ans += ALPHA * pow(x[i + 1] - x[i] * x[i], 2);
        ans += pow(1 - x[i], 2);
    }
    return ans;
}

VEC g(const VEC &x) {
    VEC res(x.len());
    res.clear();
    for (int i = 0; i < N; i += 2) {
        res[i + 1] += 2 * ALPHA * (x[i + 1] - x[i] * x[i]);
        res[i] += 2 * ALPHA * (x[i + 1] - x[i] * x[i]) * (-2 * x[i]);
        res[i] += 2 * (1 - x[i]) * -1;
    }
    return res;
}

MAT H(const VEC &x) {
    MAT res(x.len());
    res.clear();
    for (int i = 0; i < N; i += 2) {
        res[i + 1][i + 1] = 2 * ALPHA;
        res[i + 1][i] = res[i][i + 1] = -4 * ALPHA * x[i];
        res[i][i] = ALPHA * (12 * x[i] * x[i] - 4 * x[i + 1]) + 2;
    }
    return res;
}

void test_fminunc(void) {
    VEC x0(N);
    x0.fill(-1);
    Problem ptag = {
        .f = f, .g = g, .H = H, .tolerance = TOL, .max_iteration = 10000};
    Result restag;
    std::vector<std::string> algorithm = {"Steepest_Descent", "Newton"};
    for (auto &algo_name : algorithm) {
        ptag.algorithm = algo_name;
        std::cout << "----------" << algo_name << "----------" << std::endl;
        VEC ans(fminunc(x0, &ptag, &restag));
        display_result(ans, &ptag, &restag);
    }
}

int main(void) {
    using std::cout;
    using std::endl;
    VEC x0(N);
    x0.fill(-1);
    Problem ptag = {
        .f = f, .g = g, .H = H, .tolerance = TOL, .max_iteration = 10000};
    Result restag;

    cout << "----------Objective function------------" << endl;
    cout << "Rosenbrock's function" << endl;
    cout << "alpha = " << ALPHA << endl;

    // cout << "-----------Steepest_Descent-------------" << endl;
    // VEC ans_Steepest_Descent(N);
    // ans_Steepest_Descent = Steepest_Descent(x0, &ptag, &restag);
    // display_result(ans_Steepest_Descent, &ptag, &restag);
    //
    // cout << "----------------Newton------------------" << endl;
    // VEC ans_Newton(N);
    // ans_Newton = Newton(x0, &ptag, &restag);
    // display_result(ans_Newton, &ptag, &restag);

    cout << "----------------LBFGS-------------------" << endl;
    VEC ans_LBFGS(N);
    ans_LBFGS = LBFGS(x0, M, &ptag, &restag);
    cout << "M = " << M << endl;
    display_result(ans_LBFGS, &ptag, &restag);

    cout << "-------------test fminunc---------------" << endl;
    test_fminunc();
    return 0;
}
