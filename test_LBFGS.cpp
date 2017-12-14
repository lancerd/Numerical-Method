#include "VEC.h"
#include "line_search.h"
#include <cmath>
#include <iostream>

static const int N = 10;
static const int M = 3;
static const double TOL = 1e-5;
static const double ALPHA = 100;

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

int main(void) {
    using std::cout;
    using std::endl;
    VEC x0(N);
    for (int i = 0; i < N; ++i) {
        x0[i] = -1;
    }
    VEC ans(N);
    Problem ptag = {
        .f = f, .g = g, .showIterationNumber = true, .tolerance = TOL};

    ans = LBFGS(x0, M, &ptag);

    cout << endl;
    cout << "m = " << M << endl;
    cout << "alpha = " << ALPHA << endl;
    cout << "x = " << ans << endl;
    cout << "error = " << f(ans) << endl;
    cout << "gradient norm = " << g(ans).norm() << endl;
    return 0;
}
