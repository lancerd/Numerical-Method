#include "line_search.h"
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

static const int MAXITERATION = 1e7;

static bool check_ptag(const Problem *ptag) {
    // check whether ptag is valid
    if (!ptag) {
        std::cout << "problem tag not set" << std::endl;
        return false;
    }
    if (ptag->f == NULL) {
        std::cout << "function not set" << std::endl;
        return false;
    }
    if (ptag->g == NULL) {
        std::cout << "gradient not set" << std::endl;
        return false;
    }
    if (ptag->algorithm == "Newton" && ptag->H == NULL) {
        std::cout << "Hessian not set" << std::endl;
        return false;
    }
    if (ptag->tolerance <= 0) {
        std::cout << "tolerance not set" << std::endl;
        return false;
    }
    if (ptag->max_iteration <= 0) {
        std::cout << "max iteration not set" << std::endl;
        return false;
    }
    return true;
}

template <class T>
static void display_result_help(const char *name, const T &value) {
    std::cout << name << std::endl << std::setw(5) << "" << value << std::endl;
}

void display_result(const VEC &x, const Problem *ptag, const Result *restag) {
    std::cout << "Result:" << std::endl;
    display_result_help("x", x);
    display_result_help("iterations", restag->iterations);
    display_result_help("function value", ptag->f(x));
    display_result_help("1st order optimality", restag->first_order_optimality);
}

double step_length(const VEC &x, const VEC &p, const Problem *ptag,
                   const double initial_guess, const double c1, const double c2,
                   const double alpha_max) {
    // select a step length satisfying the Wolfe conditions
    // use step length selection algorithm 3.5, 3.6 in
    // Nocedal, J. and S. J. Wright. Numerical Optimization, Second Edition

    // x : starting point
    // p : a decrease direction
    // initial_guess: initial guess of the step length
    //               e.g. for Newton Method, initialguess = 1
    // 0 < c1 < c2 < 1 : constant for Wolfe conditions
    double (*f)(const VEC &) = ptag->f;
    VEC (*g)(const VEC &) = ptag->g;

    // phi: one dimensional minimization problem
    auto phi = [&f, &x, &p](double alpha) { return f(x + alpha * p); };
    auto phip = [&g, &x, &p](double alpha) { return g(x + alpha * p) * p; };

    double phi_0 = phi(0), phi_p0 = phip(0);

    if (initial_guess > 0 &&
        phi(initial_guess) - phi_0 <= c1 * initial_guess * phi_p0 &&
        phip(initial_guess) >= c2 * phi_p0)
        return initial_guess;

    auto zoom = [&x, &p, &phi, &phip, phi_0, phi_p0, c1,
                 c2](double alpha_lo, double alpha_hi) -> double {
        while (1) {
            double alpha_j = (alpha_lo + alpha_hi) / 2;
            double val = phi(alpha_j);
            if (val - phi_0 > c1 * alpha_j * phi_p0 || val >= phi(alpha_lo)) {
                alpha_hi = alpha_j;
            } else {
                double valp = phip(alpha_j);
                if (fabs(valp) <= -c2 * phi_p0) {
                    return alpha_j;
                }
                if (valp * (alpha_hi - alpha_lo) >= 0) {
                    alpha_hi = alpha_lo;
                }
                alpha_lo = alpha_j;
            }
        }
    };

    double alpha_i = alpha_max / 2, alpha_i1 = 0;
    for (int i = 1;; ++i) {
        double val = phi(alpha_i), valp = phip(alpha_i);
        if (val - phi_0 > c1 * alpha_i * phi_p0 ||
            (val >= phi(alpha_i1) && i > 1)) {
            return zoom(alpha_i1, alpha_i);
        }
        if (fabs(valp) <= -c2 * phi_p0) {
            return alpha_i;
        }
        if (valp >= 0) {
            return zoom(alpha_i, alpha_i1);
        }
        alpha_i1 = alpha_i;
        alpha_i = (alpha_i + alpha_max) / 2;
    }
}

static void LBFGS_setH(MAT &H0,
                       const std::vector<std::pair<VEC, VEC>> &tmp_sy) {
    double gamma;
    if (!tmp_sy.empty()) {
        const VEC &tmps = tmp_sy.back().first;
        const VEC &tmpy = tmp_sy.back().second;
        gamma = (tmps * tmpy) / (tmpy * tmpy);
    } else {
        gamma = 1;
    }
    H0 = MAT::identity(H0.dim(), gamma);
}

static VEC LBFGS_solve(const MAT &H0, const VEC &g,
                       const std::vector<std::pair<VEC, VEC>> &tmp_sy) {
    VEC q(g);
    std::vector<double> alpha(tmp_sy.size());
    // rho_k = 1 / (y_k * s_k)
    for (int i = (int)tmp_sy.size() - 1; i >= 0; --i) {
        const VEC &s_i = tmp_sy[i].first, &y_i = tmp_sy[i].second;
        const double rho_i = 1 / (y_i * s_i);
        alpha[i] = rho_i * s_i * q;
        q = q - alpha[i] * y_i;
    }
    VEC r(H0 * q);
    for (int i = 0; i < (int)tmp_sy.size(); ++i) {
        const VEC &s_i = tmp_sy[i].first, &y_i = tmp_sy[i].second;
        const double rho_i = 1 / (y_i * s_i);
        double beta = rho_i * y_i * r;
        r = r + s_i * (alpha[i] - beta);
    }
    return r;
};

VEC LBFGS(const VEC &x0, const int m, const Problem *ptag, Result *restag) {
    if (!check_ptag(ptag)) {
        return x0;
    }
    VEC (*g)(const VEC &) = ptag->g;

    VEC x(x0);
    int iter = 0;
    std::vector<std::pair<VEC, VEC>> tmp_sy;
    while (g(x).norm() > ptag->tolerance) {
        MAT H0(x0.len());
        LBFGS_setH(H0, tmp_sy);
        VEC p(LBFGS_solve(H0, -g(x), tmp_sy));
        double alpha = step_length(x, p, ptag, 1);
        VEC xp = x + p * alpha;
        // keep only the m most recent (s_k, y_k)
        if (iter > m) {
            tmp_sy.erase(tmp_sy.begin());
        }
        tmp_sy.push_back(std::make_pair(xp - x, g(xp) - g(x)));
        x = xp;
        ++iter;
        if ((ptag->max_iteration >= 0 && iter == ptag->max_iteration) ||
            iter == MAXITERATION) {
            break;
        }
    }
    if (restag) {
        restag->iterations = iter;
        restag->first_order_optimality = g(x).norm();
    }
    return x;
}

static VEC line_search_Steepest_Descent_update(const VEC &x,
                                               const Problem *ptag) {
    VEC p(-ptag->g(x));
    double alpha = step_length(x, p, ptag);
    return x + p * alpha;
}

static VEC line_search_Newton_update(const VEC &x, const Problem *ptag) {
    VEC p(-luSolve(ptag->H(x), ptag->g(x)));
    double alpha = step_length(x, p, ptag, 1);
    return x + p * alpha;
}

static VEC line_search(const VEC &x0, const Problem *ptag, Result *restag,
                       VEC (*update_strategy)(const VEC &, const Problem *)) {
    if (update_strategy == NULL) {
        return x0;
    }
    VEC (*g)(const VEC &) = ptag->g;

    VEC x(x0);
    int iter = 0;
    for (iter = 0; g(x).norm() > ptag->tolerance && iter <= ptag->max_iteration;
         ++iter) {
        x = update_strategy(x, ptag);
    }
    if (restag) {
        restag->iterations = iter;
        restag->first_order_optimality = g(x).norm();
    }
    return x;
}

VEC fminunc(const VEC &x0, const Problem *ptag, Result *restag) {
    if (!check_ptag(ptag)) {
        return x0;
    }

    VEC (*update_strategy)(const VEC &, const Problem *) = NULL;
    if (ptag->algorithm == "Steepest_Descent") {
        update_strategy = line_search_Steepest_Descent_update;
    } else if (ptag->algorithm == "Newton") {
        update_strategy = line_search_Newton_update;
    } else {
        std::cout << "Please specify the algorithm in fminunc" << std::endl;
    }
    return line_search(x0, ptag, restag, update_strategy);
}
