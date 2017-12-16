#include "line_search.h"
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

static const int MAXITERATION = 1e7;

bool check_ptag(const Problem *ptag) {
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
    if (ptag->tolerance <= 0) {
        std::cout << "tolerance not set" << std::endl;
        return false;
    }
    return true;
}

template <class T>
static void display_result_help(const char *name, const T &value) {
    std::cout << name << std::endl << std::setw(5) << "" << value << std::endl;
}

void display_result(const VEC &x, const Problem *ptag) {
    std::cout << "Result:" << std::endl;
    display_result_help("x", x);
    display_result_help("error", ptag->f(x));
    display_result_help("gradient", ptag->g(x).norm());
}

double step_length(const VEC &x, const VEC &p, const Problem *ptag,
                   const double initialguess, const double c1, const double c2,
                   const double alpha_max) {
    // return step length satisfying the Wolfe conditions
    double (*f)(const VEC &x) = ptag->f;
    VEC (*g)(const VEC &x) = ptag->g;

    auto phi = [&f, &x, &p](double alpha) { return f(x + alpha * p); };
    auto phip = [&g, &x, &p](double alpha) { return g(x + alpha * p) * p; };

    double phi_0 = phi(0), phi_p0 = phip(0);

    if (initialguess > 0 &&
        phi(initialguess) - phi_0 <= c1 * initialguess * phi_p0 &&
        phip(initialguess) >= c2 * phi_p0)
        return initialguess;

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

VEC Steepest_Descent(const VEC &x0, const Problem *ptag) {
    if (!check_ptag(ptag)) {
        return x0;
    }
    VEC (*g)(const VEC &x) = ptag->g;

    VEC x(x0);
    int iter = 0;
    while (g(x).norm() > ptag->tolerance) {
        VEC p(-g(x));
        double alpha = step_length(x, p, ptag);
        x += p * alpha;
        ++iter;
        if ((ptag->MaxIter >= 0 && iter == ptag->MaxIter) ||
            iter == MAXITERATION) {
            break;
        }
    }
    if (ptag->showIterationNumber) {
        std::cout << "Number of iterations = " << iter << std::endl;
    }
    return x;
}

VEC LBFGS(const VEC &x0, const int m, const Problem *ptag) {
    typedef std::vector<std::pair<VEC, VEC>> vecVV;
    auto setH = [&x0](MAT &H0, const vecVV &tmp_sy) -> void {
        double gamma;
        if (!tmp_sy.empty()) {
            const VEC &tmps = tmp_sy.back().first;
            const VEC &tmpy = tmp_sy.back().second;
            gamma = (tmps * tmpy) / (tmpy * tmpy);
        } else {
            gamma = 1;
        }
        H0 = MAT::identity(H0.dim(), gamma);
    };

    auto LBFGS_solve = [](const MAT &H0, const VEC &gx,
                          const vecVV &tmp_sy) -> VEC {
        VEC q(gx);
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

    if (!check_ptag(ptag)) {
        return x0;
    }
    VEC (*g)(const VEC &x) = ptag->g;

    VEC x(x0);
    int iter = 0;
    vecVV tmp_sy;
    while (g(x).norm() > ptag->tolerance) {
        MAT H0(x0.len());
        setH(H0, tmp_sy);
        VEC p(LBFGS_solve(H0, -g(x), tmp_sy));
        double alpha = step_length(x, p, ptag, 1);
        VEC xp = x + p * alpha;
        // keep the m most recent (s_k, y_k)
        if (iter > m) {
            tmp_sy.erase(tmp_sy.begin());
        }
        tmp_sy.push_back(std::make_pair(xp - x, g(xp) - g(x)));
        x = xp;
        ++iter;
        if ((ptag->MaxIter >= 0 && iter == ptag->MaxIter) ||
            iter == MAXITERATION) {
            break;
        }
    }
    if (ptag->showIterationNumber) {
        std::cout << "Number of iterations = " << iter << std::endl;
    }
    return x;
}
