// Name: Wang Wei Luan
// Sutdent ID: 102022191
#include "MAT.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

MAT::MAT(int dim) {
    n = dim;
    va = (VEC **)malloc(n * sizeof(VEC *));
    for (int i = 0; i < n; ++i)
        va[i] = newVEC(n);
}
MAT::MAT(const MAT &m) {
    n = m.n;
    va = (VEC **)malloc(n * sizeof(VEC *));
    for (int i = 0; i < n; ++i) {
        va[i] = newVEC(n);
        (*va[i]) = (*m.va[i]);
    }
}
MAT::MAT(int dim, double *v) {
    n = dim;
    va = (VEC **)malloc(n * sizeof(VEC *));
    for (int i = 0; i < n; ++i) {
        va[i] = newVEC(n);
        for (int j = 0; j < n; ++j)
            (*va[i])[j] = *(v++);
    }
}
MAT::~MAT() {
    for (int i = n - 1; i >= 0; --i)
        free(va[i]);
    free(va);
}
int MAT::dim() const { return n; }
MAT MAT::tpose() const {
    MAT mt(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            mt[i][j] = (*va[j])[i];
    return mt;
}
MAT MAT::operator-() const {
    MAT mt(n);
    for (int i = 0; i < n; ++i)
        (*mt.va[i]) = -(*va[i]);
    return mt;
}
MAT &MAT::operator=(const MAT m) {
    for (int i = 0; i < n; ++i)
        (*va[i]) = (*m.va[i]);
    return *this;
}
MAT &MAT::operator+=(const MAT m) {
    for (int i = 0; i < n; ++i)
        (*va[i]) += (*m.va[i]);
    return *this;
}
MAT &MAT::operator-=(const MAT m) {
    for (int i = 0; i < n; ++i)
        (*va[i]) -= (*m.va[i]);
    return *this;
}
MAT &MAT::operator*=(double a) {
    for (int i = 0; i < n; ++i)
        (*va[i]) *= a;
    return *this;
}
MAT &MAT::operator/=(double a) {
    for (int i = 0; i < n; ++i)
        (*va[i]) /= a;
    return *this;
}
MAT MAT::operator+(const MAT m) const {
    MAT mt(*this);
    mt += m;
    return mt;
}
MAT MAT::operator-(const MAT m) const {
    MAT mt(*this);
    mt -= m;
    return mt;
}
MAT MAT::operator*(const MAT m) const {
    MAT mt(n);
    MAT mtpose = m.tpose();
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            mt[i][j] = (*va[i]) * mtpose[j];
    return mt;
}
VEC &MAT::operator[](int ind) {
#ifdef DEBUG
    if (ind < 0) {
        ind = 0;
        std::cout << "index < 0" << std::endl;
    } else if (ind >= n) {
        ind = n - 1;
        std::cout << "index >= dim" << std::endl;
    }
#endif
    return *va[ind];
}
const VEC &MAT::operator[](int ind) const {
#ifdef DEBUG
    if (ind < 0) {
        ind = 0;
        std::cout << "index < 0" << std::endl;
    } else if (ind >= n) {
        ind = n - 1;
        std::cout << "index >= dim" << std::endl;
    }
#endif
    return *va[ind];
}
VEC MAT::operator*(const VEC v) const {
    VEC vt(n);
    for (int i = 0; i < n; ++i)
        vt[i] = (*va[i]) * v;
    return vt;
}
MAT MAT::operator*(double a) const {
    MAT mt(*this);
    mt *= a;
    return mt;
}
MAT MAT::operator/(double a) const {
    MAT mt(*this);
    mt /= a;
    return mt;
}
void MAT::clear() {
    for (int i = 0; i < n; ++i)
        va[i]->clear();
}
MAT MAT::zero(int dim) {
    MAT mt(dim);
    mt.clear();
    return mt;
}
MAT MAT::identity(int dim, double r) {
    MAT mt(dim);
    mt.clear();
    for (int i = 0; i < dim; ++i)
        mt[i][i] = r;
    return mt;
}
MAT operator*(double a, MAT m) {
    m *= a;
    return m;
}
VEC operator*(const VEC v, const MAT m) {
    MAT mt = m.tpose();
    return mt * v;
}

MAT &luFact(MAT &m) {
    int n = m.dim();
    for (int i = 0; i < n; ++i) {
        // copy no need due to in-place LU
        for (int j = i + 1; j < n; ++j) { // form l[j][i]
            m[j][i] /= m[i][i];
        }
        for (int j = i + 1; j < n; ++j) { // update submatrix
            for (int k = i + 1; k < n; ++k) {
                m[j][k] -= m[j][i] * m[i][k];
            }
        }
    }
    return m;
}

VEC fwdSubs(const MAT &m, VEC b) {
    VEC y(b);
    int n = m.dim();
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            y[j] -= m[j][i] * y[i];
    return y;
}

VEC bckSubs(const MAT &m, VEC b) {
    VEC x(b);
    int n = m.dim();
    for (int i = n - 1; i >= 0; --i) {
        x[i] /= m[i][i];
        for (int j = i - 1; j >= 0; --j)
            x[j] -= m[j][i] * x[i];
    }
    return x;
}

VEC luSolve(const MAT &m, VEC b) {
    MAT mt(m);
    VEC y(b.len()), x(b.len());
    luFact(mt);
    y = fwdSubs(mt, b);
    x = bckSubs(mt, y);
    return x;
}

MAT &cholesky(MAT &m) {
    int n = m.dim();
    for (int i = 0; i < n; ++i) {
        m[i][i] = sqrt(m[i][i]);
        for (int j = i + 1; j < n; ++j)
            m[j][i] /= m[i][i];
        for (int j = i + 1; j < n; ++j)
            for (int k = i + 1; k <= j; ++k)
                m[j][k] -= m[j][i] * m[k][i];
    }
    return m;
}

VEC choSolve(const MAT &L, VEC b) {
    int n = L.dim();
    VEC x(b);
    for (int i = 0; i < n; ++i) {
        x[i] /= L[i][i];
        for (int j = i + 1; j < n; ++j)
            x[j] -= x[i] * L[j][i];
    }
    for (int i = n - 1; i >= 0; --i) {
        x[i] /= L[i][i];
        for (int j = i - 1; j >= 0; --j)
            x[j] -= L[i][j] * x[i];
    }
    return x;
}

int jacobi(const MAT &A, const VEC &b, VEC &x, int maxIter, double tol) {
    int n = A.dim();
    int iter_cnt;
    VEC prev(n);
    for (iter_cnt = 1; iter_cnt <= maxIter; ++iter_cnt) {
        prev = x;
        x = b;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j)
                    x[i] -= A[i][j] * prev[j];
        for (int i = 0; i < n; ++i)
            x[i] /= A[i][i];
        // stopping criteria
        if ((prev - x).norm() < tol)
            break;
    }
    return iter_cnt;
}

int gaussSeidel(const MAT &A, const VEC &b, VEC &x, int maxIter, double tol) {
    int n = A.dim();
    int iter_cnt;
    VEC prev(n);
    for (iter_cnt = 1; iter_cnt <= maxIter; ++iter_cnt) {
        prev = x;
        for (int i = 0; i < n; ++i) {
            x[i] = b[i];
            for (int j = 0; j < n; ++j)
                if (i != j)
                    x[i] -= A[i][j] * x[j];
            x[i] /= A[i][i];
        }
        // stopping criteria
        if ((prev - x).norm() < tol)
            break;
    }
    return iter_cnt;
}

int sgs(const MAT &A, const VEC &b, VEC &x, int maxIter, double tol) {
    int n = A.dim();
    int iter_cnt;
    VEC prev(n);
    for (iter_cnt = 1; iter_cnt <= maxIter; ++iter_cnt) {
        prev = x;
        for (int i = 0; i < n; ++i) {
            x[i] = b[i];
            for (int j = 0; j < n; ++j)
                if (i != j)
                    x[i] -= A[i][j] * x[j];
            x[i] /= A[i][i];
        }
        for (int i = n - 1; i >= 0; --i) {
            x[i] = b[i];
            for (int j = 0; j < n; ++j)
                if (i != j)
                    x[i] -= A[i][j] * x[j];
            x[i] /= A[i][i];
        }
        // stopping criteria
        if ((prev - x).norm() < tol)
            break;
    }
    return iter_cnt;
}

int cg(const MAT &A, const VEC &b, VEC &x, int maxIter, double tol) {
    // Conjugate Gradient
    int n = A.dim();
    int iter_cnt;
    VEC p(n), r(n), Ap(n);
    double r2; // r * r
    p = r = b - (A * x);
    r2 = (r * r);
    for (iter_cnt = 1; iter_cnt <= maxIter; ++iter_cnt) {
        Ap = A * p;
        double alpha = r2 / (p * Ap);
        x += p * alpha;
        r -= Ap * alpha;
        double r2old = r2;
        double beta = (r2 = (r * r)) / r2old;
        p = r + (p * beta);
        // stopping criteria
        double err = sqrt(r2 / n);
        if (err < tol)
            break;
    }
    return iter_cnt;
}

int power_method(const MAT &A, VEC &x, double &lambda, int maxIter,
                 double tol) {
    // Power Method
    int n = A.dim();
    x /= x.norm();

    int iter_cnt;
    VEC r(n), Ax(n);
    for (iter_cnt = 1; iter_cnt <= maxIter; ++iter_cnt) {
        x = A * x;
        x /= x.norm();
        Ax = A * x;
        lambda = x * Ax;
        // stopping criteria
        r = Ax - x * lambda;
        if (r.norm() < tol)
            break;
    }
    // std::cout << "res = " << r.norm() << std::endl;
    return iter_cnt;
}

int inv_power_method(const MAT &A, VEC &x, double &lambda, double shift,
                     int maxIter, double tol) {
    // Inverse Power Method with Shifting
    MAT At(A);
    int n = A.dim();
    x /= x.norm();

    // shifting
    for (int i = 0; i < n; ++i)
        At[i][i] -= shift;
    // for inverse iteration
    luFact(At);

    int iter_cnt;
    VEC z(n), r(n), Ax(n);
    for (iter_cnt = 1; iter_cnt <= maxIter; ++iter_cnt) {
        z = fwdSubs(At, x);
        z = bckSubs(At, z);
        x = z / z.norm();
        Ax = A * x;
        lambda = x * Ax;
        // stopping criteria
        r = Ax - x * lambda;
        if (r.norm() < tol)
            break;
    }
    // std::cout << "res = " << r.norm() << std::endl;
    return iter_cnt;
}

void QRDecomp(const MAT &A, MAT &Q, MAT &R) {
    // QR Decomposition
    MAT At = A.tpose();
    R.clear();
    int n = A.dim();
    for (int i = 0; i < n; ++i) {
        Q[i] = At[i];
        for (int j = 0; j < i; ++j) {
            R[j][i] = (Q[j] * Q[i]);
            Q[i] -= Q[j] * R[j][i];
        }
        R[i][i] = Q[i].norm();
        Q[i] /= R[i][i];
    }
    Q = Q.tpose();
}

int EVqr(MAT &A, double tol, int maxiter) {
    // QR Iteration
    int n = A.dim();
    int iter;
    MAT Q(n), R(n);
    for (iter = 1; iter <= maxiter; ++iter) {
        // similar transformation
        QRDecomp(A, Q, R);
        A = R * Q;

        // stopping criteria
        double error = -1;
        for (int i = 1; i < n; ++i)
            error = fmax(error, fabs(A[i][i - 1]));
        if (error < tol)
            break;
    }
    return iter;
}

int EVqrShifted(MAT &A, double mu, double tol, int maxiter) {
    // Shifted QR Iteration
    int n = A.dim();
    int iter;
    MAT Q(n), R(n);
    for (iter = 1; iter <= maxiter; ++iter) {
        // shifting
        for (int i = 0; i < n; ++i)
            A[i][i] -= mu;
        // similar transformation
        QRDecomp(A, Q, R);
        A = R * Q;
        // shifting back
        for (int i = 0; i < n; ++i)
            A[i][i] += mu;

        // stopping criteria
        double error = -1;
        for (int i = 1; i < n; ++i)
            error = fmax(error, fabs(A[i][i - 1]));
        if (error < tol)
            break;
    }
    return iter;
}

double Lagrange(double x, const VEC &XDATA, const VEC &YDATA) {
    // Lagrange interpolation
    double f = 0;
    int n = XDATA.len();
    for (int i = 0; i < n; ++i) {
        double tmp = 1;
        for (int j = 0; j < n; ++j)
            if (j != i)
                tmp *= (x - XDATA[j]) / (XDATA[i] - XDATA[j]);
        f += YDATA[i] * tmp;
    }
    return f;
}

void splineM(int N, const VEC &X, const VEC &Y, VEC &M) {
    // generate spline momentum M
    // check X sorted
    for (int i = 1; i < N; ++i)
        if (X[i - 1] >= X[i]) {
            puts("ERROR: spline X NOT sorted");
            break;
        }
    VEC h(N), lambda(N), mu(N), diag(N);
    // set parameter
    for (int i = 1; i < N; ++i)
        h[i] = X[i] - X[i - 1];
    for (int i = 0; i < N; ++i)
        diag[i] = 2;
    for (int i = 1; i < N - 1; ++i)
        M[i] = 6.0 * ((Y[i + 1] - Y[i]) / h[i + 1] - (Y[i] - Y[i - 1]) / h[i]) /
               (h[i] + h[i + 1]);
    for (int i = 1; i < N - 1; ++i)
        lambda[i] = h[i + 1] / (h[i] + h[i + 1]);
    for (int i = 1; i < N - 1; ++i)
        mu[i] = h[i] / (h[i] + h[i + 1]);
    lambda[0] = M[0] = M[N - 1] = mu[N - 1] = 0;
    // solve linear system
    for (int i = 1; i < N; ++i) {
        double tmp = mu[i] / diag[i - 1];
        M[i] -= M[i - 1] * tmp;
        diag[i] -= lambda[i - 1] * tmp;
    }
    for (int i = N - 2; i >= 0; --i) {
        double tmp = lambda[i] / diag[i + 1];
        M[i] -= M[i + 1] * tmp;
    }
    for (int i = 0; i < N; ++i)
        M[i] /= diag[i];
}
double spline(double x, int N, const VEC &X, const VEC &Y, const VEC &M) {
    // spline interp at x
    if (x < X[0] || x > X[N - 1]) {
        puts("ERROR: spline interpolation x NOT in range");
        return 0;
    }
    // binary search to find X[ind-1] <= x <= X[ind]
    int lower = 1, upper = N - 1;
    while (upper > lower) {
        int mid = (lower + upper) / 2;
        if (x <= X[mid])
            upper = mid;
        else
            lower = mid + 1;
    }
    int ind = lower;
    double hi, alpha, beta, gamma, delta, dx, ans;
    hi = X[ind] - X[ind - 1];
    alpha = Y[ind - 1];
    beta = (Y[ind] - Y[ind - 1]) / hi - hi * (M[ind] + 2 * M[ind - 1]) / 6.0;
    gamma = M[ind - 1] / 2.0;
    delta = (M[ind] - M[ind - 1]) / (6.0 * hi);
    dx = x - X[ind - 1];
    ans = alpha + dx * (beta + dx * (gamma + dx * delta));
    return ans;
}
