#include "VEC.h"
#include <cmath>
#include <cstdlib>
#include <cstring>

#ifdef DEBUG
#include <iostream>
#endif

VEC::VEC(int n) {
    dim = n;
    val = (double *)malloc(n * sizeof(double));
}
VEC::VEC(const VEC &v) {
    dim = v.dim;
    val = (double *)malloc(dim * sizeof(double));
    memcpy(val, v.val, dim * sizeof(double));
}
VEC::VEC(int n, double *v) {
    dim = n;
    val = (double *)malloc(dim * sizeof(double));
    memcpy(val, v, dim * sizeof(double));
}
VEC::~VEC() { free(val); }
int VEC::len() const { return dim; }
VEC &VEC::operator-() {
    for (int i = 0; i < dim; ++i)
        val[i] = -val[i];
    return *this;
}
VEC &VEC::operator=(const VEC v) {
    dim = v.dim;
    memcpy(val, v.val, dim * sizeof(double));
    return *this;
}
VEC &VEC::operator+=(const VEC v) {
    for (int i = 0; i < dim; ++i)
        val[i] += v.val[i];
    return *this;
}
VEC &VEC::operator-=(const VEC v) {
    for (int i = 0; i < dim; ++i)
        val[i] -= v.val[i];
    return *this;
}
VEC &VEC::operator*=(double a) {
    for (int i = 0; i < dim; ++i)
        val[i] *= a;
    return *this;
}
VEC &VEC::operator/=(double a) {
#ifdef DEBUG
    if (a == 0)
        std::cout << "divider == 0" << std::endl;
#endif
    for (int i = 0; i < dim; ++i)
        val[i] /= a;
    return *this;
}
VEC VEC::operator+(const VEC v) const {
    VEC s(*this);
    s += v;
    return s;
}
VEC VEC::operator-(const VEC v) const {
    VEC s(*this);
    s -= v;
    return s;
}
double VEC::operator*(const VEC v) const {
    double s = 0;
    for (int i = 0; i < dim; ++i)
        s += val[i] * v.val[i];
    return s;
}
VEC VEC::operator*(double a) const {
    VEC s(*this);
    for (int i = 0; i < dim; ++i)
        s.val[i] *= a;
    return s;
}
VEC VEC::operator/(double a) const {
#ifdef DEBUG
    if (a == 0)
        std::cout << "divider == 0" << std::endl;
#endif
    VEC s(*this);
    s /= a;
    return s;
}
double &VEC::operator[](int n) {
#ifdef DEBUG
    if (n < 0) {
        n = 0;
        std::cout << "index < 0" << std::endl;
    } else if (n >= dim) {
        n = dim - 1;
        std::cout << "index >= dim" << std::endl;
    }
#endif
    return val[n];
}
const double &VEC::operator[](int n) const {
#ifdef DEBUG
    if (n < 0) {
        n = 0;
        std::cout << "index < 0" << std::endl;
    } else if (n >= dim) {
        n = dim - 1;
        std::cout << "index >= dim" << std::endl;
    }
#endif
    return val[n];
}
void VEC::clear() { memset(val, 0, dim * sizeof(double)); }
double VEC::norm() const {
    double ans = 0;
    for (int i = 0; i < dim; ++i)
        ans += val[i] * val[i];
    return sqrt(ans);
}
VEC operator*(double a, VEC v) {
    v *= a;
    return v;
}
VEC *newVEC(int n) {
    VEC *vptr;
    vptr = (VEC *)malloc(sizeof(VEC));
    vptr->dim = n;
    vptr->val = (double *)malloc(n * sizeof(double));
    return vptr;
}
std::ostream &operator<<(std::ostream &out, const VEC &v) {
    out << "[ ";
    for (int i = 0; i < v.len(); ++i)
        out << v[i] << ", ";
    out << "\b\b ]";
    return out;
}
