// Name: Wang Wei Luan
// Sutdent ID: 102022191
#ifndef MAT_H
#define MAT_H
#include "VEC.h"
class MAT {
  private:
    int n;
    VEC **va;

  public:
    MAT(int);                          // construct with dimension
    MAT(const MAT &m2);                // copy constructor
    MAT(int, double *);                // construct with array
    ~MAT();                            // destructor
    int dim() const;                   // dimension
    MAT tpose() const;                 // transpose
    MAT operator-() const;             // negation
    MAT &operator=(const MAT m2);      // m = m2, assignment
    MAT &operator+=(const MAT m2);     // m += m2
    MAT &operator-=(const MAT m2);     // m -= m2
    MAT &operator*=(double dbl);       // m *= dbl
    MAT &operator/=(double dbl);       // m /= dbl
    MAT operator+(const MAT m2) const; // m + m2
    MAT operator-(const MAT m2) const; // m - m2
    MAT operator*(const MAT m2) const; // matrx-matrix multiplication
    VEC &operator[](int);              // indexing
    const VEC &operator[](int) const;  // indexing
    VEC operator*(const VEC) const;    // vector-matrix multiplication
    MAT operator*(double dbl) const;   // m * dbl
    MAT operator/(double dbl) const;   // m / dbl
    void clear();                      // reset
    static MAT zero(int dim);          // zero matrix of dim
    static MAT identity(int dim,
                        double r);           // r * identity matrix of dim
    friend MAT operator*(double dbl, MAT m); // dbl * m
    friend VEC operator*(const VEC v,
                         const MAT m); // vector-matrix multiplication
};
MAT operator*(double, MAT);
VEC operator*(const VEC, const MAT);
MAT &luFact(MAT &m);               // LU decomposition
VEC fwdSubs(const MAT &m, VEC b);  // forward substitution
VEC bckSubs(const MAT &m, VEC b);  // backward substitution
VEC luSolve(const MAT &m, VEC b);  // LU solve
MAT &cholesky(MAT &m);             // Cholesky decomposition
VEC choSolve(const MAT &L, VEC b); // forward and backward substitutions
int jacobi(const MAT &A, const VEC &b, VEC &x, int maxIter, double tol);
int gaussSeidel(const MAT &A, const VEC &b, VEC &x, int maxIter, double tol);
int sgs(const MAT &A, const VEC &b, VEC &x, int maxIter, double tol);
int cg(const MAT &A, const VEC &b, VEC &x, int maxIter,
       double tol); // Conjugate Gradient
int power_method(const MAT &A, VEC &x, double &lambda, int maxIter,
                 double tol); // Power Method
int inv_power_method(const MAT &A, VEC &x, double &lambda, double shift,
                     int maxIter,
                     double tol); // Inverse Power Method with Shifting
void QRDecomp(const MAT &A, MAT &Q, MAT &R); // QR Decomposition
int EVqr(MAT &A, double tol, int maxiter);   // QR Iteration
int EVqrShifted(MAT &A, double mu, double tol,
                int maxiter); // Shifted QR Iteration
double Lagrange(double x, const VEC &XDATA,
                const VEC &YDATA); // Lagrange interpolation
void splineM(int N, const VEC &X, const VEC &Y,
             VEC &M); // generate spline momentum M
double spline(double x, int N, const VEC &X, const VEC &Y,
              const VEC &M); // spline interp at x
#endif
