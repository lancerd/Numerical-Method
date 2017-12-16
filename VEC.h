#ifndef VEC_H
#define VEC_H
#include <iostream>
#include <vector>
class VEC {
  private:
    int dim;
    double *val;

  public:
    VEC(int);                             // construct with dimension
    VEC(const VEC &v2);                   // copy constructor
    VEC(int, double *);                   // construct with array
    ~VEC();                               // destructor
    int len() const;                      // dimension
    VEC &operator-();                     // negation
    VEC &operator=(const VEC v2);         // v = v2, assignment
    VEC &operator+=(const VEC v2);        // v += v2
    VEC &operator-=(const VEC v2);        // v -= v2
    VEC &operator*=(double dbl);          // v *= dbl
    VEC &operator/=(double dbl);          // v /= dbl
    VEC operator+(const VEC v2) const;    // v + v2
    VEC operator-(const VEC v2) const;    // v - v2
    double operator*(const VEC v2) const; // inner product
    VEC operator*(double dbl) const;      // v * dbl
    VEC operator/(double dbl) const;      // v / dbl
    double &operator[](int);              // indexing
    const double &operator[](int) const;  // indexing
    void clear();                         // reset
    void fill(double dbl);                // fill all element with dbl
    double norm() const;
    friend VEC operator*(double dbl, VEC v); // dbl * v
    friend VEC *newVEC(int);                 // aloocate new vec
};
VEC operator*(double, VEC);
VEC *newVEC(int);
std::ostream &operator<<(std::ostream &out, const VEC &v);
#endif
