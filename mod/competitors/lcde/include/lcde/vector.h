#ifndef LCDE_VECTOR_H_
#define LCDE_VECTOR_H_

#include <iostream>
#include <cmath>
#include <type_traits>
// #include <boost/multiprecision/eigen.hpp>
// #include <boost/multiprecision/gmp.hpp>
#include <Eigen/Dense>
#include <Eigen/Core>

namespace lcde {

// using mpf = boost::multiprecision::mpf_float_50;
// using mpf = long double;
using mpf = double;

// Basic typedefs. 
// Eigen Matrix and Vector types with unspecified type names. 
template <typename T>
using MatrixT = Eigen::Matrix<T, 3, Eigen::Dynamic>;

// Eigen matrix with dynamic row, col and type T. 
template <typename T>
using DynamicMatrixT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
using VectorT = Eigen::Matrix<T, 1, Eigen::Dynamic>;

template <typename T>
static void print_vector(const VectorT<T>& v) {
  std::cout << std::fixed;
  std::cout.precision(7);
  std::cout << v << std::endl;
}

// Performs element-wise multiplication of two different types
// and returns a VectorT of the type of the first vector
// Caution: the multiplication of two vectors of different types will cause
// the resulting vector to be of type T, the type of the lhs vector. 
template <typename T, typename U>
VectorT<T> operator*(const VectorT<T>& lhs, const VectorT<U>& rhs) {
  VectorT<T> r(lhs.cols());
  for (int i = 0; i < lhs.cols(); ++i) {
    r(i) = lhs(i) * rhs(i);
  }
  return r;
}

}   // namespace lcde

#endif    // LCDE_VECTOR_H_