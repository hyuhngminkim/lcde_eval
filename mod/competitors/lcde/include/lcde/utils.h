#ifndef LCDE_UTILS_H_
#define LCDE_UTILS_H_

#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <functional>
// https://en.cppreference.com/w/cpp/types/numeric_limits/infinity
#include <limits>
#include <iostream>
#include <type_traits> // for type checking
#include <iomanip>

#include "vector.h"

// This header file implements basic / utility functions
// required in the code. 

namespace lcde {

// https://www.appsloveworld.com/cplus/100/92/c-cast-template
// Usage :
// U some_u;
// T some_t = cast<T>(some_u);
template <typename T, typename U>
static T inline cast(U& u) {
  return static_cast<T>(u);
}

template <typename T>
void print_vector(const std::vector<T> v) {
  for (auto i : v) std::cout << i << ", ";
  std::cout << std::endl;
  return;
}

// Concatenates two vectors
// Requires the two vectors to be of same type
template <typename T>
VectorT<T> inline concatenate(const VectorT<T>& a, const VectorT<T>& b) {
  VectorT<T> r = a;
  r.conservativeResize(1, a.cols() + b.cols());
  r.block(0, a.cols(), 1, b.cols()) = b;
  return r;
}

// Two versions of concatenating a VectorT and an element
// 1. element + VectorT
// 2. VectorT + element
template <typename T, typename U>
VectorT<T> inline concatenate(const U& val, const VectorT<T>& vec) {
  VectorT<T> r(vec.cols() + 1);
  r(0) = val;
  r.block(0, 1, 1, vec.cols()) = vec;
  return r;
}

template <typename T, typename U>
VectorT<T> inline concatenate(const VectorT<T>& vec, const U& val) {
  VectorT<T> r(vec.cols() + 1);
  r.block(0, 0, 1, vec.cols()) = vec;
  r(vec.cols()) = val;
  return r;
}

// Prepend n zeros at the front of vector
template <typename T>
VectorT<T> inline prependN(int n, const VectorT<T>& vec) {
  VectorT<T> zeros = VectorT<T>::Zero(n);
  return concatenate(zeros, vec);
}

// Append n zeros at the end of vector
template <typename T>
VectorT<T> inline appendN(const VectorT<T>& vec, int n) {
  VectorT<T> zeros = VectorT<T>::Zero(n);
  return concatenate(vec, zeros);
}


// appends element val at the end of the vector
template <typename T>
VectorT<T> inline append(const VectorT<T>& v, const T& val) {
  VectorT<T> r(v.cols() + 1);
  r << v, val;
  return r;
}

// must be used only when the vector is mpf and the appended value can be 
// automatically converted into mpf
template <typename T, typename U>
VectorT<T> inline append(const VectorT<T>& v, const U& val) {
  VectorT<T> r(v.cols() + 1);
  r << v, val;
  return r;
}

// prepends element val at the end of the vector
template <typename T>
VectorT<T> inline prepend(const VectorT<T>& v, const T& val) {
  VectorT<T> r(v.cols() + 1);
  r << val, v;
  return r;
}

// cbind functions
// make a column-wise concatenation between scalars <-> matrices
// cbind does not perform type casting - the input scalar values must be cast into 
// whatever type the target matrix takes. 
template <typename T, typename U>
MatrixT<T> inline cbind(const U& val, const MatrixT<T>& m) {
  MatrixT<T> r(3, m.cols() + 1);
  VectorT<T> aux = VectorT<T>::Constant(3, val);
  r << aux.transpose(), m;
  return r;
}

template <typename T, typename U>
MatrixT<T> inline cbind(const MatrixT<T>& m, const U& val) {
  MatrixT<T> r(3, m.cols() + 1);
  VectorT<T> aux = VectorT<T>::Constant(3, val);
  r << m, aux.transpose();
  return r;
}

template <typename T, typename U>
DynamicMatrixT<T> inline cbind(const U& val, const DynamicMatrixT<T>& m) {
  DynamicMatrixT<T> r(m.rows(), m.cols() + 1);
  VectorT<T> aux = VectorT<T>::Constant(m.rows(), val);
  r << aux.transpose(), m;
  return r;
}

template <typename T, typename U>
DynamicMatrixT<T> inline cbind(const DynamicMatrixT<T>& m, const U& val) {
  DynamicMatrixT<T> r(m.rows(), m.cols() + 1);
  VectorT<T> aux = VectorT<T>::Constant(m.rows(), val);
  r << m, aux.transpose();
  return r;
}

// The following cbind functions take Eigen::Blocks 
// instead of matrices as parameters
template <typename T, typename U>
DynamicMatrixT<T> inline cbind(const U& val, 
                               const Eigen::Block<const MatrixT<T>>& m) {
  DynamicMatrixT<T> r(m.rows(), m.cols() + 1);
  VectorT<T> aux = VectorT<T>::Constant(m.rows(), val);
  r << aux.transpose(), m;
  return r;
}

template <typename T, typename U>
DynamicMatrixT<T> inline cbind(const Eigen::Block<const MatrixT<T>>& m, 
                               const U& val) {
  DynamicMatrixT<T> r(m.rows(), m.cols() + 1);
  VectorT<T> aux = VectorT<T>::Constant(m.rows(), val);
  r << m, aux.transpose();
  return r;
}

template <typename T, typename U>
DynamicMatrixT<T> inline cbind(const U& val, 
                               const Eigen::Block<const DynamicMatrixT<T>>& m) {
  DynamicMatrixT<T> r(m.rows(), m.cols() + 1);
  VectorT<T> aux = VectorT<T>::Constant(m.rows(), val);
  r << aux.transpose(), m;
  return r;
}

template <typename T, typename U>
DynamicMatrixT<T> inline cbind(const Eigen::Block<const DynamicMatrixT<T>>& m, 
                               const U& val) {
  DynamicMatrixT<T> r(m.rows(), m.cols() + 1);
  VectorT<T> aux = VectorT<T>::Constant(m.rows(), val);
  r << m, aux.transpose();
  return r;
}

// Makes a slice from the start_index to the last element
// Simply an alias of seq(start_index, last)
template <typename T>
VectorT<T> inline sliceToLast(const VectorT<T>& v, const int& start_index) {
  assert(v.cols() > start_index);
  return v(Eigen::seq(start_index, Eigen::last));
}

// Makes a slice from the first element to finish_index
// Simply an alias of seq(0, finish_index)
template <typename T>
VectorT<T> inline sliceFromFirst(const VectorT<T>& v, const int& finish_index) {
  assert(finish_index < v.cols());
  return v(Eigen::seq(0, finish_index));
}

// Makes a slice from start_index to finish_index
// Simply an alias of seq(start_index, finish_index)
template <typename T>
VectorT<T> inline sliceRange(const VectorT<T>& v, const int& s, const int& f) {
  assert(v.cols() > s && f < v.cols() && s <= f);
  return v(Eigen::seq(s, f));
}

template <typename T>
VectorT<T> inline dropElementByIndex(const VectorT<T>& v, const int& idx) {
  int cols = v.cols();
  assert(idx >= 0 && idx < cols);
  VectorT<T> res(cols - 1);
  int res_idx = 0;
  for (int i = 0; i < cols; ++i) {
    if (i == idx) continue;
    *(res.data() + (res_idx++)) = *(v.data() + i);
  }
  return res;
}

template <typename T>
MatrixT<T> inline dropElementByIndex(const MatrixT<T>& m, const int& idx) {
  assert(idx >= 0 && idx < m.cols());
  MatrixT<T> res(3, m.cols() - 1);
  res << m.block(0, 0, 3, idx), m.block(0, idx + 1, 3, m.cols() - idx - 1);
  return res;
}

// Requires:
// 1. x and v must be sorted
// 2. x is contained within v (i.e. v(0) <= x(0) && x(x.cols() - 1) <= v(v.cols() - 1))
// 3. v.cols() == target.cols()
template <typename T, typename U>
VectorT<T> inline indx(const VectorT<U>& target, const VectorT<T>& x, 
                       const VectorT<T>& v) {
  assert(target.cols() == v.cols());
  int idx = 0, cols = x.cols();
  VectorT<U> r(cols);
  for (long int i = 0; i < cols; ++i) {
    assert(*(x.data() + i) >= v(0) && *(x.data() + i) <= v(v.cols() - 1));
    while (*(x.data() + i) >= *(v.data() + idx)) ++idx;
    *(r.data() + i) = *(target.data() + idx - 1);
  }
  return r;
}

template <typename T, typename U>
VectorT<T> inline unsorted_indx(const VectorT<U>& target, const VectorT<T>& x, 
                                const VectorT<T>& v) {
  assert(target.cols() == v.cols());
  VectorT<U> r(x.cols());
  int i, j;
  auto head = v(0);
  auto tail = v(v.cols() - 1);
  for (i = 0; i < x.cols(); ++i) {
    if (x(i) <= head) r(i) = target(0);
    else if (x(i) >= tail) r(i) = target(v.cols() - 1);
    else {
      for (j = 0; j < v.cols(); ++j) {
        if (x(i) < v(j)) {
          r(i) = target(j - 1);
          break;
        }
      }
    }
  }
  return r;
}

// Returns the index of x with respect to v, i.e. for every x_i of x, 
// getIndexOnly() returns j such that v_j <= x_i < v_(j+1)
template <typename T>
std::vector<int> getIndexOnly(const VectorT<T>& x, const VectorT<T>& v) {
  std::vector<int> r(x.cols());
  int vcols = v.cols();
  int i, j;
  auto head = v(0);
  auto tail = v(vcols - 1);
  for (i = 0; i < x.cols(); ++i) {
    if (x(i) <= head) r[i] = 0;
    else if (x(i) >= tail) r[i] = vcols - 1;
    else {
      for (j = 0; j < vcols; ++j) {
        if (x(i) < v(j)) {
          r[i] = j - 1;
          break;
        }
      }
    }
  }
  return r;
}

// Same function as getIndexOnly but assumes that `v` is sorted, so performs a 
// single iteration over `v`. 
template <typename T>
std::vector<int> getIndexOnlyFromSorted(const VectorT<T>& x, 
                                        const VectorT<T>& v) {
  std::vector<int> r(x.cols());
  int vcols = v.cols();
  int i, idx = 0;
  auto head = v(0);
  auto tail = v(vcols - 1);
  for (i = 0; i < x.cols(); ++i) {
    if (x(i) <= head) r[i] = 0;
    else if(x(i) >= tail) r[i] = vcols - 1;
    else {
      for (; idx < vcols; ++idx) {
        if (x(i) < v(idx)) {
          r[i] = idx - 1;
          break;
        }
      }
    }
  }
  return r;
}

// Auxiliary function for getIndexOnlyFromInterval
// Receives a seq() as input (i.e. an Eigen::Block) 
// and returns the index of the max coef
template <typename T>
int inline getMaxIndexFromRange(const VectorT<T>& v, int start_index, 
                                int end_index) {
  int max_idx = start_index;
  T max_coef = v(start_index);
  for (int i = start_index; i <= end_index; ++i) {
    if (*(v.data() + i) > max_coef) {
      max_idx = i;
      max_coef = *(v.data() + i);
    }
  }
  return max_idx;
}

// The vector of interval indices is assumed to cover the target vector to the end. 
// That is, the vector `idx` must contain the last element of `vec` as its last index. 
// For example, if vec = {1, 2, 3, 4, 5}, then for whatever interval `idx` represents, 
// `idx` must contain index 0 and 4. 
// getIndexOnlyFromIntervals returns the index that holds the largest coefficient for each
// interval provided by `idx`. Therefore, the return vector has one less element than `idx`. 
template <typename T>
std::vector<int> getIndexOnlyFromIntervals(const VectorT<T>& vec, 
                                           const std::vector<int> idx) {
  size_t idx_size = idx.size();
  std::vector<int> r(idx_size - 1);
  for (size_t i = 0; i < idx_size - 1; ++i) {
    r[i] = getMaxIndexFromRange(vec, idx[i], idx[i + 1]);
  }
  return r;
}

// A sly modification of getIndexOnlyFromIntervals - takes a boolean parameter `isInterleave`
// which indicates that elements in idx must be preserved in building the result vector
template <typename T>
std::vector<int> getIndexOnlyFromIntervals(const VectorT<T>& vec, 
                                           const std::vector<int> idx,
                                           bool isInterleave,
                                           bool noSides) {
  if (isInterleave) {
    std::vector<int> r;
    size_t idx_size = idx.size();
    int temp;
    for (size_t i = 0; i < idx_size - 1; ++i) {
      r.push_back(idx[i]);
      temp = getMaxIndexFromRange(vec, idx[i], idx[i + 1]);
      if (temp > idx[i] && temp < idx[i + 1]) r.push_back(temp);
    }
    if (noSides) {
      r.erase(r.begin());
    } else {
      r.push_back(idx[idx_size - 1]);
    }
    return r;
  } else {
    return getIndexOnlyFromIntervals(vec, idx);
  }
}

// Returns the elements of the `target` VectorT from the indices retrieved from
// getIndexOnly(). Best if the input `idx` is created from getIndexOnly()
template <typename T>
VectorT<T> inline getElementFromIndex(const VectorT<T>& target, 
                                      const std::vector<int>& idx) {
  return target(0, idx);
}

template <typename T>
MatrixT<T> inline getElementFromIndex(const MatrixT<T>& target, 
                                      const std::vector<int>& idx) {
  return target(Eigen::all, idx);
}

template <typename T>
DynamicMatrixT<T> inline getElementFromIndex(const DynamicMatrixT<T>& target, 
                                             const std::vector<int>& idx) {
  return target(Eigen::all, idx);
}

// Function designed specifically for indexing from cpk of lcd
template <typename T>
DynamicMatrixT<T> 
inline getElementFromIndex(const Eigen::Block<MatrixT<T>>& target, 
                           const std::vector<int>& idx) {
  DynamicMatrixT<T> res(target.rows(), idx.size());
  res << target(Eigen::all, idx);
  return res;
}

// cumulative sum
// If zero = true, appends a zero at the front of the result
template <typename T>
VectorT<T> inline cumsum(const VectorT<T>& v, bool zero) {
  if (zero) {
    VectorT<T> r = VectorT<T>::Zero(v.cols() + 1);
    std::partial_sum(v.begin(), v.end(), r.begin() + 1, std::plus<T>());
    return r;
  } else {
    VectorT<T> r = VectorT<T>::Zero(v.cols());
    std::partial_sum(v.begin(), v.end(), r.begin(), std::plus<T>());
    return r;
  }
}

template <typename T>
VectorT<T> inline cumsum(const VectorT<T>& v, bool zero, bool neg) {
  int cols = v.cols();
  VectorT<T> r(zero ? cols + 1 : cols);
  r(0) = zero ? 0 : (neg ? -1 * v(0) : v(0));
  int step = zero ? 1 : 0;
  int limit = zero ? cols + 1 : cols;
  if (neg) {
    for (int i = 1; i < limit; ++i) {
      r(i) = r(i - 1) - v(i - step);
    }
  } else { 
    for (int i = 1; i < limit; ++i) {
      r(i) = r(i - 1) + v(i - step);
    }
  }
  return r;
}

template <typename T>
VectorT<T> inline revCumsum(const VectorT<T>& v) {
  int cols = v.cols();
  if (cols <= 1) return v;
  VectorT<T> r(cols);
  r(cols - 1) = v(cols - 1);
  for (int i = cols - 2; i >= 0; --i) {
    r(i) = *(v.data() + i) + r(i + 1);
  }
  return r;
}

// rev(cumsum(rev(w))) * x - rev(cumsum(rev(x * w)))
// vector<int> x = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
// vector<int> w = {1, 2, 2, 1, 4, 6, 1, 2, 7, 10};
// C++ : 240ns 
// R   : 5860ns
template <typename T>
void inline xsq(VectorT<T>& target, const VectorT<T>& x, 
                const VectorT<int>& w) {
  size_t cols = x.cols();
  T lhs = *(w.data() + cols - 1);
  T rhs = *(x.data() + cols - 1) * *(w.data() + cols - 1);
  for (int i = cols - 2; i >= 0; --i) {
    lhs += *(w.data() + i);
    rhs += *(x.data() + i) * *(w.data() + i);
    *(target.data() + i) = (lhs * *(x.data() + i)) - rhs;
  }
  return ;
}

// Requires the vectors x and w to be allocated memory outside of this function
template <typename T>
void inline weight(const std::vector<T>& input, VectorT<T>& x, 
                   VectorT<int>& w) {
  const size_t input_size = input.size();

  size_t idx = 0;
  T max_val = std::numeric_limits<T>::lowest();
  for (size_t i = 0; i < input_size; ++i) {
    if (input[i] > max_val) {
      max_val = input[i];
      *(x.data() + idx) = max_val;
      *(w.data() + idx) = 1;
      ++idx;
    } else {
      ++(*(w.data() + idx - 1));
    }
  }

  // Must use conservativeResize() since for size reduction, resize() removes
  // all existing data whereas conservativeResize() retains data
  x.conservativeResize(1, idx);
  w.conservativeResize(1, idx);
  return;
}

// Requires the vectors x and w to be allocated memory outside of this function
template <typename T, typename U>
void inline weight(const std::vector<T>& input, VectorT<U>& x, 
                   VectorT<int>& w) {
  const size_t input_size = input.size();

  size_t idx = 0;
  U max_val = std::numeric_limits<U>::lowest();
  for (size_t i = 0; i < input_size; ++i) {
    if (input[i].x > max_val) {
      max_val = input[i].x;
      *(x.data() + idx) = max_val;
      *(w.data() + idx) = 1;
      ++idx;
    } else {
      ++(*(w.data() + idx - 1));
    }
  }

  // Must use conservativeResize() since for size reduction, resize() removes
  // all existing data whereas conservativeResize() retains data
  x.conservativeResize(1, idx);
  w.conservativeResize(1, idx);
  return;
}

// Requires the vectors x and w to be allocated memory outside of this function
void inline weight(std::vector<mpf>::iterator begin, std::vector<mpf>::iterator end,
                   VectorT<mpf>& x, VectorT<int>& w) {
  size_t idx = 0;
  mpf max_val = std::numeric_limits<mpf>::lowest();

  for (auto i = begin; i < end; ++i) {
    if (*i > max_val) {
      max_val = *i;
      *(x.data() + idx) = max_val;
      *(w.data() + idx) = 1;
      ++idx;
    } else {
      ++(*(w.data() + idx - 1));
    }
  }

  // Must use conservativeResize() since for size reduction, resize() removes
  // all existing data whereas conservativeResize() retains data
  x.conservativeResize(1, idx);
  w.conservativeResize(1, idx);
  return;
}

// Apply exponential to all elements in vector
template <typename T>
VectorT<T> inline exp(const VectorT<T> &v) {
  int cols = v.cols();
  VectorT<T> r(cols);
  for (int i = 0; i < cols; ++i) r(i) = exp(*(v.data() + i));
  return r;
}

template <typename T>
DynamicMatrixT<T> tcrossprod(const VectorT<T>& v) {
  return v.transpose() * v;
}

// Repeats the vector v n times. 
// ex) v = {1, 2, 3}, n = 2 => {1, 2, 3, 1, 2, 3}
template <typename T>
VectorT<T> inline repeatN(const VectorT<T>& v, int n) {
  return v.replicate(1, n);
}

// Repeats each element in v n times. 
// ex) v = {1, 2, 3}, n = 2 => {1, 1, 2, 2, 3, 3}
template <typename T>
VectorT<T> repeatElementwiseN(const VectorT<T>& v, int n) {
  int v_cols = v.cols();
  VectorT<T> r(v_cols * n);
  for (int i = 0; i < v_cols; ++i) {
    for (int j = 0; j < n; ++j) {
      r(0, i * n + j) = *(v.data() + i);
    }
  }
  return r;
}

template <typename T>
DynamicMatrixT<T> inline columnWiseMult(const DynamicMatrixT<T>& d, 
                                        const VectorT<T>& v) {
  // int cols = d.cols(), rows = d.rows();
  // assert(rows == v.cols());
  // DynamicMatrixT<T> r(d.rows(), d.cols());
  // for (int i = 0; i < cols; ++i) {
  //   for (int j = 0; j < rows; ++j) {
  //     r(j, i) = d(j, i) * v(j);
  //   }
  // }
  return d.array() * v.transpose().replicate(1, d.cols()).array();
  // return r;
}

// Same function but takes std::vector as input
template <typename T>
DynamicMatrixT<T> columnWiseMult(const DynamicMatrixT<T>& d, 
                                 const std::vector<T>& v) {
  unsigned int cols = d.cols(), rows = d.rows();
  assert(rows == v.size());
  DynamicMatrixT<T> r(d.rows(), d.cols());
  for (unsigned int i = 0; i < cols; ++i) {
    for (unsigned int j = 0; j < rows; ++j) {
      r(j, i) = d(j, i) * v[j];
    }
  }
  return r;
}

template <typename T, typename U>
DynamicMatrixT<T> columnWiseMult(const U& d, const std::vector<T>& v) {
  unsigned int cols = d.cols(), rows = d.rows();
  assert(rows == v.size());
  DynamicMatrixT<T> r(d.rows(), d.cols());
  for (unsigned int i = 0; i < cols; ++i) {
    for (unsigned int j = 0; j < rows; ++j) {
      r(j, i) = d(j, i) * v[j];
    }
  }
  return r;
}

// auxiliary function for matrix division
// the left hand side argument is assumed to be a VectorT
template <typename T, typename U>
VectorT<T> auxiliaryDivision(const U& lhs, const std::vector<T>& rhs) {
  unsigned int cols = lhs.cols();
  assert(cols == rhs.size());
  VectorT<T> r(cols);
  for (unsigned int i = 0; i < cols; ++i) {
    r(i) = lhs(i) / rhs[i];
  }
  return r;
}

template <typename T>
DynamicMatrixT<T> auxiliaryMatrixAddition(const DynamicMatrixT<T>& d,
                                             const VectorT<T>& v) {
  int cols = d.cols(), rows = d.rows();
  assert(rows * rows == v.cols());
  DynamicMatrixT<T> r(rows, cols);
  for (int i = 0; i < cols; ++i) {
    for (int j = 0, v_idx = rows * i; j < rows; ++j) {
      r(j, i) = d(j, i) + v(v_idx + j);
    }
  }
  return r;
}

template <typename T>
DynamicMatrixT<T> auxiliaryMatrixSubtraction(const DynamicMatrixT<T>& d,
                                             const VectorT<T>& v) {
  int cols = d.cols(), rows = d.rows();
  assert(rows * rows == v.cols());
  DynamicMatrixT<T> r(rows, cols);
  for (int i = 0; i < cols; ++i) {
    for (int j = 0, v_idx = rows * i; j < rows; ++j) {
      r(j, i) = d(j, i) - v(v_idx + j);
    }
  }
  return r;
}

// static double** cast(const DynamicMatrixT<mpf>& d) {
//   int cols = d.cols(), rows = d.rows();
//   double** res = new double*[rows];
//   for (int i = 0; i < rows; ++i) {
//     res[i] = new double[cols];
//     for (int j = 0; j < cols; ++j) {
//       // res[i][j] = d(i, j).convert_to<double>();
//       res[i][j] = static_cast<double>(d(i,j));
//     }
//   }
//   return res;
// }

// static double* cast(const VectorT<mpf>& v) {
//   int cols = v.cols();
//   double* res = new double[cols];
//   for (int i = 0; i < cols; ++i) {
//     // res[i] = v(i).convert_to<double>();
//     res[i] = static_cast<double>(*(v.data() + i));
//   }
//   return res;
// }

// Fast approximation of the exponential function
// https://github.com/simonpf/fastexp
template <typename T, size_t degree, size_t i>
struct ExpRecursion {
  static T evaluate(T x) {
    x = ExpRecursion<T, degree, i + 1>::evaluate(x);
    return x * x;
  }
};

template <typename T, size_t degree>
struct ExpRecursion<T, degree, degree> {
  static T evaluate(T x) {
    constexpr T c = 1.0 / static_cast<T>(1u << degree);
    x = 1.0 + c * x;
    return x;
  }
};

template <typename T>
T expApprox(T x) {
  return ExpRecursion<T, 30, 0>::evaluate(x);
}

bool inline vectorIsInvalid(const VectorT<mpf>& vec) {
  for (auto& el : vec) {
    if (!std::isfinite(el) || el == 0) return true;
  }
  return false;
}

}         // namespace lcde

#endif    // LCDE_UTILS_H_