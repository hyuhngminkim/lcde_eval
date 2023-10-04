#ifndef LCDE_PNNLS_H_
#define LCDE_PNNLS_H_

#include <iostream>
#include <cassert>

#include "vector.h"
#include "utils.h"

namespace lcde {

extern "C" {
  void pnnls_(double** fa_, int* fmda_, int* fm_, int* fn_, double* fb_, 
              double* fx_, double* frnorm_, double* fw_, double* fzz_,
              int* findex_, int* fmode_, int* fk_);
}

/*!
 * @name        pnnls
 * @brief       Computes the non-negative least squares
 * @details     A C wrapper function for the pnnls fortran subroutine by 
 *              Yong Wang, Charles L. Lawson, and Richard J. Hanson. 
 * @param   a   The M by N matrix
 * @param   b   The M vector
 * @param   k   The first k variables are not NN-restricted
 */
static VectorT<double> pnnls(DynamicMatrixT<double> a, VectorT<double> b, int k) {
  int m = a.rows();
  int n = a.cols();
  assert(m == b.cols());

  // https://stackoverflow.com/questions/54282732/c-convert-array-of-elements-into-2-d-matrix
  double (*fa)[m] = (double(*)[m])(a.data());
  double* fb = b.data();

  double x[n];            // only for output
  double rnorm[1];        // only for output

  double w[n];            // n-vector of working space
  double zz[m];           // m-vector of working space
  int index[n];           // n-vector index, only for output
  int mode[1];            // success/failure flag, 1 = success

  pnnls_((double**)fa, &m, &m, &n, fb, x, rnorm, w, zz, index, mode, &k);
  
  Eigen::Map<VectorT<double>> res(x, n);
  return res;
}

}           // namespace lcde

#endif      // LCDE_PNNLS_H_