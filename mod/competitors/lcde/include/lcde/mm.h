#ifndef LCDE_MM_H_
#define LCDE_MM_H_

#include <iostream>
#include <cassert>

#include "utils.h"

namespace lcde {

static DynamicMatrixT<mpf> inline computeDiff(const VectorT<mpf>& knots, 
                              const MatrixT<mpf>& cpkr, int nk) {
  VectorT<mpf> temp_vec = cpkr.row(2).replicate(1, nk).array() 
              - (knots.replicate(1, nk) + repeatElementwiseN(knots, nk)).array() 
              * cpkr.row(1).replicate(1, nk).array();

  return auxiliaryMatrixAddition(
          (columnWiseMult(tcrossprod(knots), cpkr.row(0).eval())), 
          temp_vec);
}

// Copies the lower diagonal half of a square-shaped matrix into the other
// top diagonal half
// ex) x << 1, 2, 3,
//          4, 5, 6,
//          7, 8, 9
// 
//     ===> 1, 4, 7,
//          4, 5, 8, 
//          7, 8, 9
static void inline copyDiagonal(DynamicMatrixT<mpf>& d) {
  int cols = d.cols(), rows = d.rows();
  assert(cols == rows);
  for (int i = 0; i < cols; ++i) {
    for (int j = i + 1; j < rows; ++j) {
      d(i, j) = d(j, i);
    }
  }
}

}         // namespace lcde

#endif    // LCDE_MM_H_