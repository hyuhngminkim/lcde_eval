#ifndef LCDE_GRADIENT_H_
#define LCDE_GRADIENT_H_

#include <iostream>
#include <vector>
#include <limits>

#include "lcd.h"

namespace lcde {

#define ISX     false
#define ISKNOTS true

static DynamicMatrixT<mpf> cpx(const LCD& lcd, const VectorT<mpf>& x, 
                               int order = 0) {
  const VectorT<mpf>& knots = lcd.knots;
  
  std::vector<int> idx = getIndexOnly(x, knots);
  VectorT<mpf> a = getElementFromIndex(knots, idx);
  VectorT<mpf> intercept = getElementFromIndex(lcd.intercept, idx);
  VectorT<mpf> slope = getElementFromIndex(lcd.slope, idx);
  VectorT<mpf> fkk = getElementFromIndex(lcd.fk, idx);

  DynamicMatrixT<mpf> dpx = DynamicMatrixT<mpf>::Zero(order + 1, x.cols());
  mpf logC = std::log(lcd.C);

  for (int i = 0, cols = x.cols(); i < cols; ++i) {
    if (*(a.data() + i) == *(x.data() + i)) continue;
    else if (*(slope.data() + i) == 0) { // horizontal segment
      dpx(0, i) = (*(x.data() + i) - *(a.data() + i)) * *(fkk.data() + i);
      if (order >= 1) {
        dpx(1, i) = dpx(0, i) * (*(x.data() + i) + *(a.data() + i)) * 0.5;
        if (order >= 2) {
          dpx(2, i) = 
          dpx(0, i) * (pow(*(x.data() + i), 2) + *(x.data() + i) * 
                                 *(a.data() + i) + pow(*(a.data() + i), 2)) / 3;
        }
      }
    } else {                  // non-horizontal segment
      mpf x1 = *(fkk.data() + i) / *(slope.data() + i);
      mpf y1 = std::exp(*(intercept.data() + i) + *(slope.data() + i) * 
              *(x.data() + i) - logC) / *(slope.data() + i);
      dpx(0, i) = y1 - x1;
      if (order >= 1) {
        dpx(1, i) = *(x.data() + i) * y1 - *(a.data() + i) * x1 - dpx(0, i) / *(slope.data() + i);
        if (order >= 2) {
          dpx(2, i) = 
          pow(*(x.data() + i), 2) * y1 - pow(*(a.data() + i), 2) * x1 - dpx(1, i) * 2 / *(slope.data() + i);
        }
      }
    }
  }
  if (knots.cols() > 1) {
    return dpx + 
           getElementFromIndex(cbind(0.0, lcd.cpk.topRows(order + 1)), idx);
  } else return dpx;
}

// Gradient function
// if t = 0 : theta = x
// if t = 1 : theta = lcd.knots
static VectorT<mpf> gradient_lcd(const LCD& lcd, const VectorT<mpf>& x, 
                        const VectorT<int>& w, const VectorT<mpf>& xx, bool t) {
  const VectorT<mpf>& knots = lcd.knots;
  int nk = knots.cols() - 1;
  if (t) {  // theta = lcd.knots
    VectorT<mpf> xxt;
    if (knots.cols() == 2) xxt = unsorted_indx(xx, knots, x);
    else xxt = indx(xx, knots, x);
    VectorT<mpf> px1 = 
    concatenate(0, lcd.cpk.row(0)(Eigen::seq(0, nk - 1)).eval());
    VectorT<mpf> px2 = 
    concatenate(0, lcd.cpk.row(1)(Eigen::seq(0, nk - 1)).eval());
    return ((lcd.cpk(1, nk) - px2.array() 
           - (1 - px1.array()) * knots.array()) * w.sum()) + xxt.array();
  } else {  // theta = x
    const DynamicMatrixT<mpf>& cpx_matrix = cpx(lcd, x, 1);
    return xx.array() + ((lcd.cpk(1, nk) - cpx_matrix.row(1).array() 
                                 - (1 - cpx_matrix.row(0).array()) * x.array()) * w.sum());
  }
}

static VectorT<mpf>
maxima_gradient(LCD& lcd, const VectorT<mpf>& x, const VectorT<int>& w, 
      const VectorT<mpf>& xx, double tol = std::numeric_limits<double>::min()) {
  VectorT<mpf> knots = concatenate(lcd.knots, lcd.upper);
  VectorT<mpf> grad = gradient_lcd(lcd, x, w, xx, ISX);
  std::vector<int> idx = getIndexOnlyFromSorted(knots, x); // index of intervals
  std::vector<int> max_grad_index = 
  getIndexOnlyFromIntervals(grad, idx, true, true);
  VectorT<mpf> theta = getElementFromIndex(x, max_grad_index);
  return theta;
}

}           // namespace lcde

#endif      // LCDE_GRADIENT_H_