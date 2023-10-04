#ifndef LCDE_LOGLIK_H_
#define LCDE_LOGLIK_H_

#include <algorithm>
#include <numeric>
#include <math.h>

#include "lcd.h"
#include "utils.h"


namespace lcde {
  
// x must be a sorted array
static VectorT<mpf> dlcd(LCD& lcd, VectorT<mpf>& x, bool islog = false) {
  long int x_cols = x.cols();
  long int k_cols = lcd.knots.cols() - 1;
  VectorT<mpf> logd(x_cols);

  int j = 0;
  for (int i = 0; i < x_cols; ++i) {
    // while (j < k_cols && lcd.knots(j + 1) <= x(i)) ++j;
    if (j < k_cols && lcd.knots(j + 1) == x(i)) ++j;
    logd(i) = lcd.intercept(j) + lcd.slope(j) * (x(i)) - std::log(lcd.C); 
  }
  if (islog) {
    return logd;
  } else {
    return logd.array().exp();
  }
}

static void loglik(LCD& lcd, VectorT<mpf> &x, const VectorT<int>& w) {
  VectorT<mpf> r = dlcd(lcd, x, true);
  r = r * w;
  lcd.ll = r.sum();
  return;
}

}         // namespace lcde

#endif    // LCDE_LOGLIK_H_