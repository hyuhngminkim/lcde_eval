#ifndef LCDE_LINE_LCD_H_
#define LCDE_LINE_LCD_H_

#include <iostream>
#include <vector>

#include "utils.h"
#include "lcd.h"
#include "gradient.h"
#include "loglik.h"

namespace lcde {

/*!
 * @name          line_lcd
 * @brief         Performs line search using Armijo's rule
 * @details       The line search is complemented by Armijo's rule to ensure a
 *                monotone yet sufficiently large increase of the log-likelihood
 *                after each iteration. The NNLS operation is performed inside
 *                this function, since performing NNLS in the main cnmlcd would
 *                
 * @param   lcd   Original lcd object formulated during the iteration
 * @param   x     Target dataset
 * @param   w     Weight of the target dataset
 * @param   xx    Cumulative cardinality
 */
static void line_lcd(LCD& lcd, VectorT<mpf>& x, const VectorT<int>& w, 
                     const VectorT<mpf>& xx, const VectorT<double>& nnls_d, 
                     const mpf& ll0) {
  VectorT<mpf> grad = gradient_lcd(lcd, x, w, xx, ISKNOTS);
  VectorT<mpf> nnls = nnls_d.cast<mpf>();
  grad(0) = -grad(0);
  mpf delta = 
  (grad.array() * (nnls - concatenate(lcd.alpha, lcd.pi)).array()).sum() * 0.333;
  mpf alpha = 1;
  mpf tol = 1e-10;
  auto nnls_len = nnls.cols();
  mpf nnls_alpha = nnls(0);
  VectorT<mpf> nnls_pi = nnls.tail(nnls_len - 1);
  mpf prev_ll = lcd.ll;
  while (true) {
    VectorT<mpf> npa = nnls_pi * alpha;
    VectorT<mpf> lpa = (1 - alpha) * lcd.pi;
    VectorT<mpf> qwerty = lpa + npa;

    LCD temp_lcd = LCD((1 - alpha) * lcd.alpha + alpha * nnls_alpha,
                        lcd.lower, 
                        lcd.upper, 
                        lcd.theta, 
                        qwerty);
    // ERROR
    if (!std::isfinite(temp_lcd.C)) return;
    loglik(temp_lcd, x, w);
    if (prev_ll >= temp_lcd.ll) {
      lcd = temp_lcd;
      return;
    }
    if (temp_lcd.ll >= ll0 + alpha * delta) {
      lcd = temp_lcd;
      return;
    }
    if (alpha < tol) {
      lcd = temp_lcd;
      return;
    }
    alpha *= 0.5;
    prev_ll = temp_lcd.ll;
  }
}

}         // namespace lcde

#endif    // LCDE_LINE_H_