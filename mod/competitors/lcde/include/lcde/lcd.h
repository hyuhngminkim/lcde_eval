#ifndef LCDE_LCD_H_
#define LCDE_LCD_H_

#include <vector>
#include <math.h>

#include "utils.h"

namespace lcde {

class LCD {
  public:
  mpf alpha;                // initial slope on [x_1, x_2]
  mpf C;                    // normalizing constant
  mpf ll;                   // log likelihood of this lcd object
  mpf lower;                // lower boundary (= x_1)
  mpf upper;                // upper boundary (= x_n)
  VectorT<mpf> theta;       // interior knots (knots other than x_1 and x_n)
  VectorT<mpf> knots;       // lower + theta // = VectorT<T>(2);
  VectorT<mpf> pi;          // changes of slope at theta
  VectorT<mpf> intercept;   // intercepts over all segments
  VectorT<mpf> slope;       // slopes over all segments
  VectorT<mpf> fk;          // density at knots [lower, theta]
  MatrixT<mpf> dpk;         // Int x^o f(x) dx over each segment between knots 
                            // for o = 0, 1, 2
  MatrixT<mpf> cpk;         // cumulative sum of dpk for all segments

  LCD() = default;

  LCD& operator=(const LCD& rhs) {
    if (this != &rhs) {
      alpha = rhs.alpha;
      C = rhs.C;
      ll = rhs.ll;
      lower = rhs.lower;
      upper = rhs.upper;
      theta = rhs.theta;
      knots = rhs.knots;
      pi = rhs.pi;
      intercept = rhs.intercept;
      slope = rhs.slope;
      fk = rhs.fk;
      dpk = rhs.dpk;
      cpk = rhs.cpk;
    }
    return *this;
  }

  LCD(mpf a /* alpha */, mpf l /* lower */, mpf u /* upper */, 
    VectorT<mpf>& t /* theta */, VectorT<mpf>& p /* pi */) {
    if (t.cols() > 0) {
      if (p.cols() < t.cols()) p.resize(t.cols());
    } 
    alpha = a;
    lower = l;
    upper = u;
    theta = t;
    pi = p;
    knots = prepend(theta, lower);

    intercept = cumsum(theta * pi, true).array() - (alpha * lower); // intercepts
    slope = cumsum(pi, true, true).array() + a;                     // slopes

    VectorT<mpf> knots2 = append(theta, upper);

    VectorT<mpf> dk = knots2 - knots;
    int nk = dk.cols();
    VectorT<mpf> fk1 = ((knots * slope) + intercept).array().exp();
    VectorT<mpf> fk2 = append(dropElementByIndex(fk1, 0), 
     std::exp(*(intercept.data() + nk - 1) + *(slope.data() + nk - 1) * upper));

    dpk.conservativeResize(3, nk);
    for (size_t i = 0, cols = slope.cols(); i < cols; ++i) {
      if (*(slope.data() + i) == 0) {     // horizontal segments
        dpk(0, i) = fk1(i) * dk(i);
        dpk(1, i) = dpk(0, i) * (knots2(i) + knots(i)) * 0.5;
        dpk(2, i) = dpk(0, i) * 
           (pow(knots2(i), 2) + knots2(i) * knots(i) + pow(knots(i), 2)) / 3;
      } else {                            // nonhorizontal segments
        mpf x1 = fk1(i) / slope(i);
        mpf y1 = fk2(i) / slope(i);
        dpk(0, i) = y1 - x1;
        dpk(1, i) = knots2(i) * y1 - knots(i) * x1 - dpk(0, i) / slope(i);
        dpk(2, i) = pow(knots2(i), 2) * y1
                  - pow(knots(i), 2) * x1 - 2 * dpk(1, i) / slope(i);
      }
    }

    C = dpk.row(0).sum();
    fk = fk1 / C;

    cpk = dpk = dpk / C;

    for (size_t i = 1, cpk_cols = cpk.cols(); i < cpk_cols; ++i) {
      cpk(0, i) += cpk(0, i - 1);
      cpk(1, i) += cpk(1, i - 1);
      cpk(2, i) += cpk(2, i - 1);
    }
  }

  LCD(bool dummy) {

  }

  void simplify() {
    std::vector<int> nzi_base;        // non zero index
    std::vector<int> nzi_pk;          // non zero index for dpk and cpk
    std::vector<int> nzi_coef(1, 0);  // non zero index for intercept and slope

    for (int i = 0, cols = pi.cols(); i < cols; ++i) {
      if (pi(i) != 0) nzi_base.push_back(i);
    }

    pi = getElementFromIndex(pi, nzi_base);
    theta = getElementFromIndex(theta, nzi_base);

    nzi_pk = nzi_base;
    nzi_pk.push_back(theta.cols() + 1);
    dpk = getElementFromIndex(dpk, nzi_pk);
    cpk = getElementFromIndex(cpk, nzi_pk);

    std::for_each(nzi_base.begin(), nzi_base.end(), [](int &n){ n+=1; });
    nzi_coef.insert(nzi_coef.end(), nzi_base.begin(), nzi_base.end());

    intercept = getElementFromIndex(intercept, nzi_coef);
    slope = getElementFromIndex(slope, nzi_coef);
    fk = getElementFromIndex(fk, nzi_coef);

    knots = prepend(theta, lower);
  }

  void print_lcd() const {
    std::cout << "\n\033[1;91;47m Printing out class elements for testing purposes \033[m\n";
    std::cout << std::setprecision(9);
    std::cout << "\033[1;31mLOWER   \033[0m" << lower << std::endl;
    std::cout << "\033[1;31mUPPER   \033[0m" << upper << std::endl;
    std::cout << "\033[1;31mC \033[0m\033[31mnormalising constant   \033[0m " << C << std::endl;
    std::cout << "\033[1;31mTHETA\033[0m \033[31minternal points at which the slopes change\033[0m\n";
    for (auto el : theta) std::cout << el << " ";
    std::cout << "\n\033[1;31mINTERCEPTS\033[0m \033[31mintercept values of piecewise polynomials\033[0m\n";
    for (auto el : intercept) std::cout << el << " ";
    std::cout << "\n\033[1;31mSLOPES\033[0m \033[31mslope values of piecewise polynomials\033[0m\n";
    for (auto el : slope) std::cout << el << " ";
    // obsolete printing code
    // std::cout << "\n\033[1;31malpha\033[0m\n" << alpha << std::endl;
    std::cout << "\n\033[1;31mlog likelihood\033[0m\n" << ll << std::endl;
    // std::cout << "\033[1;31mpi\033[0m\n" << pi << std::endl;
    // std::cout << "\n\033[1;31mdpk\033[0m\n" << dpk << std::endl;
    std::cout << "\033[1;31mcpk\033[0m\n" << cpk.row(0) << std::endl;
    // std::cout << "\033[1;31mfk\033[0m\n" << fk << std::endl;
    std::cout << "\n\033[1;91;47m                                                  \033[m\n\n";
  }

  // Calculate the approximate memory consumption of lcd. 
  void calculate_memory() const {
    size_t approx_memory = 0;
    size_t sizeof_T = sizeof(mpf);
    approx_memory += (sizeof_T                                // alpha
                    + sizeof_T                                // C
                    + sizeof_T                                // ll
                    + sizeof_T                                // lower
                    + sizeof_T                                // upper
                    + sizeof_T * theta.cols()                 // theta
                    + sizeof_T * knots.cols()                 // knots
                    + sizeof_T * pi.cols()                    // pi
                    + sizeof_T * intercept.cols()             // intercept
                    + sizeof_T * slope.cols()                 // slope
                    + sizeof_T * fk.cols()                    // fk
                    + sizeof_T * dpk.cols() * dpk.rows()      // dpk
                    + sizeof_T * cpk.cols() * cpk.cols());    // cpk

    std::cout << approx_memory << " bytes\n";
    return;
  }
};        // class LCD

}         // namespace lcde

#endif    // LCDE_LCD_H_