#ifndef LCDE_KNOT_OBJECT_H_
#define LCDE_KNOT_OBJECT_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <map>
#include <utility>
#include <iterator>
#include <cmath>
#include <Eigen/Eigenvalues>

#include "utils.h"
#include "lcd.h"

namespace lcde {

class Knot {
  public:
  mpf slope;
  mpf intercept;
  mpf theta;
  mpf addend;
  long error;

  Knot() = default;
  Knot(mpf s, mpf i, mpf t, mpf a, long e)
   : slope(s), intercept(i), theta(t), addend(a), error(e) {}

  void printKnotObject() const {
    std::cout << std::fixed << std::setprecision(30);
    std::cout << "\n\033[1;91;47m Printing contents of knotObject \033[m\n";
    std::cout << "\033[1;31mSLOPE   \033[0m" << slope << std::endl;
    std::cout << "\033[1;31mINTER   \033[0m" << intercept << std::endl;
    std::cout << "\033[1;31mTHETA   \033[0m" << theta << std::endl;
    std::cout << "\033[1;31mADDEND  \033[0m" << addend << std::endl;
    std::cout << "\033[1;31mERROR   \033[0m" << error << std::endl;
    std::cout << "\033[1;91;47m                                 \033[m\n\n";
  } 

  static size_t getSize() {
    return sizeof(mpf) * 4 + sizeof(long);
  }
};      // class Knot

template <typename KeyType>
class KnotObject {
  public:
  KnotObject() = default;
  KnotObject(mpf slope, mpf intercept, long data_size, long fanout)
   : linear_layer_slope(slope), linear_layer_intercept(intercept), data_size(data_size), fanout(fanout) {}

  mpf linear_layer_slope;
  mpf linear_layer_intercept;
  long data_size;
  long fanout;
  
  std::vector<std::vector<Knot>> knots;
  
  void setParameters(mpf _slope, mpf _intercept, long _data_size, long _fanout) {
    linear_layer_slope = _slope;
    linear_layer_intercept = _intercept;
    data_size = _data_size;
    fanout = _fanout;
  }

  void append(std::vector<Knot> k) {
    knots.push_back(k);
  }

  inline mpf exp1a(mpf x) const {
    x = 1.0 + x / 64.0;
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x; x *= x;
    return x;
  }

  // Initial implementation of find function
  // Currently the fastest implementation
  std::pair<size_t, size_t> findOne(const KeyType& key) const {
    long rank = static_cast<long>(std::fma(linear_layer_slope, key, linear_layer_intercept));
    rank = std::max(0L, std::min(fanout - 1, rank));

    const std::vector<Knot>& knot = knots[rank];

    // MODIFIED ////////////////////////////////////////////////////////////////
    // const mpf tmp_key = std::max(knot[0].theta, std::min(knot[knot.size() - 1].theta, static_cast<mpf>(key)));
    // auto iter = std::upper_bound(knot.begin(), knot.end(), tmp_key, 
    // [](const mpf& k, const Knot& ko) {
    //   return k < ko.theta;
    // });
    // iter = std::max(knot.begin(), std::min(iter - 1, knot.end() - 2));
    // const mpf& a = iter->addend;
    // const mpf& s = iter->slope;
    // const mpf& i = iter->intercept;

    // long search_result = data_size * (
    //   s == 0 ? (std::fma(i, tmp_key, a)) 
    //         //  : (s > 0 ? (a + exp1(std::fma(s, tmp_key, i))) 
    //         //           : (a - exp1(std::fma(s, tmp_key, i))))
    //         //  : (s > 0 ? (a + expApprox<mpf>(std::fma(s, tmp_key, i))) 
    //         //           : (a - expApprox<mpf>(std::fma(s, tmp_key, i))))
    //           : (s > 0 ? (a + std::exp(std::fma(s, tmp_key, i))) 
    //                    : (a - std::exp(std::fma(s, tmp_key, i))))
    // );

    // const uint& error = iter->error;
    ////////////////////////////////////////////////////////////////////////////

    // ORIGINAL ////////////////////////////////////////////////////////////////
    long search_result;
    long error;
    if (knot.size() == 1) {
      search_result = knot[0].addend;
      error = knot[0].error;
    } else {
      // const mpf tmp_key = std::max(knot[0].theta, std::min(knot[knot.size() - 1].theta, static_cast<mpf>(key)));
      const mpf tmp_key = static_cast<mpf>(key);
      auto iter = std::upper_bound(knot.begin(), knot.end(), tmp_key, 
      [](const mpf& k, const Knot& ko) {
        return k < ko.theta;
      });
      iter = std::min(iter - 1, knot.end() - 2);
      // iter = std::max(knot.begin(), std::min(iter - 1, knot.end() - 2));
      const mpf& a = iter->addend;
      const mpf& s = iter->slope;
      const mpf& i = iter->intercept;
      error = iter->error;

      // search_result = data_size * (
      //   s == 0 ? (std::fma(i, tmp_key, a)) 
      //         //  : (s > 0 ? (a + exp1(std::fma(s, tmp_key, i))) 
      //         //           : (a - exp1(std::fma(s, tmp_key, i))))
      //          : (s > 0 ? (a + expApprox<mpf>(std::fma(s, tmp_key, i))) 
      //                   : (a - expApprox<mpf>(std::fma(s, tmp_key, i))))
      //         //  : (s > 0 ? (a + std::exp(std::fma(s, tmp_key, i))) 
      //         //           : (a - std::exp(std::fma(s, tmp_key, i))))
      // );

      if (s > 0) {
        search_result = static_cast<long>(std::fma(data_size, std::exp(std::fma(s, tmp_key, i)), a));
        // search_result = static_cast<long>(std::fma(data_size, exp1(std::fma(s, tmp_key, i)), a));
        // search_result = static_cast<long>(std::fma(data_size, expApprox<mpf>(std::fma(s, tmp_key, i)), a));
      } else if (s < 0) {
        search_result = static_cast<long>(std::fma(data_size, -std::exp(std::fma(s, tmp_key, i)), a));
        // search_result = static_cast<long>(std::fma(data_size, -exp1(std::fma(s, tmp_key, i)), a));
        // search_result = static_cast<long>(std::fma(data_size, -expApprox<mpf>(std::fma(s, tmp_key, i)), a));
      } else {
        search_result = static_cast<long>(std::fma(i, key, a));
      }
      // return std::make_pair(static_cast<size_t>(std::max(search_result - error, 0L)),
      //                       static_cast<size_t>(std::min(data_size, search_result + error)));
    }
    ////////////////////////////////////////////////////////////////////////////
    

    return std::make_pair(static_cast<size_t>(std::max(search_result - error, 0L)),
                          static_cast<size_t>(std::min(data_size, search_result + error)));
  }

  size_t getSize() const {
    size_t model_size = sizeof(mpf) * 2 + sizeof(long);
    size_t knot_size = Knot::getSize();
    for (const auto& knot_vector : knots) {
      model_size += knot_size * knot_vector.size();
    }
    return model_size;
  }

  void printKnots() const {
    std::cout << "Size of global knot vector : ";
    std::cout << knots.size() << std::endl;
    std::cout << "Size of individual knot vectors : \n";
    for (auto knotv : knots) {
      for (auto k : knotv) {
        std::cout << k.error << " ";
        // k.printKnotObject();
      }
      std::cout << "\n";
    }
    return;
  }

  long getError(const KeyType& key) const {
    long rank = static_cast<long>(std::fma(linear_layer_slope, key, linear_layer_intercept));
    rank = std::max(0L, std::min(fanout - 1, rank));

    const std::vector<Knot>& knot = knots[rank];

    long error;
    if (knot.size() == 1) {
      error = knot[0].error;
    } else {
      const mpf tmp_key = static_cast<mpf>(key);
      auto iter = std::upper_bound(knot.begin(), knot.end(), tmp_key, 
      [](const mpf& k, const Knot& ko) {
        return k < ko.theta;
      });
      iter = std::min(iter - 1, knot.end() - 2);
      error = iter->error;
    }

    return error;
  }

};     // class KnotObject

}      // namespace lcde

#endif // LCDE_KNOT_OBJECT_H_