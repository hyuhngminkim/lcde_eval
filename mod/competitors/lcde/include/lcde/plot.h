#ifndef LCDE_TEST_UNIT_H_
#define LCDE_TEST_UNIT_H_

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

#include <chrono>

#include "utils.h"
#include "lcd.h"
#include "loglik.h"
#include "gradient.h"
#include "mm.h"
#include "pnnls.h"
#include "line.h"

namespace lcde {

// The purpose of class TestUnit is twofold:
// 1. Serve as a wrapper class for the builder such that we can plot the cdf of
//    the resulting LCDE object. 
// 2. Serve as a testing unit that we can test individual operations of LCDE.
//    Current LCDE implementation is built for the SOSD benchmark, and therefore
//    lacks the function of testing individual LCDE objects. 
// As such, class TestUnit has many similarities with its original class Builder.
template <typename KeyType>
class TestUnit {
  public:
  TestUnit() {
    testObject.setOuter(*this);
  }

  private:
  // The sampling ratio
  double sampling_rate;
  // Minimum number of elements for an estimation => Under test
  unsigned long min_size = 3;
  // The number of second layer lcds => Under test
  long fanout;
  // Slope and intercept for the linear regression of the first layer
  mpf slope;
  mpf intercept;
  // Auxiliary data structure for storing sampled data and its sample density
  struct point {
    mpf x;
    mpf y;

    point(mpf lhs, mpf rhs) : x(lhs), y(rhs) {}
  };
  // For storing sampled data and their sample density
  std::vector<point> sample_data;

  mpf current_base = 0;
  mpf current_ratio = 0;
  long data_size;

  // Parameters for calculating index
  // The original equation for calculating the CDF is
  /*
    (base + ratio * (cpk - (fk / slope))) + (ratio / (C * slope)) * e^{slope * x + intercept}
  = (base + ratio * (cpk - (fk / slope))) + e^{slope * x + intercept + ln(ratio / (C * slope))}
  */
  // which requires 6 parameters. 
  // We reduce the number of parameters into 3, by precomputing each term.
  // The resulting equation is
  /*
    addend + e^{slope * x + newIntercept}
  */
  struct knotObject {
    mpf slope;
    mpf intercept;
    mpf theta;
    mpf addend;
    long error;

    public:
    static size_t getSize() {
      return sizeof(mpf) * 4 + sizeof(long);
    }
  };

  // The block size of the permanent storage
  // Subject to change according to the input parameter
  int block_size = 4 * 1024;
  // Size of a single "knot" of a Knot
  const size_t knot_size = knotObject::getSize();
  // Size of level1 parameters
  const size_t level1_size = sizeof(mpf) * 2 + sizeof(long);

  // Struct cdfPoint is used to reduce the memory consumption.
  // However, memory consumption is not a matter in this testing environment.
  // We therefore define a new struct that stores all necessary elements for
  // plotting. 
  class {
    private:
    // Reference to the outer class TestUnit
    // Required for storing current_base and current_ratio
    TestUnit* tu = nullptr;

    public:
    // Constructor for testObject
    void setOuter(TestUnit& testUnitInstance) {
      tu = &testUnitInstance;
    }

    std::vector<std::vector<knotObject>> knots;

    void append(const LCD& lcd) {
      // C
      const mpf t_C = lcd.C;
      // slopes
      const auto& t_slope = lcd.slope;
      // intercepts
      const auto& t_intercept = lcd.intercept;
      // cpk
      VectorT<mpf> cpk_first_row = lcd.cpk.row(0);
      std::vector<mpf> t_cpk = convertToSTL(cpk_first_row);
      t_cpk.insert(t_cpk.begin(), 0.);
      // fk
      const auto& t_fk = lcd.fk;
      // Append theta
      std::vector<mpf> t_theta = concatenate3(lcd.lower, lcd.theta, lcd.upper);

      std::vector<mpf> t_addend;
      std::vector<mpf> t_newIntercept;
      const mpf& t_base = tu->current_base;
      const mpf& t_ratio = tu->current_ratio;

      std::vector<knotObject> t_ko;

      const size_t interval_size = t_slope.size();
      for (size_t i = 0; i < interval_size; ++i) {
        knotObject ko;
        mpf addend;
        mpf newIntercept;
        if (t_slope[i] != 0) {
          addend = (t_base + t_ratio * (t_cpk[i] - (t_fk[i] / t_slope[i])));
          mpf t_ln;
          if (t_slope[i] < 0)  {
            t_ln = std::log(t_ratio / (t_C * -t_slope[i]));
          } else {
            t_ln = std::log(t_ratio / (t_C * t_slope[i]));
          }
          newIntercept = t_intercept[i] + t_ln;
        } else {
          // When slope = 0, we cannot divide the integral by slope. 
          // We calclate the cdf by addend + intercept * x for two reasons:
          // 1. we must keep the slope variable to check if it is zero
          // 2. we do not use the intercept variable
          // We therefore store the value multiplied to x in intercept not slope
          addend = t_base + t_ratio * (t_cpk[i] - (std::exp(t_intercept[i]) * t_theta[i]) / t_C);
          newIntercept = (t_ratio * std::exp(t_intercept[i])) / t_C;
        } 
        ko.slope = t_slope[i];
        ko.intercept = newIntercept;
        ko.theta = t_theta[i];
        ko.addend = addend;
        t_ko.push_back(ko);
      }

      // Last element for theta[-1]
      knotObject ko = t_ko[interval_size - 1];
      ko.theta = t_theta[interval_size];
      t_ko.push_back(ko);

      knots.push_back(t_ko);
    }

    void append() {
      std::vector<knotObject> t_ko;
      knotObject ko;
      // Set slope = 0, intercept = 0, addend = base. 
      // theta has no meaning in calculation, and we set 0 as dummy data.
      ko.slope = 0;
      ko.intercept = 0;
      ko.theta = 0;
      ko.addend = tu->current_base;

      t_ko.push_back(ko);
      knots.push_back(t_ko);
    }

    void getSize() const {
      size_t model_size = tu->level1_size;
      size_t knot_size = tu->knot_size;
      for (const auto& knot_vector : knots) {
        model_size += knot_size * knot_vector.size();
      }
      std::cout << "\nSize of index : " << model_size << std::endl;
    }

    void printJSON(const std::string& dataset) const {
      std::ofstream file;
      std::string filename = "parameters/" + dataset + "_plot_parameters.json";
      file.open(filename);

      file << std::fixed << std::setprecision(30);
      
      file << "{\"parameters\":[";

      for (long i = 0; i < tu->fanout; ++i) {
        size_t knot_size = knots[i].size();
        if (i > 0) file << ",";
        file << "{";
        // print slopes
        file << "\"slope\":[";
        for (size_t j = 0; j < knot_size - 1; ++j) {
          if (j > 0) file << ",";
          file << knots[i][j].slope;
        }
        file << "],";
        // print thetas
        file << "\"theta\":[";
        for (size_t j = 0; j < knot_size; ++j) {
          if (j > 0) file << ",";
          file << knots[i][j].theta;
        }
        file << "],";
        // print addends
        file << "\"addend\":[";
        for (size_t j = 0; j < knot_size - 1; ++j) {
          if (j > 0) file << ",";
          file << knots[i][j].addend;
        }
        file << "],";
        // print intercepts
        file << "\"intercept\":[";
        for (size_t j = 0; j < knot_size - 1; ++j) {
          if (j > 0) file << ",";
          file << knots[i][j].intercept;
        }
        file << "]";
        file << "}";
      }

      file << "]}";

      file.close();

      return;
    }

    private:
    inline std::vector<mpf> convertToSTL(const VectorT<mpf>& v) {
      std::vector<mpf> res(v.data(), v.data() + v.size());
      return res;
    }

    // Auxiliary function to concatenate lower, theta, and upper
    inline std::vector<mpf> concatenate3(mpf lower, const VectorT<mpf> t, mpf upper) {
      VectorT<mpf> temp = concatenate(concatenate(lower, t), upper);
      return convertToSTL(temp);
    }

  } testObject;

  public:
  // Build LCDE into a single object. 
  size_t buildSingle(const std::vector<point>& data, size_t budget) {
    int max_iter = 100;
    mpf tol = 1e-6;
    auto data_size = data.size();
    VectorT<mpf> x(data_size);
    VectorT<int> w(data_size);

    weight(data, x, w);
    double lambda = 1e-15;
    int n = data_size;
    int nx = x.cols();

    if (nx <= 2) {
      testObject.append();
      // MOD
      return budget - 1;
    }

    mpf lower = x(0);
    mpf upper = x(nx - 1);
    VectorT<mpf> theta;
    VectorT<mpf> pi;

    LCD lcd_ = LCD(0, lower, upper, theta, pi); 
    if (lower < lcd_.lower) lcd_.lower = lower;
    if (upper > lcd_.upper) lcd_.upper = upper;

    VectorT<mpf> xx(x.cols());
    xsq(xx, x, w);
    loglik(lcd_, x, w);
    mpf ll_old = std::numeric_limits<mpf>::lowest();

    // main loop
    for (int i = 0; i < max_iter; ++i) {
      // MOD
      if (lcd_.ll <= ll_old + tol || lcd_.knot_count() >= budget) break;
      // if (lcd_.ll <= ll_old + tol) break;
      ll_old = lcd_.ll;
      VectorT<mpf> g_theta = maxima_gradient(lcd_, x, w, xx);

      if (g_theta.cols() >= 1) {
        VectorT<mpf> g_pi = VectorT<mpf>::Zero(g_theta.cols());
        int oi = 0;     // index of theta of original lcd
        int ni = 0;     // index of new theta
        while (oi < lcd_.theta.cols()) {
          if (lcd_.theta(oi) == g_theta(ni)) g_pi(ni) = lcd_.pi(oi++);
          ++ni;
        }
        lcd_ = LCD(lcd_.alpha, lcd_.lower, lcd_.upper, g_theta, g_pi);
      }

      VectorT<mpf>& knots = lcd_.knots;
      MatrixT<mpf>& cpk = lcd_.cpk;

      int nk = knots.cols();

      MatrixT<mpf> cpkr = getElementFromIndex(cpk, {nk - 1}).replicate(1, nk)
                          - cbind(0.0, dropElementByIndex(cpk, nk - 1));

      // E{(X - theta)_+}
      VectorT<mpf> mu = cpkr.row(1).array() - knots.array() * cpkr.row(0).array();
      // index of knots against x
      std::vector<int> knots_idx = getIndexOnly(knots, x);
      // gradient vector
      VectorT<mpf> grad = (mu * n) + getElementFromIndex(xx, knots_idx);
      // E{(X - theta_j)_+ (X - theta_k)_+}
      DynamicMatrixT<mpf> mm = computeDiff(knots, cpkr, nk);
      // negative Hessian matrix
      mm.triangularView<Eigen::Upper>() = mm.transpose();
      DynamicMatrixT<mpf> H = mm - tcrossprod(mu);

      // Catalogue: 
      // https://eigen.tuxfamily.org/dox/group__TopicLinearAlgebraDecompositions.html
      // https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html
      Eigen::SelfAdjointEigenSolver<DynamicMatrixT<mpf>> es(H);

      std::vector<mpf> v2;
      std::vector<int> R_idx;
      v2.reserve(nk);
      R_idx.reserve(nk);
      // Note that the eigen() function of R sorts the eigenvalues in decreasing
      // order, but SelfAdjointEigenSolver of Eigen sorts the eigenvalues in
      // increasing order.
      mpf lower_bound = es.eigenvalues()(nk - 1) * lambda;
      for (auto i = 0; i < nk; ++i) {
        auto tmp_ev = std::abs(es.eigenvalues()(i));
        if (tmp_ev != 0 && tmp_ev >= lower_bound) {
          R_idx.push_back(i);
          v2.push_back(sqrt(tmp_ev));
        }
      }

      size_t kr = v2.size();
      DynamicMatrixT<mpf> first_kr_transposed = 
                es.eigenvectors()(Eigen::all, Eigen::seqN(0, kr))
                  .transpose();
      DynamicMatrixT<mpf> R = columnWiseMult(first_kr_transposed, v2);
      VectorT<mpf> p = grad / n + concatenate(-lcd_.alpha, lcd_.pi) * H;
      VectorT<mpf> b = 
      auxiliaryDivision((first_kr_transposed * p.transpose()).transpose(), v2);

      VectorT<double> nnls = pnnls(R.cast<double>(), b.cast<double>(), 1);
      nnls(0) = -nnls(0);

      if (vectorIsInvalid(nnls)) break;

      // perform line search
      // MOD
      line_lcd(lcd_, x, w, xx, nnls, ll_old, budget);
      // remove zero slope changes if exists
      if ((lcd_.pi.array() == 0).any()) lcd_.simplify();
    }

    testObject.append(lcd_);

    // MOD
    return (budget >= lcd_.knot_count() ? budget - lcd_.knot_count() : 0);
  }
  // void buildSingle(const std::vector<point>& data) {
  // size_t buildSingle(const std::vector<point>& data, size_t budget) {
  //   int max_iter = 100;
  //   mpf tol = 1e-6;      // experimental changes
  //   auto data_size = data.size();
  //   VectorT<mpf> x(data_size);
  //   VectorT<int> w(data_size);

  //   weight(data, x, w);
  //   double lambda = 1e-15;
  //   int n = data_size;
  //   int nx = x.cols();

  //   if (nx <= 2) {
  //     testObject.append();
  //     return budget - 1;
  //   }

  //   mpf lower = x(0);
  //   mpf upper = x(nx - 1);
  //   VectorT<mpf> theta;
  //   VectorT<mpf> pi;

  //   LCD lcd_ = LCD(0, lower, upper, theta, pi); 
  //   if (lower < lcd_.lower) lcd_.lower = lower;
  //   if (upper > lcd_.upper) lcd_.upper = upper;

  //   VectorT<mpf> xx(x.cols());
  //   xsq(xx, x, w);
  //   loglik(lcd_, x, w);
  //   mpf ll_old = std::numeric_limits<mpf>::lowest();

  //   // main loop
  //   for (int i = 0; i < max_iter; ++i) {
  //     if (lcd_.ll <= ll_old + tol || lcd_.knot_count() >= budget) break;
  //     // if (lcd_.ll <= ll_old + tol) break;
  //     ll_old = lcd_.ll;
  //     VectorT<mpf> g_theta = maxima_gradient(lcd_, x, w, xx);

  //     if (g_theta.cols() >= 1) {
  //       VectorT<mpf> g_pi = VectorT<mpf>::Zero(g_theta.cols());
  //       int oi = 0;     // index of theta of original lcd
  //       int ni = 0;     // index of new theta
  //       while (oi < lcd_.theta.cols()) {
  //         if (lcd_.theta(oi) == g_theta(ni)) g_pi(ni) = lcd_.pi(oi++);
  //         ++ni;
  //       }
  //       lcd_ = LCD(lcd_.alpha, lcd_.lower, lcd_.upper, g_theta, g_pi);
  //     }

  //     VectorT<mpf>& knots = lcd_.knots;
  //     MatrixT<mpf>& cpk = lcd_.cpk;

  //     int nk = knots.cols();

  //     MatrixT<mpf> cpkr = getElementFromIndex(cpk, {nk - 1}).replicate(1, nk)
  //                         - cbind(0.0, dropElementByIndex(cpk, nk - 1));

  //     // E{(X - theta)_+}
  //     VectorT<mpf> mu = cpkr.row(1).array() - knots.array() * cpkr.row(0).array();
  //     // index of knots against x
  //     std::vector<int> knots_idx = getIndexOnly(knots, x);
  //     // gradient vector
  //     VectorT<mpf> grad = (mu * n) + getElementFromIndex(xx, knots_idx);
  //     // E{(X - theta_j)_+ (X - theta_k)_+}
  //     DynamicMatrixT<mpf> mm = computeDiff(knots, cpkr, nk);
  //     // negative Hessian matrix
  //     mm.triangularView<Eigen::Upper>() = mm.transpose();
  //     DynamicMatrixT<mpf> H = mm - tcrossprod(mu);

  //     // Catalogue: 
  //     // https://eigen.tuxfamily.org/dox/group__TopicLinearAlgebraDecompositions.html
  //     // https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html
  //     Eigen::SelfAdjointEigenSolver<DynamicMatrixT<mpf>> es(H);

  //     std::vector<mpf> v2;
  //     std::vector<int> R_idx;
  //     v2.reserve(nk);
  //     R_idx.reserve(nk);
  //     // Note that the eigen() function of R sorts the eigenvalues in decreasing
  //     // order, but SelfAdjointEigenSolver of Eigen sorts the eigenvalues in
  //     // increasing order.
  //     mpf lower_bound = es.eigenvalues()(nk - 1) * lambda;
  //     for (auto i = 0; i < nk; ++i) {
  //       auto tmp_ev = std::abs(es.eigenvalues()(i));
  //       if (tmp_ev != 0 && tmp_ev >= lower_bound) {
  //         R_idx.push_back(i);
  //         v2.push_back(sqrt(tmp_ev));
  //       }
  //     }

  //     size_t kr = v2.size();
  //     DynamicMatrixT<mpf> first_kr_transposed = 
  //               es.eigenvectors()(Eigen::all, Eigen::seqN(0, kr))
  //                 .transpose();
  //     DynamicMatrixT<mpf> R = columnWiseMult(first_kr_transposed, v2);
  //     VectorT<mpf> p = grad / n + concatenate(-lcd_.alpha, lcd_.pi) * H;
  //     VectorT<mpf> b = 
  //     auxiliaryDivision((first_kr_transposed * p.transpose()).transpose(), v2);

  //     VectorT<double> nnls = pnnls(R.cast<double>(), b.cast<double>(), 1);
  //     nnls(0) = -nnls(0);

  //     if (vectorIsInvalid(nnls)) break;

  //     // perform line search
  //     line_lcd(lcd_, x, w, xx, nnls, ll_old, budget);
  //     // remove zero slope changes if exists
  //     if ((lcd_.pi.array() == 0).any()) lcd_.simplify();
  //   }

  //   testObject.append(lcd_);
  //   // return;
  //   return (budget >= lcd_.knot_count() ? budget - lcd_.knot_count() : 0);
  // }

  // Build LCDE using the standard linear model - fanout structure. 
  // Parameters are brought in as a vector for scalability and testing purposes.
  void build(const std::vector<KeyType>& data, std::vector<double> params) {
    // Set parameters
    // params[0]: sampling rate
    sampling_rate = params[0];
    // params[1]: fanout
    fanout = static_cast<long>(params[1]);

    // MOD: strict budget cap
    size_t size_in_b = (size_t)(params[2] * 1024);
    if (size_in_b <= level1_size) {
      std::cerr << "Index size insufficiently small. Aborting...\n";
      exit(EXIT_FAILURE);
    }
    const size_t average = 3;
    size_t budget = 0;
    fanout = (size_t)((size_in_b - level1_size)/ (average * knot_size));
    std::cout << "Fanout : " << fanout << std::endl;

    data_size = data.size();

    const long input_size = data.size();
    const long sample_size = std::min<long>(
      input_size, std::max<long>(sampling_rate * input_size, min_size)
    );

    slope = 1. / static_cast<mpf>(data.back() - data.front());
    intercept = -slope * static_cast<mpf>(data.front());

    slope *= fanout;
    intercept *= fanout;

    sample_data.reserve(sample_size);

    long offset = static_cast<long>(1. * input_size / sample_size);
    std::vector<std::vector<point>> training_data(fanout);

    for (long i = 0; i < input_size; i += offset) {
      mpf cast_data = static_cast<mpf>(data[i]);
      long rank = static_cast<long>(slope * cast_data + intercept);
      rank = std::max(0L, std::min(fanout - 1, rank));
      training_data[rank].push_back({cast_data,
                                     static_cast<mpf>(1. * (i) / input_size)});
    }

    // Build first layer 
    point min = training_data.front().front();
    point max = training_data.back().back();

    // Train each subdata 
    for (long model_idx = 0; model_idx < fanout; ++model_idx) {
      budget += average;
      std::vector<point>& current_training_data = training_data[model_idx];
      size_t current_training_data_size = current_training_data.size();
      // The case for when the current_training_data.size() < 0 is a 
      // safety measure for when a slot of training_data has no element.
      // Bring the largest element of the previous training_data that has at least
      // a single element. If the first training_data slot has no element, insert
      // 0 as its element.
      if (model_idx == 0) {                      // First model
        if (current_training_data_size < min_size) {
          current_training_data.push_back(point(0, 0));
          testObject.append();
          --budget;
          if (current_training_data_size != 0) {
            current_base += current_training_data[current_training_data_size - 1].y;
          }
        } else {
          max = current_training_data.back();

          current_ratio = max.y;
          budget = buildSingle(current_training_data, budget);
          current_base += current_ratio;
        }
      } else if (model_idx == fanout - 1) {      // Last model
        if (current_training_data_size < min_size) {
          testObject.append();
          --budget;
        } else {
          min = training_data[model_idx - 1].back();

          current_ratio = 1 - min.y;
          budget = buildSingle(current_training_data, budget);
        }
      } else {                                    // Intermediate models
        if (current_training_data_size == 0) {
          current_training_data.push_back(training_data[model_idx - 1].back());
          testObject.append();
          --budget;
        } else {
          min = training_data[model_idx - 1].back();
          max = current_training_data.back();
          current_ratio = max.y - min.y;
          if (current_training_data_size < min_size) {
            testObject.append();
            --budget;
          } else {
            budget = buildSingle(current_training_data, budget);
          }
          current_base += current_ratio;
        }
      }
    }
  }

  void getSize() const {
    testObject.getSize();
  }

  void printJSON(const std::string& dataset) const {
    testObject.printJSON(dataset);
  }
};

}

#endif // LCDE_TEST_UNIT_H_