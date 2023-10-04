#ifndef LEVELDB_LEARNED_LINEAR_H_
#define LEVELDB_LEARNED_LINEAR_H_

#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <utility>
#include <vector>

namespace linear {

class LinearModel {
  public:
  double slope;
  double intercept;
  uint64_t data_size;
  uint64_t error;

  LinearModel() = default;
  LinearModel(double s, double i, uint64_t sz, uint64_t err)
   : slope(s), intercept(i), data_size(sz), error(err) {}

  // Fits a single linear segment over the sorted data
  void Learn(const std::vector<uint64_t> data) {
    uint64_t min_key = data.front();
    uint64_t max_key = data.back();
    slope = 1. / (max_key - min_key);
    intercept = -slope * min_key;

    data_size = static_cast<uint64_t>(data.size());
    uint64_t max_error = 0;
    for (uint64_t i = 0; i < data_size; ++i) {
      uint64_t prediction = static_cast<uint64_t>(data[i] * slope + intercept);
      uint64_t local_error = prediction > i ? prediction - i : i - prediction;
      if (local_error > max_error) max_error = local_error;
    }

    error = max_error;
  }

  std::pair<uint64_t, uint64_t> Get(uint64_t key) const {
    int64_t prediction = static_cast<int64_t>(slope * key + intercept);
    uint64_t lhs = static_cast<uint64_t>(std::max(prediction, 0L));
    lhs = lhs > error ? lhs - error : 0;
    uint64_t rhs = std::min(data_size, static_cast<uint64_t>(prediction) + error);
    return std::make_pair(lhs, rhs);
  }

  size_t Size() const {
    return sizeof(double) * 2 + sizeof(uint64_t) * 2;
  }

};     // class LinearModel

}

#endif // LEVELDB_LEARNED_LINEAR_H_