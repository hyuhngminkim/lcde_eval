//
// Created by daiyi on 2020/02/02.
//

#include "learned_index.h"

#include "db/version_set.h"
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <utility>

#include "util/mutexlock.h"

#include "util.h"

namespace adgMod {

std::pair<uint64_t, uint64_t> LearnedIndexData::GetPosition(
    const Slice& target_x) const {
  // assert(string_segments.size() > 1);
  ++served;

  // check if the key is within the model bounds
  uint64_t target_int = SliceToInteger(target_x);
  if (target_int > max_key) return std::make_pair(size, size);
  if (target_int < min_key) return std::make_pair(0, 0);

#ifdef LCDE
  return ko.findOne(target_int);
#endif

#ifdef BOURBON
  // binary search between segments
  uint32_t left = 0, right = (uint32_t)string_segments.size() - 1;
  while (left != right - 1) {
    uint32_t mid = (right + left) / 2;
    if (target_int < string_segments[mid].x)
      right = mid;
    else
      left = mid;
  }

  // calculate the interval according to the selected segment
  double result =
      target_int * string_segments[left].k + string_segments[left].b;
  result = is_level ? result / 2 : result;
  uint64_t lower =
      result - error > 0 ? (uint64_t)std::floor(result - error) : 0;
  uint64_t upper = (uint64_t)std::ceil(result + error);
  if (lower >= size) return std::make_pair(size, size);
  upper = upper < size ? upper : size - 1;
  return std::make_pair(lower, upper);
#endif

#ifdef RS
  rs::SearchBound sb = rs_.GetSearchBound(target_int);
  return std::make_pair(sb.begin, sb.end);
#endif

#ifdef PGM
  auto approx_range = pgm_.find_approximate_position(target_int);
  return std::make_pair(approx_range.lo, approx_range.hi + 1);
#endif

#ifdef CHT
  cht::SearchBound sb = cht_.GetSearchBound(target_int);
  return std::make_pair(sb.begin, sb.end);
#endif

#ifdef LINEAR
  return lm_.Get(target_int);
#endif

  // default return: search entire key set
  return std::make_pair(0, size);
}

uint64_t LearnedIndexData::MaxPosition() const { return size - 1; }

double LearnedIndexData::GetError() const { return error; }

// Actual function doing learning
bool LearnedIndexData::Learn() {

  // check if data if filled
  if (string_keys.empty()) assert(false);

  // fill in some bounds for the model
  uint64_t temp = string_keys.back();
  min_key = string_keys.front();
  max_key = string_keys.back();
  size = string_keys.size();

  // actual training
#ifdef LCDE
  lcde::BuilderObject<uint64_t> bo = lcde::BuilderObject<uint64_t>(ko);
  std::vector<double> params = {1, 100};
  bo.build(string_keys, params);
  learned.store(true);

  return true;
#endif

#ifdef BOURBON
  PLR plr = PLR(error);
  // PLR plr = PLR(8);
  std::vector<Segment> segs = plr.train(string_keys, !is_level);
  if (segs.empty()) return false;
  // fill in a dummy last segment (used in segment binary search)
  segs.push_back((Segment){temp, 0, 0});
  string_segments = std::move(segs);

  learned.store(true);
  return true;
#endif

#ifdef RS
  rs::Builder<uint64_t> rsb(min_key, max_key, 
                            rs_configs[rs_variant].first, rs_configs[rs_variant].second);
  for (const auto& key : string_keys) rsb.AddKey(key);
  rs_ = rsb.Finalize();

  learned.store(true);
  return true;
#endif

#ifdef PGM
  pgm_ = decltype(pgm_)(string_keys.begin(), string_keys.end());

  learned.store(true);
  return true;
#endif

#ifdef CHT
  size_t num_bins_ = cht_configs[cht_variant].first;
  size_t max_error_ = cht_configs[cht_variant].second;
  cht::Builder<uint64_t> chtb(min_key, max_key, num_bins_, max_error_, false, false);
  for (const auto& key : string_keys) chtb.AddKey(key);
  cht_ = chtb.Finalize();

  learned.store(true);
  return true;
#endif

#ifdef LINEAR
  lm_.Learn(string_keys);

  learned.store(true);
  return true;
#endif

  return false;
}

// static learning function to be used with LevelDB background scheduling
// file learning
uint64_t LearnedIndexData::FileLearn(void* arg) {
  Stats* instance = Stats::GetInstance();
  bool entered = false;
  instance->StartTimer(11);

  MetaAndSelf* mas = reinterpret_cast<MetaAndSelf*>(arg);
  LearnedIndexData* self = mas->self;
  self->level = mas->level;

  Version* c = db->GetCurrentVersion();
  if (self->FillData(c, mas->meta)) {
    self->Learn();
    entered = true;
  } else {
    self->learning.store(false);
  }
  adgMod::db->ReturnCurrentVersion(c);

  auto time = instance->PauseTimer(11, true);

  if (entered) {
    // count how many file learning are done.
    self->cost = time.second - time.first;
    learn_counter_mutex.Lock();
    events[1].push_back(new LearnEvent(time, 1, self->level, true));
    levelled_counters[11].Increment(mas->level, time.second - time.first);
    learn_counter_mutex.Unlock();
  }

  if (!fresh_write) delete mas->meta;
  delete mas;
  return entered ? time.second - time.first : 0;
}

// general model checker
bool LearnedIndexData::Learned() {
  if (learned_not_atomic)
    return true;
  else if (learned.load()) {
    learned_not_atomic = true;
    return true;
  } else
    return false;
}

// level model checker, used to be also learning trigger
bool LearnedIndexData::Learned(Version* version, int v_count, int level) {
  if (learned_not_atomic)
    return true;
  else if (learned.load()) {
    learned_not_atomic = true;
    return true;
  }
  return false;
}

// file model checker, used to be also learning trigger
bool LearnedIndexData::Learned(Version* version, int v_count,
                               FileMetaData* meta, int level) {
  if (learned_not_atomic)
    return true;
  else if (learned.load()) {
    learned_not_atomic = true;
    return true;
  } else
    return false;
}

bool LearnedIndexData::FillData(Version* version, FileMetaData* meta) {
  // if (filled) return true;

  if (version->FillData(adgMod::read_options, meta, this)) {
    // filled = true;
    return true;
  }
  return false;
}

void LearnedIndexData::WriteModel(const string& filename) {
  if (!learned.load()) return;

  std::ofstream output_file(filename);
  output_file.precision(15);
  output_file << adgMod::block_num_entries << " " << adgMod::block_size << " "
              << adgMod::entry_size << "\n";
  // for (Segment& item : string_segments) {
  //   output_file << item.x << " " << item.k << " " << item.b << "\n";
  // }
  output_file << "StartAcc"
              << " " << min_key << " " << max_key << " " << size << " " << level
              << " " << cost << "\n";
  // for (auto& pair : num_entries_accumulated.array) {
  //   output_file << pair.first << " " << pair.second << "\n";
  // }
}

void LearnedIndexData::ReadModel(const string& filename) {
  std::ifstream input_file(filename);

  if (!input_file.good()) return;
  input_file >> adgMod::block_num_entries >> adgMod::block_size >>
      adgMod::entry_size;
  while (true) {
    string x;
    double k, b;
    input_file >> x;
    if (x == "StartAcc") break;
    input_file >> k >> b;
    // string_segments.emplace_back(atoll(x.c_str()), k, b);
  }
  input_file >> min_key >> max_key >> size >> level >> cost;
  while (true) {
    uint64_t first;
    string second;
    if (!(input_file >> first >> second)) break;
    // num_entries_accumulated.Add(first, std::move(second));
  }

  learned.store(true);
}

// MOD: Write model into binary
void LearnedIndexData::WriteModelBinary(const string& filename) {
  if (!learned.load()) return;

  std::ofstream output_file(filename, std::ios::binary);
  if (!output_file.is_open()) {
    std::cerr << "Unable to open file " << filename << std::endl;
    return;
  }

  // std::cout << "Size of data being learned : " << string_keys.size() << std::endl;

  output_file.write(reinterpret_cast<char*>(&adgMod::block_num_entries), sizeof(uint64_t));
  output_file.write(reinterpret_cast<char*>(&adgMod::block_size), sizeof(uint64_t));
  output_file.write(reinterpret_cast<char*>(&adgMod::entry_size), sizeof(uint64_t));
  output_file.write(reinterpret_cast<char*>(&min_key), sizeof(uint64_t));
  output_file.write(reinterpret_cast<char*>(&max_key), sizeof(uint64_t));
  output_file.write(reinterpret_cast<char*>(&size), sizeof(uint64_t));
  output_file.write(reinterpret_cast<char*>(&level), sizeof(int));
  output_file.write(reinterpret_cast<char*>(&cost), sizeof(uint64_t));

#ifdef LCDE
  output_file.write(reinterpret_cast<char*>(&(ko.linear_layer_slope)), sizeof(double));
  output_file.write(reinterpret_cast<char*>(&(ko.linear_layer_intercept)), sizeof(double));
  output_file.write(reinterpret_cast<char*>(&(ko.data_size)), sizeof(long));
  output_file.write(reinterpret_cast<char*>(&(ko.fanout)), sizeof(long));
  size_t knots_size = ko.knots.size();
  output_file.write(reinterpret_cast<char*>(&knots_size), sizeof(size_t));
  for (size_t i = 0; i < knots_size; ++i) {
    size_t internal_knots_size = ko.knots[i].size();
    output_file.write(reinterpret_cast<char*>(&internal_knots_size), sizeof(size_t));
    for (size_t j = 0; j < internal_knots_size; ++j) {
      output_file.write(reinterpret_cast<char*>(&(ko.knots[i][j].slope)), sizeof(double));
      output_file.write(reinterpret_cast<char*>(&(ko.knots[i][j].intercept)), sizeof(double));
      output_file.write(reinterpret_cast<char*>(&(ko.knots[i][j].theta)), sizeof(double));
      output_file.write(reinterpret_cast<char*>(&(ko.knots[i][j].addend)), sizeof(double));
      output_file.write(reinterpret_cast<char*>(&(ko.knots[i][j].error)), sizeof(long));
    }
  }
#endif

#ifdef BOURBON
  size_t b_segments_size = string_segments.size();
  output_file.write(reinterpret_cast<char*>(&b_segments_size), sizeof(size_t));
  for (Segment& item : string_segments) {
    output_file.write(reinterpret_cast<char*>(&item.x), sizeof(uint64_t));
    output_file.write(reinterpret_cast<char*>(&item.k), sizeof(double));
    output_file.write(reinterpret_cast<char*>(&item.b), sizeof(double));
  }
#endif

#ifdef RS
  output_file.write(reinterpret_cast<char*>(&rs_.num_radix_bits_), sizeof(size_t));
  output_file.write(reinterpret_cast<char*>(&rs_.num_shift_bits_), sizeof(size_t));
  output_file.write(reinterpret_cast<char*>(&rs_.max_error_), sizeof(size_t));

  size_t radix_table_size = rs_.radix_table_.size();
  output_file.write(reinterpret_cast<char*>(&radix_table_size), sizeof(size_t));
  for (size_t i = 0; i < radix_table_size; ++i) {
    output_file.write(reinterpret_cast<char*>(&(rs_.radix_table_[i])), sizeof(uint32_t));
  }

  size_t splie_points_size = rs_.spline_points_.size();
  output_file.write(reinterpret_cast<char*>(&splie_points_size), sizeof(size_t));
  for (size_t i = 0; i < splie_points_size; ++i) {
    output_file.write(reinterpret_cast<char*>(&(rs_.spline_points_[i].x)), sizeof(uint64_t));
    output_file.write(reinterpret_cast<char*>(&(rs_.spline_points_[i].y)), sizeof(double));
  }
#endif

#ifdef PGM
  output_file.write(reinterpret_cast<char*>(&(pgm_.n)), sizeof(size_t));
  output_file.write(reinterpret_cast<char*>(&(pgm_.first_key)), sizeof(uint64_t));

  size_t segments_size = pgm_.segments.size();
  size_t levels_sizes_size = pgm_.levels_sizes.size();
  size_t levels_offsets_size = pgm_.levels_offsets.size();
  output_file.write(reinterpret_cast<char*>(&segments_size), sizeof(size_t));
  output_file.write(reinterpret_cast<char*>(&levels_sizes_size), sizeof(size_t));
  output_file.write(reinterpret_cast<char*>(&levels_offsets_size), sizeof(size_t));

  for (size_t i = 0; i < segments_size; ++i) {
    output_file.write(reinterpret_cast<char*>(&(pgm_.segments[i].key)), sizeof(uint64_t));
    output_file.write(reinterpret_cast<char*>(&(pgm_.segments[i].slope)), sizeof(double));
    output_file.write(reinterpret_cast<char*>(&(pgm_.segments[i].intercept)), sizeof(int32_t));
  }

  for (size_t i = 0; i < levels_sizes_size; ++i) {
    output_file.write(reinterpret_cast<char*>(&(pgm_.levels_sizes[i])), sizeof(size_t));
  }

  for (size_t i = 0; i < levels_offsets_size; ++i) {
    output_file.write(reinterpret_cast<char*>(&(pgm_.levels_offsets[i])), sizeof(size_t));
  }

#endif

#ifdef CHT
  output_file.write(reinterpret_cast<char*>(&(cht_.min_key_)), sizeof(uint64_t));
  output_file.write(reinterpret_cast<char*>(&(cht_.max_key_)), sizeof(uint64_t));
  output_file.write(reinterpret_cast<char*>(&(cht_.num_keys_)), sizeof(size_t));
  output_file.write(reinterpret_cast<char*>(&(cht_.num_bins_)), sizeof(size_t));
  output_file.write(reinterpret_cast<char*>(&(cht_.log_num_bins_)), sizeof(size_t));
  output_file.write(reinterpret_cast<char*>(&(cht_.max_error_)), sizeof(size_t));
  output_file.write(reinterpret_cast<char*>(&(cht_.shift_)), sizeof(size_t));

  size_t c_table_size = cht_.table_.size();
  output_file.write(reinterpret_cast<char*>(&c_table_size), sizeof(size_t));

  for (size_t i = 0; i < c_table_size; ++i) {
    output_file.write(reinterpret_cast<char*>(&(cht_.table_[i])), sizeof(unsigned));
  }
#endif

#ifdef LINEAR
  output_file.write(reinterpret_cast<char*>(&lm_.slope), sizeof(double));
  output_file.write(reinterpret_cast<char*>(&lm_.intercept), sizeof(double));
  output_file.write(reinterpret_cast<char*>(&lm_.error), sizeof(uint64_t));
#endif

  output_file.close();
}

// MOD: Read model from binary
void LearnedIndexData::ReadModelBinary(const string& filename) {
  std::ifstream input_file(filename, std::ios::binary);
  if (!input_file.good()) return;

  input_file.read(reinterpret_cast<char*>(&adgMod::block_num_entries), sizeof(uint64_t));
  input_file.read(reinterpret_cast<char*>(&adgMod::block_size), sizeof(uint64_t));
  input_file.read(reinterpret_cast<char*>(&adgMod::entry_size), sizeof(uint64_t));
  input_file.read(reinterpret_cast<char*>(&min_key), sizeof(uint64_t));
  input_file.read(reinterpret_cast<char*>(&max_key), sizeof(uint64_t));
  input_file.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
  input_file.read(reinterpret_cast<char*>(&level), sizeof(int));
  input_file.read(reinterpret_cast<char*>(&cost), sizeof(uint64_t));

#ifdef LCDE
  double lls, lli;
  long ds, fo;
  input_file.read(reinterpret_cast<char*>(&lls), sizeof(double));
  input_file.read(reinterpret_cast<char*>(&lli), sizeof(double));
  input_file.read(reinterpret_cast<char*>(&ds), sizeof(long));
  input_file.read(reinterpret_cast<char*>(&fo), sizeof(long));
  ko = lcde::KnotObject<uint64_t>(lls, lli, ds, fo);

  size_t knots_size, internal_knots_size;
  double ts, ti, tt, ta;
  long te;
  input_file.read(reinterpret_cast<char*>(&knots_size), sizeof(size_t));
  for (size_t i = 0; i < knots_size; ++i) {
    input_file.read(reinterpret_cast<char*>(&internal_knots_size), sizeof(size_t));
    std::vector<lcde::Knot> knot;
    for (size_t j = 0; j < internal_knots_size; ++j) {
      input_file.read(reinterpret_cast<char*>(&ts), sizeof(double));
      input_file.read(reinterpret_cast<char*>(&ti), sizeof(double));
      input_file.read(reinterpret_cast<char*>(&tt), sizeof(double));
      input_file.read(reinterpret_cast<char*>(&ta), sizeof(double));
      input_file.read(reinterpret_cast<char*>(&te), sizeof(long));
      knot.emplace_back(ts, ti, tt, ta, te);
    }
    ko.knots.push_back(knot);
  }
#endif

#ifdef BOURBON
  size_t b_segments_size;
  uint64_t tmp_x;
  double tmp_k;
  double tmp_b;

  // string_segments.clear();

  input_file.read(reinterpret_cast<char*>(&b_segments_size), sizeof(size_t));
  string_segments.reserve(b_segments_size);

  for (uint64_t i = 0; i < b_segments_size; ++i) {
    input_file.read(reinterpret_cast<char*>(&tmp_x), sizeof(uint64_t));
    input_file.read(reinterpret_cast<char*>(&tmp_k), sizeof(uint64_t));
    input_file.read(reinterpret_cast<char*>(&tmp_b), sizeof(uint64_t));
    Segment seg = Segment(tmp_x, tmp_k, tmp_b);
    string_segments.push_back(seg);
  }
#endif

#ifdef RS
  size_t num_radix_bits_, num_shift_bits_, max_error_;
  size_t radix_table_size, spline_points_size;
  std::vector<uint32_t> radix_table;
  std::vector<rs::Coord<uint64_t>> spline_points;

  input_file.read(reinterpret_cast<char*>(&num_radix_bits_), sizeof(size_t));
  input_file.read(reinterpret_cast<char*>(&num_shift_bits_), sizeof(size_t));
  input_file.read(reinterpret_cast<char*>(&max_error_), sizeof(size_t));
  
  input_file.read(reinterpret_cast<char*>(&radix_table_size), sizeof(size_t));
  for (size_t i = 0; i < radix_table_size; ++i) {
    uint32_t val;
    input_file.read(reinterpret_cast<char*>(&val), sizeof(uint32_t));
    radix_table.emplace_back(val);
  }

  input_file.read(reinterpret_cast<char*>(&spline_points_size), sizeof(size_t));
  for (size_t i = 0; i < spline_points_size; ++i) {
    uint64_t x;
    double y;
    input_file.read(reinterpret_cast<char*>(&x), sizeof(uint64_t));
    input_file.read(reinterpret_cast<char*>(&y), sizeof(double));
    rs::Coord<uint64_t> c(x, y);
    spline_points.push_back(c);
  }
  rs_ = rs::RadixSpline<uint64_t>(min_key, max_key, size, num_radix_bits_, num_shift_bits_, max_error_,
                                  radix_table, spline_points);
#endif

#ifdef PGM
  input_file.read(reinterpret_cast<char*>(&(pgm_.n)), sizeof(size_t));
  input_file.read(reinterpret_cast<char*>(&(pgm_.first_key)), sizeof(uint64_t));

  size_t segments_size, levels_sizes_size, levels_offsets_size;
  input_file.read(reinterpret_cast<char*>(&segments_size), sizeof(size_t));
  input_file.read(reinterpret_cast<char*>(&levels_sizes_size), sizeof(size_t));
  input_file.read(reinterpret_cast<char*>(&levels_offsets_size), sizeof(size_t));

  for (size_t i = 0; i < segments_size; ++i) {
    uint64_t p_key;
    double p_slope;
    int32_t p_intercept;
    input_file.read(reinterpret_cast<char*>(&p_key), sizeof(uint64_t));
    input_file.read(reinterpret_cast<char*>(&p_slope), sizeof(double));
    input_file.read(reinterpret_cast<char*>(&p_intercept), sizeof(int32_t));
    pgm_.segments.emplace_back(p_key, p_slope, static_cast<double>(p_intercept));
  }

  for (size_t i = 0; i < levels_sizes_size; ++i) {
    size_t p_levels_size;
    input_file.read(reinterpret_cast<char*>(&p_levels_size), sizeof(size_t));
    pgm_.levels_sizes.push_back(p_levels_size);
  }

  for (size_t i = 0; i < levels_offsets_size; ++i) {
    size_t p_levels_offset_size;
    input_file.read(reinterpret_cast<char*>(&p_levels_offset_size), sizeof(size_t));
    pgm_.levels_offsets.push_back(p_levels_offset_size);
  }
#endif

#ifdef CHT
  uint64_t c_min_key, c_max_key;
  size_t c_num_keys, c_num_bins, c_log_num_bins, c_max_error, c_shift;
  size_t c_table_size;
  input_file.read(reinterpret_cast<char*>(&c_min_key), sizeof(uint64_t));
  input_file.read(reinterpret_cast<char*>(&c_max_key), sizeof(uint64_t));
  input_file.read(reinterpret_cast<char*>(&c_num_keys), sizeof(size_t));
  input_file.read(reinterpret_cast<char*>(&c_num_bins), sizeof(size_t));
  input_file.read(reinterpret_cast<char*>(&c_log_num_bins), sizeof(size_t));
  input_file.read(reinterpret_cast<char*>(&c_max_error), sizeof(size_t));
  input_file.read(reinterpret_cast<char*>(&c_shift), sizeof(size_t));
  input_file.read(reinterpret_cast<char*>(&c_table_size), sizeof(size_t));

  std::vector<unsigned> c_table;
  c_table.reserve(c_table_size);

  for (size_t i = 0; i < c_table_size; ++i) {
    unsigned c_table_value;
    input_file.read(reinterpret_cast<char*>(&c_table_value), sizeof(unsigned));
    c_table.push_back(c_table_value);
  }

  cht_ = cht::CompactHistTree<uint64_t>(c_min_key, c_max_key, c_num_keys, c_num_bins, c_log_num_bins, c_max_error, c_shift, c_table);
#endif

#ifdef LINEAR
  double l_slope, l_intercept;
  uint64_t l_error;
  input_file.read(reinterpret_cast<char*>(&l_slope), sizeof(double));
  input_file.read(reinterpret_cast<char*>(&l_intercept), sizeof(double));
  input_file.read(reinterpret_cast<char*>(&l_error), sizeof(uint64_t));

  lm_ = linear::LinearModel(l_slope, l_intercept, size, l_error);
#endif

  input_file.close();

  learned.store(true);
}

// Reports the size of the model in bytes. 
// Note that this does not include the size indicators of vectors. That is, 
// when saving a vector as binary, the size of the vector as well as its contents
// are saved. In this function, the size of the size variable is not reported. 
size_t LearnedIndexData::ReportStats() {
  size_t model_size = 0;
#ifdef LCDE
  model_size = ko.getSize();
#endif

#ifdef BOURBON
  model_size = segment_size * string_segments.size();
#endif

#ifdef RS
  model_size = rs_.GetSize();
#endif

#ifdef PGM
  model_size = pgm_.size_in_bytes();
#endif

#ifdef CHT
  model_size = cht_.GetSize();
#endif

#ifdef LINEAR
  model_size = lm_.Size();
#endif

  model_size += fixed_header_size;
  // printf("%d %d %lu %lu %lu\n", level, served, model_size, cost, size);
  return model_size;
}

LearnedIndexData* FileLearnedIndexData::GetModel(int number) {
  leveldb::MutexLock l(&mutex);
  if (file_learned_index_data.size() <= number)
    file_learned_index_data.resize(number + 1, nullptr);
  if (file_learned_index_data[number] == nullptr)
    file_learned_index_data[number] = new LearnedIndexData(file_allowed_seek, false);
  return file_learned_index_data[number];
}

bool FileLearnedIndexData::FillData(Version* version, FileMetaData* meta) {
  LearnedIndexData* model = GetModel(meta->number);
  return model->FillData(version, meta);
}

std::vector<uint64_t>& FileLearnedIndexData::GetData(FileMetaData* meta) {
  auto* model = GetModel(meta->number);
  return model->string_keys;
}

bool FileLearnedIndexData::Learned(Version* version, FileMetaData* meta,
                                   int level) {
  LearnedIndexData* model = GetModel(meta->number);
  return model->Learned(version, db->version_count, meta, level);
}

std::pair<uint64_t, uint64_t> FileLearnedIndexData::GetPosition(
    const Slice& key, int file_num) {
  return file_learned_index_data[file_num]->GetPosition(key);
}

FileLearnedIndexData::~FileLearnedIndexData() {
  leveldb::MutexLock l(&mutex);
  for (auto pointer : file_learned_index_data) {
    delete pointer;
  }
}

void FileLearnedIndexData::Report() {
  leveldb::MutexLock l(&mutex);

  std::set<uint64_t> live_files;
  adgMod::db->versions_->AddLiveFiles(&live_files);

  size_t accumulated_file_model_size = 0;

  for (size_t i = 0; i < file_learned_index_data.size(); ++i) {
    auto pointer = file_learned_index_data[i];
    if (pointer != nullptr && pointer->cost != 0) {
      printf("FileModel %lu %d \n", i, i > watermark);
      accumulated_file_model_size += pointer->ReportStats();
    }
  }
  printf("AccumulatedFileModelSize:%zu\n", accumulated_file_model_size);
}

}  // namespace adgMod
