
#include <cassert>
#include <chrono>
#include <iostream>
#include "leveldb/db.h"
#include "leveldb/comparator.h"
#include "util.h"
#include "stats.h"
#include "learned_index.h"
#include <cstring>
#include <sstream>
#include "cxxopts.hpp"
#include <unistd.h>
#include <fstream>
#include "../db/version_set.h"
#include <cmath>
#include <random>

// For YCSB
#include "ycsb_utils.h"
#include "generators/generator.h"
#include "generators/acknowledged_counter_generator.h"
#include "generators/const_generator.h"
#include "generators/counter_generator.h"
#include "generators/discrete_generator.h"
#include "generators/random_byte_generator.h"
#include "generators/scrambled_zipfian_generator.h"
#include "generators/skewed_latest_generator.h"
#include "generators/uniform_generator.h"
#include "generators/zipfian_generator.h"

using namespace leveldb;
using namespace adgMod;
using std::string;
using std::cout;
using std::endl;
using std::to_string;
using std::vector;
using std::map;
using std::ifstream;
using std::ofstream;
using std::string;

int num_pairs_base = 1000;
int mix_base = 20;

enum LoadType {
  Ordered = 0,
  Reversed = 1,
  ReversedChunk = 2,
  Random = 3,
  RandomChunk = 4
};

enum Operation {
  INSERT = 0,
  READ,
  UPDATE,
  SCAN,
  READMODIFYWRITE,
  DELETE,
  INSERT_FAILED,
  READ_FAILED,
  UPDATE_FAILED,
  SCAN_FAILED,
  READMODIFYWRITE_FAILED,
  DELETE_FAILED,
  MAXOPTYPE
};

map<string, string> parse_workload(ifstream& input) {
  map<string, string> parse_result;
  if (!input.is_open()) {
    std::cerr << "Workload file is not open. Aborting...\n";
    exit(1);
  }
  while (!input.eof() && !input.bad()) {
    string line;
    getline(input, line);
    if (line[0] == '#') continue; // pass comment lines
    size_t pos = line.find_first_of('=');
    if (pos == string::npos) continue;
    parse_result[ycsbc::utils::Trim(line.substr(0, pos))] = ycsbc::utils::Trim(line.substr(pos + 1));
  }
  return parse_result;
}

int main(int argc, char *argv[]) {
  int rc;
  int num_operations, num_iteration, num_mix;
  float test_num_segments_base;
  float num_pair_step;
  string db_location, profiler_out, input_filename, distribution_filename, ycsb_filename;
  bool print_single_timing, print_file_info, evict, unlimit_fd, use_distribution = false, pause, use_ycsb = false;
  bool change_level_load, change_file_load, change_level_learning, change_file_learning;
  int load_type, insert_bound;
  string db_location_copy;
  string request_dist;

  bool print_timestamps;

  // YCSB workload type
  // workloads/workloada, workloadb, ... , workloadf
  string ycsb_workload_name;

  string output;

  cxxopts::Options commandline_options("leveldb read test", "Testing leveldb read performance.");
    commandline_options.add_options()
    ("s,step", "the step of the loop of the size of db", cxxopts::value<float>(num_pair_step)->default_value("1"))
    ("h,help", "print help message", cxxopts::value<bool>()->default_value("false"))
    ("single_timing", "print the time of every single get", cxxopts::value<bool>(print_single_timing)->default_value("false"))
    ("file_info", "print the file structure info", cxxopts::value<bool>(print_file_info)->default_value("false"))
    ("test_num_segments", "test: number of segments per level", cxxopts::value<float>(test_num_segments_base)->default_value("1"))
    ("string_mode", "test: use string or int in model", cxxopts::value<bool>(adgMod::string_mode)->default_value("false"))
    ("file_model_error", "error in file model", cxxopts::value<uint32_t>(adgMod::file_model_error)->default_value("18"))
    ("level_model_error", "error in level model", cxxopts::value<uint32_t>(adgMod::level_model_error)->default_value("1"))
    ("f,input_file", "the filename of input file", cxxopts::value<string>(input_filename)->default_value(""))
    ("multiple", "test: use larger keys", cxxopts::value<uint64_t>(adgMod::key_multiple)->default_value("1"))
    ("c,uncache", "evict cache", cxxopts::value<bool>(evict)->default_value("false"))
    ("u,unlimit_fd", "unlimit fd", cxxopts::value<bool>(unlimit_fd)->default_value("false"))
    ("l,load_type", "load type", cxxopts::value<int>(load_type)->default_value("3"))
    ("filter", "use filter", cxxopts::value<bool>(adgMod::use_filter)->default_value("false"))
    ("change_level_load", "load level model", cxxopts::value<bool>(change_level_load)->default_value("false"))
    ("change_file_load", "enable level learning", cxxopts::value<bool>(change_file_load)->default_value("false"))
    ("p,pause", "pause between operation", cxxopts::value<bool>(pause)->default_value("false"))
    ("policy", "learn policy", cxxopts::value<int>(adgMod::policy)->default_value("0"))
    // Modify from here
    // ("name", "name of the index used", cxxopts::value<string>(adgMod::index_name)->default_value("wisckey"))
    ("m,modification", "if set, run our modified version, 7 file-model bourbon, 8 wiskey, 9 ours", cxxopts::value<int>(adgMod::MOD)->default_value("10"))
    ("print", "print all recorded timestamps", cxxopts::value<bool>(print_timestamps)->default_value("false")) 
    ("w,write", "writedb", cxxopts::value<bool>(fresh_write)->default_value("false"))
    ("d,directory", "the directory of db", cxxopts::value<string>(db_location)->default_value("/tmp/wisckey_ycsb"))
    ("k,key_size", "the size of key", cxxopts::value<int>(adgMod::key_size)->default_value("16"))
    ("v,value_size", "the size of value", cxxopts::value<int>(adgMod::value_size)->default_value("64"))
    ("dataset", "name of the dataset", cxxopts::value<string>(adgMod::dataset_name)->default_value("osm_cellids_10M_uint64"))
    ("i,iteration", "the number of iterations of a same size", cxxopts::value<int>(num_iteration)->default_value("3"))
    ("n,num_operations", "the number of operations", cxxopts::value<int>(num_operations)->default_value("100000"))
    ("r,req_distribution", "the distribution of operations", cxxopts::value<string>(request_dist)->default_value("uniform"))
    ("workload", "YCSB workload", cxxopts::value<string>(ycsb_workload_name)->default_value("workloadc"));

  auto result = commandline_options.parse(argc, argv);
  if (result.count("help")) {
    std::cout << commandline_options.help() << std::endl;
    exit(0);
  }

  // // set name to lowercase
  // std::transform(index_name.begin(), index_name.end(), index_name.begin(), tolower);
#ifdef WISCKEY
  cout << "index:WISCKEY" << endl;
  index_name = "wisckey";
  adgMod::MOD = 8;
#endif
#ifdef LCDE
  cout << "index:LCDE" << endl;
  index_name = "lcde";
  adgMod::MOD = 7;
#endif
#ifdef BOURBON
  cout << "index:BOURBON" << endl;
  index_name = "bourbon";
  adgMod::MOD = 7;
#endif
#ifdef RS
  cout << "index:RS" << endl;
  index_name = "rs";
  adgMod::MOD = 7;
#endif
#ifdef PGM
  cout << "index:PGM" << endl;
  index_name = "pgm";
  adgMod::MOD = 7;
#endif
#ifdef CHT
  cout << "index:CHT" << endl;
  index_name = "cht";
  adgMod::MOD = 7;
#endif
#ifdef LINEAR
  cout << "index:LINEAR" << endl;
  index_name = "linear";
  adgMod::MOD = 7;
#endif
  if (adgMod::MOD == 10) {
    std::cerr << "Index is invalid. Aborting...\n";
    exit(EXIT_FAILURE);
  }
  
  // configure db location
  int index_variant;
#ifdef INDEX_VARIANT
  index_variant = INDEX_VARIANT;
#else
  index_variant = 0;
#endif
  db_location = "/tmp/" + index_name + "_" + to_string(index_variant) + "_ycsb";
  
  // configure ycsb workload name
  ycsb_workload_name = "../workloads/" + ycsb_workload_name;
  
  std::default_random_engine e1(0), e2(255), e3(0);
  srand(0);

  adgMod::fd_limit = unlimit_fd ? 1024 * 1024 : 1024;
  adgMod::restart_read = true;
  if (adgMod::MOD == 7) {        // Bourbon; only allow file learning
    adgMod::file_learning_enabled = true;
    adgMod::load_file_model = true;
  } else if (adgMod::MOD == 9) { // Modified Bourbon; only allow level learning.
                                 // this configuration is not used.
    adgMod::level_learning_enabled = true;
    adgMod::load_level_model = true;
  }

  // vector storing all keys
  vector<string> keys;
  // if the workload is fully loaded, |distribution| == |ycsb_is_write|
  vector<uint64_t> distribution;
  vector<bool> ycsb_is_write;

  // print benchmark informations
  cout << "location:" << db_location << endl;
  cout << "dataset:" << dataset_name << endl;
  cout << "workload:" << ycsb_workload_name << endl;
  cout << "operations:" <<num_operations << endl;
  // load keys from dataset
  vector<uint64_t> sosd_keys = ycsbc::utils::load_data<uint64_t>(dataset_name);
  // for calculating minimum distance between keys
  uint64_t cur_key = 0;
  // min_distance is used to generate a non-existing key for negative queries
  // or fresh insert transactions
  uint64_t min_distance = (uint64_t)1 << 63;
  for (uint64_t sosd_key : sosd_keys) {
    sosd_key >>= 14;
    if (sosd_key - cur_key < min_distance) {
      min_distance = sosd_key - cur_key;
    }
    cur_key = sosd_key;
    string raw_key_string = to_string(sosd_key);
    string the_key = generate_key2(raw_key_string);
    keys.push_back(the_key);
  }
  
  const bool copy_out = true;

  adgMod::Stats* instance = adgMod::Stats::GetInstance();
  vector<vector<size_t>> times(20);
  string values(1024 * 1024, '0');
  
  if (copy_out) {
    // Flush cache
    rc = system("sync; echo 3 | sudo tee /proc/sys/vm/drop_caches");
  }

  // Random generator for value
  std::uniform_int_distribution<uint64_t > uniform_dist_value(0, (uint64_t) values.size() - adgMod::value_size - 1);

  //////////////////////////////////////////////////////////////////////////////
  // Declare DB and Options
  //////////////////////////////////////////////////////////////////////////////
  // DB* db;
  // Options options;
  // ReadOptions& read_options = adgMod::read_options;
  // WriteOptions& write_options = adgMod::write_options;
  // Status status;

  // options.create_if_missing = true;
  // write_options.sync = false;
  // instance->ResetAll();

  //////////////////////////////////////////////////////////////////////////////
  // Load Dataset
  //////////////////////////////////////////////////////////////////////////////
  if (fresh_write) {
    ////////////////////////////////////////////////////////////////////////////
    // Declare DB and Options
    ////////////////////////////////////////////////////////////////////////////
    DB* db;
    Options options;
    ReadOptions& read_options = adgMod::read_options;
    WriteOptions& write_options = adgMod::write_options;
    Status status;

    options.create_if_missing = true;
    write_options.sync = false;
    instance->ResetAll();

    // Remove DB if exists
    string command = "rm -rf " + db_location;
    rc = system(command.c_str());

    // Open DB
    status = DB::Open(options, db_location, &db);
    assert(status.ok() && "Open Error");

    instance->StartTimer(9);
    switch (load_type) {
      case Ordered: {
        for (int i = 0; i < (int)keys.size(); ++i) {
          status = db->Put(write_options, keys[i], generate_value(uniform_dist_value(e2)));
          if (!status.ok()) cout << status.ToString() << endl;
          assert(status.ok() && "File Put Error");
        }
        break;
      }
      case ReversedChunk: {
        for (int i = (int)keys.size(); i >= 0; --i) {
          status = db->Put(write_options, keys[i], generate_value(uniform_dist_value(e2)));
          assert(status.ok() && "File Put Error");
        }
        break;
      }
      case Random: {
        std::random_shuffle(keys.begin(), keys.end());
        for (int i = 0; i < (int)keys.size(); ++i) {
          if (keys[i] =="0224184788660671") cout << "It's in bro...\n";
          status = db->Put(write_options, keys[i], generate_value(uniform_dist_value(e2)));
          assert(status.ok() && "File Put Error");
        }
        break;
      }
      case RandomChunk: {
        break;
      }
      default: assert(false && "Unsupported load type.");
    }

    adgMod::db->vlog->Sync();
    instance->PauseTimer(9, true);
    instance->ReportTime();
    cout << "Put Complete" << endl;

    if (print_file_info) db->PrintFileInfo();

    if (print_timestamps) {
      // cout << "\n===== Reporting events =====\n";
      // for (auto& event_array : events) {
      //   for (Event* e : event_array) e->Report();
      // }

      // cout << "\n===== Reporting levelled counters =====\n";
      // for (Counter& c : levelled_counters) c.Report();

      cout << "\n===== Reporting file data =====\n";
      file_data->Report();

      // cout << "\n===== Level model stats =====" << endl;
      // Version* current = adgMod::db->versions_->current();
      // for (int i = 1; i < config::kNumLevels; ++i) {
      //   current->learned_index_data_[i]->ReportStats();
      // }

      // cout << "\n===== Report file stats =====" << endl;
      // for (auto it : file_stats) {
      //   printf("FileStats %d %d %lu %lu %u %u %lu %d\n", it.first, it.second.level, it.second.start,
      //     it.second.end, it.second.num_lookup_pos, it.second.num_lookup_neg, it.second.size, it.first < file_data->watermark ? 0 : 1);
      // }

      // cout << "\n===== Report learn CB model =====" << endl;
      // adgMod::learn_cb_model->Report();
    }


    ////////////////////////////////////////////////////////////////////////////
    // Offline Learning Phase
    ////////////////////////////////////////////////////////////////////////////
    // adgMod::db->WaitForBackground();
    // delete db;
    // status = DB::Open(options, db_location, &db);
    // adgMod::db->WaitForBackground();

    // // Begin offline learning
    // if (adgMod::MOD == 6 || adgMod::MOD == 7 || adgMod::MOD == 9) {
    //   Version* current = adgMod::db->versions_->current();

    //   // for (int i = 1; i < config::kNumLevels; ++i) {
    //   //   LearnedIndexData::LevelLearn(new VersionAndSelf{current, adgMod::db->version_count, current->learned_index_data_[i].get(), i});
    //   // }

    //   current->FileLearn();
    // }

    // status = DB::Open(options, db_location, &db);
    // adgMod::db->WaitForBackground();

    // if (print_timestamps) {
    //   cout << "\n===== Reporting events =====\n";
    //   for (auto& event_array : events) {
    //     for (Event* e : event_array) e->Report();
    //   }

    //   cout << "\n===== Reporting levelled counters =====\n";
    //   for (Counter& c : levelled_counters) c.Report();

    //   cout << "\n===== Reporting file data =====\n";
    //   file_data->Report();

    //   cout << "\n===== Level model stats =====" << endl;
    //   Version* current = adgMod::db->versions_->current();
    //   for (int i = 1; i < config::kNumLevels; ++i) {
    //     current->learned_index_data_[i]->ReportStats();
    //   }

    //   cout << "\n===== Report file stats =====" << endl;
    //   for (auto it : file_stats) {
    //     printf("FileStats %d %d %lu %lu %u %u %lu %d\n", it.first, it.second.level, it.second.start,
    //       it.second.end, it.second.num_lookup_pos, it.second.num_lookup_neg, it.second.size, it.first < file_data->watermark ? 0 : 1);
    //   }

    //   cout << "\n===== Report learn CB model =====" << endl;
    //   adgMod::learn_cb_model->Report();
    // }

    // // cout << "Offline learning complete... turning off DB" << endl;
    adgMod::db->WaitForBackground();
    delete db;

    return 0;
  }

  // for (int iter = 0; iter < num_iteration; ++iter) {
  ////////////////////////////////////////////////////////////////////////////
  // Declare DB and Options
  ////////////////////////////////////////////////////////////////////////////
  DB* db;
  Options options;
  ReadOptions& read_options = adgMod::read_options;
  WriteOptions& write_options = adgMod::write_options;
  Status status;

  options.create_if_missing = true;
  write_options.sync = false;
  instance->ResetAll();

  ////////////////////////////////////////////////////////////////////////////
  // Open DB for Transactions
  ////////////////////////////////////////////////////////////////////////////
  cout << "Starting up DB for transactions" << endl;
  status = DB::Open(options, db_location, &db);
  adgMod::db->WaitForBackground();
  assert(status.ok() && "Open Error");

  // testing something wrong...
  // string op_result;
  // status = db->Get(read_options, keys[1404998], &op_result);
  // if (!status.ok()) {
  //   cout << "Could not find " << keys[1404998] << endl;
  // }
  // exit(EXIT_SUCCESS);

  ////////////////////////////////////////////////////////////////////////////
  // Load YCSB Workload Distribution
  ////////////////////////////////////////////////////////////////////////////
  ifstream ycsb_workload_file(ycsb_workload_name);
  cout << ycsb_workload_name << endl;
  map<string, string> parsed_workload = parse_workload(ycsb_workload_file);

  double read_proportion, update_proportion, insert_proportion, scan_proportion,
        readmodifywrite_proportion;
  // int op_count;

  read_proportion = stod(parsed_workload["readproportion"] == "" ? "0" : parsed_workload["readproportion"]);
  update_proportion = stod(parsed_workload["updateproportion"] == "" ? "0" : parsed_workload["updateproportion"]);
  scan_proportion = stod(parsed_workload["scanproportion"] == "" ? "0" : parsed_workload["scanproportion"]);
  insert_proportion = stod(parsed_workload["insertproportion"] == "" ? "0" : parsed_workload["insertproportion"]);
  readmodifywrite_proportion = stod(parsed_workload["readmodifywriteproportion"] == "" ? "0" : parsed_workload["readmodifywriteproportion"]);
  // op_count = stoi(parsed_workload["operationcount"] == "" ? "0" : parsed_workload["operationcount"]);
  request_dist = parsed_workload["requestdistribution"] == "" ? "zipfian" : parsed_workload["requestdistribution"];

  ////////////////////////////////////////////////////////////////////////////
  // Declare Variables for Transaction
  ////////////////////////////////////////////////////////////////////////////
  uint64_t last_read = 0, last_write = 0;
  int last_level = 0, last_file = 0, last_baseline = 0, last_succeeded = 0, last_false = 0, last_compaction = 0, last_learn = 0;
  // std::vector<uint64_t> detailed_times;
  bool start_new_event = true;

  ////////////////////////////////////////////////////////////////////////////
  // Declare Generators
  ////////////////////////////////////////////////////////////////////////////
  ycsbc::Generator<uint64_t> *field_len_generator_;
  ycsbc::DiscreteGenerator<Operation> op_chooser_;
  ycsbc::Generator<uint64_t> *key_chooser_; // transaction key gen
  ycsbc::Generator<uint64_t> *field_chooser_;
  ycsbc::Generator<uint64_t> *scan_len_chooser_;
  // ycsbc::CounterGenerator *insert_key_sequence_; // load insert key gen
  ycsbc::AcknowledgedCounterGenerator *transaction_insert_key_sequence_; // transaction insert key gen

  ////////////////////////////////////////////////////////////////////////////
  // Initialize Generators
  ////////////////////////////////////////////////////////////////////////////
  // 1. Initialize op chooser
  //    = choose the next operation given a workload
  if (read_proportion > 0) {
    op_chooser_.AddValue(READ, read_proportion);
  }
  if (update_proportion > 0) {
    op_chooser_.AddValue(UPDATE, update_proportion);
  }
  if (insert_proportion > 0) {
    op_chooser_.AddValue(INSERT, insert_proportion);
  }
  if (scan_proportion > 0) {
    op_chooser_.AddValue(SCAN, scan_proportion);
  }
  if (readmodifywrite_proportion > 0) {
    op_chooser_.AddValue(READMODIFYWRITE, readmodifywrite_proportion);
  }
  // 2. Initialize transaction insert key generator
  //    = choose the next key to be inserted
  transaction_insert_key_sequence_ = new ycsbc::AcknowledgedCounterGenerator(keys.size());
  // 3. Initialize request distribution generator
  //    = defines the underlying distribution by which the next key is chosen
  if (request_dist == "uniform") {
    key_chooser_ = new ycsbc::UniformGenerator(0, keys.size() - 1);
  } else if (request_dist == "zipfian") {
    // use the default zipfian const = 0.99
    key_chooser_ = new ycsbc::ScrambledZipfianGenerator(keys.size() - 1);
  } else if (request_dist == "latest") {
    key_chooser_ = new ycsbc::SkewedLatestGenerator(*transaction_insert_key_sequence_);
  } else {
    std::cerr << "Unknown request distribution: " << request_dist << endl;
    exit(1);
  }
  scan_len_chooser_ = new ycsbc::UniformGenerator(0, 100);

  // Generate dummy keys for inserts
  uint64_t max_uint_key = std::stoull(keys[keys.size() - 1]) + 1;

  

  uint64_t write_i = 0;
  uint64_t key_idx;         // index of the key for reading or inserting
  uint64_t tmp_key_idx;     // index of the key for generating dummy insert key
                            // we make a dummy key for inserting new, previously
                            // not inserted keys
  string op_key;            // actual key value to be queried
  string op_value;          // placeholder for read or inserted value
  int scan_len = 100; // scan 100 entries
  uint64_t wrong_ops = 0;
  for (int i = 0; i < num_operations; ++i) {
    // if (start_new_event) {
    //   detailed_times.push_back(instance->GetTime());
    //   start_new_event = false;
    // }
    do {
      key_idx = key_chooser_->Next();
    } while (key_idx > transaction_insert_key_sequence_->Last()); 
    op_key = keys[key_idx];

    // cout << op_key;
    scan_len = scan_len_chooser_->Next();

    instance->StartTimer(13);

    switch (op_chooser_.Next()) {
      case READ: {
        instance->StartTimer(4);
        status = db->Get(read_options, op_key, &op_value);
        instance->PauseTimer(4);
        if (!status.ok()) {
          // cout << "KeyNotFound(READ)  : iteration " << i << ", " \
          //     << "key index " << key_idx << endl;
          ++wrong_ops;
        }
        break;
      }
      case UPDATE: {
        op_value = generate_value(uniform_dist_value(e2));
        instance->StartTimer(10);
        status = db->Put(write_options, op_key, op_value);
        instance->PauseTimer(10);
        if (!status.ok()) {
          // cout << "UpdateFail(UPDATE) : iteration " << i << ", " \
          //     << "key index " << key_idx << endl;
          ++wrong_ops;
        }
        break;
      }
      case INSERT: {
        // cout << " insert ";
        // tmp_key_idx = key_idx + 1;
        // op_key = ycsbc::utils::average_keys(op_key, keys[tmp_key_idx]);
        op_key = to_string(max_uint_key++);
        op_key = generate_key2(op_key);
        // cout << " insert ";
        op_value = generate_value(uniform_dist_value(e2));
        // cout << " insert ";
        instance->StartTimer(10);
        status = db->Put(write_options, op_key, op_value);
        // cout << " insert ";
        instance->PauseTimer(10);
        if (!status.ok()) {
          // cout << "InsertFail(INSERT) : iteration " << i << ", " \
          //     << "key index " << key_idx << endl;
          ++wrong_ops;
        }
        break;
      }
      case SCAN: {
        instance->StartTimer(17);
        leveldb::Iterator *db_iter = db->NewIterator(read_options);
        vector<string> tmp_result_vector;
        db_iter->Seek(op_key);
        for (int j = 0; db_iter->Valid() && j < scan_len; ++j) {
          tmp_result_vector.push_back(db_iter->value().ToString());
          db_iter->Next();
        }
        delete db_iter;
        instance->PauseTimer(17);
        break;
      }
      case READMODIFYWRITE: {
        // read
        instance->StartTimer(4);
        status = db->Get(read_options, op_key, &op_value);
        instance->PauseTimer(4);
        if (!status.ok()) {
          // cout << "KeyNotFound(RMW)   : iteration " << i << ", " \
          //     << "key index " << key_idx << endl;
          ++wrong_ops;
          break;
        }
        // modify
        op_value = generate_value(uniform_dist_value(e2));
        // write
        instance->StartTimer(10);
        status = db->Put(write_options, op_key, op_value);
        instance->PauseTimer(10);
        if (!status.ok()) {
          // cout << "InsertFail(RMW)    : iteration " << i << ", " \
          //     << "key index " << key_idx << endl;
          ++wrong_ops;
        }
        break;
      }
      default: {
        std::cerr << "Operation request is not recognized. Aborting...\n";
        exit(EXIT_FAILURE);
      }
    }
    instance->PauseTimer(13, true);
  }

  cout << "Correct : " << num_operations << " Wrong : " << wrong_ops << " Wrong ratio : " << (double)wrong_ops / (double)num_operations << endl;


  // report various data after the run
  // if (print_timestamps) {
    instance->ReportTime();
  // }
  for (int s = 0; s < times.size(); ++s) {
    times[s].push_back(instance->ReportTime(s));
  }
  adgMod::db->WaitForBackground();
  sleep(2);

  if (print_timestamps) {
    cout << "\n===== Reporting events =====\n";
    for (auto& event_array : events) {
      for (Event* e : event_array) e->Report();
    }

    cout << "\n===== Reporting levelled counters =====\n";
    for (Counter& c : levelled_counters) c.Report();

    cout << "\n===== Reporting file data =====\n";
    file_data->Report();

    cout << "\n===== Level model stats =====" << endl;
    Version* current = adgMod::db->versions_->current();
    for (int i = 1; i < config::kNumLevels; ++i) {
      current->learned_index_data_[i]->ReportStats();
    }

    cout << "\n===== Report file stats =====" << endl;
    for (auto it : file_stats) {
      printf("FileStats %d %d %lu %lu %u %u %lu %d\n", it.first, it.second.level, it.second.start,
        it.second.end, it.second.num_lookup_pos, it.second.num_lookup_neg, it.second.size, it.first < file_data->watermark ? 0 : 1);
    }

    cout << "\n===== Report learn CB model =====" << endl;
    adgMod::learn_cb_model->Report();
  }

  adgMod::db->WaitForBackground();
  delete db;
  // }

  // for (int s = 0; s < times.size(); ++s) {
  //   vector<uint64_t>& time = times[s];
  //   vector<double> diff(time.size());
  //   if (time.empty()) continue;

  //   double sum = std::accumulate(time.begin(), time.end(), 0.0);
  //   double mean = sum / time.size();
  //   std::transform(time.begin(), time.end(), diff.begin(), [mean] (double x) { return x - mean; });
  //   double stdev = std::sqrt(std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0) / time.size());

  //   printf("Timer %d MEAN: %lu, STDDEV: %f\n", s, (uint64_t) mean, stdev);
  // }

  return 0;
}
