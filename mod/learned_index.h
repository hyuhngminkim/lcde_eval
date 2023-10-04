//
// Created by daiyi on 2020/02/02.
//

#ifndef LEVELDB_LEARNED_INDEX_H
#define LEVELDB_LEARNED_INDEX_H


#include <vector>
#include <cstring>
#include "util.h"
#include <atomic>
#include <iostream>

#ifdef LCDE
#include "./competitors/lcde/include/lcde/builderObject.h"
#include "./competitors/lcde/include/lcde/knotObject.h"
#endif

#ifdef BOURBON
#include "plr.h"
#endif

#ifdef RS
// #include "./competitors/RadixSpline/include/rs/common.h"
#include "./competitors/rs/include/rs/builder.h"
#include "./competitors/rs/include/rs/radix_spline.h"
#endif

#ifdef PGM
#include "./competitors/PGM-index/include/pgm_index.hpp"
#endif

#ifdef CHT
#include "./competitors/CHT/include/cht/builder.h"
#include "./competitors/CHT/include/cht/cht.h"
#endif

#ifdef LINEAR
#include "./competitors/linear/linear.h"
#endif

using std::string;
using leveldb::Slice;
using leveldb::Version;
using leveldb::FileMetaData;



namespace adgMod {

    class LearnedIndexData;

    class VersionAndSelf {
    public:
        Version* version;
        int v_count;
        LearnedIndexData* self;
        int level;
    };

    class MetaAndSelf {
    public:
        Version* version;
        int v_count;
        FileMetaData* meta;
        LearnedIndexData* self;
        int level;
    };

    // The structure for learned index. Could be a file model or a level model
    class LearnedIndexData {
        friend class leveldb::Version;
        friend class leveldb::VersionSet;
    private:
        // predefined model error
        double error;
        // some flags used in online learning to control the state of the model
        std::atomic<bool> learned;
        std::atomic<bool> aborted;
        bool learned_not_atomic;
        std::atomic<bool> learning;
        // some params for level triggering policy, deprecated
        int allowed_seek;
        int current_seek;
    public:
        // is the data of this model filled (ready for learning)
        bool filled;
        // is this a level model
        bool is_level;

        // Learned linear segments and some other data needed
        // std::vector<Segment> string_segments;

#ifdef LCDE
        lcde::KnotObject<uint64_t> ko;
#endif

#ifdef BOURBON
        // int bb_variant = 7;
        int bb_variant = INDEX_VARIANT;
        std::vector<double> bb_configs = 
        {64, 32, 28, 24, 20, 16, 12, 8, 4, 2};
        std::vector<Segment> string_segments;
        const size_t segment_size = sizeof(uint64_t) + 2 * sizeof(double);
#endif

#ifdef RS
        // int rs_variant = 2;
        int rs_variant = INDEX_VARIANT;
        rs::RadixSpline<uint64_t> rs_;
        std::vector<std::pair<size_t, size_t>> rs_configs = 
        {{4, 8192}, {8, 4096}, {10, 2048}, {14, 512}, {20, 320}, {20, 160}, {24, 40}, {24, 20}, {26, 8}, {26, 3}};
#endif

#ifdef PGM
        std::vector<size_t> pgm_max_error = 
        {4, 8, 16, 32, 64, 256, 1024, 2048, 4096, 8192};
        PGMIndex<uint64_t, 4, 4> pgm_;
#endif

#ifdef CHT
        int cht_variant = INDEX_VARIANT;
        // int cht_variant = 4;
        cht::CompactHistTree<uint64_t> cht_;
        std::vector<std::pair<size_t, size_t>> cht_configs = 
        {{8, 2048}, {16, 2048}, {16, 1024}, {32, 1024}, {32, 512}, {64, 256}, {64, 128}, {64, 32}, {64, 16}, {1024, 16}};
#endif

#ifdef LINEAR
        linear::LinearModel lm_;
#endif

        // accumulated header size in the written binary file
        // block_num_entries, block_size, entry_size, min_key, max_key, size,
        // level, cost
        const size_t fixed_header_size = sizeof(uint64_t) * 7 + sizeof(int);
        uint64_t min_key;
        uint64_t max_key;
        uint64_t size;



    public:
        // all keys in the file/level to be leraned from
        std::vector<uint64_t> string_keys;

        int level;
        mutable int served;
        uint64_t cost;

        explicit LearnedIndexData(int allowed_seek, bool level_model) : error(level_model?level_model_error:file_model_error), learned(false), aborted(false), learning(false),
            learned_not_atomic(false), allowed_seek(allowed_seek), current_seek(0), filled(false), is_level(level_model), level(0), served(0), cost(0) {};
        LearnedIndexData(const LearnedIndexData& other) = delete;

        // Inference function. Return the predicted interval.
        // If the key is in the training set, the output interval guarantees to include the key
        // otherwise, the output is undefined!
        // If the output lower bound is larger than MaxPosition(), the target key is not in the file
        std::pair<uint64_t, uint64_t> GetPosition(const Slice& key) const;
        uint64_t MaxPosition() const;
        double GetError() const;
        
        // Learning function and checker (check if this model is available)
        bool Learn();
        bool Learned();
        bool Learned(Version* version, int v_count, int level);
        bool Learned(Version* version, int v_count, FileMetaData* meta, int level);
        // static void LevelLearn(void* arg, bool no_lock=false);
        static uint64_t FileLearn(void* arg);

        // Load all the keys in the file/level
        bool FillData(Version* version, FileMetaData* meta);

        // writing this model to disk and load this model from disk
        void WriteModel(const string& filename);
        void ReadModel(const string& filename);

        // writing model and loading disk from disk (binary)
        void WriteModelBinary(const string& filename);
        void ReadModelBinary(const string& filename);
        
        // print model stats
        // MOD: modified to return file model size so that it can be accumulated
        size_t ReportStats();

        bool Learn(bool file);
    };

    // an array storing all file models and provide similar access interface with multithread protection
    class FileLearnedIndexData {
    private:
        leveldb::port::Mutex mutex;
        std::vector<LearnedIndexData*> file_learned_index_data;
    public:
        uint64_t watermark;


        bool Learned(Version* version, FileMetaData* meta, int level);
        bool FillData(Version* version, FileMetaData* meta);
        std::vector<uint64_t>& GetData(FileMetaData* meta);
        std::pair<uint64_t, uint64_t> GetPosition(const Slice& key, int file_num);
        LearnedIndexData* GetModel(int number);
        void Report();
        ~FileLearnedIndexData();
    };

    class LevelLearnedIndexData {
     private:
      leveldb::port::Mutex mutex;
      std::vector<LearnedIndexData*> level_learned_index_data;
     public:

    };


}

#endif //LEVELDB_LEARNED_INDEX_H
