#ifndef BESTPARTH
#define BESTPARTH

#include <math.h>
#include <omp.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

#include "Bloom.h"
#include "ExponentialBloom.h"
#include "best.h"
#include "utils.h"

using namespace std;
using namespace filesystem;

template <class T>
class BestPart {
   public:
    vector<Best<T>*> buckets;
    omp_lock_t* mutex_array;
    uint K;
    uint64_t offsetUpdatekmer;
    uint64_t offsetUpdateminimizer;
    uint64_t size;
    uint64_t leaf_number;
    uint64_t core_number;
    uint64_t number_bit_set;
    uint64_t disk_space_used;
    uint64_t number_bit_set_abt;
    uint number_hash_function;
    uint small_minimizer_size;
    uint large_minimizer_size;
    uint bucket_number;
    uint large_minimizer_number;
    bool filter;
    bool use_double_index;
    string w_dir;

    BestPart(const uint64_t Isize, const uint Inumber_hash_function, const uint Ik, bool Ifilter, const string Iwdir, bool Idouble, uint bucketing, uint Icore_number) {
        core_number = Icore_number;
        K = Ik;
        use_double_index = Idouble;
        w_dir = (Iwdir);
        path p{w_dir};
        if (exists(p)) {
        } else {
            create_directory(p);
        }
        filter = Ifilter;
        leaf_number = 0;
        small_minimizer_size = bucketing;
        bucket_number = 1 << (2 * small_minimizer_size);
        large_minimizer_size = small_minimizer_size + 3;
        large_minimizer_number = 1 << (2 * large_minimizer_size);
        size = Isize;
        number_hash_function = Inumber_hash_function;
        buckets.resize(bucket_number, NULL);
        mutex_array = new omp_lock_t[bucket_number];
        for (uint32_t i = 0; i < bucket_number; ++i) {
            buckets[i] = new Best<T>(size / bucket_number, number_hash_function, K, w_dir + "/P" + to_string(i));
            omp_init_lock(&mutex_array[i]);
        }
        offsetUpdatekmer = 1;
        offsetUpdatekmer <<= 2 * K;
        offsetUpdateminimizer = 1;
        offsetUpdateminimizer <<= (2 * large_minimizer_size);
        disk_space_used = 0;
        number_bit_set_abt = 0;
    }

    BestPart(const string& existing_index, uint Icore) {
        core_number = Icore;
        cout << "Loading " << existing_index << "..." << endl;
        path initial_path = current_path();
        w_dir = existing_index;
        path p{existing_index};
        if (exists(p)) {
        } else {
            cout << "I cannot found you index sorry ..." << endl;
            return;
        }
        current_path(p);
        string filename("MainIndex");
        zstr::ifstream in(filename);
        // VARIOUS INTEGERS
        in.read(reinterpret_cast<char*>(&K), sizeof(K));
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        in.read(reinterpret_cast<char*>(&number_hash_function), sizeof(number_hash_function));
        in.read(reinterpret_cast<char*>(&small_minimizer_size), sizeof(small_minimizer_size));
        in.read(reinterpret_cast<char*>(&large_minimizer_size), sizeof(large_minimizer_size));
        in.read(reinterpret_cast<char*>(&leaf_number), sizeof(leaf_number));
        in.read(reinterpret_cast<char*>(&use_double_index), sizeof(use_double_index));
        bucket_number = 1 << (2 * small_minimizer_size);
        large_minimizer_number = 1 << (2 * large_minimizer_size);
        mutex_array = new omp_lock_t[bucket_number];
        offsetUpdatekmer = 1;
        offsetUpdatekmer <<= 2 * K;
        offsetUpdateminimizer = 1;
        offsetUpdateminimizer <<= (2 * large_minimizer_size);
        for (uint32_t i = 0; i < bucket_number; ++i) {
            omp_init_lock(&mutex_array[i]);
        }
        buckets.resize(bucket_number, NULL);
        // buckets
        for (size_t i = 0; i < bucket_number; ++i) {
            buckets[i] = new Best<T>(size / bucket_number, number_hash_function, K, existing_index + "/P" + to_string(i), false);
            // buckets[i]->load(leaf_number,use_double_index);
        }
        current_path(initial_path);
        cout << "The index  use " << K << "mers with Bloom filters of size " << intToString(size) << " with " << number_hash_function << " hash functions  using " << intToString(1 << (2 * small_minimizer_size)) << " partitions " << endl;
    }

    ~BestPart() {
        for (uint32_t i = 0; i < bucket_number; ++i) {
            if (buckets[i] != NULL) {
                delete buckets[i];
            }
        }
        delete[] mutex_array;
    }

    // LOW LEVEL FUNCTIONS
    void check_key_leaf(const uint64_t key, const uint level) const;
    uint64_t rcb(uint64_t min, uint64_t n) const;
    void updateK(uint64_t& min, char nuc) const;
    void updateM(uint64_t& min, char nuc) const;
    void updateRCK(uint64_t& min, char nuc) const;
    void updateRCM(uint64_t& min, char nuc) const;
    void insert_last_leaf_trunk();
    uint64_t regular_minimizer_pos(uint64_t seq, uint64_t& position);
    uint64_t canonize(uint64_t x, uint64_t n);
    void insert_keys(const vector<uint64_t>& key, uint minimizer, uint level, Bloom<T>* unique_filter);
    vector<T> query_keys(const vector<uint64_t>& key, uint minimizer);
    vector<pair<vector<uint64_t>, uint64_t> > get_super_kmers(const string& ref);
    // HIGH LEVEL FUNCTIONS
    void insert_sequence(const string& reference, uint level, Bloom<T>* unique_filter);
    void insert_file(const string& filename, uint level, uint32_t indice_bloom);
    void insert_file_of_file(const string& filename);
    vector<uint> query_key(const uint64_t key);
    vector<uint32_t> query_sequence(const string& reference);
    void load_super_kmer(vector<vector<pair<uint64_t, uint32_t> > >& colored_kmer_per_bucket, uint32_t query_id, const string& reference);
    uint64_t query_bucket(const vector<pair<uint64_t, uint32_t> >& colored_kmer, vector<pair<string, vector<uint32_t> > >& result, uint bucket_number);
    uint64_t query_bucket(const vector<pair<uint64_t, uint32_t> >& colored_kmer, vector<pair<string, vector<uint32_t> > >& result, uint bucket_number, vector<bool>* bouly);
    uint64_t query_bucket2(const vector<pair<uint64_t, uint32_t> >& colored_kmer, vector<pair<string, vector<uint32_t> > >& result, uint bucket_number);
    void insert_previous_index(const string& filename);

    void get_stats() const;
    void optimize();
    void query_file(const string& reference, const string& output);
    void serialize() const;
    void load(const string& existing_index);
    void index();
    void double_index();
    void add_leaf();
    uint64_t nb_bit_set() const;
    uint64_t nb_insertions() const;
};

template class BestPart<uint8_t>;
template class BestPart<uint16_t>;
template class BestPart<uint32_t>;

#endif
