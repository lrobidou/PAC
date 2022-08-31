#ifndef PACH
#define PACH

#include <libpac/Bloom.h>
#include <libpac/utils.h>
#include <math.h>
#include <omp.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <libpac/Bucket.hpp>
#include <vector>

template <typename T>
class Pac {
   private:
    uint64_t _nbCells;
    uint64_t _nbBitsPerCell;
    uint _k;
    uint _s;
    uint _number_hash_function;
    uint _small_minimizer_size;
    uint _large_minimizer_size;
    uint _large_minimizer_number;
    uint _bucket_number;

    uint64_t _core_number;
    std::string _w_dir;
    uint64_t _leaf_number;
    uint64_t _number_bit_set;
    uint64_t _number_bit_set_abt;
    uint64_t _disk_space_used;

    omp_lock_t* _mutex_array;
    uint64_t _offsetUpdatekmer;
    uint64_t _offsetUpdateminimizer;

    bool _use_double_index;

   public:
    std::vector<Bucket<T>*> _buckets;  // TODO
    Pac(const uint64_t nbCells, const uint64_t nbBitsPerCell, const uint number_hash_function, const uint k, const uint z, const std::string wdir, bool use_double_index, uint bucketing, uint core_number);
    // Pac(const std::string& existing_index, uint Icore);
    ~Pac();

    // LOW LEVEL FUNCTIONS
    void check_key_leaf(const uint64_t key, const uint level) const;
    // uint64_t rcb(uint64_t min, uint64_t n) const;
    void updateK(uint64_t& min, char nuc) const;
    void updateM(uint64_t& min, char nuc) const;
    void updateRCK(uint64_t& min, char nuc) const;
    void updateRCM(uint64_t& min, char nuc) const;
    void insert_last_leaf_trunk();
    uint64_t regular_minimizer_pos(uint64_t seq, uint64_t& position) const;
    // uint64_t canonize(uint64_t x, uint64_t n);
    void insert_keys(const std::vector<uint64_t>& keys, uint minimizer, uint level, uint64_t value);
    std::vector<T> query_keys(const std::vector<uint64_t>& keys, uint minimizer);
    std::vector<std::pair<std::vector<uint64_t>, uint64_t>> get_super_kmers(const std::string& ref) const;

    // HIGH LEVEL FUNCTIONS
    void insert_sequence(const std::string& reference, uint level, uint64_t value);
    void insert_file(const std::string& filename, uint level, uint32_t indice_bloom);
    void insert_file_of_file(const std::string& filename);
    std::vector<uint> query_key(const uint64_t key);
    std::vector<uint32_t> query_sequence(const std::string& reference);
    void load_super_kmer(std::vector<std::vector<std::pair<uint64_t, uint32_t>>>& colored_kmer_per_bucket, uint32_t query_id, const std::string& reference);
    uint64_t query_bucket(const std::vector<std::pair<uint64_t, uint32_t>>& colored_kmer, std::vector<std::pair<std::string, std::vector<uint32_t>>>& result, uint bucket_number);
    uint64_t query_bucket(const std::vector<std::pair<uint64_t, uint32_t>>& colored_kmer, std::vector<std::pair<std::string, std::vector<uint32_t>>>& result, uint bucket_number, std::vector<bool>* bouly);
    // uint64_t query_bucket2(const std::vector<std::pair<uint64_t, uint32_t> >& colored_kmer, std::vector<std::pair<std::string, std::vector<uint32_t> > >& result, uint bucket_number);
    void insert_previous_index(const std::string& filename);
    std::vector<std::vector<uint64_t>> query(const std::string& reference) const;
    std::vector<std::vector<uint64_t>> fimpera(const std::string& reference) const;
    void optimize();
    void query_file(const std::string& reference, const std::string& output);
    void serialize() const;
    void load(const std::string& existing_index);
    void index();
    void double_index();
    void add_leaf();
    uint64_t nb_bit_set() const;
    uint64_t nb_insertions() const;

    uint64_t get_leaf_numbers() const;
};

template class Pac<uint8_t>;
template class Pac<uint16_t>;
template class Pac<uint32_t>;

// template class BestPart<uint8_t>;
// template class BestPart<uint16_t>;
// template class BestPart<uint32_t>;
#endif
