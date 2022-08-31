#ifndef BUCKETH
#define BUCKETH

#include <bm.h>        // BitMagic
#include <bmserial.h>  // BitMagic
#include <bmundef.h>   // BitMagic
#include <libpac/AggregatedBloom.h>
#include <libpac/utils.h>
#include <math.h>

#include <cmath>
#include <iostream>
#include <libpac/CBF.hpp>
#include <vector>

template <typename T>
class Bucket {
   private:
    uint64_t _nbCells;  // TODO use nbCells
    uint _k;
    uint64_t _nbBitsPerCell;  // TODO initialize
    uint _number_hash_function;
    std::string _prefix;
    bool _write;
    uint64_t _number_bit_set;
    uint64_t _number_bit_set_abt;

    std::vector<CBF*> _leaf_filters;

    bool _filter;
    AggregatedBloom<T>* _trunk;
    AggregatedBloom<T>* _reverse_trunk;

    uint64_t _disk_space_used;
    uint64_t _leaf_number;
    zstr::ofstream* _out;

   public:
    Bucket(const uint64_t nbCells, const uint number_hash_function, const uint k, const std::string prefix, uint64_t nbBitsPerCell, bool write = true);
    Bucket(const Bucket<T>& that);
    ~Bucket();

   private:
    // LOW LEVEL FUNCTIONS

    void check_key_leaf(const uint64_t key, const uint level) const;
    uint64_t rcb(uint64_t min) const;
    void update_kmer(uint64_t& min, char nuc) const;
    void update_kmer_RC(uint64_t& min, char nuc) const;
    void insert_last_leaf_trunk(uint level, AggregatedBloom<T>* EB);
    void insert_leaf_trunk(uint indice, uint level, AggregatedBloom<T>* EB);

    // HIGH LEVEL FUNCTIONS
    // void insert_sequence(const std::string& reference);
    void insert_file(const std::string filename);

   public:
    // LOW LEVEL FUNCTIONS
    // insert value for given key in filter of corresponding level
    void insert_key(const uint64_t key, uint level, uint64_t value);
    void load(uint64_t leaf_number, bool double_index);

    // HIGH LEVEL FUNCTIONS
    std::vector<std::tuple<T, uint64_t>> query_key(const uint64_t key) const;
    void construct_reverse_trunk();
    uint64_t construct_trunk();
    void serialize();
    void optimize(uint i);
    void dump(uint i, bm::serializer<bm::bvector<>>& bvs);
    void add_leaf();
    void free_ram();
    void load_bf(uint64_t leaf_number);
    void load(uint64_t leaf_number, bool double_index, std::vector<bool>* BBV);
    std::tuple<bool, T> query_key_max(const uint64_t key) const;
    std::tuple<bool, T> query_key_min(const uint64_t key) const;

    // GETTERS
    std::vector<CBF*> get_leaf_filters();  // TODO remove ?
    // AggregatedBloom<T>* get_trunk();       // TODO remove
    uint64_t get_nbCells();
    uint get_number_hash_function();
    uint64_t get_number_bit_set_abt();
    uint64_t get_disk_space_used();
    uint64_t get_leaf_number();
    zstr::ofstream* get_out();
    bool operator==(const Bucket<T>& that) const;
    bool operator!=(const Bucket<T>& that) const;
};

template class Bucket<uint8_t>;
template class Bucket<uint16_t>;
template class Bucket<uint32_t>;

#endif
