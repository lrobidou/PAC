#ifndef BESTH
#define BESTH

#include <bm.h>        // BitMagic
#include <bmserial.h>  // BitMagic
#include <bmundef.h>   // BitMagic
#include <libpac/Bloom.h>
#include <libpac/ExponentialBloom.h>
#include <libpac/utils.h>
#include <math.h>

#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

template <class T>
class Bloom;

template <class T>
class Best {
   private:
    vector<Bloom<T>*> _leaf_filters;
    string _prefix;
    bool _filter;
    ExponentialBloom<T>* _trunk;
    ExponentialBloom<T>* _reverse_trunk;
    uint _K;
    uint64_t _size;
    uint _number_hash_function;
    uint64_t _number_bit_set;
    uint64_t _number_bit_set_abt;
    uint64_t _disk_space_used;
    uint64_t _leaf_number;
    zstr::ofstream* _out;
    bool _write;

   public:
    Best(const uint64_t Isize, const uint Inumber_hash_function, const uint Ik, const string Iprefix, bool Iwrite = true) {
        _K = Ik;
        _prefix = Iprefix;
        _write = Iwrite;
        if (_write) {
            _out = new zstr::ofstream(Iprefix, ios::trunc);
        }
        _size = Isize - 1;
        _number_hash_function = Inumber_hash_function;
        _trunk = NULL;
        _reverse_trunk = NULL;
        _number_bit_set_abt = _number_bit_set = _disk_space_used = 0;
    }

    ~Best() {
        if (_trunk != NULL) {
            delete _trunk;
        }
        if (_reverse_trunk != NULL) {
            delete _reverse_trunk;
        }
        for (uint i = 0; i < _leaf_filters.size(); ++i) {
            delete _leaf_filters[i];
        }
        if (write) {
            delete _out;
        } else {
        }
    }

    // LOW LEVEL FUNCTIONS
    void insert_key(const uint64_t key, uint level);
    void check_key_leaf(const uint64_t key, const uint level) const;
    uint64_t rcb(uint64_t min) const;
    void update_kmer(uint64_t& min, char nuc) const;
    void update_kmer_RC(uint64_t& min, char nuc) const;
    void insert_last_leaf_trunk(uint level, ExponentialBloom<T>* EB);
    void insert_leaf_trunk(uint indice, uint level, ExponentialBloom<T>* EB);
    void load(uint64_t leaf_number, bool double_index);

    // HIGH LEVEL FUNCTIONS
    void insert_sequence(const string& reference);
    void insert_file(const string filename);
    void insert_file_of_file(const string filename);
    vector<T> query_key(const uint64_t key);
    vector<uint> query_sequence(const string& reference);
    void construct_reverse_trunk();
    uint64_t construct_trunk();
    void get_stats() const;
    void serialize();
    void optimize();
    void optimize(uint i);
    void dump(uint i, bm::serializer<bm::bvector<> >& bvs);
    void add_leaf();
    void free_ram();
    void load_bf(uint64_t leaf_number);
    void load(uint64_t leaf_number, bool double_index, vector<bool>* BBV);
    T query_key_max(const uint64_t key);
    T query_key_min(const uint64_t key);

    // GETTERS
    vector<Bloom<T>*> get_leaf_filters();
    uint64_t get_size();
    uint get_number_hash_function();
    uint64_t get_number_bit_set_abt();
    uint64_t get_disk_space_used();
    uint64_t get_leaf_number();
    zstr::ofstream* get_out();

    // SETTERS
    void set_leaf_number(uint64_t leaf_number);
};

template class Best<uint8_t>;
template class Best<uint16_t>;
template class Best<uint32_t>;

#endif
