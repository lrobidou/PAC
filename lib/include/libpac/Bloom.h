#ifndef BLOOMH
#define BLOOMH
#include <bm.h>        // BitMagic
#include <bmserial.h>  // BitMagic
#include <bmundef.h>   // BitMagic
#include <libpac/utils.h>
#include <math.h>
#include <omp.h>

#include <cmath>
#include <iostream>
#include <vector>

class Bloom {
   private:
    uint64_t _size;
    uint _number_hash_function;

   public:
    bm::bvector<> BV;
    void insert_key(uint64_t key);
    bool check_key(uint64_t key);
    const bm::bvector<>::enumerator get_first() const;
    const bm::bvector<>::enumerator get_end() const;
    void optimize();
    uint64_t get_cardinality() const;
    uint64_t dump_disk(bm::serializer<bm::bvector<> >& bvs, zstr::ofstream* out);
    void load_disk(zstr::ifstream* in);
    void free_ram();

    // constructor
    Bloom(uint64_t size, uint number_hash_function);

    // operator
    bool operator==(const Bloom& that) const;
    bool operator!=(const Bloom& that) const;

    // getter
    uint64_t get_filter_size();
};

#endif
