#pragma once

#include <bm.h>        // BitMagic
#include <bmserial.h>  // BitMagic
#include <bmundef.h>   // BitMagic

#include <zstr.hpp>

class CBF {
   private:
    // number of "cell" int he filter (each cell containing an integer)
    uint64_t _nbCells;
    // size of each cell in bit
    uint64_t _nbBitsPerCell;
    // max value in a bucket, "cached" to prevent recomputing it
    uint64_t _limitValueInBucket;
    // TODO
    uint64_t _mask;
    // number of hash function to apply on key
    uint64_t _number_hash_function;
    // retrieve value from a hash
    uint64_t get_from_hash(const uint64_t& hash) const;

   public:
    // bit vector (data, i.e. the filter itself)
    bm::bvector<> _bits;  // TODO move
    // constructor

    // create a new counting Bloom filter
    CBF(uint64_t nbCells, uint64_t nbBitsPerCell, uint64_t number_hash_function);

    // methods

    // insert a value in the filter given a key
    uint64_t insert_key(uint64_t key, uint64_t value);
    // retrieve the value stored in the filter given a key
    uint64_t get_key(const uint64_t& index) const;
    // get an enumerator pointing on the first non-zero bit
    const bm::bvector<>::enumerator get_first() const;
    // get an enumerator pointing on the next bit after the last
    const bm::bvector<>::enumerator get_end() const;
    // TODO doc
    void optimize();
    // save the counting Bloom filter on disk
    uint64_t dump_disk(bm::serializer<bm::bvector<> >& bvs, zstr::ofstream* out);
    // load the Bloom filter from disk
    void load_disk(zstr::ifstream* in);  // TODO: use a new constructor instead ?
    // release RAM
    void free_ram();  // TODO: does not work

    // operator

    // check equality between two counting Bloom filter
    bool operator==(const CBF& that) const;
    // check inequality between two counting Bloom filter
    bool operator!=(const CBF& that) const;

    // getter

    // get the size of the inner counting Bloom filter
    uint64_t get_filter_size();
};
