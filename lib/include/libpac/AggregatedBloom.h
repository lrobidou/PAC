#ifndef AGGREGATEDBLOOMH
#define AGGREGATEDBLOOMH

#include <bm.h>
#include <bmserial.h>  // BitMagic
#include <bmundef.h>   // BitMagic
#include <math.h>
#include <omp.h>

#include <cmath>
#include <iostream>
#include <libpac/CBF.hpp>  // TODO remove
#include <vector>
#include <zstr.hpp>

template <class T>
class AggregatedBloom {
   private:
    uint64_t _size;
    uint _number_hash_functions;
    std::vector<T> _filter;  // TODO change

   public:
    // void insert_key(uint64_t key, uint level);
    T check_key(uint64_t key) const;
    uint64_t serialize(zstr::ofstream& out) const;
    // bool is_empty(uint64_t key) const;

    // "stack" a filter
    std::tuple<uint64_t, uint64_t> aggregateFilter(uint indice, const CBF* const filterOfThatLevel, const uint64_t nbBitsPerCell = 1);

    AggregatedBloom(const uint64_t size, const uint number_hash_functions);
    AggregatedBloom(const uint64_t size, const uint number_hash_functions, zstr::ifstream& in);
    AggregatedBloom(const AggregatedBloom<T>& that);
    bool operator==(const AggregatedBloom<T>& that) const;
    bool operator!=(const AggregatedBloom<T>& that) const;
};

template class AggregatedBloom<uint8_t>;
template class AggregatedBloom<uint16_t>;
template class AggregatedBloom<uint32_t>;

#endif
