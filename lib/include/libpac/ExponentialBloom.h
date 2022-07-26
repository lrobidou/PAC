#ifndef EXPONENTIALBLOOMH
#define EXPONENTIALBLOOMH

#include <math.h>
#include <omp.h>

#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

template <class T>
class ExponentialBloom {
   public:
    uint64_t size;
    uint number_hash_functions;
    vector<T>* filter;
    void insert_key(uint64_t key, uint level);
    T check_key(uint64_t key) const;
    bool check_key(uint64_t key, uint level) const;

    ExponentialBloom(const uint64_t Isize, const uint Inumber_hash_functions) {
        size = Isize;
        number_hash_functions = Inumber_hash_functions;
        filter = new vector<T>(size + 1, 0);
    }

    ~ExponentialBloom() { delete filter; }
};

template class ExponentialBloom<uint8_t>;
template class ExponentialBloom<uint16_t>;
template class ExponentialBloom<uint32_t>;

#endif
