#include <libpac/AggregatedBloom.h>
#include <libpac/Bloom.h>  //TODO remove
#include <libpac/utils.h>

// TODO change Isize to size
// TODO why +1
template <class T>
AggregatedBloom<T>::AggregatedBloom(const uint64_t size, const uint number_hash_functions) : _size(size), _number_hash_functions(number_hash_functions), _filter(size + 1, 0) {
}

template <class T>
AggregatedBloom<T>::AggregatedBloom(const uint64_t size, const uint number_hash_functions, zstr::ifstream& in) : _size(size), _number_hash_functions(number_hash_functions), _filter(size + 1, 0) {
    void* point = _filter.data();
    in.read((char*)point, sizeof(T) * _size);
}

template <class T>
AggregatedBloom<T>::AggregatedBloom(const AggregatedBloom<T>& that) : _size(that._size), _number_hash_functions(that._number_hash_functions), _filter(that._filter) {
}

// template <class T>
// bool AggregatedBloom<T>::is_empty(uint64_t key) const {
//     for (uint64_t i = 0; i < number_hash_functions; ++i) {
//         uint64_t h = hash_family(key, i) & size;
//         if (filter[h] == 0) {
//             return false;
//         }
//     }
//     return true;
// }

template <class T>
std::tuple<uint64_t, uint64_t> AggregatedBloom<T>::aggregateFilter(uint indice, const CBF* const filterOfThatLevel, const uint64_t nbBitsPerCell) {
    const bm::bvector<>::enumerator& startFilter = filterOfThatLevel->get_first();
    const bm::bvector<>::enumerator& endFilter = filterOfThatLevel->get_end();
    uint64_t number_bit_set_abt = 0;
    uint64_t number_bit_set = 0;
    uint64_t last_bit_set = 0;
    bool first = true;
    for (auto it = startFilter; it != endFilter; it++) {
        bm::bvector<>::size_type value = *it / nbBitsPerCell;  // TODO check
        if (first || (last_bit_set != value)) {
            last_bit_set = value;

            if (_filter[value] == 0) {
                _filter[value] = indice;  // indice = nb filter having a one
                number_bit_set_abt++;
            }
            number_bit_set++;
        }
        if (first) {
            first = false;
        }
    }
    return {number_bit_set, number_bit_set_abt};
}

template <class T>
T AggregatedBloom<T>::check_key(uint64_t key) const {
    T result = -1;  // warning: not actually -1: this template accepts only unsigned, so it is the max that can ba stored in this type
    for (uint64_t i = 0; i < _number_hash_functions; ++i) {
        uint64_t h = hash_family(key, i, _size);
        double value = _filter[h];
        if (value < result) {
            result = value;
        }
    }
    return result;
}

// TODO ifthere is a serialization bug, maybe it is here (got rid of size parameter)
template <class T>
uint64_t AggregatedBloom<T>::serialize(zstr::ofstream& out) const {
    const void* point = _filter.data();
    out.write((char*)point, sizeof(T) * _size);
    return sizeof(T) * _size;
}

template <typename T>
bool AggregatedBloom<T>::operator==(const AggregatedBloom<T>& that) const {
    return (
        (this->_size == that._size) &&
        (this->_number_hash_functions == that._number_hash_functions) &&
        (this->_filter == that._filter));
}

template <typename T>
bool AggregatedBloom<T>::operator!=(const AggregatedBloom<T>& that) const {
    return (!((*this) == that));
}