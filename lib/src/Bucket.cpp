#include <bm.h>
#include <bmserial.h>  // BitMagic
#include <bmundef.h>   // BitMagic
#include <libpac/utils.h>
#include <math.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <libpac/Bucket.hpp>
#include <libpac/CBF.hpp>
#include <vector>

template <typename T>
Bucket<T>::Bucket(const uint64_t nbCells, const uint number_hash_function, const uint k, const std::string prefix, uint64_t nbBitsPerCell, bool write) : _nbCells(nbCells), _k(k), _nbBitsPerCell(nbBitsPerCell), _number_hash_function(number_hash_function), _prefix(prefix), _write(write), _trunk(nullptr), _reverse_trunk(nullptr), _number_bit_set(0), _number_bit_set_abt(0), _disk_space_used(0), _leaf_number(0) {
    // TODO why leaf number was not initialized ?
    if (_write) {
        _out = new zstr::ofstream(_prefix, std::ios::trunc);
    } else {
        _out = nullptr;
    }
}

template <typename T>
Bucket<T>::Bucket(const Bucket<T>& that) : _nbCells(that._nbCells), _k(that._k), _nbBitsPerCell(that._nbBitsPerCell), _number_hash_function(that._number_hash_function), _prefix(that._prefix), _write(that._write), _number_bit_set(that._number_bit_set), _number_bit_set_abt(that._number_bit_set_abt), _leaf_filters(that._leaf_filters), _filter(that._filter), _disk_space_used(that._disk_space_used), _leaf_number(that._leaf_number), _out(nullptr) {
    if (that._trunk != nullptr) {
        _trunk = new AggregatedBloom<T>(*that._trunk);
    } else {
        _trunk = nullptr;
    }
    if (that._reverse_trunk != nullptr) {
        _reverse_trunk = new AggregatedBloom<T>(*that._reverse_trunk);
    } else {
        _reverse_trunk = nullptr;
    }
}

template <typename T>
Bucket<T>::~Bucket() {
    if (_trunk != nullptr) {
        delete _trunk;
    }
    if (_reverse_trunk != nullptr) {
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

template <typename T>
void Bucket<T>::construct_reverse_trunk() {
    _reverse_trunk = new AggregatedBloom<T>(_nbCells, _number_hash_function);
    uint64_t lfs(_leaf_filters.size());
    for (int i = _leaf_filters.size() - 1; i >= 0; i--) {
        insert_leaf_trunk(i, i, _reverse_trunk);
    }
}

template <typename T>
uint64_t Bucket<T>::construct_trunk() {
    // TODO lrobidou: memory leaks ?
    _trunk = new AggregatedBloom<T>(_nbCells, _number_hash_function);
    uint64_t lfs(_leaf_filters.size());
    for (uint i = 0; i < lfs; i++) {
        insert_leaf_trunk(i, lfs - i, _trunk);
    }
    return _number_bit_set;
}

template <typename T>
void Bucket<T>::insert_key(const uint64_t key, uint level, uint64_t value) {
    _leaf_filters[level]->insert_key(key, value);
}

template <typename T>
void Bucket<T>::add_leaf() {
    _leaf_filters.push_back(new CBF(_nbCells, _nbBitsPerCell, _number_hash_function));
    _leaf_number++;
}

template <typename T>
void Bucket<T>::optimize(uint i) {
    if (i < _leaf_filters.size()) {
        _leaf_filters[i]->optimize();
    }
    // TODO exception?
}

template <typename T>
void Bucket<T>::dump(uint i, bm::serializer<bm::bvector<>>& bvs) {
    _disk_space_used += _leaf_filters[i]->dump_disk(bvs, _out);
    _leaf_filters[i]->free_ram();
}

template <typename T>
void Bucket<T>::free_ram() {
    for (uint i = 0; i < _leaf_filters.size(); ++i) {
        _leaf_filters[i]->free_ram();
    }
    if (_trunk != NULL) {
        delete _trunk;
        _trunk = NULL;
    }
    if (_reverse_trunk != NULL) {
        delete _reverse_trunk;
        _reverse_trunk = NULL;
    }
}

template <typename T>
uint64_t Bucket<T>::rcb(uint64_t min) const {
    uint64_t res(0);
    uint64_t offset(1);
    offset <<= (2 * _k - 2);
    for (uint i(0); i < _k; ++i) {
        res += (3 - (min % 4)) * offset;
        min >>= 2;
        offset >>= 2;
    }
    return res;
}

template <typename T>
void Bucket<T>::insert_leaf_trunk(uint level, uint indice, AggregatedBloom<T>* EB) {
    CBF* filterOfThatLevel(_leaf_filters[level]);
    const auto& [number_bit_set_tmp, number_bit_set_abt_tmp] = EB->aggregateFilter(indice, filterOfThatLevel, _nbBitsPerCell);
    _number_bit_set_abt += number_bit_set_abt_tmp;
    _number_bit_set += number_bit_set_tmp;
}

template <typename T>
std::tuple<bool, T> Bucket<T>::query_key_min(const uint64_t key) const {
    if (_trunk != NULL) {
        T value = _trunk->check_key(key);
        return {(value != 0), _leaf_number - value};
    }
    return {true, 0};
}

template <typename T>
std::tuple<bool, T> Bucket<T>::query_key_max(const uint64_t key) const {
    if (_reverse_trunk != NULL) {
        T value = _reverse_trunk->check_key(key);
        return {(value != 0), value};  // TODO
    }
    return {true, _leaf_number - 1};  // TODO
}

template <typename T>
std::vector<std::tuple<T, uint64_t>> Bucket<T>::query_key(const uint64_t key) const {
    std::vector<std::tuple<T, uint64_t>> result;
    const auto& [found_min, min_level] = query_key_min(key);
    const auto& [found_max, max_level] = query_key_max(key);
    if (!found_min || !found_max) {
        return result;
    }

    for (uint i = min_level; i <= max_level; i++) {
        uint64_t value = _leaf_filters[i]->get_key(key);
        if (value) {
            result.push_back({i, value});
        }
    }
    return result;
}

template <typename T>
void Bucket<T>::serialize() {
    if (_out == nullptr) {
        std::cerr << "error" << std::endl;  // DEBUG exception here
        return;
    }
    if (_trunk == nullptr) {
        std::cerr << "error trunk null" << std::endl;  //
    }
    uint64_t moresize = _trunk->serialize(*_out);
    _disk_space_used += moresize;
    if (_reverse_trunk != nullptr) {
        _disk_space_used += _reverse_trunk->serialize(*_out);
    }
    _out->flush();
}

template <typename T>
void Bucket<T>::load_bf(uint64_t leaf_number) {
    zstr::ifstream in(_prefix);
    _leaf_filters.clear();
    _leaf_number = leaf_number;
    for (uint i = 0; i < leaf_number; ++i) {
        _leaf_filters.push_back(new CBF(_nbCells, _nbBitsPerCell, _number_hash_function));
    }
    for (uint i = 0; i < leaf_number; ++i) {
        uint32_t indiceBloom;
        in.read(reinterpret_cast<char*>(&indiceBloom), sizeof(indiceBloom));
        _leaf_filters[indiceBloom]->load_disk(&in);
    }
}

template <typename T>
void Bucket<T>::load(uint64_t leaf_number, bool double_index) {
    zstr::ifstream in(_prefix);
    _leaf_filters.clear();
    _leaf_number = leaf_number;
    for (uint i = 0; i < leaf_number; ++i) {
        _leaf_filters.push_back(new CBF(_nbCells, _nbBitsPerCell, _number_hash_function));
    }

    for (uint i = 0; i < _leaf_number; ++i) {
        uint32_t indiceBloom;
        in.read(reinterpret_cast<char*>(&indiceBloom), sizeof(indiceBloom));
        _leaf_filters[indiceBloom]->load_disk(&in);
    }
    //
    _trunk = new AggregatedBloom<T>(_nbCells, _number_hash_function, in);
    if (double_index) {
        _reverse_trunk = new AggregatedBloom<T>(_nbCells, _number_hash_function, in);
    }
}

template <typename T>
void Bucket<T>::load(uint64_t leaf_number, bool double_index, std::vector<bool>* BBV) {
    _leaf_filters.clear();
    // TODO memory leak ?
    CBF* leon = new CBF(_nbCells, _nbBitsPerCell, _number_hash_function);
    zstr::ifstream in(_prefix);
    bm::bvector<>::enumerator en;
    bm::bvector<>::enumerator en_end;
    for (uint i = 0; i < _leaf_number; ++i) {
        uint32_t indiceBloom;
        in.read(reinterpret_cast<char*>(&indiceBloom), sizeof(indiceBloom));
        leon->load_disk(&in);
        en = leon->get_first();
        en_end = leon->get_end();
        while (en < en_end) {
            (*BBV)[(*en) * leaf_number + indiceBloom] = true;
            ++en;
        }
        // lrobidou: why ? It will be destroyed soon anyway
        // leon->BV.clear();
    }
    delete leon;
    _trunk = new AggregatedBloom<T>(_nbCells, _number_hash_function, in);
    if (double_index) {
        _reverse_trunk = new AggregatedBloom<T>(_nbCells, _number_hash_function, in);
    }
}

// TODO prevent copy ?
template <typename T>
std::vector<CBF*> Bucket<T>::get_leaf_filters() {
    return _leaf_filters;
}

template <typename T>
uint64_t Bucket<T>::get_nbCells() {
    return _nbCells;
}

template <typename T>
uint Bucket<T>::get_number_hash_function() {
    return _number_hash_function;
}

template <typename T>
uint64_t Bucket<T>::get_number_bit_set_abt() {
    return _number_bit_set_abt;
}

template <typename T>
uint64_t Bucket<T>::get_disk_space_used() {
    return _disk_space_used;
}

template <typename T>
uint64_t Bucket<T>::get_leaf_number() {
    return _leaf_number;
}

// TODO lrobidou rename, c'est confus
template <typename T>
zstr::ofstream* Bucket<T>::get_out() {
    return _out;
}

// template <typename T>
// AggregatedBloom<T>* Bucket<T>::get_trunk() {
//     return _trunk;
// }

template <typename T>
bool Bucket<T>::operator==(const Bucket<T>& that) const {
    // first, let's check for values that are quick to read
    bool quick_to_check_values_equality =
        (this->_nbCells == that._nbCells) &&
        (this->_k == that._k) &&
        (this->_nbBitsPerCell == that._nbBitsPerCell) &&
        (this->_number_hash_function == that._number_hash_function) &&
        (this->_prefix == that._prefix) &&
        (this->_write == that._write) &&
        (this->_number_bit_set == that._number_bit_set) &&
        (this->_number_bit_set_abt == that._number_bit_set_abt) &&
        // (this->_disk_space_used == that._disk_space_used) && // TODO
        (this->_filter == that._filter) &&
        (this->_leaf_number == that._leaf_number);
    if (!quick_to_check_values_equality) {
        return false;
    }
    if (this->_leaf_filters.size() == that._leaf_filters.size()) {
        if (!(this->_leaf_filters == that._leaf_filters)) {
            bool all_equal = true;
            for (uint64_t i = 0; i < this->_leaf_filters.size(); i++) {
                if (!(*(this->_leaf_filters[i]) != *(that._leaf_filters[i])))
                    all_equal = false;
            }
            if (!all_equal) {
                return false;
            }
        }
    } else {
        return false;
    }

    // let's check reverse_trunk first, as it is more likely to be nullptr

    // checking reverse_trunk (is one of them nullptr ?)
    if ((this->_reverse_trunk == nullptr) && (that._reverse_trunk != nullptr)) {
        return false;
    } else if ((this->_reverse_trunk == nullptr) && (that._reverse_trunk != nullptr)) {
        return false;
    }

    // checking trunk (is one of them nullptr ?)
    if ((this->_trunk == nullptr) && (that._trunk != nullptr)) {
        return false;
    } else if ((this->_trunk == nullptr) && (that._trunk != nullptr)) {
        return false;
    }

    // reverse trunk
    if ((this->_reverse_trunk != nullptr) && (that._reverse_trunk != nullptr)) {
        if (*(this->_reverse_trunk) != *(that._reverse_trunk)) {
            return false;
        }
    }
    if ((this->_trunk != nullptr) && (that._trunk != nullptr)) {
        if (*(this->_trunk) != *(that._trunk)) {
            return false;
        }
    }
    return true;
}

template <typename T>
bool Bucket<T>::operator!=(const Bucket<T>& that) const {
    return (!((*this) == that));
}