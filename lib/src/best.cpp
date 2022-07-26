#include <bm.h>
#include <bmserial.h>  // BitMagic
#include <bmundef.h>   // BitMagic
#include <libpac/Bloom.h>
#include <libpac/ExponentialBloom.h>
#include <libpac/best.h>
#include <libpac/utils.h>
#include <math.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

using namespace std;
using namespace filesystem;

template <class T>
void Best<T>::construct_reverse_trunk() {
    _reverse_trunk = new ExponentialBloom<T>(_size, _number_hash_function);
    uint64_t lfs(_leaf_filters.size());
    for (int i = _leaf_filters.size() - 1; i >= 0; i--) {
        insert_leaf_trunk(i, lfs - i, _reverse_trunk);
    }
}

template <class T>
uint64_t Best<T>::construct_trunk() {
    _trunk = new ExponentialBloom<T>(_size, _number_hash_function);
    uint64_t lfs(_leaf_filters.size());
    for (uint i = 0; i < lfs; i++) {
        insert_leaf_trunk(i, lfs - i, _trunk);
    }
    return _number_bit_set;
}

template <class T>
void Best<T>::insert_key(const uint64_t key, uint level) {
    _leaf_filters[level]->insert_key(key);
}

template <class T>
void Best<T>::add_leaf() {
    _leaf_filters.push_back(new Bloom<T>(_size, _number_hash_function));
}

template <class T>
void Best<T>::optimize() {
    for (uint64_t i = 0; i < _leaf_filters.size(); ++i) {
        _leaf_filters[i]->optimize();
    }
}

template <class T>
void Best<T>::optimize(uint i) {
    if (i < _leaf_filters.size()) {
        _leaf_filters[i]->optimize();
    }
    // TODO exception?
}

template <class T>
void Best<T>::dump(uint i, bm::serializer<bm::bvector<> >& bvs) {
    _disk_space_used += _leaf_filters[i]->dump_disk(bvs, _out, i);
    _leaf_filters[i]->free_ram();
}

template <class T>
void Best<T>::free_ram() {
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

template <class T>
uint64_t Best<T>::rcb(uint64_t min) const {
    uint64_t res(0);
    uint64_t offset(1);
    offset <<= (2 * _K - 2);
    for (uint i(0); i < _K; ++i) {
        res += (3 - (min % 4)) * offset;
        min >>= 2;
        offset >>= 2;
    }
    return res;
}

template <class T>
void Best<T>::insert_leaf_trunk(uint level, uint indice, ExponentialBloom<T>* EB) {
    Bloom<T>* leon(_leaf_filters[level]);
    bm::bvector<>::enumerator en = leon->BV->first();
    bm::bvector<>::enumerator en_end = leon->BV->end();
    while (en < en_end) {
        if ((*(EB->filter))[*en] == 0) {
            (*(EB->filter))[*en] = indice;
            _number_bit_set_abt++;
        }
        ++en;
        _number_bit_set++;
    }
}

template <class T>
T Best<T>::query_key_min(const uint64_t key) {
    T min_level = 0;
    if (_trunk != NULL) {
        min_level = _leaf_number - _trunk->check_key(key);
    }
    return min_level;
}

template <class T>
T Best<T>::query_key_max(const uint64_t key) {
    T max_level = _leaf_number;
    if (_reverse_trunk != NULL) {
        max_level = _leaf_number - _reverse_trunk->check_key(key);
    }
    return max_level;
}

template <class T>
vector<T> Best<T>::query_key(const uint64_t key) {
    vector<T> result;
    uint min_level = 0;
    if (_trunk != NULL) {
        min_level = _leaf_number - _trunk->check_key(key);
    }
    uint max_level = _leaf_number;
    if (_reverse_trunk != NULL) {
        max_level = _reverse_trunk->check_key(key);
    }
    for (uint i = min_level; i < max_level; i++) {
        if (_leaf_filters[i]->check_key(key)) {
            result.push_back(i);
        }
    }
    return result;
}

template <class T>
void Best<T>::serialize() {
    void* point = _trunk->filter->data();
    _out->write((char*)point, sizeof(T) * _size);
    _disk_space_used += sizeof(T) * _size;
    if (_reverse_trunk != NULL) {
        void* point = _reverse_trunk->filter->data();
        _out->write((char*)point, sizeof(T) * _size);
        _disk_space_used += sizeof(T) * _size;
    }
    _out->flush();
}

template <class T>
void Best<T>::load_bf(uint64_t leaf_number) {
    zstr::ifstream in(_prefix);
    _leaf_filters.clear();
    for (uint i = 0; i < leaf_number; ++i) {
        _leaf_filters.push_back(new Bloom<T>(this->_size, this->_number_hash_function));
    }
    for (uint i = 0; i < leaf_number; ++i) {
        uint32_t indiceBloom;
        in.read(reinterpret_cast<char*>(&indiceBloom), sizeof(indiceBloom));
        _leaf_filters[indiceBloom]->load_disk(&in);
    }
}

template <class T>
void Best<T>::load(uint64_t leaf_number, bool double_index) {
    _leaf_filters.clear();
    for (uint i = 0; i < leaf_number; ++i) {
        _leaf_filters.push_back(new Bloom<T>(_size, _number_hash_function));
    }
    zstr::ifstream in(_prefix);
    for (uint i = 0; i < _leaf_number; ++i) {
        uint32_t indiceBloom;
        in.read(reinterpret_cast<char*>(&indiceBloom), sizeof(indiceBloom));
        _leaf_filters[indiceBloom]->load_disk(&in);
    }
    _trunk = new ExponentialBloom<T>(_size, _number_hash_function);
    void* point = (_trunk->filter->data());
    in.read((char*)point, sizeof(T) * _size);
    if (double_index) {
        _reverse_trunk = new ExponentialBloom<T>(_size, _number_hash_function);
        void* point = _reverse_trunk->filter->data();
        in.read((char*)point, sizeof(T) * _size);
    }
}

template <class T>
void Best<T>::load(uint64_t leaf_number, bool double_index, vector<bool>* BBV) {
    _leaf_filters.clear();
    Bloom<T>* leon(new Bloom<T>(_size, _number_hash_function));
    zstr::ifstream in(_prefix);
    bm::bvector<>::enumerator en;
    bm::bvector<>::enumerator en_end;
    for (uint i = 0; i < _leaf_number; ++i) {
        uint32_t indiceBloom;
        in.read(reinterpret_cast<char*>(&indiceBloom), sizeof(indiceBloom));
        leon->load_disk(&in);
        en = leon->BV->first();
        en_end = leon->BV->end();
        while (en < en_end) {
            (*BBV)[(*en) * leaf_number + indiceBloom] = true;
            ++en;
        }
        leon->BV->clear();
    }
    delete leon;
    _trunk = new ExponentialBloom<T>(_size, _number_hash_function);
    void* point = (_trunk->filter->data());
    in.read((char*)point, sizeof(T) * _size);
    if (double_index) {
        _reverse_trunk = new ExponentialBloom<T>(_size, _number_hash_function);
        void* point = _reverse_trunk->filter->data();
        in.read((char*)point, sizeof(T) * _size);
    }
}

// TODO prevent copy ?
template <typename T>
vector<Bloom<T>*> Best<T>::get_leaf_filters() {
    return _leaf_filters;
}

template <typename T>
uint64_t Best<T>::get_size() {
    return _size;
}

template <typename T>
uint Best<T>::get_number_hash_function() {
    return _number_hash_function;
}

template <typename T>
uint64_t Best<T>::get_number_bit_set_abt() {
    return _number_bit_set_abt;
}

template <typename T>
uint64_t Best<T>::get_disk_space_used() {
    return _disk_space_used;
}

template <typename T>
uint64_t Best<T>::get_leaf_number() {
    return _leaf_number;
}

template <typename T>
zstr::ofstream* Best<T>::get_out() {
    return _out;
}

template <typename T>
void Best<T>::set_leaf_number(uint64_t leaf_number) {
    _leaf_number = leaf_number;
};