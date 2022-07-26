#include <bm.h>        // BitMagic
#include <bmserial.h>  // BitMagic
#include <bmundef.h>   // BitMagic
#include <libpac/Bloom.h>
#include <libpac/utils.h>

#include <filesystem>

using namespace std;
using namespace filesystem;

template <class T>
Bloom<T>::Bloom(uint64_t size, uint number_hash_function) : _size(size), _number_hash_function(number_hash_function) {
    BV = new bm::bvector<>(_size + 1, bm::BM_GAP);
}

template <class T>
void Bloom<T>::insert_key(uint64_t key) {
    for (uint64_t i = 0; i < _number_hash_function; ++i) {
        uint64_t h = (hash_family(key, i)) & (_size);
        (*BV)[h] = true;
    }
}

template <class T>
bool Bloom<T>::check_key(uint64_t key) {
    for (uint64_t i = 0; i < _number_hash_function; ++i) {
        uint64_t h = hash_family(key, i) & _size;
        if ((*BV)[h] == false) {
            return false;
        }
    }
    return true;
}

template <class T>
uint64_t Bloom<T>::dump_disk(bm::serializer<bm::bvector<> >& bvs, zstr::ofstream* out, uint32_t i) {
    bm::serializer<bm::bvector<> >::buffer sbuf;
    unsigned char* buf = 0;
    bvs.serialize(*(BV), sbuf);
    buf = sbuf.data();
    uint64_t sz = sbuf.size();
    auto point2 = &buf[0];
    out->write(reinterpret_cast<const char*>(&sz), sizeof(sz));
    out->write((char*)point2, sz);
    out->flush();
    return sz;
}

template <class T>
void Bloom<T>::load_disk(zstr::ifstream* in) {
    if (BV == NULL) {
        BV = new bm::bvector<>(_size + 1, bm::BM_GAP);
    }
    uint64_t sz;
    in->read(reinterpret_cast<char*>(&sz), sizeof(sz));
    uint8_t* buff = new uint8_t[sz];
    in->read((char*)buff, sz);
    bm::deserialize(*(BV), buff);
    delete[] buff;
}

template <class T>
void Bloom<T>::free_ram() {
    BV->reset();
}
