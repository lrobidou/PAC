#include <bm.h>        // BitMagic
#include <bmserial.h>  // BitMagic
#include <bmundef.h>   // BitMagic
#include <libpac/Bloom.h>
#include <libpac/utils.h>

#include <filesystem>

Bloom::Bloom(uint64_t size, uint number_hash_function) : _size(size), _number_hash_function(number_hash_function), BV(_size + 1, bm::BM_GAP) {
}

void Bloom::insert_key(uint64_t key) {
    for (uint64_t i = 0; i < _number_hash_function; ++i) {
        uint64_t h = hash_family(key, i, _size);
        BV[h] = true;
    }
}

bool Bloom::check_key(uint64_t key) {
    for (uint64_t i = 0; i < _number_hash_function; ++i) {
        uint64_t h = hash_family(key, i, _size);
        if (BV[h] == false) {
            return false;
        }
    }
    return true;
}

const bm::bvector<>::enumerator Bloom::get_first() const {
    return BV.first();
}

const bm::bvector<>::enumerator Bloom::get_end() const {
    return BV.end();
}

uint64_t Bloom::dump_disk(bm::serializer<bm::bvector<> >& bvs, zstr::ofstream* out) {
    bm::serializer<bm::bvector<> >::buffer sbuf;
    unsigned char* buf = 0;
    bvs.serialize(BV, sbuf);
    buf = sbuf.data();
    uint64_t sz = sbuf.size();
    auto point2 = &buf[0];
    out->write(reinterpret_cast<const char*>(&sz), sizeof(sz));
    out->write((char*)point2, sz);
    out->flush();
    return sz;
}

void Bloom::load_disk(zstr::ifstream* in) {
    BV = bm::bvector<>(_size + 1, bm::BM_GAP);
    uint64_t sz;
    in->read(reinterpret_cast<char*>(&sz), sizeof(sz));
    uint8_t* buff = new uint8_t[sz];
    in->read((char*)buff, sz);
    bm::deserialize(BV, buff);
    delete[] buff;
}

void Bloom::optimize() {
    BV.optimize(NULL, bm::bvector<>::opt_compress);
}

void Bloom::free_ram() {
    BV.reset();
}

uint64_t Bloom::get_filter_size() {
    return BV.size();
}

bool Bloom::operator==(const Bloom& that) const {
    return (
        (this->_size == that._size) &&
        (this->_number_hash_function == that._number_hash_function) &&
        (this->BV == that.BV));
}

bool Bloom::operator!=(const Bloom& that) const {
    return (!((*this) == that));
}