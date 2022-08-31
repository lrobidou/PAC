#include <libpac/utils.h>
#include <math.h>

#include <algorithm>  // std::min_element
#include <libpac/CBF.hpp>

CBF::CBF(uint64_t nbCells, uint64_t nbBitsPerCell, uint64_t number_hash_function) : _bits(nbCells * nbBitsPerCell, bm::BM_GAP), _nbCells(nbCells), _nbBitsPerCell(nbBitsPerCell), _limitValueInBucket(pow(2, _nbBitsPerCell) - 1), _mask(pow(2, _nbBitsPerCell - 1)), _number_hash_function(number_hash_function) {
}

uint64_t CBF::get_from_hash(const uint64_t& hash) const {
    uint64_t valueInFilter = 0;
    for (uint64_t i = 0; i < _nbBitsPerCell; i++) {
        valueInFilter <<= 1;
        valueInFilter += _bits[hash * _nbBitsPerCell + i];
    }
    return valueInFilter;
}

uint64_t CBF::insert_key(uint64_t key, uint64_t value) {
    // sorry, no value smaller than 0 accepted
    if (value <= 0) {
        return 0;
    }
    // can't set values higher because of space limit
    if (value > _limitValueInBucket) {
        value = _limitValueInBucket;
    }

    // value we would found if we check the filter after insertion
    uint64_t valueInserted = _limitValueInBucket;

    for (uint64_t i = 0; i < _number_hash_function; i++) {
        uint64_t hash = hash_family(key, i, _nbCells);
        uint64_t valueInFilterForThatHash = get_from_hash(hash);
        if (valueInFilterForThatHash >= value) {
            valueInserted = std::min(valueInFilterForThatHash, value);
            // inserting "value" in the filter won't change the filter's state
            // no need to continue ths iteration, let's skip to the next one
            continue;
        }

        // from here: value > valueInFilterForThatHash

        // save data
        uint64_t mask = _mask;
        // we save each bit one by one
        uint64_t start = hash * _nbBitsPerCell;
        for (uint64_t j = 0; j < _nbBitsPerCell; j++) {
            _bits[start + j] = value & mask;
            mask >>= 1;
        }
    }

    return valueInserted;
}

uint64_t CBF::get_key(const uint64_t& key) const {
    std::vector<uint64_t> valuesInFilter(_number_hash_function, 0);
    for (uint64_t i = 0; i < _number_hash_function; ++i) {
        // TODO lrobidou: why & and not %
        uint64_t hash = hash_family(key, i, _nbCells);
        valuesInFilter[i] = get_from_hash(hash);
    }
    return *std::min_element(std::begin(valuesInFilter), std::end(valuesInFilter));
}

const bm::bvector<>::enumerator CBF::get_first() const {
    return _bits.first();
}

const bm::bvector<>::enumerator CBF::get_end() const {
    return _bits.end();
}

void CBF::optimize() {
    _bits.optimize(NULL, bm::bvector<>::opt_compress);
}

uint64_t CBF::dump_disk(bm::serializer<bm::bvector<> >& bvs, zstr::ofstream* out) {
    // TODO
}

void CBF::load_disk(zstr::ifstream* in)  // TODO: use a new constructor instead ?
{
    // TODO
}

void CBF::free_ram() {
    _bits.reset();
}

bool CBF::operator==(const CBF& that) const {
    return (
        (this->_nbBitsPerCell == that._nbBitsPerCell) &&
        (this->_nbCells == that._nbCells) &&
        (this->_number_hash_function == that._number_hash_function) &&
        (this->_bits == that._bits));
}

bool CBF::operator!=(const CBF& that) const {
    return (!((*this) == that));
}

// getter

uint64_t get_filter_size() {
    // TODO
}