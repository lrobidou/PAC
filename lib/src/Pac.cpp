#include <libpac/Bloom.h>
#include <libpac/utils.h>
#include <math.h>

#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fimpera-lib/finderec.hpp>
#include <iostream>
#include <libpac/Bucket.hpp>
#include <libpac/Pac.hpp>
#include <libpac/fimpera.hpp>
#include <libpac/generators/KMCReader.hpp>
#include <string>
#include <vector>

// shift x n times
uint shift(uint x, uint n) {
    return x << n;
}

template <typename T>
Pac<T>::Pac(const uint64_t nbCells, const uint64_t nbBitsPerCell, const uint number_hash_function, const uint k, const uint z, const std::string wdir, bool use_double_index, uint bucketing, uint core_number) : _nbCells(nbCells), _nbBitsPerCell(nbBitsPerCell), _k(k), _s(k - z), _number_hash_function(number_hash_function), _small_minimizer_size(bucketing), _large_minimizer_size(_small_minimizer_size + 3), _large_minimizer_number(shift(1, 2 * _large_minimizer_size)), _bucket_number(shift(1, 2 * _small_minimizer_size)), _buckets(_bucket_number, nullptr), _core_number(core_number), _use_double_index(use_double_index), _w_dir(wdir), _leaf_number(0), _number_bit_set(0), _number_bit_set_abt(0), _disk_space_used(0) {
    createPathIfNotExists(_w_dir);

    _mutex_array = new omp_lock_t[_bucket_number];
    for (uint32_t i = 0; i < _bucket_number; ++i) {
        _buckets[i] = new Bucket<T>(_nbCells / _bucket_number, _number_hash_function, _s, _w_dir + "/P" + std::to_string(i), _nbBitsPerCell);
        omp_init_lock(&_mutex_array[i]);
    }
    // TODO smer
    // TODO initialization list
    _offsetUpdatekmer = 1;
    _offsetUpdatekmer <<= 2 * _s;
    _offsetUpdateminimizer = 1;
    _offsetUpdateminimizer <<= (2 * _large_minimizer_size);
}

template <typename T>
Pac<T>::~Pac() {
    for (uint32_t i = 0; i < _bucket_number; ++i) {
        if (_buckets[i] != nullptr) {
            delete _buckets[i];
        }
    }
    delete[] _mutex_array;
}

template <typename T>
void Pac<T>::insert_file(const std::string& filename, uint level, uint32_t indice_bloom) {
    for (const auto& [kmer, abundance] : libpac::generators::KMCReader(filename)) {
        insert_sequence(kmer, level, abundance);
    }

    // TODO lrobidou: I broke it
    // bm::serializer<bm::bvector<>> bvs;
    // bvs.byte_order_serialization(false);
    // bvs.gap_length_serialization(false);
    // for (uint i = 0; i < _buckets.size(); i++) {
    //     _buckets[i]->optimize(level);
    //     omp_set_lock(&_mutex_array[i]);  // TODO THE LOCK SHOULD NOT INCLUDE THE SERIALIZATION
    //     _buckets[i]->get_out()->write(reinterpret_cast<const char*>(&indice_bloom), sizeof(indice_bloom));
    //     _buckets[i]->dump(level, bvs);
    //     omp_unset_lock(&_mutex_array[i]);
    // }
}

template <typename T>
void Pac<T>::insert_file_of_file(const std::string& filename) {
    std::cout << "I index " << _s << "-mers with counting Bloom filters having " << intToString(_nbCells) << " cells each having " << _nbBitsPerCell << " with " << _number_hash_function << " hash functions using " << intToString(1 << (2 * _small_minimizer_size)) << " partitions." << std::endl;
    std::filesystem::path absolute_input_path(std::filesystem::absolute(filename));
    std::filesystem::path prefix = absolute_input_path.parent_path();
    zstr::ifstream in(absolute_input_path);
    std::cout << "Insert file of file " << absolute_input_path << std::endl;

    // DEBUG what if core_number < nb of file ??
    for (uint i(0); i < _core_number; ++i) {
        add_leaf();
    }

    auto start = std::chrono::system_clock::now();
    ;
#pragma omp parallel num_threads(core_number)
    {
        std::string filename;
        uint level;
        bool file_valid;
        int idt = omp_get_thread_num();
        while (not in.eof()) {
#pragma omp critical(inputfile)
            {
                std::getline(in, filename);
                filename = prefix / filename;
                if (exists_test(filename)) {
                    level = _leaf_number;
                    _leaf_number++;
                    file_valid = true;
                } else {
                    file_valid = false;
                }
            }
            if (file_valid) {
                insert_file(filename, level, idt);
                if (level % 100 == 0) {
                    std::cout << level << " files" << std::endl;
                }
            }
        }
    }
    if (_leaf_number == 0) {
        std::cerr << "No valid files found in " << absolute_input_path << "." << std::endl;
        return;
    }
    auto middle = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = middle - start;
    std::cout << "Bloom construction time: " << elapsed_seconds.count() << "s\n";
    std::cout << intToString(getMemorySelfMaxUsed()) << " KB total" << std::endl;
    std::cout << intToString(getMemorySelfMaxUsed() / (_leaf_number)) << " KB per file (" << _leaf_number << " files)" << std::endl;

    index();

    nb_bit_set();
    auto end = std::chrono::system_clock::now();
    elapsed_seconds = end - middle;
    std::cout << "Exponential Bloom construction time: " << elapsed_seconds.count() << "s\n";
    elapsed_seconds = end - start;
    std::cout << "Total Index time: " << elapsed_seconds.count() << "s\n";
    std::cout << "Index time per file: " << elapsed_seconds.count() / _leaf_number << "s\n";
    std::cout << intToString(getMemorySelfMaxUsed()) << " kB total" << std::endl;
    std::cout << intToString(getMemorySelfMaxUsed() / (_leaf_number)) << " kB per file" << std::endl;
}

template <typename T>
void Pac<T>::insert_keys(const std::vector<uint64_t>& keys, uint minimizer, uint level, uint64_t value) {
    for (uint i = 0; i < keys.size(); i++) {
        _buckets[minimizer]->insert_key(keys[i], level, value);
    }
}

template <typename T>
void Pac<T>::insert_sequence(const std::string& reference, uint level, uint64_t value) {
    for (const auto& [superkmers, minimizer] : get_super_kmers(reference)) {
        insert_keys(superkmers, minimizer, level, value);
    }
}

template <typename T>
void Pac<T>::add_leaf() {
    //~ #pragma omp parallel for
    for (uint i = 0; i < _buckets.size(); ++i) {
        _buckets[i]->add_leaf();
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// lrobidou: why not returning a tuple ? AFAIU position is not a parameter
template <typename T>
uint64_t Pac<T>::regular_minimizer_pos(uint64_t seq, uint64_t& position) const {
    uint64_t mini, mmer;
    mmer = seq % _large_minimizer_number;
    mini = mmer = canonize(mmer, _large_minimizer_size);
    uint64_t hash_mini = (unrevhash(mmer));
    position = 0;
    for (uint64_t i(1); i <= _s - _large_minimizer_size; i++) {
        seq >>= 2;
        mmer = seq % _large_minimizer_number;
        mmer = canonize(mmer, _large_minimizer_size);
        uint64_t hash = (unrevhash(mmer));
        if (hash_mini > hash) {
            position = _s - _large_minimizer_size - i;
            mini = mmer;
            hash_mini = hash;
        }
    }
    return mini;
}

template <typename T>
void Pac<T>::updateK(uint64_t& min, char nuc) const {
    min <<= 2;
    min += nuc2int(nuc);
    min %= _offsetUpdatekmer;
}

template <typename T>
void Pac<T>::updateM(uint64_t& min, char nuc) const {
    min <<= 2;
    min += nuc2int(nuc);
    min %= _offsetUpdateminimizer;
}

template <typename T>
void Pac<T>::updateRCK(uint64_t& min, char nuc) const {
    min >>= 2;
    min += (nuc2intrc(nuc) << (2 * _s - 2));
}

template <typename T>
void Pac<T>::updateRCM(uint64_t& min, char nuc) const {
    min >>= 2;
    min += (nuc2intrc(nuc) << (2 * _large_minimizer_size - 2));
}

template <typename T>
std::vector<std::pair<std::vector<uint64_t>, uint64_t>> Pac<T>::get_super_kmers(const std::string& ref) const {
    std::vector<std::pair<std::vector<uint64_t>, uint64_t>> result;
    uint64_t old_minimizer, minimizer;
    std::vector<uint64_t> superkmer;
    old_minimizer = minimizer = _large_minimizer_number;
    uint64_t last_position(0);
    // FOREACH KMER
    uint64_t seq(str2num(ref.substr(0, _s)));
    uint64_t rcseq = rcb(seq, _s);
    uint64_t canon = std::min(seq, rcseq);
    uint64_t position_min;
    uint64_t min_seq = (str2num(ref.substr(_s - _large_minimizer_size, _large_minimizer_size)));
    uint64_t min_rcseq(rcb(min_seq, _large_minimizer_size));
    uint64_t min_canon(std::min(min_seq, min_rcseq));
    minimizer = regular_minimizer_pos(seq, position_min);
    old_minimizer = minimizer;
    uint64_t hash_min = unrevhash(minimizer);
    for (uint64_t i = 0; i + _s < ref.size(); ++i) {  // TODO : _s ou _k ?
        updateK(seq, ref[i + _s]);
        updateRCK(rcseq, ref[i + _s]);
        canon = std::min(seq, rcseq);
        updateM(min_seq, ref[i + _s]);
        updateRCM(min_rcseq, ref[i + _s]);
        min_canon = (std::min(min_seq, min_rcseq));
        uint64_t new_h = unrevhash(min_canon);
        // THE NEW mmer is a MINIMIZER
        if (new_h < hash_min) {
            minimizer = (min_canon);
            hash_min = new_h;
            position_min = i + _s - _large_minimizer_size + 1;
        } else {
            // the previous minimizer is outdated
            if (i >= position_min) {
                minimizer = regular_minimizer_pos(seq, position_min);
                hash_min = unrevhash(minimizer);
                position_min += (i + 1);
            } else {
            }
        }
        // COMPUTE KMER MINIMIZER
        if (revhash(old_minimizer) % _bucket_number != revhash(minimizer) % _bucket_number) {
            old_minimizer = (revhash(old_minimizer) % _bucket_number);
            result.push_back({superkmer, old_minimizer});
            superkmer.clear();
            last_position = i + 1;
            old_minimizer = minimizer;
        }
        superkmer.push_back(canon);
    }
    if (ref.size() - last_position > _s - 1) {
        old_minimizer = (revhash(old_minimizer) % _bucket_number);
        result.push_back({superkmer, old_minimizer});
        superkmer.clear();
    }
    return result;
}

template <typename T>
void Pac<T>::query_file(const std::string& filename, const std::string& output) {
    std::cout << "Query file " << filename << " in " << output << std::endl;
    auto start = std::chrono::system_clock::now();
    zstr::ofstream out(output, std::ios::trunc);
    char type = get_data_type(filename);
    zstr::ifstream in(filename);
    std::filesystem::path old(std::filesystem::current_path());
    std::vector<uint32_t> colors_count(_leaf_number, 0);

    std::filesystem::current_path(_w_dir);
    std::vector<std::pair<std::string, std::vector<uint32_t>>> result;
    std::vector<std::vector<std::pair<uint64_t, uint32_t>>> colored_kmer_per_bucket(_bucket_number);
#pragma omp parallel num_threads(core_number)
    {
        std::string ref, header;
        uint32_t query_id;
        while (not in.eof()) {
#pragma omp critical(inputfile)
            {
                std::tie(ref, header) = Biogetline(&in, type, _s);
                if (ref.size() > _s) {  // TODO k ou s ?
                    result.push_back({header, colors_count});
                    query_id = result.size() - 1;
                }
            }
            if (ref.size() > _s) {
                load_super_kmer(colored_kmer_per_bucket, query_id, ref);
            }
        }
    }
    uint64_t sum_query = 0;
#pragma omp parallel num_threads(core_number)
    {
        std::vector<bool>* BBV = new std::vector<bool>(_leaf_number * _nbCells / _bucket_number, false);
#pragma omp for
        for (uint i = 0; i < _bucket_number; i++) {
            if (not(colored_kmer_per_bucket)[i].empty()) {
                uint64_t sq(query_bucket(colored_kmer_per_bucket[i], result, i, BBV));
#pragma omp atomic
                sum_query += sq;
            }
        }
        delete BBV;
    }
    std::cout << "query done" << std::endl;
    for (uint i = 0; i < result.size(); i++) {
        {
            out << result[i].first << "\n";
            for (uint j = 0; j < result[i].second.size(); j++) {
                out << result[i].second[j] << ' ';
            }
            out << "\n";
        }
    }
    out.close();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Query time: " << elapsed_seconds.count() << "s for " << intToString(sum_query) << " kmer query or " << intToString((double)sum_query / elapsed_seconds.count()) << " query per second  \n";
    current_path(old);
}

template <typename T>
void Pac<T>::insert_previous_index(const std::string& filename) {
    zstr::ifstream in(filename + "/MainIndex");
    uint64_t osef, leaf_number_to_insert;
    // TODO verifier ordre
    in.read(reinterpret_cast<char*>(&osef), sizeof(_s));
    in.read(reinterpret_cast<char*>(&osef), sizeof(_nbCells));
    in.read(reinterpret_cast<char*>(&osef), sizeof(_number_hash_function));
    in.read(reinterpret_cast<char*>(&osef), sizeof(_small_minimizer_size));
    in.read(reinterpret_cast<char*>(&osef), sizeof(_large_minimizer_size));
    in.read(reinterpret_cast<char*>(&leaf_number_to_insert), sizeof(leaf_number_to_insert));
    for (uint i = 0; i < _buckets.size(); ++i) {
        zstr::ifstream in(filename + "/P" + std::to_string(i));
        for (uint ii = 0; ii < leaf_number_to_insert; ++ii) {
            uint32_t ibloom;
            in.read(reinterpret_cast<char*>(&ibloom), sizeof(ibloom));
            uint64_t sz;
            in.read(reinterpret_cast<char*>(&sz), sizeof(sz));
            uint8_t* buff = new uint8_t[sz];
            in.read((char*)buff, sz);
            _buckets[i]->get_out()->write(reinterpret_cast<const char*>(&ibloom), sizeof(ibloom));
            _buckets[i]->get_out()->write(reinterpret_cast<const char*>(&sz), sizeof(sz));
            _buckets[i]->get_out()->write((char*)&buff[0], sz);
        }
        _buckets[i]->get_out()->flush();
    }
    _leaf_number += leaf_number_to_insert;
    std::cout << _leaf_number << " datasets loaded from previous index" << std::endl;
}

template <typename T>
void Pac<T>::load_super_kmer(std::vector<std::vector<std::pair<uint64_t, uint32_t>>>& colored_kmer_per_bucket, uint32_t query_id, const std::string& reference) {
    for (const auto& [superkmer, minimizer] : get_super_kmers(reference)) {
        omp_set_lock(&_mutex_array[minimizer]);
        for (const auto& kmer : superkmer) {
            colored_kmer_per_bucket[minimizer].push_back({kmer, query_id});
        }
        omp_unset_lock(&_mutex_array[minimizer]);
    }
}

template <typename T>
uint64_t Pac<T>::query_bucket(const std::vector<std::pair<uint64_t, uint32_t>>& colored_kmer, std::vector<std::pair<std::string, std::vector<uint32_t>>>& result, uint bucket_id, std::vector<bool>* BBV) {
    uint64_t sum(0);
    uint64_t size_partition = _nbCells / _bucket_number;
    // _buckets[bucket_id]->load(_leaf_number, _use_double_index, BBV);

    uint64_t mask(size_partition - 1);
    for (uint32_t i = 0; i < colored_kmer.size(); ++i) {
        uint64_t max_local(_leaf_number - 1);
        const auto& [found_min, min_local] = _buckets[bucket_id]->query_key_min(colored_kmer[i].first);
        bool found_max = true;
        if (_use_double_index) {
            const auto& [x, y] = _buckets[bucket_id]->query_key_max(colored_kmer[i].first);
            found_max = x;
            max_local = y;
        }
        if (found_min && found_max) {
            std::cout << "found" << std::endl;
            uint64_t hash = hash_family(colored_kmer[i].first, 0, mask) * _leaf_number;

            for (uint64_t l = min_local + hash; l <= max_local + hash; ++l) {
                sum++;
                if ((*BBV)[l]) {
                    result[colored_kmer[i].second].second[l - hash]++;
                }
            }
        }
    }
    BBV->clear();
    // TODO uncomment ?
    // delete _buckets[bucket_id];
    // _buckets[bucket_id] = nullptr;
    return sum;
}

template <typename T>
std::vector<T> Pac<T>::query_keys(const std::vector<uint64_t>& keys, uint minimizer) {
    std::vector<T> result;
    std::vector<T> colors;
    for (const auto& key : keys) {
        for (const auto& [level, value] : _buckets[minimizer]->query_key(key)) {
            // TODO handle value (tmp.second)
            colors.push_back(level);
        }
        // colors = buckets[minimizer]->query_key(key[i]);
        result.insert(result.end(), colors.begin(), colors.end());
        colors.clear();
    }
    return result;
}

template <typename T>
void Pac<T>::serialize() const {
    std::cout << "Dump the index in " << _w_dir << std::endl;
    std::filesystem::path initial_path = std::filesystem::current_path();
    std::filesystem::path p{_w_dir};
    if (std::filesystem::exists(p)) {
    } else {
        std::filesystem::create_directory(p);
    }
    std::filesystem::current_path(p);
    std::string filename("MainIndex");
    std::filebuf fb;
    fb.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    zstr::ostream out(&fb);
    // VARIOUS INTEGERS
    // TODO _nbBitsPerCell
    out.write(reinterpret_cast<const char*>(&_s), sizeof(_s));
    out.write(reinterpret_cast<const char*>(&_nbCells), sizeof(_nbCells));
    out.write(reinterpret_cast<const char*>(&_number_hash_function), sizeof(_number_hash_function));
    out.write(reinterpret_cast<const char*>(&_small_minimizer_size), sizeof(_small_minimizer_size));
    out.write(reinterpret_cast<const char*>(&_large_minimizer_size), sizeof(_large_minimizer_size));
    out.write(reinterpret_cast<const char*>(&_leaf_number), sizeof(_leaf_number));
    out.write(reinterpret_cast<const char*>(&_use_double_index), sizeof(_use_double_index));
    out << std::flush;
    fb.close();
    std::filesystem::current_path(initial_path);
    std::cout << "Index use " << intToString(_disk_space_used) << " Bytes on disk" << std::endl;
}

template <typename T>
void Pac<T>::index() {
#pragma omp parallel for num_threads(core_number)
    for (uint i = 0; i < _buckets.size(); ++i) {
        // _buckets[i]->load_bf(_leaf_number);
        uint64_t nbbs(_buckets[i]->construct_trunk());
        if (_use_double_index) {
            _buckets[i]->construct_reverse_trunk();
        }
#pragma omp atomic
        _number_bit_set += nbbs;
        // _buckets[i]->serialize();
        nbbs = _buckets[i]->get_disk_space_used();
#pragma omp atomic
        _disk_space_used += nbbs;
        nbbs = _buckets[i]->get_number_bit_set_abt();
#pragma omp atomic
        _number_bit_set_abt += nbbs;
        // TODO uncomment ?
        // delete _buckets[i];
        // _buckets[i] = nullptr;
    }
}

template <typename T>
uint64_t Pac<T>::nb_bit_set() const {
    std::cout << intToString(_number_bit_set) << " bits sets" << std::endl;
    // TODO verify
    std::cout << (double)100 * _number_bit_set / ((double)_leaf_number * _nbCells * _nbBitsPerCell) << "% of bits are set" << std::endl;
    std::cout << intToString(_number_bit_set_abt) << " non empty cell in ABT" << std::endl;
    std::cout << (double)100 * _number_bit_set_abt / (_nbCells * _nbBitsPerCell) << "% of cells are used" << std::endl;
    return _number_bit_set;
}

// OPTIMIZE
template <typename T>
std::vector<std::vector<uint64_t>> Pac<T>::query(const std::string& reference) const {
    std::vector<std::vector<uint64_t>> result(_leaf_number);  // one vector per file
    for (uint64_t i = 0; i < _leaf_number; i++) {
        result[i].resize(reference.size() - _s + 1);
    }
    uint64_t pos = 0;
    for (const auto& [superkmers, minimizer] : get_super_kmers(reference)) {
        for (const auto& kmer : superkmers) {
            auto levels_and_values = _buckets[minimizer]->query_key(kmer);
            for (const auto& [level, value] : levels_and_values) {
                result[level][pos] = value;
            }
            pos++;
        }
    }
    return result;
}

// OPTIMIZE
template <typename T>
std::vector<std::vector<uint64_t>> Pac<T>::fimpera(const std::string& reference) const {
    std::vector<std::vector<uint64_t>> smer_responses = query(reference);
    std::vector<std::vector<uint64_t>> kmer_responses;
    kmer_responses.reserve(smer_responses.size());

    for (const auto& smer_response : smer_responses) {
        kmer_responses.push_back(::fimpera::fimpera(smer_response, _k, _k - _s));
    }

    return kmer_responses;
}

template <typename T>
uint64_t Pac<T>::get_leaf_numbers() const {
    return _leaf_number;
}