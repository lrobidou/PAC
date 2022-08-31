#include <libpac/Bloom.h>
#include <libpac/bestpart.h>
#include <libpac/utils.h>
#include <math.h>

#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <libpac/Bucket.hpp>
#include <string>
#include <vector>

template <typename T>
void BestPart<T>::insert_keys(const std::vector<uint64_t>& key, uint minimizer, uint level, uint64_t value, Bloom* unique_filter) {
    if (filter) {
        for (uint i = 0; i < key.size(); ++i) {
            if (unique_filter->check_key(key[i])) {
                buckets[minimizer]->insert_key(key[i], level, value);
            } else {
                // TODO lrobidou: what do I do with value ?
                unique_filter->insert_key(key[i]);
            }
        }
    } else {
        for (uint i = 0; i < key.size(); ++i) {
            buckets[minimizer]->insert_key(key[i], level, value);
        }
    }
}

// TODO CAN BE IMPROVED
// rcb: reverse complement, binary
template <typename T>
uint64_t BestPart<T>::rcb(uint64_t min, uint64_t n) const {
    uint64_t res(0);
    uint64_t offset(1);  // used like a mask
    offset <<= (2 * n - 2);
    for (uint i(0); i < n; ++i) {
        res += (3 - (min % 4)) * offset;
        min >>= 2;
        offset >>= 2;
    }
    return res;
}

template <typename T>
uint64_t BestPart<T>::canonize(uint64_t x, uint64_t n) {
    return std::min(x, rcb(x, n));
}

// lrobidou: why not returning a tuple ? AFAIU position is not a parameter
template <typename T>
uint64_t BestPart<T>::regular_minimizer_pos(uint64_t seq, uint64_t& position) {
    uint64_t mini, mmer;
    mmer = seq % large_minimizer_number;
    mini = mmer = canonize(mmer, large_minimizer_size);
    uint64_t hash_mini = (unrevhash(mmer));
    position = 0;
    for (uint64_t i(1); i <= K - large_minimizer_size; i++) {
        seq >>= 2;
        mmer = seq % large_minimizer_number;
        mmer = canonize(mmer, large_minimizer_size);
        uint64_t hash = (unrevhash(mmer));
        if (hash_mini > hash) {
            position = K - large_minimizer_size - i;
            mini = mmer;
            hash_mini = hash;
        }
    }
    return mini;
}

template <typename T>
void BestPart<T>::updateK(uint64_t& min, char nuc) const {
    min <<= 2;
    min += nuc2int(nuc);
    min %= offsetUpdatekmer;
}

template <typename T>
void BestPart<T>::updateM(uint64_t& min, char nuc) const {
    min <<= 2;
    min += nuc2int(nuc);
    min %= offsetUpdateminimizer;
}

template <typename T>
void BestPart<T>::updateRCK(uint64_t& min, char nuc) const {
    min >>= 2;
    min += (nuc2intrc(nuc) << (2 * K - 2));
}

template <typename T>
void BestPart<T>::updateRCM(uint64_t& min, char nuc) const {
    min >>= 2;
    min += (nuc2intrc(nuc) << (2 * large_minimizer_size - 2));
}

template <typename T>
std::vector<std::pair<std::vector<uint64_t>, uint64_t> > BestPart<T>::get_super_kmers(const std::string& ref) {
    std::vector<std::pair<std::vector<uint64_t>, uint64_t> > result;
    uint64_t old_minimizer, minimizer;
    std::vector<uint64_t> superkmer;
    old_minimizer = minimizer = large_minimizer_number;
    uint64_t last_position(0);
    // FOREACH KMER
    uint64_t seq(str2num(ref.substr(0, K)));
    uint64_t rcseq = rcb(seq, K);
    uint64_t canon = std::min(seq, rcseq);
    uint64_t position_min;
    uint64_t min_seq = (str2num(ref.substr(K - large_minimizer_size, large_minimizer_size)));
    uint64_t min_rcseq(rcb(min_seq, large_minimizer_size));
    uint64_t min_canon(std::min(min_seq, min_rcseq));
    minimizer = regular_minimizer_pos(seq, position_min);
    old_minimizer = minimizer;
    uint64_t hash_min = unrevhash(minimizer);
    uint64_t i(0);
    for (; i + K < ref.size(); ++i) {
        updateK(seq, ref[i + K]);
        updateRCK(rcseq, ref[i + K]);
        canon = std::min(seq, rcseq);
        updateM(min_seq, ref[i + K]);
        updateRCM(min_rcseq, ref[i + K]);
        min_canon = (std::min(min_seq, min_rcseq));
        uint64_t new_h = unrevhash(min_canon);
        // THE NEW mmer is a MINIMIZER
        if (new_h < hash_min) {
            minimizer = (min_canon);
            hash_min = new_h;
            position_min = i + K - large_minimizer_size + 1;
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
        if (revhash(old_minimizer) % bucket_number != revhash(minimizer) % bucket_number) {
            old_minimizer = (revhash(old_minimizer) % bucket_number);
            result.push_back({superkmer, old_minimizer});
            superkmer.clear();
            last_position = i + 1;
            old_minimizer = minimizer;
        }
        superkmer.push_back(canon);
    }
    if (ref.size() - last_position > K - 1) {
        old_minimizer = (revhash(old_minimizer) % bucket_number);
        result.push_back({superkmer, old_minimizer});
        superkmer.clear();
    }
    return result;
}

template <typename T>
void BestPart<T>::insert_sequence(const std::string& reference, uint level, uint64_t value, Bloom* unique_filter) {
    for (const auto& [superkmers, minimizer] : get_super_kmers(reference)) {
        insert_keys(superkmers, minimizer, level, value, unique_filter);
    }
}

template <typename T>
void BestPart<T>::query_file(const std::string& filename, const std::string& output) {
    std::cout << "Query file " << filename << " in " << output << std::endl;
    auto start = std::chrono::system_clock::now();
    zstr::ofstream out(output, std::ios::trunc);
    char type = get_data_type(filename);
    zstr::ifstream in(filename);
    std::filesystem::path old(std::filesystem::current_path());
    std::vector<uint32_t> colors_count(leaf_number, 0);

    // lrobidou: why pointers ?
    std::filesystem::current_path(w_dir);
    auto result = new std::vector<std::pair<std::string, std::vector<uint32_t> > >;
    auto colored_kmer_per_bucket = new std::vector<std::vector<std::pair<uint64_t, uint32_t> > >(bucket_number);
#pragma omp parallel num_threads(core_number)
    {
        std::string ref, header;
        uint32_t query_id;
        while (not in.eof()) {
#pragma omp critical(inputfile)
            {
                // TODO uncomment
                // Biogetline(&in, ref, type, K, header);
                if (ref.size() > K) {
                    result->push_back({header, colors_count});
                    query_id = result->size() - 1;
                }
            }
            if (ref.size() > K) {
                load_super_kmer(*colored_kmer_per_bucket, query_id, ref);
            }
        }
    }
    uint64_t sum_query = 0;
#pragma omp parallel num_threads(core_number)
    {
        std::vector<bool>* BBV = new std::vector<bool>(leaf_number * size / bucket_number, false);
#pragma omp for
        for (uint i = 0; i < bucket_number; i++) {
            if (not(*colored_kmer_per_bucket)[i].empty()) {
                uint64_t sq(query_bucket((*colored_kmer_per_bucket)[i], *result, i, BBV));
#pragma omp atomic
                sum_query += sq;
            }
        }
        delete BBV;
    }
    std::cout << "query done" << std::endl;
    for (uint i = 0; i < result->size(); i++) {
        {
            out << (*result)[i].first << "\n";
            for (uint j = 0; j < (*result)[i].second.size(); j++) {
                out << (*result)[i].second[j] << ' ';
            }
            out << "\n";
        }
    }
    out.close();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Query time: " << elapsed_seconds.count() << "s for " << intToString(sum_query) << " kmer query or " << intToString((double)sum_query / elapsed_seconds.count()) << " query per second  \n";
    delete result;
    delete colored_kmer_per_bucket;
    current_path(old);
}

template <typename T>
void BestPart<T>::insert_previous_index(const std::string& filename) {
    zstr::ifstream in(filename + "/MainIndex");
    uint64_t osef, leaf_number_to_insert;
    in.read(reinterpret_cast<char*>(&osef), sizeof(K));
    in.read(reinterpret_cast<char*>(&osef), sizeof(size));
    in.read(reinterpret_cast<char*>(&osef), sizeof(number_hash_function));
    in.read(reinterpret_cast<char*>(&osef), sizeof(small_minimizer_size));
    in.read(reinterpret_cast<char*>(&osef), sizeof(large_minimizer_size));
    in.read(reinterpret_cast<char*>(&leaf_number_to_insert), sizeof(leaf_number_to_insert));
    for (uint i = 0; i < buckets.size(); ++i) {
        zstr::ifstream in(filename + "/P" + std::to_string(i));
        for (uint ii = 0; ii < leaf_number_to_insert; ++ii) {
            uint32_t ibloom;
            in.read(reinterpret_cast<char*>(&ibloom), sizeof(ibloom));
            uint64_t sz;
            in.read(reinterpret_cast<char*>(&sz), sizeof(sz));
            uint8_t* buff = new uint8_t[sz];
            in.read((char*)buff, sz);
            buckets[i]->get_out()->write(reinterpret_cast<const char*>(&ibloom), sizeof(ibloom));
            buckets[i]->get_out()->write(reinterpret_cast<const char*>(&sz), sizeof(sz));
            buckets[i]->get_out()->write((char*)&buff[0], sz);
        }
        buckets[i]->get_out()->flush();
    }
    leaf_number += leaf_number_to_insert;
    std::cout << leaf_number << " datasets loaded from previous index" << std::endl;
}

template <typename T>
void BestPart<T>::load_super_kmer(std::vector<std::vector<std::pair<uint64_t, uint32_t> > >& colored_kmer_per_bucket, uint32_t query_id, const std::string& reference) {
    auto V(get_super_kmers(reference));
    for (uint32_t i = 0; i < V.size(); ++i) {
        omp_set_lock(&mutex_array[V[i].second]);
        for (uint32_t j = 0; j < V[i].first.size(); ++j) {
            (colored_kmer_per_bucket[V[i].second]).push_back({V[i].first[j], query_id});
        }
        omp_unset_lock(&mutex_array[V[i].second]);
    }
}

// TODO remove ?
// template <typename T>
// uint64_t BestPart<T>::query_bucket2(const std::vector<std::pair<uint64_t, uint32_t> >& colored_kmer, std::vector<std::pair<std::string, std::vector<uint32_t> > >& result, uint bucket_id) {
//     uint64_t sum(0);
//     buckets[bucket_id]->load(leaf_number, use_double_index);
//     std::vector<uint32_t> minV(colored_kmer.size(), 0);
//     uint32_t minimal_value(leaf_number), maximal_value(leaf_number - 1);
//     std::vector<uint32_t> maxV(colored_kmer.size(), leaf_number - 1);
//     for (uint32_t i = 0; i < colored_kmer.size(); ++i) {
//         minV[i] = buckets[bucket_id]->query_key_min(colored_kmer[i].first);
//         if (minV[i] < minimal_value) {
//             minimal_value = minV[i];
//         }
//     }
//     if (use_double_index) {
//         maxV.resize(colored_kmer.size(), 0);
//         maximal_value = 0;
//         for (uint32_t i = 0; i < colored_kmer.size(); ++i) {
//             maxV[i] = buckets[bucket_id]->query_key_max(colored_kmer[i].first);
//             if (maxV[i] > maximal_value) {
//                 maximal_value = maxV[i];
//             }
//         }
//     }

//     for (uint32_t l = (minimal_value); l < maximal_value; ++l) {
//         for (uint32_t i = 0; i < colored_kmer.size(); ++i) {
//             if (l >= minV[i] and l <= maxV[i]) {
//                 sum++;
//                 if (buckets[bucket_id]->get_leaf_filters()[l]->get_key(colored_kmer[i].first)) {
//                     // TODO use value
// #pragma omp atomic
//                     result[colored_kmer[i].second].second[l]++;
//                 }
//             }
//         }
//     }
//     //~ std::cout<<"end read BBv"<<std::endl;
//     delete buckets[bucket_id];
//     buckets[bucket_id] = NULL;
//     return sum;
// }

template <typename T>
uint64_t BestPart<T>::query_bucket(const std::vector<std::pair<uint64_t, uint32_t> >& colored_kmer, std::vector<std::pair<std::string, std::vector<uint32_t> > >& result, uint bucket_id, std::vector<bool>* BBV) {
    uint64_t sum(0);
    uint64_t size_partition = size / bucket_number;
    buckets[bucket_id]->load(leaf_number, use_double_index, BBV);
    //~ std::vector<uint32_t> minV(colored_kmer.size(),0);
    //~ uint32_t minimal_value(leaf_number),maximal_value(leaf_number-1);
    //~ std::vector<uint32_t> maxV(colored_kmer.size(),leaf_number-1);
    //~ for(uint32_t i = 0; i < colored_kmer.size();++i){
    //~ minV[i]=buckets[bucket_id]->query_key_min(colored_kmer[i].first);
    //~ if(minV[i] < minimal_value){
    //~ minimal_value=minV[i];
    //~ }
    //~ }
    //~ if(use_double_index){
    //~ maxV.resize(colored_kmer.size(),0);
    //~ maximal_value=0;
    //~ for(uint32_t i = 0; i < colored_kmer.size();++i){
    //~ maxV[i]=buckets[bucket_id]->query_key_max(colored_kmer[i].first);
    //~ if(maxV[i]>maximal_value){
    //~ maximal_value=maxV[i];
    //~ }
    //~ }
    //~ }
    uint64_t mask(size_partition - 1);
    for (uint32_t i = 0; i < colored_kmer.size(); ++i) {
        uint64_t max_local(leaf_number - 1);
        const auto& [found_min, min_local] = buckets[bucket_id]->query_key_min(colored_kmer[i].first);
        bool found_max = true;
        if (use_double_index) {
            const auto& [x, y] = buckets[bucket_id]->query_key_max(colored_kmer[i].first);
            found_max = x;
            max_local = y;
        }
        if (found_min && found_max) {
            uint64_t hash(hash_family(colored_kmer[i].first, 0, mask) * leaf_number);

            for (uint64_t l = min_local + hash; l <= max_local + hash; ++l) {
                sum++;
                if ((*BBV)[l]) {
                    result[colored_kmer[i].second].second[l - hash]++;
                }
            }
        }
    }
    BBV->clear();
    delete buckets[bucket_id];
    buckets[bucket_id] = NULL;
    return sum;
}

template <typename T>
std::vector<T> BestPart<T>::query_keys(const std::vector<uint64_t>& key, uint minimizer) {
    std::vector<T> result;
    std::vector<T> colors;
    for (uint i = 0; i < key.size(); ++i) {
        auto tmp = buckets[minimizer]->query_key(key[i]);
        // TODO handle value (tmp.second)
        for (const auto& [level, value] : tmp) {
            colors.push_back(level);
        }
        // colors = buckets[minimizer]->query_key(key[i]);
        result.insert(result.end(), colors.begin(), colors.end());
        colors.clear();
    }
    return result;
}

template <typename T>
void BestPart<T>::insert_file(const std::string& filename, uint level, uint32_t indice_bloom) {
    char type = get_data_type(filename);
    zstr::ifstream in(filename);
    Bloom* unique_filter;
    if (filter) {
        Bucket<T>* first = buckets[0];
        unique_filter = new Bloom(first->get_nbCells(), first->get_number_hash_function());
    }
    {
        std::string ref;
        while (not in.eof()) {
            {
                // TODO uncomment
                // Biogetline(&in, ref, type, K);
            }
            if (ref.size() > K) {
                insert_sequence(ref, level, 1, unique_filter);  // TODO
            }
        }
    }
    if (filter) {
        delete unique_filter;
    }
    bm::serializer<bm::bvector<> > bvs;
    bvs.byte_order_serialization(false);
    bvs.gap_length_serialization(false);
    for (uint i = 0; i < buckets.size(); ++i) {
        buckets[i]->optimize(level);
        omp_set_lock(&mutex_array[i]);  // TODO THE LOCK SHOULD NOT INCLUDE THE SERIALIZATION
        buckets[i]->get_out()->write(reinterpret_cast<const char*>(&indice_bloom), sizeof(indice_bloom));
        buckets[i]->dump(level, bvs);
        omp_unset_lock(&mutex_array[i]);
    }
}

template <typename T>
void BestPart<T>::insert_file_of_file(const std::string& filename) {
    std::cout << "I index " << K << "mers with Bloom filters of size " << intToString(size) << " with " << number_hash_function << " hash functions  using " << intToString(1 << (2 * small_minimizer_size)) << " partitions " << std::endl;

    std::filesystem::path old(std::filesystem::current_path());
    std::filesystem::path vodka(std::filesystem::absolute(filename));
    zstr::ifstream in(vodka);
    std::cout << "Insert file of file " << filename << std::endl;
    std::filesystem::path initial_path = std::filesystem::current_path();
    for (uint i(0); i < core_number; ++i) {
        add_leaf();
    }

    auto start = std::chrono::system_clock::now();
    ;
#pragma omp parallel num_threads(core_number)
    {
        std::string ref;
        uint level;
        bool go;
        int idt = omp_get_thread_num();
        while (not in.eof()) {
#pragma omp critical(inputfile)
            {
                std::getline(in, ref);
                if (exists_test(ref)) {
                    level = leaf_number;
                    leaf_number++;
                    go = true;
                } else {
                    go = false;
                }
            }
            if (go) {
                insert_file(ref, idt, level);
                if (level % 100 == 0) {
                    std::cout << level << " files" << std::endl;
                }
            }
        }
    }
    auto middle = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = middle - start;
    std::cout << "Bloom construction time: " << elapsed_seconds.count() << "s\n";
    std::cout << intToString(getMemorySelfMaxUsed()) << " KB total" << std::endl;
    std::cout << intToString(getMemorySelfMaxUsed() / (leaf_number)) << " KB per file (" << leaf_number << " files)" << std::endl;
    index();
    nb_bit_set();
    auto end = std::chrono::system_clock::now();
    elapsed_seconds = end - middle;
    std::cout << "Exponential Bloom construction time: " << elapsed_seconds.count() << "s\n";
    elapsed_seconds = end - start;
    std::cout << "Total Index time: " << elapsed_seconds.count() << "s\n";
    std::cout << "Index time per file: " << elapsed_seconds.count() / leaf_number << "s\n";
    std::cout << intToString(getMemorySelfMaxUsed()) << " kB total" << std::endl;
    std::cout << intToString(getMemorySelfMaxUsed() / (leaf_number)) << " kB per file" << std::endl;
    std::filesystem::current_path(old);
}

template <typename T>
void BestPart<T>::serialize() const {
    std::cout << "Dump the index in " << w_dir << std::endl;
    std::filesystem::path initial_path = std::filesystem::current_path();
    std::filesystem::path p{w_dir};
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
    out.write(reinterpret_cast<const char*>(&K), sizeof(K));
    out.write(reinterpret_cast<const char*>(&size), sizeof(size));
    out.write(reinterpret_cast<const char*>(&number_hash_function), sizeof(number_hash_function));
    out.write(reinterpret_cast<const char*>(&small_minimizer_size), sizeof(small_minimizer_size));
    out.write(reinterpret_cast<const char*>(&large_minimizer_size), sizeof(large_minimizer_size));
    out.write(reinterpret_cast<const char*>(&leaf_number), sizeof(leaf_number));
    out.write(reinterpret_cast<const char*>(&use_double_index), sizeof(use_double_index));
    out << std::flush;
    fb.close();
    std::filesystem::current_path(initial_path);
    std::cout << "Index use " << intToString(disk_space_used) << " Bytes on disk" << std::endl;
}

template <typename T>
void BestPart<T>::index() {
#pragma omp parallel for num_threads(core_number)
    for (uint i = 0; i < buckets.size(); ++i) {
        buckets[i]->load_bf(leaf_number);
        uint64_t nbbs(buckets[i]->construct_trunk());
        if (use_double_index) {
            buckets[i]->construct_reverse_trunk();
        }
#pragma omp atomic
        number_bit_set += nbbs;
        buckets[i]->serialize();
        nbbs = buckets[i]->get_disk_space_used();
#pragma omp atomic
        disk_space_used += nbbs;
        nbbs = buckets[i]->get_number_bit_set_abt();
#pragma omp atomic
        number_bit_set_abt += nbbs;
        delete buckets[i];
        buckets[i] = NULL;
    }
}

template <typename T>
void BestPart<T>::add_leaf() {
    //~ #pragma omp parallel for
    for (uint i = 0; i < buckets.size(); ++i) {
        buckets[i]->add_leaf();
    }
}

template <typename T>
uint64_t BestPart<T>::nb_bit_set() const {
    std::cout << intToString(number_bit_set) << " bits sets" << std::endl;
    std::cout << (double)100 * number_bit_set / ((double)leaf_number * size) << "% of bits are set" << std::endl;
    std::cout << intToString(number_bit_set_abt) << " non empty cell in ABT" << std::endl;
    std::cout << (double)100 * number_bit_set_abt / (size) << "% of cells are used" << std::endl;
    return number_bit_set;
}
