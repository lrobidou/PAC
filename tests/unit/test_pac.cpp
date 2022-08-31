#include <gtest/gtest.h>

#include <libpac/Pac.hpp>

inline std::string random_string(size_t length) {
    auto randchar = []() -> char {
        const char charset[] = "ACTG";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[rand() % max_index];
    };
    std::string str(length, 0);
    std::generate_n(str.begin(), length, randchar);
    return str;
}

class SimplePac : public testing::Test {
   protected:
    const uint64_t nbCells;
    const uint64_t nbBitsPerCell;
    const uint number_hash_function;
    const uint k;
    const uint z;
    const std::string wdir;
    bool use_double_index;
    uint bucketing;
    uint core_number;

    Pac<uint8_t> pac;

    SimplePac() : nbCells(10000), nbBitsPerCell(5), number_hash_function(1), k(31), z(0), wdir("test_pac"), use_double_index(true), bucketing(4), core_number(4), pac(nbCells, nbBitsPerCell, number_hash_function, k, z, wdir, use_double_index, bucketing, core_number) {}
};

class TwoPac : public testing::Test {
   protected:
    const uint64_t nbCells;
    const uint64_t nbBitsPerCell;
    const uint number_hash_function;
    const uint k;
    const uint z;
    const std::string wdir;
    bool use_double_index;
    uint bucketing;
    uint core_number;

    Pac<uint8_t> pac;
    Pac<uint8_t> pac_fimpera_enabled;

    // TODO why k=30 and not 31 ??
    // DEBUG
    TwoPac() : nbCells(194957), nbBitsPerCell(5), number_hash_function(1), k(31), z(5), wdir("test_pac"), use_double_index(true), bucketing(4), core_number(4), pac(nbCells, nbBitsPerCell, number_hash_function, k - 1, 0, wdir, use_double_index, bucketing, core_number), pac_fimpera_enabled(nbCells, nbBitsPerCell, number_hash_function, k - 1, z, wdir, use_double_index, bucketing, core_number) {}
};

// shift x n times
uint inline shift(uint x, uint n) {
    return x << n;
}

double get_fpr(const Pac<uint8_t>& pac, int nbElemnTest, int k) {
    int pos = 0;
    int neg = 0;
    for (int i = 0; i < nbElemnTest; i++) {
        std::string s = random_string(k);
        std::vector<std::vector<size_t>> res = pac.fimpera(s);
        EXPECT_EQ(res.size(), pac.get_leaf_numbers());
        for (const auto res_for_a_file : res) {
            EXPECT_EQ(res_for_a_file.size(), 2);
            EXPECT_EQ(res_for_a_file[1], 0);
            int x = res_for_a_file[0];
            if (x == 0) {
                neg++;
            } else {
                pos++;
            }
        }
    }

    EXPECT_EQ(pos + neg, nbElemnTest * pac.get_leaf_numbers());

    double fpr = (double)pos / (double)(pos + neg);

    return fpr;
}

TEST_F(SimplePac, getters) {
    // TODO use getters here
}

TEST_F(SimplePac, insert_file) {
    // TODO test
    uint level = 0;
    const std::string& filename = "../tests/unit/data/KMC_kmers.txt";

    pac.add_leaf();
    pac.insert_file(filename, level, level + 1);
}

TEST_F(TwoPac, noFN) {
    const std::string& input_filename = "../tests/unit/data/10000_31-mers_fof.txt";
    const std::string& query_filename = "../tests/unit/data/query.fasta";

    pac.insert_file_of_file(input_filename);
    pac_fimpera_enabled.insert_file_of_file(input_filename);

    std::ifstream my_file_fof(input_filename);
    std::string filename;
    uint64_t n_file = 0;
    while (std::getline(my_file_fof, filename)) {  // \n are removed
        std::ifstream my_file(filename);
        std::string line;
        while (std::getline(my_file, filename)) {  // \n are removed
            std::string kmer;
            std::string abundanceStr;
            std::stringstream linestream(line);
            std::getline(linestream, kmer, '\t');
            std::getline(linestream, abundanceStr, '\t');
            uint64_t abundance = std::stoll(abundanceStr);
            auto query_result = pac.query(kmer);
            EXPECT_EQ(query_result.size(), 3);
            EXPECT_EQ(query_result[0].size(), 1);
            EXPECT_GE(query_result[n_file][0], abundance);
        }
        n_file++;
    }
}

TEST_F(TwoPac, FPRate) {
    const int nbElemnTest = 100000;
    const std::string& input_filename = "../tests/unit/data/10000_31-mers_fof.txt";
    const std::string& query_filename = "../tests/unit/data/query.fasta";

    pac.insert_file_of_file(input_filename);
    pac_fimpera_enabled.insert_file_of_file(input_filename);

    double fpr_pac = get_fpr(pac, nbElemnTest, k);
    EXPECT_GE(fpr_pac, 0.048);
    EXPECT_LE(fpr_pac, 0.052);

    double fpr_pac_fimpera = get_fpr(pac_fimpera_enabled, nbElemnTest, k);
    EXPECT_GE(fpr_pac_fimpera, 0.0029);
    EXPECT_LE(fpr_pac_fimpera, 0.0031);
}