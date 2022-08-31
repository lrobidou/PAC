#include <gtest/gtest.h>

#include <libpac/CBF.hpp>

TEST(pac_cbf, simple_retrieval_test) {
    // basic test to check if the Bloom filter works
    uint number_hash_function = 1;
    uint64_t nbCell = 10000;
    uint64_t nbBitsPerCell = 10;

    CBF bloom = CBF(nbCell, nbBitsPerCell, number_hash_function);

    bloom.insert_key(0, 4);
    bloom.insert_key(1, 4);
    bloom.insert_key(3, 4);
    bloom.insert_key(4, 4);

    EXPECT_EQ(bloom.get_key(0), 4);
    EXPECT_EQ(bloom.get_key(1), 4);
    EXPECT_EQ(bloom.get_key(2), 0);
    EXPECT_EQ(bloom.get_key(3), 4);
    EXPECT_EQ(bloom.get_key(4), 4);
}

TEST(pac_cbf, one_set) {
    uint number_hash_function = 1;
    uint64_t nbCell = 10000;
    uint64_t nbBitsPerCell = 10;

    uint64_t key = 5;
    uint64_t value = 1;

    CBF cbf = CBF(nbCell, nbBitsPerCell, number_hash_function);
    cbf.insert_key(key, value);
    EXPECT_EQ(cbf.get_key(key), value);
}

inline uint64_t random_key(uint64_t max) {
    return rand() % max;
}

TEST(pac_cbf, multiple_set) {
    uint number_hash_function = 1;
    uint64_t nbCell = 10000;
    uint64_t nbBitsPerCell = 10;

    uint64_t key = 5;
    uint64_t value = 1;

    CBF cbf = CBF(nbCell, nbBitsPerCell, number_hash_function);

    cbf.insert_key(key, value);
    cbf.insert_key(key, value);
    cbf.insert_key(key, value);
    cbf.insert_key(key, value);
    cbf.insert_key(key, value);
    EXPECT_EQ(cbf.get_key(key), value);
}

TEST(pac_cbf, mutilple_set_check_max) {
    uint number_hash_function = 1;
    uint64_t nbCell = 10000;
    uint64_t nbBitsPerCell = 10;

    uint64_t key = 5;

    CBF cbf = CBF(nbCell, nbBitsPerCell, number_hash_function);

    uint64_t expected = 8;

    cbf.insert_key(key, 1);
    cbf.insert_key(key, 2);
    cbf.insert_key(key, 8);
    cbf.insert_key(key, 3);
    cbf.insert_key(key, 1);
    EXPECT_EQ(cbf.get_key(key), expected);
}

TEST(pac_cbf, mutilple_set_check_max_limit) {
    uint number_hash_function = 1;
    uint64_t nbCell = 2000 / 5;
    uint64_t nbBitsPerCell = 5;

    uint64_t key = 5;

    CBF cbf = CBF(nbCell, nbBitsPerCell, number_hash_function);

    uint64_t maxLimit = 31;
    cbf.insert_key(key, 1);
    cbf.insert_key(key, 2);
    cbf.insert_key(key, 8);
    cbf.insert_key(key, 100);
    cbf.insert_key(key, 1);
    EXPECT_EQ(cbf.get_key(key), maxLimit);
}

// TODO write a function to get fpr (and test it)
// formula to get fpr (1-(1-1/m)**(n*k))**k
// here:
// n = nbElemnInsertion
// m = size
// k = 1 (not the size of kmers, the number of hash functions!)
TEST(pac_cbf, fpr_1_hash_1) {
    uint number_hash_function = 1;
    uint64_t nbCell = 194957;
    uint64_t nbBitsPerCell = 1;

    CBF cbf = CBF(nbCell, nbBitsPerCell, number_hash_function);

    int nbElemnInsertion = 10000;
    int nbElemnTest = 100000;
    int key_max = 1000000000;

    for (int i = 0; i < nbElemnInsertion; i++) {
        uint64_t key = random_key(key_max);
        cbf.insert_key(key, 1);
    }

    int pos = 0;
    int neg = 0;
    for (int i = 0; i < nbElemnTest; i++) {
        uint64_t key = random_key(key_max);
        if (cbf.get_key(key) == 0) {
            neg++;
        } else {
            pos++;
        }
    }

    double fpr = (double)pos / (double)(pos + neg);
    EXPECT_EQ(pos + neg, nbElemnTest);
    EXPECT_GE(fpr, 0.048);
    EXPECT_LE(fpr, 0.052);
}

TEST(pac_cbf, fpr_1_hash_4) {
    uint number_hash_function = 4;
    uint64_t nbCell = 194957;
    uint64_t nbBitsPerCell = 1;

    CBF cbf = CBF(nbCell, nbBitsPerCell, number_hash_function);

    int nbElemnInsertion = 10000;
    int nbElemnTest = 1000000;
    int key_max = 1000000000;

    for (int i = 0; i < nbElemnInsertion; i++) {
        uint64_t key = random_key(key_max);
        cbf.insert_key(key, 1);
    }

    int pos = 0;
    int neg = 0;
    for (int i = 0; i < nbElemnTest; i++) {
        uint64_t key = random_key(key_max);
        if (cbf.get_key(key) == 0) {
            neg++;
        } else {
            pos++;
        }
    }

    double fpr = (double)pos / (double)(pos + neg);
    EXPECT_EQ(pos + neg, nbElemnTest);
    // in theory, FPR is about 0.00118
    EXPECT_GE(fpr, 0.00105);
    EXPECT_LE(fpr, 0.00131);
}

TEST(pac_cbf, fpr_1_hash_15) {
    uint number_hash_function = 15;
    uint64_t nbCell = 194957;
    uint64_t nbBitsPerCell = 1;

    CBF cbf = CBF(nbCell, nbBitsPerCell, number_hash_function);

    int nbElemnInsertion = 10000;
    int nbElemnTest = 100000;
    int key_max = 1400000000;

    for (int i = 0; i < nbElemnInsertion; i++) {
        uint64_t key = random_key(key_max);
        cbf.insert_key(key, 1);
    }

    int pos = 0;
    int neg = 0;
    for (int i = 0; i < nbElemnTest; i++) {
        uint64_t key = random_key(key_max);
        if (cbf.get_key(key) == 0) {
            neg++;
        } else {
            pos++;
        }
    }

    double fpr = (double)pos / (double)(pos + neg);
    EXPECT_EQ(pos + neg, nbElemnTest);
    // in theory, FPR is about 8.832775369727666e-05
    EXPECT_GE(fpr, 0);
    EXPECT_LE(fpr, 0.0001);
}

TEST(pac_cbf, fpr_n_hash_1) {
    uint number_hash_function = 1;
    uint64_t nbCell = 194957;
    uint64_t nbBitsPerCell = 5;

    CBF cbf = CBF(nbCell, nbBitsPerCell, number_hash_function);

    int nbElemnInsertion = 10000;
    int nbElemnTest = 100000;
    int key_max = 100000000;

    for (int i = 0; i < nbElemnInsertion; i++) {
        uint64_t key = random_key(key_max);
        cbf.insert_key(key, 1);
    }

    int pos = 0;
    int neg = 0;
    for (int i = 0; i < nbElemnTest; i++) {
        uint64_t key = random_key(key_max);
        if (cbf.get_key(key) == 0) {
            neg++;
        } else {
            pos++;
        }
    }

    double fpr = (double)pos / (double)(pos + neg);
    EXPECT_EQ(pos + neg, nbElemnTest);
    EXPECT_GE(fpr, 0.048);
    EXPECT_LE(fpr, 0.052);
}

TEST(pac_cbf, fpr_n_hash_4) {
    uint number_hash_function = 4;
    uint64_t nbCell = 194957;
    uint64_t nbBitsPerCell = 5;

    CBF cbf = CBF(nbCell, nbBitsPerCell, number_hash_function);

    int nbElemnInsertion = 10000;
    int nbElemnTest = 1000000;
    int key_max = 1500000000;

    for (int i = 0; i < nbElemnInsertion; i++) {
        uint64_t key = random_key(key_max);
        cbf.insert_key(key, 1);
    }

    int pos = 0;
    int neg = 0;
    for (int i = 0; i < nbElemnTest; i++) {
        uint64_t key = random_key(key_max);
        if (cbf.get_key(key) == 0) {
            neg++;
        } else {
            pos++;
        }
    }

    double fpr = (double)pos / (double)(pos + neg);
    EXPECT_EQ(pos + neg, nbElemnTest);
    EXPECT_GE(fpr, 0.00114);
    EXPECT_LE(fpr, 0.00122);
}