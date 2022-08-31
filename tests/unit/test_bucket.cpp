#include <gtest/gtest.h>

#include <libpac/Bucket.hpp>

class SimpleBucket : public testing::Test {
   protected:
    uint64_t nbCells;
    uint number_hash_function;
    uint k;
    std::string prefix;
    uint64_t nbBitsPerCell;
    bool write;
    Bucket<uint8_t> bucket;

    SimpleBucket() : nbCells(10000), number_hash_function(1), k(31), prefix("prefix"), nbBitsPerCell(5), write(true), bucket(Bucket<uint8_t>(nbCells, number_hash_function, k, prefix, nbBitsPerCell, write)) {
    }
};

class BucketOneLeaf : public testing::Test {
   protected:
    uint64_t nbCells;
    uint number_hash_function;
    uint k;
    std::string prefix;
    uint64_t nbBitsPerCell;
    bool write;
    Bucket<uint8_t> bucket;
    uint64_t key0;

    BucketOneLeaf() : nbCells(10000), number_hash_function(1), k(31), prefix("prefix"), nbBitsPerCell(5), write(true), bucket(Bucket<uint8_t>(nbCells, number_hash_function, k, prefix, nbBitsPerCell, write)), key0(0) {
        uint8_t level = 0;

        key0 = 0;
        uint64_t key1 = 5;
        uint64_t key2 = 23;

        bucket.add_leaf();
        bucket.insert_key(key0, level, 15);
        bucket.insert_key(key1, level, 14);
        bucket.insert_key(key2, level, 15);
    }
};

class BucketFiveLeafs : public testing::Test {
   protected:
    uint64_t nbCells;
    uint number_hash_function;
    uint k;
    std::string prefix;
    uint64_t nbBitsPerCell;
    bool write;
    Bucket<uint8_t> bucket;
    uint64_t key00;
    uint64_t key08;
    std::vector<uint64_t> keys;

    BucketFiveLeafs() : nbCells(10000), number_hash_function(1), k(31), prefix("prefix"), nbBitsPerCell(5), write(true), bucket(Bucket<uint8_t>(nbCells, number_hash_function, k, prefix, nbBitsPerCell, write)) {
        uint8_t level = 0;

        key00 = 0;
        uint64_t key01 = 5;
        uint64_t key02 = 23;
        uint64_t key03 = 85;
        uint64_t key04 = 100;
        uint64_t key05 = 50000;
        uint64_t key06 = 60000;
        uint64_t key07 = 455225;
        key08 = 465465;
        uint64_t key09 = 9654755;
        uint64_t key10 = 13541354;
        uint64_t key11 = 2451;

        keys = {
            key00,
            key01,
            key02,
            key03,
            key04,
            key05,
            key06,
            key07,
            key08,
            key09,
            key10,
            key11};

        bucket.add_leaf();
        bucket.add_leaf();
        bucket.add_leaf();
        bucket.add_leaf();
        bucket.add_leaf();

        bucket.insert_key(key00, level + 0, 15);
        bucket.insert_key(key01, level + 0, 14);
        bucket.insert_key(key02, level + 0, 15);

        bucket.insert_key(key00, level + 1, 25);
        bucket.insert_key(key03, level + 1, 29);
        bucket.insert_key(key04, level + 1, 12);
        bucket.insert_key(key05, level + 1, 20);

        bucket.insert_key(key00, level + 2, 12);
        bucket.insert_key(key02, level + 2, 15);
        bucket.insert_key(key06, level + 2, 12);
        bucket.insert_key(key07, level + 2, 20);

        bucket.insert_key(key00, level + 3, 11);
        bucket.insert_key(key01, level + 3, 14);
        bucket.insert_key(key02, level + 3, 15);

        bucket.insert_key(key00, level + 4, 1);
        bucket.insert_key(key03, level + 4, 29);
        bucket.insert_key(key04, level + 4, 12);
        bucket.insert_key(key05, level + 4, 20);
        bucket.insert_key(key00, level + 4, 12);
        bucket.insert_key(key02, level + 4, 15);
        bucket.insert_key(key06, level + 4, 12);
        bucket.insert_key(key07, level + 4, 20);
    }
};

template <typename U, typename V>
std::tuple<U, V> t(U a, V b) {
    return {a, b};
}

inline void check_keys_min_max_before_trunk_construction(const std::vector<uint64_t> keys, const Bucket<uint8_t>& bucket) {
    // change from 5 to 4 as max
    EXPECT_EQ(bucket.query_key_min(keys[0]), t(1, 0));
    EXPECT_EQ(bucket.query_key_max(keys[0]), t(1, 4));
    EXPECT_EQ(bucket.query_key_min(keys[1]), t(1, 0));
    EXPECT_EQ(bucket.query_key_max(keys[1]), t(1, 4));
    EXPECT_EQ(bucket.query_key_min(keys[2]), t(1, 0));
    EXPECT_EQ(bucket.query_key_max(keys[2]), t(1, 4));
    EXPECT_EQ(bucket.query_key_min(keys[3]), t(1, 0));
    EXPECT_EQ(bucket.query_key_max(keys[3]), t(1, 4));
    EXPECT_EQ(bucket.query_key_min(keys[4]), t(1, 0));
    EXPECT_EQ(bucket.query_key_max(keys[4]), t(1, 4));
    EXPECT_EQ(bucket.query_key_min(keys[5]), t(1, 0));
    EXPECT_EQ(bucket.query_key_max(keys[5]), t(1, 4));
    EXPECT_EQ(bucket.query_key_min(keys[6]), t(1, 0));
    EXPECT_EQ(bucket.query_key_max(keys[6]), t(1, 4));
    EXPECT_EQ(bucket.query_key_min(keys[7]), t(1, 0));
    EXPECT_EQ(bucket.query_key_max(keys[7]), t(1, 4));
    EXPECT_EQ(bucket.query_key_min(keys[8]), t(1, 0));
    EXPECT_EQ(bucket.query_key_max(keys[8]), t(1, 4));
    EXPECT_EQ(bucket.query_key_min(keys[9]), t(1, 0));
    EXPECT_EQ(bucket.query_key_max(keys[9]), t(1, 4));
    EXPECT_EQ(bucket.query_key_min(keys[10]), t(1, 0));
    EXPECT_EQ(bucket.query_key_max(keys[10]), t(1, 4));
    EXPECT_EQ(bucket.query_key_min(keys[11]), t(1, 0));
    EXPECT_EQ(bucket.query_key_max(keys[11]), t(1, 4));
}

inline void check_keys_min_max_after_trunk_construction(const std::vector<uint64_t> keys, const Bucket<uint8_t>& bucket) {
    EXPECT_EQ(bucket.query_key_min(keys[0]), t(true, 0));
    EXPECT_EQ(bucket.query_key_max(keys[0]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[1]), t(true, 0));
    EXPECT_EQ(bucket.query_key_max(keys[1]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[2]), t(true, 0));
    EXPECT_EQ(bucket.query_key_max(keys[2]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[3]), t(true, 1));
    EXPECT_EQ(bucket.query_key_max(keys[3]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[4]), t(true, 1));
    EXPECT_EQ(bucket.query_key_max(keys[4]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[5]), t(true, 1));
    EXPECT_EQ(bucket.query_key_max(keys[5]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[6]), t(true, 2));
    EXPECT_EQ(bucket.query_key_max(keys[6]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[7]), t(true, 2));
    EXPECT_EQ(bucket.query_key_max(keys[7]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[8]), t(false, 5));
    EXPECT_EQ(bucket.query_key_max(keys[8]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[9]), t(false, 5));
    EXPECT_EQ(bucket.query_key_max(keys[9]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[10]), t(false, 5));
    EXPECT_EQ(bucket.query_key_max(keys[10]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[11]), t(false, 5));
    EXPECT_EQ(bucket.query_key_max(keys[11]), t(true, 4));
}

inline void check_keys_min_max_after_both_trunk_construction(const std::vector<uint64_t> keys, const Bucket<uint8_t>& bucket) {
    EXPECT_EQ(bucket.query_key_min(keys[0]), t(true, 0));
    EXPECT_EQ(bucket.query_key_max(keys[0]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[1]), t(true, 0));
    EXPECT_EQ(bucket.query_key_max(keys[1]), t(true, 3));
    EXPECT_EQ(bucket.query_key_min(keys[2]), t(true, 0));
    EXPECT_EQ(bucket.query_key_max(keys[2]), t(true, 4));

    EXPECT_EQ(bucket.query_key_min(keys[3]), t(true, 1));
    EXPECT_EQ(bucket.query_key_max(keys[3]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[4]), t(true, 1));
    EXPECT_EQ(bucket.query_key_max(keys[4]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[5]), t(true, 1));
    EXPECT_EQ(bucket.query_key_max(keys[5]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[6]), t(true, 2));
    EXPECT_EQ(bucket.query_key_max(keys[6]), t(true, 4));
    EXPECT_EQ(bucket.query_key_min(keys[7]), t(true, 2));
    EXPECT_EQ(bucket.query_key_max(keys[7]), t(true, 4));

    EXPECT_EQ(bucket.query_key_min(keys[8]), t(false, 5));
    EXPECT_EQ(bucket.query_key_max(keys[8]), t(false, 0));
    EXPECT_EQ(bucket.query_key_min(keys[9]), t(false, 5));
    EXPECT_EQ(bucket.query_key_max(keys[9]), t(false, 0));
    EXPECT_EQ(bucket.query_key_min(keys[10]), t(false, 5));
    EXPECT_EQ(bucket.query_key_max(keys[10]), t(false, 0));
    EXPECT_EQ(bucket.query_key_min(keys[11]), t(false, 5));
    EXPECT_EQ(bucket.query_key_max(keys[11]), t(false, 0));
}

TEST_F(SimpleBucket, getters) {
    EXPECT_EQ(bucket.get_leaf_filters(), std::vector<CBF*>());
    EXPECT_EQ(bucket.get_nbCells(), nbCells);
    EXPECT_EQ(bucket.get_number_hash_function(), number_hash_function);
    EXPECT_EQ(bucket.get_number_bit_set_abt(), 0);
    EXPECT_EQ(bucket.get_disk_space_used(), 0);
    EXPECT_EQ(bucket.get_leaf_number(), 0);
}

TEST_F(BucketOneLeaf, one_cbf) {
    // EXPECT_EQ(bucket.get_leaf_filters(), std::vector<CBF*>());
    EXPECT_EQ(bucket.get_nbCells(), nbCells);
    EXPECT_EQ(bucket.get_number_hash_function(), number_hash_function);
    EXPECT_EQ(bucket.get_number_bit_set_abt(), 0);
    EXPECT_EQ(bucket.get_disk_space_used(), 0);
    EXPECT_EQ(bucket.get_leaf_number(), 1);

    std::vector<std::tuple<uint8_t, uint64_t> > ids_in_which_is_key_0 = bucket.query_key(key0);
    std::vector<std::tuple<uint8_t, uint64_t> > expected;
    expected.push_back({0, 15});

    EXPECT_EQ(ids_in_which_is_key_0, expected);
    EXPECT_EQ(ids_in_which_is_key_0.size(), expected.size());
    for (size_t i = 0; i < ids_in_which_is_key_0.size(); i++) {
        EXPECT_EQ(ids_in_which_is_key_0[i], expected[i]);
    }
}

TEST(pac_bucket, three_cbf) {
    const uint64_t nbCells = 10000;
    const uint number_hash_function = 1;
    const uint k = 31;
    const std::string prefix = "prefix";  // TODO
    const uint64_t nbBitsPerCell = 5;
    const bool write = true;

    Bucket<uint8_t> bucket(nbCells, number_hash_function, k, prefix, nbBitsPerCell, write);

    uint8_t level = 0;

    uint64_t key0 = 0;
    uint64_t key1 = 5;
    uint64_t key2 = 23;
    uint64_t key3 = 85;
    uint64_t key4 = 100;
    uint64_t key5 = 50000;
    uint64_t key6 = 60000;
    uint64_t key7 = 455225;

    bucket.add_leaf();
    bucket.add_leaf();
    bucket.add_leaf();

    bucket.insert_key(key0, level + 0, 15);
    bucket.insert_key(key1, level + 0, 14);
    bucket.insert_key(key2, level + 0, 15);
    bucket.insert_key(key0, level + 1, 25);
    bucket.insert_key(key3, level + 1, 29);
    bucket.insert_key(key4, level + 1, 12);
    bucket.insert_key(key5, level + 1, 20);
    bucket.insert_key(key0, level + 2, 12);
    bucket.insert_key(key2, level + 2, 15);
    bucket.insert_key(key6, level + 2, 12);
    bucket.insert_key(key7, level + 2, 20);

    std::vector<std::tuple<uint8_t, uint64_t> > ids_in_which_is_key_0 = bucket.query_key(key0);
    std::vector<std::tuple<uint8_t, uint64_t> > expected;

    EXPECT_EQ(bucket.get_nbCells(), nbCells);
    EXPECT_EQ(bucket.get_number_hash_function(), number_hash_function);
    EXPECT_EQ(bucket.get_number_bit_set_abt(), 0);
    EXPECT_EQ(bucket.get_disk_space_used(), 0);
    EXPECT_EQ(bucket.get_leaf_number(), 3);

    ids_in_which_is_key_0 = bucket.query_key(key0);
    expected.push_back({0, 15});
    expected.push_back({1, 25});
    expected.push_back({2, 12});
    EXPECT_EQ(ids_in_which_is_key_0.size(), expected.size());
    for (size_t i = 0; i < ids_in_which_is_key_0.size(); i++) {
        EXPECT_EQ(ids_in_which_is_key_0[i], expected[i]);
    }
}

TEST_F(BucketFiveLeafs, five_cbf) {
    // EXPECT_EQ(bucket.get_leaf_filters(), std::vector<CBF*>());
    EXPECT_EQ(bucket.get_nbCells(), nbCells);
    EXPECT_EQ(bucket.get_number_hash_function(), number_hash_function);
    EXPECT_EQ(bucket.get_number_bit_set_abt(), 0);
    EXPECT_EQ(bucket.get_disk_space_used(), 0);
    EXPECT_EQ(bucket.get_leaf_number(), 5);

    std::vector<std::tuple<uint8_t, uint64_t> > ids_in_which_is_key_00 = bucket.query_key(key00);
    std::vector<std::tuple<uint8_t, uint64_t> > expected_key00;
    expected_key00.push_back({0, 15});
    expected_key00.push_back({1, 25});
    expected_key00.push_back({2, 12});
    expected_key00.push_back({3, 11});
    expected_key00.push_back({4, 12});

    EXPECT_EQ(ids_in_which_is_key_00, expected_key00);
    EXPECT_EQ(ids_in_which_is_key_00.size(), expected_key00.size());
    for (size_t i = 0; i < ids_in_which_is_key_00.size(); i++) {
        EXPECT_EQ(ids_in_which_is_key_00[i], expected_key00[i]);
    }

    std::vector<std::tuple<uint8_t, uint64_t> > ids_in_which_is_key_08 = bucket.query_key(key08);
    std::vector<std::tuple<uint8_t, uint64_t> > expected_key08;

    EXPECT_EQ(ids_in_which_is_key_08, expected_key08);
    EXPECT_EQ(ids_in_which_is_key_08.size(), expected_key08.size());
    for (size_t i = 0; i < ids_in_which_is_key_08.size(); i++) {
        EXPECT_EQ(ids_in_which_is_key_08[i], expected_key08[i]);
    }

    check_keys_min_max_before_trunk_construction(keys, bucket);

    bucket.construct_trunk();
    check_keys_min_max_after_trunk_construction(keys, bucket);
    ids_in_which_is_key_00 = bucket.query_key(key00);
    EXPECT_EQ(ids_in_which_is_key_00, expected_key00);
    EXPECT_EQ(ids_in_which_is_key_00.size(), expected_key00.size());
    for (size_t i = 0; i < ids_in_which_is_key_00.size(); i++) {
        EXPECT_EQ(ids_in_which_is_key_00[i], expected_key00[i]);
    }
    ids_in_which_is_key_08 = bucket.query_key(key08);
    EXPECT_EQ(ids_in_which_is_key_08, expected_key08);
    EXPECT_EQ(ids_in_which_is_key_08.size(), expected_key08.size());
    for (size_t i = 0; i < ids_in_which_is_key_08.size(); i++) {
        EXPECT_EQ(ids_in_which_is_key_08[i], expected_key08[i]);
    }

    bucket.construct_reverse_trunk();
    check_keys_min_max_after_both_trunk_construction(keys, bucket);
    ids_in_which_is_key_00 = bucket.query_key(key00);
    EXPECT_EQ(ids_in_which_is_key_00, expected_key00);
    EXPECT_EQ(ids_in_which_is_key_00.size(), expected_key00.size());
    for (size_t i = 0; i < ids_in_which_is_key_00.size(); i++) {
        EXPECT_EQ(ids_in_which_is_key_00[i], expected_key00[i]);
    }
    ids_in_which_is_key_08 = bucket.query_key(key08);
    EXPECT_EQ(ids_in_which_is_key_08, expected_key08);
    EXPECT_EQ(ids_in_which_is_key_08.size(), expected_key08.size());
    for (size_t i = 0; i < ids_in_which_is_key_08.size(); i++) {
        EXPECT_EQ(ids_in_which_is_key_08[i], expected_key08[i]);
    }
}

// small test to flex how beautiful everything is
TEST(pac_bucket, show_time) {
    const uint64_t nbCells = 10000;
    const uint number_hash_function = 1;
    const uint k = 31;
    const std::string prefix = "prefix";  // TODO
    const uint64_t nbBitsPerCell = 5;
    const bool write = true;

    Bucket<uint8_t> bucket(nbCells, number_hash_function, k, prefix, nbBitsPerCell, write);

    uint64_t kmer00 = 0;
    uint64_t kmer01 = 5;
    uint64_t kmer02 = 23;
    uint64_t kmer03 = 85;
    uint64_t kmer04 = 100;
    uint64_t kmer05 = 50000;
    uint64_t kmer06 = 60000;
    uint64_t kmer07 = 455225;
    uint64_t kmer08 = 465465;
    uint64_t kmer09 = 9654755;
    uint64_t kmer10 = 13541354;
    uint64_t kmer11 = 2451;

    const std::vector<uint64_t> kmers = {
        kmer00,
        kmer01,
        kmer02,
        kmer03,
        kmer04,
        kmer05,
        kmer06,
        kmer07,
        kmer08,
        kmer09,
        kmer10,
        kmer11};

    // insert 3 leafs
    bucket.add_leaf();
    bucket.add_leaf();
    bucket.add_leaf();

    // usage: bucket.insert_key(kmer, level, abundance);
    bucket.insert_key(kmers[0], 0, 12);
    bucket.insert_key(kmers[0], 1, 25);
    // bucket.insert_key(kmers[0], 2, 12);  // pas dans 2

    std::vector<std::tuple<uint8_t, uint64_t> > ids_in_which_is_key_00 = bucket.query_key(kmers[0]);
    std::vector<std::tuple<uint8_t, uint64_t> > expected_key00;
    // expecting to find kmers[0] at level 0 with abundance 12
    expected_key00.push_back({0, 12});
    // expecting to find kmers[0] at level 1 with abundance 25
    expected_key00.push_back({1, 25});
    EXPECT_EQ(ids_in_which_is_key_00, expected_key00);

    // kmers[0] is found in filters from 0 to 2 included
    EXPECT_EQ(bucket.query_key_min(kmers[0]), t(true, 0));
    EXPECT_EQ(bucket.query_key_max(kmers[0]), t(true, 2));

    // construct the "PAC" itself for the bucket
    bucket.construct_trunk();
    bucket.construct_reverse_trunk();

    // and now it is found from 0 to 1 included, allowing to skip queries on level 2
    EXPECT_EQ(bucket.query_key_min(kmers[0]), t(true, 0));
    EXPECT_EQ(bucket.query_key_max(kmers[0]), t(true, 1));
}

TEST_F(SimpleBucket, self_equality) {
    EXPECT_EQ(bucket, bucket);
    std::cout << "eq done" << std::endl;
}

TEST(pac_bucket, copy) {
    const uint64_t nbCells = 10000;
    const uint number_hash_function = 1;
    const uint k = 31;
    const std::string prefix = "prefix";  // TODO
    const uint64_t nbBitsPerCell = 5;
    const bool write = true;

    Bucket<uint8_t> bucket(nbCells, number_hash_function, k, prefix, nbBitsPerCell, write);
    Bucket<uint8_t> bucket_copy = bucket;
    EXPECT_EQ(bucket, bucket_copy);
}

// TEST_F(BucketFiveLeafs, serialize_and_load) {
//     bucket.construct_trunk();
//     Bucket<uint8_t> bucket_copy = bucket;

//     EXPECT_EQ(bucket, bucket_copy);

//     bucket.serialize();
//     EXPECT_EQ(bucket, bucket_copy);
//     bucket.load(5, true);
//     EXPECT_NE(bucket, bucket_copy);
//     bucket.construct_trunk();
//     // EXPECT_EQ(bucket, bucket_copy); // TODO
// }