#include <gtest/gtest.h>
#include <libpac/Bloom.h>

TEST(pac_bloom, simple_retrieval_test) {
    // basic test to check if the Bloom filter works
    uint number_hash_function = 1;
    uint64_t size = 1000;

    Bloom bloom = Bloom(size, number_hash_function);

    bloom.insert_key(0);
    bloom.insert_key(1);
    bloom.insert_key(3);
    bloom.insert_key(4);

    EXPECT_TRUE(bloom.check_key(0));
    EXPECT_TRUE(bloom.check_key(1));
    EXPECT_FALSE(bloom.check_key(2));
    EXPECT_TRUE(bloom.check_key(3));
    EXPECT_TRUE(bloom.check_key(4));
}

TEST(pac_bloom, simple_retrieval_test_multiple_hash_function) {
    // basic test to check if the Bloom filter works with mutiple hash functions
    uint number_hash_function = 6;
    uint64_t size = 1000;

    Bloom bloom = Bloom(size, number_hash_function);

    bloom.insert_key(0);
    bloom.insert_key(1);
    bloom.insert_key(3);
    bloom.insert_key(4);

    EXPECT_TRUE(bloom.check_key(0));
    EXPECT_TRUE(bloom.check_key(1));
    EXPECT_FALSE(bloom.check_key(2));
    EXPECT_TRUE(bloom.check_key(3));
    EXPECT_TRUE(bloom.check_key(4));
}

TEST(pac_bloom, equality_same) {
    // check if a Bloom is equal to itself
    uint number_hash_function = 1;
    uint64_t size = 1000;

    Bloom bloom = Bloom(size, number_hash_function);
    bloom.insert_key(0);
    bloom.insert_key(1);
    bloom.insert_key(3);
    bloom.insert_key(4);

    EXPECT_EQ(bloom, bloom);
}

TEST(pac_bloom, equality_copy) {
    // copy a filter and checks wether the filter is equal to its copy
    uint number_hash_function = 1;
    uint64_t size = 1000;

    Bloom bloom = Bloom(size, number_hash_function);

    bloom.insert_key(0);
    bloom.insert_key(1);
    bloom.insert_key(3);
    bloom.insert_key(4);

    Bloom bloom2(bloom);

    EXPECT_EQ(bloom, bloom2);
}

TEST(pac_bloom, equality_same_data) {
    uint number_hash_function = 1;
    uint64_t size = 1000;

    Bloom bloom = Bloom(size, number_hash_function);
    bloom.insert_key(0);
    bloom.insert_key(1);
    bloom.insert_key(3);
    bloom.insert_key(4);

    Bloom bloom2 = Bloom(size, number_hash_function);
    bloom2.insert_key(0);
    bloom2.insert_key(1);
    bloom2.insert_key(3);
    bloom2.insert_key(4);

    EXPECT_EQ(bloom, bloom2);
}

TEST(pac_bloom, equality_different_data) {
    uint number_hash_function = 1;
    uint64_t size = 1000;

    Bloom bloom = Bloom(size, number_hash_function);
    bloom.insert_key(0);
    bloom.insert_key(1);

    Bloom bloom2 = Bloom(size, number_hash_function);
    bloom2.insert_key(3);
    bloom2.insert_key(4);

    EXPECT_NE(bloom, bloom2);
}

TEST(pac_bloom, equality_different_number_hash_function) {
    uint number_hash_function = 1;
    uint64_t size = 1000;

    Bloom bloom = Bloom(size, number_hash_function);
    bloom.insert_key(0);
    bloom.insert_key(1);
    bloom.insert_key(3);
    bloom.insert_key(4);

    Bloom bloom2 = Bloom(size, number_hash_function + 1);
    bloom2.insert_key(0);
    bloom2.insert_key(1);
    bloom2.insert_key(3);
    bloom2.insert_key(4);

    EXPECT_NE(bloom, bloom2);
}

TEST(pac_bloom, equality_different_size) {
    uint number_hash_function = 1;
    uint64_t size = 1000;

    Bloom bloom = Bloom(size, number_hash_function);
    bloom.insert_key(0);
    bloom.insert_key(1);
    bloom.insert_key(3);
    bloom.insert_key(4);

    Bloom bloom2 = Bloom(size + 1, number_hash_function);
    bloom2.insert_key(0);
    bloom2.insert_key(1);
    bloom2.insert_key(3);
    bloom2.insert_key(4);

    EXPECT_NE(bloom, bloom2);
}

TEST(pac_bloom, dump_load) {
    // create a filter, fill it, dump it and load it
    uint number_hash_function = 1;
    uint64_t size = 1000;

    zstr::ofstream out("test", std::ios::trunc);
    bm::serializer<bm::bvector<> > bvs;
    bvs.byte_order_serialization(false);
    bvs.gap_length_serialization(false);

    Bloom bloom = Bloom(size, number_hash_function);
    bloom.insert_key(0);
    bloom.insert_key(1);
    bloom.insert_key(3);
    bloom.insert_key(4);
    bloom.dump_disk(bvs, &out);

    zstr::ifstream in("test");

    Bloom bloom2 = Bloom(size, number_hash_function);
    bloom2.load_disk(&in);

    EXPECT_EQ(bloom, bloom2);
}

TEST(pac_bloom, get_size) {
    // create a filter, fill it, dump it and load it
    uint number_hash_function = 1;
    uint64_t size = 1000;

    Bloom bloom = Bloom(size, number_hash_function);

    EXPECT_EQ(bloom.get_filter_size(), size + 1);  // lrobidou: why +1 ?
}

// TODO: lrobidou: broken
// TEST(pac_bloom, free_ram) {
//     // create a filter, fill it, dump it and load it
//     uint number_hash_function = 1;
//     uint64_t size = 1000;

//     Bloom<uint8_t> bloom = Bloom<uint8_t>(size, number_hash_function);
//     bloom.free_ram();
//     EXPECT_EQ(bloom.get_filter_size(), 0);
// }