#include <gtest/gtest.h>
#include <libpac/AggregatedBloom.h>

#include <libpac/CBF.hpp>

TEST(pac_aggregate_cbf, single_empty_aggregate) {
    // test that the aggregate filter is correctly set for one single empty Bloom filter
    uint number_hash_function = 1;
    uint64_t nbCell = 10000;
    uint64_t nbBitsPerCell = 10;

    CBF filter(nbCell, nbBitsPerCell, number_hash_function);

    AggregatedBloom<uint8_t> agg(nbCell, number_hash_function);
    const auto& [number_bit_set, number_bit_set_abt] = agg.aggregateFilter(5, &filter);

    // no bit should be set if we aggregate an empty filter
    EXPECT_EQ(number_bit_set, 0);
    EXPECT_EQ(number_bit_set_abt, 0);
}

TEST(pac_aggregate_cbf, single_aggregate) {
    // test that the aggregate filter is correctly set for one single (not empty) Bloom filter
    uint number_hash_function = 1;
    uint64_t nbCell = 10000;
    uint64_t nbBitsPerCell = 10;

    CBF filter(nbCell, nbBitsPerCell, number_hash_function);

    filter.insert_key(0, 1);
    filter.insert_key(1, 1);
    filter.insert_key(3, 1);
    filter.insert_key(4, 1);

    AggregatedBloom<uint8_t> agg(nbCell, number_hash_function);
    const auto& [number_bit_set, number_bit_set_abt] = agg.aggregateFilter(5, &filter, nbBitsPerCell);

    // 4 insertions => 4 bit changes (assuming no collision)
    EXPECT_EQ(number_bit_set, 4);
    EXPECT_EQ(number_bit_set_abt, 4);
}

TEST(pac_aggregate_cbf, two_aggregate) {
    // test that the aggregate filter  is correctly set for one single Bloom filter
    uint number_hash_function = 1;
    uint64_t nbCell = 100000;
    uint64_t nbBitsPerCell = 10;

    CBF filter0(nbCell, nbBitsPerCell, number_hash_function);
    filter0.insert_key(0, 2);
    filter0.insert_key(1, 2);
    filter0.insert_key(3, 5);
    filter0.insert_key(4, 1);

    CBF filter1(nbCell, nbBitsPerCell, number_hash_function);
    filter1.insert_key(0, 8);
    filter1.insert_key(1, 1);
    filter1.insert_key(5, 8);  // not 3, 5 here
    filter1.insert_key(4, 1);

    AggregatedBloom<uint8_t> agg(nbCell, number_hash_function);
    const auto& [number_bit_set0, number_bit_set_abt0] = agg.aggregateFilter(5, &filter0, nbBitsPerCell);
    const auto& [number_bit_set1, number_bit_set_abt1] = agg.aggregateFilter(4, &filter1, nbBitsPerCell);

    // 4 insertions => 4 bit changes (assuming no collision)
    EXPECT_EQ(number_bit_set0, 4);
    EXPECT_EQ(number_bit_set_abt0, 4);
    // 4 insertions too
    EXPECT_EQ(number_bit_set1, 4);
    // but only one change this round
    EXPECT_EQ(number_bit_set_abt1, 1);

    // check values in the aggregate filter
    // 0, 1, 3, 4, 5 should be at 5
    // 5 should be at 4
    // 2 should still be at 0
    EXPECT_EQ(agg.check_key(0), 5);
    EXPECT_EQ(agg.check_key(1), 5);
    EXPECT_EQ(agg.check_key(2), 0);
    EXPECT_EQ(agg.check_key(3), 5);
    EXPECT_EQ(agg.check_key(4), 5);
    EXPECT_EQ(agg.check_key(5), 4);
}