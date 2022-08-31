/**
 * Is this file even needed now ?
 * A CBF can be just like a bloom filter
 * Is there a need to keep the AggregatedBloom compatible with Bloom filters ?
 */
// #include <gtest/gtest.h>
// #include <libpac/AggregatedBloom.h>
// #include <libpac/Bloom.h>

// TEST(pac_aggregate, single_empty_aggregate) {
//     // test that the aggregate filter is correctly set for one single empty Bloom filter
//     uint64_t size = 1000;
//     uint number_hash_function = 1;
//     Bloom filter(size, number_hash_function);
//     AggregatedBloom<uint8_t> agg(size, number_hash_function);
//     const auto& [number_bit_set, number_bit_set_abt] = agg.aggregateFilter(5, &filter);

//     // no bit should be set if we aggregate an empty filter
//     EXPECT_EQ(number_bit_set, 0);
//     EXPECT_EQ(number_bit_set_abt, 0);
// }

// TEST(pac_aggregate, single_aggregate) {
//     // test that the aggregate filter is correctly set for one single (not empty) Bloom filter
//     uint64_t size = 1000;
//     uint number_hash_function = 1;
//     Bloom filter(size, number_hash_function);
//     filter.insert_key(0);
//     filter.insert_key(1);
//     filter.insert_key(3);
//     filter.insert_key(4);

//     AggregatedBloom<uint8_t> agg(size, number_hash_function);
//     const auto& [number_bit_set, number_bit_set_abt] = agg.aggregateFilter(5, &filter);

//     // 4 insertions => 4 bit changes (assuming no collision)
//     EXPECT_EQ(number_bit_set, 4);
//     EXPECT_EQ(number_bit_set_abt, 4);
// }

// TEST(pac_aggregate, two_aggregate) {
//     // test that the aggregate filter  is correctly set for one single Bloom filter
//     uint64_t size = 1000;
//     uint number_hash_function = 1;
//     Bloom filter0(size, number_hash_function);
//     filter0.insert_key(0);
//     filter0.insert_key(1);
//     filter0.insert_key(3);
//     filter0.insert_key(4);

//     Bloom filter1(size, number_hash_function);
//     filter1.insert_key(0);
//     filter1.insert_key(1);
//     filter1.insert_key(5);  // not 3, 5 here
//     filter1.insert_key(4);

//     AggregatedBloom<uint8_t> agg(size, number_hash_function);
//     const auto& [number_bit_set0, number_bit_set_abt0] = agg.aggregateFilter(5, &filter0);
//     const auto& [number_bit_set1, number_bit_set_abt1] = agg.aggregateFilter(4, &filter1);

//     // 4 insertions => 4 bit changes (assuming no collision)
//     EXPECT_EQ(number_bit_set0, 4);
//     EXPECT_EQ(number_bit_set_abt0, 4);
//     // 4 insertions too
//     EXPECT_EQ(number_bit_set1, 4);
//     // but only one change this round
//     EXPECT_EQ(number_bit_set_abt1, 1);

//     // check values in the aggregate filter
//     // 0, 1, 3, 4, 5 should be at 5
//     // 5 should be at 4
//     // 2 should still be at 0
//     EXPECT_EQ(agg.check_key(0), 5);
//     EXPECT_EQ(agg.check_key(1), 5);
//     EXPECT_EQ(agg.check_key(2), 0);
//     EXPECT_EQ(agg.check_key(3), 5);
//     EXPECT_EQ(agg.check_key(4), 5);
//     EXPECT_EQ(agg.check_key(5), 4);
// }