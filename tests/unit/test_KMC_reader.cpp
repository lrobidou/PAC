#include <gtest/gtest.h>

#include <libpac/generators/KMCReader.hpp>
inline bool fileExists(const std::string& name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

TEST(test_pac, simple_read) {
    const std::string filename = "../tests/unit/data/KMC.txt";
    EXPECT_EQ(fileExists(filename), true);
    std::vector<std::string> expected_kmers = {
        "dfgdenfeqvbr",
        "sdfghjklmùdfrgfedz",
        "ergfvegrfezgrfer''",
        "All in all it's just another brick in the wall",
        "All in all you're just another brick in the wall",
        "it",
        "should ",
        "work",
        "with",
        "utf8",
        "(＃`Д´)",
        "(`皿´＃)",
        "( ` ω ´ )",
        "ヽ( `д´*)ノ",
        "(・`ω´・)",
        "(`ー´)",
        "ヽ(`⌒´メ)ノ",
        "凸(`△´＃)",
        "( `ε´ )",
        "ψ( ` ∇ ´ )ψ",
        "ヾ(`ヘ´)ﾉﾞ",
        "ヽ(‵﹏´)ノ",
        "(ﾒ` ﾛ ´)",
        "(╬`益´)",
        "┌∩┐(◣_◢)┌∩┐",
        "凸( ` ﾛ ´ )凸",
        "Σ(▼□▼メ)",
        "(°ㅂ°╬)",
        "ψ(▼へ▼メ)～→",
        "(ノ°益°)ノ",
        "(҂ `з´ )",
        "(‡▼益▼)",
        "(҂` ﾛ ´)凸",
        "((╬◣﹏◢))",
        "٩(╬ʘ益ʘ╬)۶",
        "(╬ Ò﹏Ó)",
        "＼＼٩(๑`^´๑)۶／／",
        "(凸ಠ益ಠ)凸",
        "↑_(ΦwΦ)Ψ",
        "←~(Ψ▼ｰ▼)∈",
        "୧((#Φ益Φ#))୨",
        "٩(ఠ益ఠ)۶",
        "(ﾉಥ益ಥ)ﾉ",
        "(≖､≖╬)"};
    std::vector<uint64_t> expected_abundance = {
        0,
        0,
        0,
        5,
        5,
        7,
        9,
        9,
        9,
        9,
        0,
        9,
        2,
        9,
        163543,
        1,
        654,
        65341,
        987,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        10000000,
        0,
        0,
        0,
        0,
        0};
    std::vector<std::string> res_kmers;
    std::vector<uint64_t> res_abundance;
    res_kmers.reserve(expected_kmers.size());
    res_abundance.reserve(expected_abundance.size());
    for (const auto [kmer, abundance] : libpac::generators::KMCReader(filename)) {
        // std::cout << abundance << std::endl;
        res_kmers.push_back(kmer);
        res_abundance.push_back(abundance);
    }
    EXPECT_EQ(res_kmers, expected_kmers);
    EXPECT_EQ(res_abundance, expected_abundance);
}