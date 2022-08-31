#pragma once

#include <string>
#include <tuple>
#include <zstr.hpp>

namespace libpac {
namespace generators {
namespace internal {
class End_iterator {};

class Iterator {
   public:
    Iterator(const std::string& filename);
    const std::tuple<const std::string&, const uint64_t&> operator*() const;
    Iterator& operator++();
    bool operator!=(End_iterator) const;

   private:
    const std::string& _filename;
    std::ifstream _myFileGz;
    zstr::istream _myFile;
    bool _over;
    std::string _line;
    std::string _kmer;
    uint64_t _abundance;
    void next();
};
}  // namespace internal

class KMCReader {
   private:
    const std::string& _filename;

   public:
    explicit KMCReader(const std::string& filename);
    libpac::generators::internal::Iterator begin();
    libpac::generators::internal::End_iterator end();
};
}  // namespace generators
}  // namespace libpac
