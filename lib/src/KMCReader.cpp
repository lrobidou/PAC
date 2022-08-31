#include <libpac/generators/KMCReader.hpp>

libpac::generators::internal::Iterator::Iterator(const std::string& filename) : _filename(filename), _myFileGz(filename), _myFile(_myFileGz), _over(false) {
    next();
}

const std::tuple<const std::string&, const uint64_t&> libpac::generators::internal::Iterator::operator*() const {
    return {_kmer, _abundance};
}

libpac::generators::internal::Iterator& libpac::generators::internal::Iterator::operator++() {
    next();
    return *this;
}

bool libpac::generators::internal::Iterator::operator!=(End_iterator) const {
    return !_over;
}

void libpac::generators::internal::Iterator::next() {
    if (std::getline(_myFile, _line)) {
        std::string abundanceStr;
        std::string tmp_kmer;
        std::stringstream linestream(_line);
        std::getline(linestream, tmp_kmer, '\t');
        std::getline(linestream, abundanceStr, '\t');

        try {
            _abundance = std::stoll(abundanceStr);
            _kmer = tmp_kmer;
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument when parsing line :\"" << _line << "\"" << std::endl;
        }
    } else {
        _over = true;
    }
}

libpac::generators::KMCReader::KMCReader(const std::string& filename) : _filename(filename) {
}

libpac::generators::internal::Iterator libpac::generators::KMCReader::begin() {
    return internal::Iterator{_filename};
}

libpac::generators::internal::End_iterator libpac::generators::KMCReader::end() {
    return libpac::generators::internal::End_iterator{};
}