#include <iostream>  // TODO remove
#include <stdexcept>
#include <vector>

namespace fimpera {

// search the minimum in every sliding windows from vector ARR.
// **Erases the ARR vector.**
// see fimpera's paper for explanation on how it works
std::vector<uint64_t> sliding_window_minimum(std::vector<uint64_t>& ARR, uint64_t w);

// apply fimpera on a vector of response from an AMQ
std::vector<uint64_t> fimpera(const std::vector<uint64_t>& input, const unsigned int& k, const unsigned int& z);
}  // namespace fimpera
