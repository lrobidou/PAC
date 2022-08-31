#include <libpac/fimpera.hpp>

std::vector<uint64_t> fimpera::sliding_window_minimum(std::vector<uint64_t>& ARR, uint64_t w) {
    // basically, the minimum in the sliding window is the minimum of two "offseted" vectors (min_left and min_right)
    // but:
    // min_left is computed on the fly (we only store its required element)
    // min_right is stored in the vector passed as parameter
    // the response is computed in the ARR vector directly to avoid memory allocation.
    // OPTIMIZE throw errors instead of returning empty response
    if (ARR.size() < w) {
        return {};
    }
    if (w == 0) {
        return ARR;
    }
    if (w == 1) {
        return ARR;
    }

    uint64_t nbWin = ARR.size() / w;
    uint64_t nb_elem_last_window = ARR.size() % w;

    uint64_t min_left = ARR[0];
    for (uint64_t i = 0; i < w; i++) {
        min_left = std::min(min_left, ARR[i]);
    }

    for (uint64_t i = 0; i < nbWin - 1; i++) {
        uint64_t start_window = i * w;
        for (uint64_t indice = start_window + w - 2; indice >= start_window; indice--) {
            // we compute "min_right" here, directly in ARR vector
            ARR[indice] = std::min(ARR[indice + 1], ARR[indice]);
        }

        for (uint64_t j = 0; j < w; j++) {
            ARR[start_window + j] = std::min(ARR[start_window + j], min_left);
            min_left = (j == 0) ? ARR[start_window + w + j] : std::min(min_left, ARR[start_window + w + j]);
        }
    }

    // last window

    // compute min_right for last window
    uint64_t start_window = (nbWin - 1) * w;

    // This code is buggy (consider the case start_window = 0)
    // for (uint64_t indice = start_window + w - 2; indice >= start_window; indice--) {
    //     ARR[indice] = std::min(ARR[indice + 1], ARR[indice]);
    // }
    // so let's do it like that
    for (uint64_t offset = w - 2 + 1; offset > 0; offset--) {
        uint64_t indice = start_window + offset - 1;
        ARR[indice] = std::min(ARR[indice + 1], ARR[indice]);
    }

    // compute the min for the last window
    for (uint64_t j = 0; j < nb_elem_last_window; j++) {
        ARR[start_window + j] = std::min(ARR[start_window + j], min_left);
        min_left = (j == 0) ? ARR[start_window + w + j] : std::min(min_left, ARR[start_window + w + j]);
    }

    ARR[start_window + nb_elem_last_window] = std::min(ARR[start_window + nb_elem_last_window], min_left);
    for (std::size_t i = 0; i < w - 1; i++) {
        ARR.pop_back();
    }

    return ARR;
}

std::vector<uint64_t> fimpera::fimpera(const std::vector<uint64_t>& input, const unsigned int& k, const unsigned int& z) {
    // if (input.size() < 1) {
    // throw std::range_error("`input.length()` (which is " + std::to_string(input.size()) + ") < 1(which is " + std::to_string(k) + ")");
    // }

    const unsigned int s = k - z;
    unsigned long long size = input.size();
    std::vector<uint64_t> response(size - z, 0);
    unsigned long long stretchLength = 0;  // number of consecutive positives kmers
    unsigned long long j = 0;              // index of the input vector
    bool extending_stretch = true;
    uint64_t amq_answer = 0;
    std::vector<uint64_t> previous_answers;
    previous_answers.reserve(response.size());
    while (j < size) {
        if (stretchLength != previous_answers.size()) {
            // TODO : placed here for debug. Remove once we are sure everything is OK
            throw std::domain_error("`stretchLength != previous_answers.size()`. This indicates a logic violation in fimpera. Can not recover from that: contact the autor.");
        }
        if (amq_answer = input[j]) {
            if (extending_stretch) {
                previous_answers.push_back(amq_answer);
                stretchLength++;
                j++;
            } else {
                extending_stretch = true;
                j = j - z;
            }
        } else {
            if (stretchLength > z) {
                unsigned long long t = j - stretchLength;
                for (const auto& minimum : sliding_window_minimum(previous_answers, z + 1)) {
                    response[t++] = minimum;
                }
            }
            previous_answers.clear();
            stretchLength = 0;
            extending_stretch = false;
            j = j + z + 1;
        }
    }

    // Last values:
    if (stretchLength > z) {
        unsigned long long t = size - s + 1 - stretchLength;
        for (const auto& minimum : sliding_window_minimum(previous_answers, z + 1)) {
            response[t++] = minimum;
        }
    }

    return response;
}
