#include "maxElementWithNan.h"
#include <algorithm>
#include <limits>
#include <cmath>

// Implement the function in the source file.
double maxElementWithNan(const std::vector<double>& vec) {
    for (auto elem : vec) {
        if (std::isnan(elem)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    return *std::max_element(vec.begin(), vec.end());
}
