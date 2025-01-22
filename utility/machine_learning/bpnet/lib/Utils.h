/**
 * @author  Gavin
 * @date    2022/2/10
 * @Email   gavinsun0921@foxmail.com
 */

#ifndef OFEC_UTILS_H
#define OFEC_UTILS_H

#include <cmath>
#include "Net.h"


namespace ofec {
    namespace Utils {
        static Real sigmoid(Real x) {
            return 1.0 / (1.0 + std::exp(-x));
        }

        std::vector<Real> getFileData(const std::string& filename);

        std::vector<Sample> getTrainData(const std::string& filename, size_t feature_dim, size_t label_dim);

        std::vector<Sample> getTestData(const std::string& filename, size_t feature_dim);
    }
}

#endif // !OFEC_UTILS_H

