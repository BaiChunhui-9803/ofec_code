/**
 * @author  Gavin
 * @date    2022/2/10
 * @Email   gavinsun0921@foxmail.com
 */

#include <iostream>
#include <fstream>
#include "Utils.h"
//#include "Config.h"

#if defined(WIN64) || defined(_WIN64) || defined(WIN32) || defined(_WIN32)
    #include <direct.h>
#else
    #include <unistd.h>
#endif

namespace ofec {
    std::vector<Real> Utils::getFileData(const std::string& filename) {
        std::vector<Real> res;

        std::ifstream in(filename);
        if (in.is_open()) {
            while (!in.eof()) {
                Real val;
                in >> val;
                res.push_back(val);
            }
            in.close();
        }
        else {
            // Data file not found
            std::cout << "[ERROR], " << filename.c_str() << "not found." << std::endl;;
            // Prompt the correct storage path of data file
            char path[256];
            getcwd(path, sizeof(path));
            std::cout << "Please check the path " << filename.c_str() << " is relative to " << path << std::endl;;
            exit(1);
        }

        return res;
    }

    std::vector<Sample> Utils::getTrainData(const std::string& filename, size_t feature_dim, size_t label_dim) {
        std::vector<Sample> trainDataSet;

        std::vector<Real> buffer = getFileData(filename);

        for (size_t i = 0; i < buffer.size(); i += (feature_dim + label_dim)) {
            
            std::vector<Real> feature;
            std::vector<Real> label;
            // Read in training sample 'feature'
            for (size_t j = 0; j < feature_dim; ++j)
                feature.push_back(buffer[i + j]);
            // Read in training sample 'label'
            for (size_t k = 0; k < label_dim; ++k)
                label.push_back(buffer[i + feature_dim + k]);
            // Add samples to the 'trainDataSet'
            Sample trainSample(feature,label);
            trainDataSet.emplace_back(trainSample);
        }
        return trainDataSet;
    }

    std::vector<Sample> Utils::getTestData(const std::string& filename, size_t feature_dim) {
        std::vector<Sample> testDataSet;

        std::vector<Real> buffer = getFileData(filename);

        for (size_t i = 0; i < buffer.size(); i += feature_dim) {
            Sample testSample;
            // Read in test sample 'feature'
            std::vector<Real> feature;
            for (size_t j = 0; j < feature_dim; ++j)
                feature.push_back(buffer[i + j]);
            // Add samples to the 'testDataSet'
            testSample.setSampleFeatures(feature);
            testDataSet.push_back(testSample);
        }

        return testDataSet;
    }
}

