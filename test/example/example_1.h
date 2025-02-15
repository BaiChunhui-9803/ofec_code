//
// Created by Administrator on 25-2-14.
//

#ifndef EXAMPLE_1_H
#define EXAMPLE_1_H

#include "../../core/global.h"
#include "../../run/interface.h"

struct ProblemConfig {
    std::string name;
    std::map<std::string, std::string> parameters;

    void print() const {
        std::cout << "Problem: " << name << std::endl;
        std::cout << "Parameters:" << std::endl;
        for (const auto& [key, value] : parameters) {
            std::cout << "  " << key << ": " << value << std::endl;
        }
    }
};

struct AlgorithmConfig {
    std::string name;
    std::map<std::string, std::string> parameters;

    void print() const {
        std::cout << "Algorithm: " << name << std::endl;
        std::cout << "Parameters:" << std::endl;
        for (const auto& [key, value] : parameters) {
            std::cout << "  " << key << ": " << value << std::endl;
        }
    }
};

struct Configuration {
    std::vector<ProblemConfig> problems;
    AlgorithmConfig algorithm;

    Configuration(std::vector<ProblemConfig> problems, AlgorithmConfig algorithm)
        : problems(std::move(problems)), algorithm(std::move(algorithm)) {}

    void print() const {
        std::cout << "Configuration:" << std::endl;
        for (const auto& problem : problems) {
            problem.print();
            std::cout << std::endl;
        }
        algorithm.print();
    }

    size_t getProblemCount() const {
        return problems.size();
    }
};

namespace ofec {
    void run(int argc, char *argv[]);
}

#endif //EXAMPLE_1_H
