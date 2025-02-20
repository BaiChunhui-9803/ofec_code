/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
******************************************************************************************
*  Paper:
*******************************************************************************************/


#ifndef EXAMPLE_1_H
#define EXAMPLE_1_H

#include "../../core/global.h"
#include "../../run/interface.h"

struct ProblemConfig {
    std::string name;
    std::map<std::string, std::string> parameters;
    ofec::Real seed = 0;

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
    ofec::Real seed;

    void print() const {
        std::cout << "Algorithm: " << name << std::endl;
        std::cout << "Parameters:" << std::endl;
        for (const auto& [key, value] : parameters) {
            std::cout << "  " << key << ": " << value << std::endl;
        }
    }
};

struct Configuration {
    std::vector<ProblemConfig> problemVec;
    std::vector<AlgorithmConfig> algorithmVec;
    std::map<std::string, ProblemConfig> problemMap;
    std::map<std::string, AlgorithmConfig> algorithmMap;


    Configuration(std::vector<ProblemConfig> problems, std::vector<AlgorithmConfig> algorithms)
        : problemVec(std::move(problems)), algorithmVec(std::move(algorithms)) {
        for (auto& problem : problemVec) {
            problemMap[problem.name] = problem;
        }
        for (auto& algorithm : algorithmVec) {
            algorithmMap[algorithm.name] = algorithm;
        }
    }

    void print() const {
        std::cout << "Configuration:" << std::endl;
        for (const auto& problem : problemVec) {
            problem.print();
            std::cout << std::endl;
        }
        for (const auto& algorithm : algorithmVec) {
            algorithm.print();
            std::cout << std::endl;
        }
    }

    size_t getProblemCount() const {
        return problemVec.size();
    }

    size_t getAlgorithmCount() const {
        return algorithmVec.size();
    }
};

namespace ofec {

    void run_example_1(int argc, char *argv[]);
}

#endif //EXAMPLE_1_H
