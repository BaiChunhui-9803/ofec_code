//
// Created by Administrator on 25-2-14.
//

#include "example_1.h"

#include "../../core/parameter/variants_to_stream.h"
#include "../../instance/problem/continuous/single_objective/global/cec2005/f1_shifted_sphere.h"
#include "../../instance/problem/continuous/single_objective/global/cec2005/f2_shifted_schwefel_1_2.h"
#include "../../instance/problem/combination/travelling_salesman/travelling_salesman.h"
#include "../../instance/algorithm/continuous/single_objective/global/canonical_de/canonical_de.h"
#include "../../instance/algorithm/combination/eax_tsp/eax_tsp_origin/eax_tsp_alg.h"


namespace ofec {

    void print_(const std::list<std::pair<ParameterVariant, int>>& param_stream) {
        if (param_stream.empty()) {
            std::cout << "null" << std::endl;
        } else {
            for (const auto& pair : param_stream) {
                if (std::holds_alternative<bool>(pair.first)) {
                    std::cout << std::get<bool>(pair.first);
                } else if (std::holds_alternative<int>(pair.first)) {
                    std::cout << std::get<int>(pair.first);
                } else if (std::holds_alternative<size_t>(pair.first)) {
                    std::cout << std::get<size_t>(pair.first);
                } else if (std::holds_alternative<char>(pair.first)) {
                    std::cout << std::get<char>(pair.first);
                } else if (std::holds_alternative<Real>(pair.first)) {
                    std::cout << std::get<Real>(pair.first);
                } else if (std::holds_alternative<std::vector<bool>>(pair.first)) {
                    auto& param_vec = std::get<std::vector<bool>>(pair.first);
                    std::cout << "\t" << param_vec.size();
                    for (const auto& it : param_vec) {
                        std::cout << "\t" << it;
                    }
                } else if (std::holds_alternative<std::vector<int>>(pair.first)) {
                    auto& param_vec = std::get<std::vector<int>>(pair.first);
                    std::cout << "\t" << param_vec.size();
                    for (const auto& it : param_vec) {
                        std::cout << "\t" << it;
                    }
                } else if (std::holds_alternative<std::vector<size_t>>(pair.first)) {
                    auto& param_vec = std::get<std::vector<size_t>>(pair.first);
                    std::cout << "\t" << param_vec.size();
                    for (const auto& it : param_vec) {
                        std::cout << "\t" << it;
                    }
                } else if (std::holds_alternative<std::string>(pair.first)) {
                    std::cout << std::get<std::string>(pair.first);
                } else if (std::holds_alternative<std::vector<Real>>(pair.first)) {
                    auto& param_vec = std::get<std::vector<Real>>(pair.first);
                    std::cout << "\t" << param_vec.size();
                    for (const auto& it : param_vec) {
                        std::cout << "\t" << it;
                    }
                } else if (std::holds_alternative<unsigned long>(pair.first)) {
                    std::cout << std::get<unsigned long>(pair.first);
                }
            }
        }
    }

    void printProblemArchivedParams(const Problem* problem) {
        ofec::ParameterVariantStream paramsStream;
        ofec::ParameterMap param_map = problem->archivedParameters();
        for (const auto& pair : param_map) {
            paramsStream << "Key: " << pair.first << ", Value: " << pair.second << "\n";
        }
        print_(paramsStream.getStream());
        // ofec::variants_stream::outputToFile(paramsStream, "test_1.txt");
    }

    void printAlgorithmArchivedParams(const Algorithm* algorithm) {
        ofec::ParameterVariantStream paramsStream;
        ofec::ParameterMap param_map = algorithm->archivedParameters();
        for (const auto& pair : param_map) {
            paramsStream << "Key: " << pair.first << ", Value: " << pair.second << "\n";
        }
        print_(paramsStream.getStream());
        // ofec::variants_stream::outputToFile(paramsStream, "test_2.txt");
    }

    void run(int argc, char *argv[]) {
        ofec::g_working_directory = "/data/Share/Student/2018/DiaoYiya/ofec_data";
        std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";

        ProblemConfig problem1 = {
            "GOP_CEC2005_F01",
            {
                {"number of variables", "2"},
            }
        };
        ProblemConfig problem2 = {
            "GOP_CEC2005_F02",
            {
                {"number of variables", "10"},
            }
        };
        ProblemConfig problem3 = {
            "TSP",
            {
                {"dataFile1", "A280.tsp"},
                {"dataDirectory", dir2},
            },
            0.3
        };
        AlgorithmConfig algorithm1 = {
            "Canonical-DE",
            {
                {"population size", "10"},
                {"maximum evaluations", "10000"}
            },
            0.5
        };
        AlgorithmConfig algorithm2 = {
            "EAX-TSP",
            {
                {"maximum evaluations", "5000"}},
            0.4
        };
        const Configuration config = {
            {problem1, problem2, problem3},
            {algorithm1, algorithm2}
        };

        registerInstance();

        // Initialize envs manually
        std::vector<std::shared_ptr<Environment>> envs;
        {
            for (const auto& prob_conf: config.problemVec) {
                std::shared_ptr<Environment> env(generateEnvironmentByFactory(prob_conf.name));
                env->recordInputParameters();
                env->initialize();
                env->setProblem(generateProblemByFactory(prob_conf.name));
                for (const auto& prob_param_conf: prob_conf.parameters) {
                    env->problem()->inputParameters().at(prob_param_conf.first)->setValue(prob_param_conf.second);
                }

                std::cout << "[Print archived_parameters before recording input_parameters for problem " << prob_conf.name << "]" << std::endl;
                printProblemArchivedParams(env->problem());
                env->problem()->recordInputParameters();
                std::cout << "[Print archived_parameters after recording input_parameters for problem " << prob_conf.name << std::endl;
                printProblemArchivedParams(env->problem());
                env->initializeProblem(prob_conf.seed);

                AlgorithmConfig alg_config = (prob_conf.name == "TSP") ? config.algorithmVec[1] : config.algorithmVec[0];

                if (checkValidation(prob_conf.name, alg_config.name)) {
                    env->setAlgorithm(generateAlgorithmByFactory(alg_config.name));
                    for (const auto& alg_param_conf: alg_config.parameters) {
                        env->algorithm()->inputParameters().at(alg_param_conf.first)->setValue(alg_param_conf.second);
                    }

                    std::cout << "[Print archived_parameters before recording input_parameters for algorithm " << alg_config.name << "]" << std::endl;
                    printAlgorithmArchivedParams(env->algorithm());
                    env->algorithm()->recordInputParameters();
                    std::cout << "[Print archived_parameters after recording input_parameters for algorithm " << alg_config.name << "]" << std::endl;
                    printAlgorithmArchivedParams(env->algorithm());
                    env->initializeAlgorithm(alg_config.seed);
                }
                envs.push_back(env);
            }
        }

        envs.clear();
        std::cout << std::endl;

        // Initialize envs from g_algorithm_for_environment_problem - the mapping of algorithms to environment-problem pairs
        {
            for (const auto& env_pair : g_algorithm_for_environment_problem) {
                std::string prob_name = env_pair.first;
                std::cout << "Problem: " << prob_name << std::endl;
                std::cout << "Algorithm: ";
                for (const auto& alg_name : env_pair.second) {
                    std::cout << alg_name << std::endl;

                    std::shared_ptr<Environment> env(generateEnvironmentByFactory(prob_name));
                    env->recordInputParameters();
                    env->initialize();
                    env->setProblem(generateProblemByFactory(prob_name));
                    for (const auto& prob_param_conf: config.problemMap.at(prob_name).parameters) {
                        env->problem()->inputParameters().at(prob_param_conf.first)->setValue(prob_param_conf.second);
                    }
                    env->problem()->recordInputParameters();
                    printProblemArchivedParams(env->problem());
                    env->initializeProblem(config.problemMap.at(prob_name).seed);

                    env->setAlgorithm(generateAlgorithmByFactory(alg_name));
                    for (const auto& alg_param_conf: config.algorithmMap.at(alg_name).parameters) {
                        env->algorithm()->inputParameters().at(alg_param_conf.first)->setValue(alg_param_conf.second);
                    }
                    env->algorithm()->recordInputParameters();
                    printAlgorithmArchivedParams(env->algorithm());
                    env->initializeAlgorithm(config.algorithmMap.at(alg_name).seed);
                    envs.push_back(env);
                }
            }
        }

        std::cout << std::endl;
        envs[1]->algorithm()->inputParameters().at("population size")->setValue("20");
        envs[1]->algorithm()->recordInputParameters();
        std::cout << "[Print archived_parameters for envs[2] after updating param \"population size\" from 10 to 20 ]" << std::endl;
        printAlgorithmArchivedParams(envs[1]->algorithm());

        for (const auto& env:envs) {
            env->runAlgorithm();
        }


    }
}