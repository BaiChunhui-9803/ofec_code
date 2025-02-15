//
// Created by Administrator on 25-2-14.
//

#include "example_1.h"
#include "../../instance/problem/continuous/single_objective/global/cec2005/f1_shifted_sphere.h"
#include "../../instance/problem/continuous/single_objective/global/cec2005/f2_shifted_schwefel_1_2.h"
#include "../../instance/algorithm/continuous/single_objective/global/canonical_de/canonical_de.h"


namespace ofec {
    void run(int argc, char *argv[]) {
        ofec::g_working_directory = "/data/Share/Student/2018/DiaoYiya/ofec_data";

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
        AlgorithmConfig algorithm = {
            "CanonicalDE",
            {
                {"population size", "10"},
                {"maximum evaluations", "100000"}
            }
        };
        const Configuration config = {{problem1, problem2}, algorithm};

        std::vector<std::shared_ptr<Environment>> envs;
        for (const auto& conf: config.problems) {
            std::shared_ptr<Environment> env(Environment::create());
            env->recordInputParameters();
            env->initialize();

            if (conf.name == "GOP_CEC2005_F01") {
                env->setProblem(GOP_CEC2005_F01::create());
            } else if (conf.name == "GOP_CEC2005_F02") {
                env->setProblem(GOP_CEC2005_F02::create());
            } else {
                throw std::invalid_argument("Unknown problem name: " + conf.name);
            }
            for (const auto& param: conf.parameters) {
                env->problem()->inputParameters().at(param.first)->setValue(param.second);
            }
            env->problem()->recordInputParameters();
            env->initializeProblem();

            if (config.algorithm.name == "CanonicalDE") {
                env->setAlgorithm(CanonicalDE:: create());
            }
            for (const auto& param: config.algorithm.parameters) {
                env->algorithm()->inputParameters().at(param.first)->setValue(param.second);
            }
            env->algorithm()->recordInputParameters();
            env->initializeAlgorithm(0.5);
            envs.push_back(env);
        }

        for (const auto& env:envs) {
            env->runAlgorithm();
        }


    }
}