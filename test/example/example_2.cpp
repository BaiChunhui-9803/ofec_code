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


#include "example_2.h"
#include <filesystem>

#include "../../core/parameter/variants_to_stream.h"

namespace ofec {

    void run_example_2(int argc, char *argv[]) {
        ofec::g_working_directory = "/data/Share/Student/2018/DiaoYiya/ofec_data";
        std::string offline_data_root_path = "../instance/problem/realworld/game";

        ProblemConfig problem = {
            "Game",
            {
                        {"rootPath", offline_data_root_path},
                        {"sceGroup", "hrl_imcbs"},
                        {"sceName", "sce_2m"},
                        {"sceType", "train"},
                    },
                0.5
        };

        registerInstance();

        std::shared_ptr<Environment> env(GameEnvironment::create());
        env->inputParameters().at("changeFre")->setValue("500");
        env->recordInputParameters();
        env->initialize();
        env->setProblem(generateProblemByFactory(problem.name));
        for (const auto& prob_param_conf: problem.parameters) {
            env->problem()->inputParameters().at(prob_param_conf.first)->setValue(prob_param_conf.second);
        }
        env->problem()->recordInputParameters();
        env->initializeProblem(problem.seed);

        if (auto game_prob = dynamic_cast<Game*>(env->problem())) {
            std::cout << "[Before problem change]" << std::endl;
            std::cout << " - number of solution records: " << game_prob->solutionRecords().size() << std::endl;
            std::cout << " - number of optima: " << game_prob->optima()->numberObjectives() << std::endl;
            std::cout << " - optima value: " << game_prob->optima()->objective(0).at(0) << std::endl;
            game_prob->change(new Random(0.201));
            std::cout << "[Problem changed, before optima update]" << std::endl;
            std::cout << " - number of solution records: " << game_prob->solutionRecords().size() << std::endl;
            std::cout << " - number of optima: " << game_prob->optima()->numberObjectives() << std::endl;
            std::cout << " - optima value: " << game_prob->optima()->objective(0).at(0) << std::endl;
            game_prob->updateOptEnv(env.get());
            std::cout << "[Problem changed and optima updated]" << std::endl;
            std::cout << " - number of solution records: " << game_prob->solutionRecords().size() << std::endl;
            std::cout << " - number of optima: " << game_prob->optima()->numberObjectives() << std::endl;
            std::cout << " - optima value: " << game_prob->optima()->objective(0).at(0) << std::endl;

            ofec::ParameterVariantStream paramsStream;

            {
                paramsStream << game_prob->optima()->numberSolutions();
                for (size_t i = 0; i < game_prob->optima()->numberSolutions(); ++i) {
                    auto sol = game_prob->optima()->solution(i);
                    game_prob->solutionToParameterVariants(sol, paramsStream);
                }
                ofec::variants_stream::outputToFile(paramsStream, offline_data_root_path + "/test" + "_optima.txt");
            }

            {
                ofec::variants_stream::inputFromFile(paramsStream, offline_data_root_path + "/test" + "_optima.txt");
                size_t n_opt;
                paramsStream >> n_opt;
                for (size_t i = 0; i < game_prob->optima()->numberSolutions(); ++i) {
                    auto sol = game_prob->optima()->solution(i);
                    game_prob->parameterVariantsToSolution(paramsStream, sol);
                    if (game_prob->findBySequence(sol.variable()) == game_prob->findBySequence(game_prob->optima()->solution(i).variable())) {
                        std::cout << "solution " << game_prob->findBySequence(sol.variable()) << " is in optima " << i << std::endl;
                    }
                }
            }
        }


    }
}
