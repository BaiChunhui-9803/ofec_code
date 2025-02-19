//
// Created by Administrator on 25-2-17.
//

#include "example_2.h"
#include <filesystem>

namespace ofec {

    void run2(int argc, char *argv[]) {
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

        if (auto dp = dynamic_cast<DynamicProblem*>(env->problem())) {
            Random rnd(0.5);
            dp->change(&rnd);
        }


    }
}