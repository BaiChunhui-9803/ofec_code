//
// Created by Administrator on 25-2-17.
//

#ifndef EXAMPLE_2_H
#define EXAMPLE_2_H

#include "../../core/global.h"
#include "../../run/interface.h"
#include "../../instance/environment/example/game/game_environment.h"
#include "../../instance/environment/template/uncertainty/dynamic_problem.h"

namespace ofec {
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

    void run2(int argc, char *argv[]);
}


#endif //EXAMPLE_2_H
