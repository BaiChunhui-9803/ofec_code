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


#ifndef EXAMPLE_2_H
#define EXAMPLE_2_H

#include "../../core/global.h"
#include "../../run/interface.h"
#include "../../instance/environment/example/game/game_environment.h"
#include "../../instance/environment/template/uncertainty/dynamic_problem.h"
#include "../../instance/problem/realworld/game/game.h"

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

    void run_example_2(int argc, char *argv[]);
}


#endif //EXAMPLE_2_H
