/********* Begin Register Information **********
{
	"name": "DMOP_T3",
	"identifier": "DMOP_T3",
	"problem tags": [ "DMOP", "ConOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 31 December 2014
// Modified: 29 Mar 2018 by Junchen Wang (wangjunchen@cug.edu.cn)

#ifndef T3_H
#define T3_H

#include "../DMOPs.h"

namespace ofec {
    namespace DMOP {
        class T3 final : public DMOPs {
        public:
            T3(const ParameterMap &v);
            T3(const std::string &name, size_t size_var);
            ~T3() {}
            void initialize();
            void generateAdLoadPF();
        private:
            void evaluateObjective(Real *x, std::vector<Real> &obj);
        };
    }
    using DMOP_T3 = DMOP::T3;
}

#endif //T3_H


