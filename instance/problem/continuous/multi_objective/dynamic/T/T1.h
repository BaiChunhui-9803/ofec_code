/********* Begin Register Information **********
{
	"name": "DMOP_T1",
	"identifier": "DMOP_T1",
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

#ifndef T1_H
#define T1_H


#include "../DMOPs.h"

namespace ofec {
    namespace DMOP {
        class T1 final : public DMOPs {
        public:
            T1(const ParameterMap &v);

            T1(const std::string &name, size_t size_var);

            ~T1() {}

            void initialize();

            void generateAdLoadPF();

        private:
            void evaluateObjective(Real *x, std::vector<Real> &obj);
        };
    }
    using DMOP_T1 = DMOP::T1;
}

#endif //T1_H
