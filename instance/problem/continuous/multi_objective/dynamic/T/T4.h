/********* Begin Register Information **********
{
	"name": "DMOP_T4",
	"identifier": "DMOP_T4",
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

#ifndef T4_H
#define T4_H


#include "../DMOPs.h"

namespace ofec {
    namespace DMOP {
        class T4 final : public DMOPs {
        public:
            T4(const ParameterMap &v);
            T4(const std::string &name, size_t size_var);
            ~T4() {}
            void initialize();
            void generateAdLoadPF();
        private:
            void evaluateObjective(Real *x, std::vector<Real> &obj);
        };
    }
    using DMOP_T4 = DMOP::T4;
}

#endif //T4_H



