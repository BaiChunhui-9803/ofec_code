/********* Begin Register Information **********
{
	"name": "DMOP_F08",
	"identifier": "DMOP_F08",
	"problem tags": [ "DMOP", "ConOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Qingshan Tan
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

* see https://github.com/Changhe160/OFEC for more information
*************************************************************************/

/************************************************************************
Aimin Zhou, Yaochu Jin, and Qingfu Zhang.
A population prediction strategy for evolutionary dynamic multiobjective optimization.
IEEE transactions on cybernetics, 44(1):40-53, 2013
************************************************************************/

// Created: 5 August 2019 by Qingshan Tan
// Last modified at 5 August 2019 by Qingshan Tan

#ifndef F8_H
#define F8_H


#include "../DMOPs.h"
#include"../metrics_dmop.h"

namespace ofec {
    namespace DMOP {
        class F8 : public DMOPs, public MetricsDMOP {
        public:
            void initialize_() override;
            void generateAdLoadPF() override;
        private:
            int updateEvaluationTag(SolutionBase& s, Algorithm *alg) override;
            void evaluateObjective(Real* x, std::vector<Real>& obj);
        };
    }
    using DMOP_F08 = DMOP::F8;
}

#endif //F8_H
