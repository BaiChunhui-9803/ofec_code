/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@gmail.com
* Language: C++
* -----------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/

/*************************************************************************
* C. Li,S. Yang, and T. T. Nguyen.
* A self-learning particle swarm optimizer for global optimization problems.
* IEEE Transactions on Systems, Man, and Cybernetics Part B: Cybernetics, vol. 42, no. 3, pp. 627¨C646, 2011
* -----------------------------------------------------------------------
* Created: 19 Jan 2015
* Last modified: 8 Apr 2022 by Junchen Wang
*************************************************************************/

#ifndef OFEC_SL_PROGRESS_H
#define OFEC_SL_PROGRESS_H

#include <vector>
#include "../../../../../../core/random/newran.h"

namespace ofec::slpso {
    class Progress {
    public:
        int m_numSelected;
        int m_numSuccess;
        double m_rewards;
        double m_ratio;
        double m_minRatio;
    public:
        Progress() {}
        virtual ~Progress() {}
        void initialize(const int num_ls_operators);
        static void updateProgress(Random *rnd, std::vector<Progress> &prog, const int num_ls_operators, double p_factor = 0.9);
        static size_t getAction(Random *rnd, const size_t num_actions, const std::vector<Progress> &p);
    };
}

#endif // !OFEC_SL_PROGRESS_H
