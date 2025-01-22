/********* Begin Register Information **********
{
    "name": "SLPSO",
    "identifier": "SLPSO",
    "tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

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

#ifndef OFEC_SLPSO_H
#define OFEC_SLPSO_H

#include "sl_swarm.h"
#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
    class SLPSO : virtual public Algorithm {
        OFEC_CONCRETE_INSTANCE(SLPSO)
    protected:
        size_t m_swarm_size;
        Real m_maxW, m_minW, m_accelerator1, m_accelerator2;
        SwarmSL m_swarm;

        void addInputParameters();
        void initialize_(Environment *env) override;
        void run_(Environment *env) override;
    };
}

#endif //!OFEC_SLPSO_H