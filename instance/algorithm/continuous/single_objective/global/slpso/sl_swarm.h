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

#ifndef OFEC_SL_SWARM_H
#define OFEC_SL_SWARM_H

#include "sl_particle.h"
#include "../../../../template/classic/particle_swarm_optimization/swarm.h"

namespace ofec {
    class SwarmSL : public Swarm<ParticleSL> {
    public:
        int m_numLearnToGbest;
        static int ms_updateFre;
        static double ms_learnRatio;
        static float ms_ratioLearnToGbest;
        std::vector<int> mv_sel;
    public:
        SwarmSL() {}
        SwarmSL(size_t size_pop, Environment *env);
        ~SwarmSL() {}
        void resize(size_t size_pop, Environment* env);
        void updateLearnToGbest(Random *rnd);
        int evolve(Environment *env, Random *rnd) override;
        void setParameters(Random *rnd);
        void updateParameters(Random *rnd);
        int getNumLearnTogbest();
        void MRandToBest(int num, Random *rnd);
        int updateBestByRatio(const int p, double ratio, Environment *env, Random *rnd);
        void calculateNumLearning(const int sfes, Algorithm *alg);
        int sel(size_t i) const { return mv_sel[i]; }
    };
}

#endif //!OFEC_SL_SWARM_H