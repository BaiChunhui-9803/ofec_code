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

#ifndef OFEC_SL_PARTICLE_H
#define OFEC_SL_PARTICLE_H

#include "../../../../template/classic/particle_swarm_optimization/particle.h"
#include "sl_progress.h"

namespace ofec {
    class SwarmSL;
    class ParticleSL : public Particle {
        friend class SwarmSL;

    protected:
        //bool m_type; inherited from Solution , flag of Learning to gbest;
        bool m_type;
        size_t m_itersUnimpr;
        size_t m_updateFre;
        double m_learnRatio;
        static const int ms_numOperators;
        std::vector<slpso::Progress> mv_prog, mv_monitor;

    public:
        virtual ~ParticleSL() {}
        ParticleSL();
        ParticleSL(size_t num_obj, size_t num_con, size_t size_var);
        ParticleSL(const ParticleSL &rhs);
        ParticleSL(ParticleSL &&rhs) noexcept;
        ParticleSL& operator=(const ParticleSL &other);
        ParticleSL& operator=(ParticleSL &&other) noexcept;

        void setType(bool type) { m_type = type; }
        bool type() const { return m_type; }
        void setSelRatio();
        void nonLearnToLearn();
        void learnToNonLearn();
        void updateSelectionRatioMonitor(Random *rnd);
        void updateSelectionRatioProg(Random *rnd);
        int selectOperator(Random *rnd);
        void nextVelocity(const Solution<> *lbest, Real w, Real c1, Real c2, Random *rnd) override;
        void normalMutation(double *avg_v, Random *rnd, Environment *env);
    };
}

#endif //!OFEC_SL_PARTICLE_H