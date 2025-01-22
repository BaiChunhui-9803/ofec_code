//
// Created by Dell on 2022/11/17.
//

#include "dfdl_particle.h"
#include "../../../../../../core/problem/continuous/continuous.h"

namespace ofec {
    void DFDLParticle::initialize(int id, Problem *pro, Random *rnd) {
        Solution::initialize(id, pro, rnd);
        initVelocity(pro, rnd);
        initVelocityMax(pro, rnd);
    }

    void DFDLParticle::initVelocityMax(Problem *pro, Random *rnd) {
        for (size_t j = 0; j < variable().size(); j++) {
            auto &range = CAST_CONOP(pro)->range(j);
            m_vel_max[j].max = (range.second - range.first) / 10;
            m_vel_max[j].min = -m_vel_max[j].max;
        }
    }

    void DFDLParticle::initVelocity(Problem *pro, Random *rnd) {
        for (size_t j = 0; j < variable().size(); j++) {
            auto &range = CAST_CONOP(pro)->initialDomain()[j].limited ?
                          CAST_CONOP(pro)->initialRange(j) : CAST_CONOP(pro)->range(j);
            m_vel[j] = (range.second - range.first) * (-0.5 + rnd->uniform.next()) / 10;
        }
    }

    void DFDLParticle::initNichingVelocity(Problem *pro, Random *rnd, Real radius) {
        for (size_t j = 0; j < variable().size(); j++) {
            auto &range = CAST_CONOP(pro)->initialDomain()[j].limited ?
                          CAST_CONOP(pro)->initialRange(j) : CAST_CONOP(pro)->range(j);
            m_vel[j] = (2 * radius) * (-0.5 + rnd->uniform.next());
        }
    }

    void DFDLParticle::brownianMove(Problem *pro, Random *rnd, double radius) {
        for (size_t j = 0; j < this->variable().size(); j++) {
            this->variable()[j] += rnd->normal.nextNonStd(0, radius);
        }
        CAST_CONOP(pro)->validateSolution(*this, Validation::kSetToBound, rnd);
    }

    int DFDLParticle::cauchyMove(Problem *pro, Algorithm *alg, Random *rnd, double radius) {
        int rf;
        auto& domain = CAST_CONOP(pro)->domain();
        for (size_t i = 0; i < this->variable().size(); i++) {
            if (radius < 0) {
                this->variable()[i] += rnd->cauchy.nextNonStd(0, (domain[i].limit.second - domain[i].limit.first) / 2);
            }
            else {
                this->variable()[i] += rnd->cauchy.nextNonStd(0, radius);
            }
        }
        CAST_CONOP(pro)->validateSolution(*this, Validation::kSetToBound, rnd);
        // pbest move: WARNING - if pbest not move, the swarm will be also stagnanted.
        rf = this->evaluate(pro, alg, true);
        this->pbest() = *this;
        return rf;
    }
}