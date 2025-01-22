#ifndef SPSOSWARM_H
#define SPSOSWARM_H

#include "spso_particle.h"
#include "../../../../template/classic/pso/swarm.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include <memory>

namespace ofec {
    class SPSOSwarm : public Swarm<SPSOParticle> {
    public:
        SPSOSwarm(size_t pop, Problem *pro);
        void initialize(Problem *pro, Random *rnd);
        void findpopSeed();
        int evolve(Problem *pro, Algorithm *alg, Random *rnd) override;
        int best_idx(Problem *pro);
        void setNeighborhood(Random *rnd) override;
        //SPSOParticle getpopSeed() { return m_seed; }
    protected:
        //SPSOParticle m_seed;
    };
}

#endif
