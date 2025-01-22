#ifndef CPSOSWARM_H
#define CPSOSWARM_H

#include "cpso_particle.h"
#include "../../../../template/classic/pso/swarm.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../problem/continuous/single_objective/dynamic/uncertianty_continuous.h"
#include <memory>

namespace ofec {
    using templateParticle = CPSOParticle;
    class CPSOSwarm : public Swarm<templateParticle> {
    public:
        CPSOSwarm(size_t pop, Problem *pro);
        void initializeParameters(Problem *pro, Algorithm *alg, Random *rnd);
        //CPSOSwarm()
        //void initialize(Problem *pro, Random *rnd) override;
        void initializeRadius();

        int localSearch(Problem *pro, Algorithm *alg, Random *rnd, int pop_size);

        void setNeighborhood(Random *rnd) override;
        void setConverge(bool flag) { m_converge_flag = flag; }
        bool getConverge() const { return m_converge_flag; }

        int learnGbest(const std::unique_ptr<templateParticle>& pi, std::unique_ptr<templateParticle>& gbest, Problem *pro, Algorithm *alg) const;

        int best_idx(Problem *pro);

        void computeCenter(Problem *pro);

        void updateCurRadius(Problem *pro);

        void checkOverCrowd(size_t max_sub_size, Problem *pro, Random *rnd);

        void sortPbest(Problem *pro);

        void merge(CPSOSwarm& pop, Problem *pro, Random *rnd);

        Real getInitCurRadius() const { return m_initradius; }
        Real getCurRadius() const { return m_radius; }

        std::unique_ptr<Solution<>> & center() { return m_center; }

    protected:
        double m_weight, m_maxW, m_minW;                // inertia weight
        size_t m_evals, m_u;
        Real m_radius;
        Real m_initradius;

        bool m_converge_flag;
        std::unique_ptr<Solution<>> m_center;
    };
}

#endif
