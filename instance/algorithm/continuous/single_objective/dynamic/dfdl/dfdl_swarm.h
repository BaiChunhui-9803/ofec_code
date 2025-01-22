#ifndef OFEC_DFDL_SWARM_H
#define OFEC_DFDL_SWARM_H

#include "../../../../template/classic/pso/swarm.h"
#include "dfdl_particle.h"
#include"../../../../../../core/problem/continuous/continuous.h"
#include <iostream>

namespace ofec {
    class DFDLSwarm : public Swarm<DFDLParticle> {
    public:
        enum class PopulationType { Explore, Exploit };
        int m_stagnanted_count = 0;
    private:
        Real m_initRadius = 0.;
        Real m_curRadius = 0.;
        Real m_curVelocity = 0.;
        bool m_convergeFlag = false;
        bool m_seedFlag = false;
        bool m_stagnantedFlag = false;
        size_t m_max_sub_size;
        size_t m_hill_located = -1;
        Real m_converge_radius = 10e-3;
        Real m_seed_radius = 10e-2;
        std::unique_ptr<Solution<>> m_center;
    public:
        PopulationType m_pop_type;
    public:
        DFDLSwarm(size_t pop, Problem *pro, PopulationType pop_type);
        void initializeParameters(Problem *pro, Random *rnd);
        void initializeNichingPops(Problem *pro, Random *rnd);
        size_t bestIndex(Problem *pro);
        void setNeighborhood(Random *rnd) override;
        void initializeRadius();
        void computeCenter(Problem *pro);
        void updateCurRadius(Problem *pro);
        void updateCurVelocity(Problem *pro);
        Real getInitRadius() const { return m_initRadius; }
        Real getCurRadius() const { return m_curRadius; }
        Real getCurVelocity() const { return m_curVelocity; }
        void merge(DFDLSwarm& pop, Problem *pro, Random *rnd);
        void checkOverCrowd(Problem *pro, Random *rnd);
        void checkConverged(Problem *pro);
        void checkStagnanted(Real avg_radius);
        void setConvergeFlag(bool flag) { m_convergeFlag = flag; }
        void setSeedFlag(bool flag) { m_seedFlag = flag; }
        bool getSeedFlag() { return m_seedFlag; }
        void setStagnantedFlag(bool flag) { m_stagnantedFlag = flag; }
        void setHillLocated(size_t hill) { m_hill_located = hill; }
        size_t getHillLocated() { return m_hill_located; }
        bool getConvergeFlag() { return m_convergeFlag; }
        bool getStagnantedFlag() { return m_stagnantedFlag; }
        int evolve(Problem *pro, Algorithm *alg, Random *rnd) override;

        friend std::ostream & operator<<(std::ostream &out, DFDLSwarm &A);
    protected:
        void sortPbest(Problem *pro);
    };
} // ofec

#endif //OFEC_DFDL_SWARM_H
