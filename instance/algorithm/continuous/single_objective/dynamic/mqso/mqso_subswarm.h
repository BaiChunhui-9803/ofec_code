#ifndef MQSO_SUBSWARM_H
#define MQSO_SUBSWARM_H

#include "mqso_particle.h"
#include "../../../../template/classic/pso/swarm.h"

namespace ofec {
	class MQSOSwarm : public Swarm<MQSOParticle> {
	public:
		enum class SwarmType { SWARM_CONVERGED, SWARM_CONVERGING, SWARM_FREE };
		int evolve(Problem *pro, Algorithm *alg, Random *rnd);
	public:
		MQSOSwarm(size_t pop, Problem *pro);
		~MQSOSwarm() = default;
		void initialize(Algorithm *alg, Problem *pro, Random *rnd);
		int reInitialize(Algorithm *alg, Problem *pro, Random *rnd);
		bool& getworstFlag() {return m_worstFlag;};
		void computeCenter(Problem *pro);
		void updateCurRadius(Problem *pro);
		void setswarmType(SwarmType);
		SwarmType getswarmType() { return m_swarmType; }
		Real getCurRadius() { return m_radius; }
		void setNeighborhood(Random *rnd) override;
		int bestIdx(Problem *pro);
		void setBestSubPop(bool flag) { m_best_subpop = flag; }
		void updateBestIdx(int idx) { m_best_idx = idx; }
		std::unique_ptr<Solution<>> m_center;
	protected:
		void assignType(Algorithm *alg);
		void computeRepulsion(const int idx, Problem *pro);

	protected:
		int m_Nplus;    // the number of charged particles of each swarm
		int m_Nq;       // the number of quantum particles of each swarm
		int m_N;        // the number of neutral particles of each swarm
		double m_Q;     // charge of particles
		double m_Rcloud; // radius of quantum swarms
		SwarmType m_swarmType;		//0 for converged swarm, 1 for converging swarm, 2 for free swarm

		int m_best_idx;
		bool m_worstFlag = false;
		bool m_best_subpop = false;
		Real m_radius;
	};
}


#endif