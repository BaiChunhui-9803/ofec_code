#ifndef OFEC_LIPS_POP_H
#define OFEC_LIPS_POP_H

#include "../../../../template/classic/particle_swarm_optimization/particle.h"
#include "../../../../template/classic/particle_swarm_optimization/swarm.h"
#include <list>

namespace ofec {
	class SwarmLIP;
	class LI_particle : public Particle {
		friend SwarmLIP;
	private:
		std::vector<Real> m_pos;           // the mean position of neighborhood
		std::vector<Real> m_rdm;
		std::vector<int> m_nbr;
	public:
		LI_particle(size_t num_obj, size_t num_con, size_t size_var);
		void nextVelocityByWeight(Real w);
	};

	class SwarmLIP : public Swarm<LI_particle> {
	protected:
		int m_M;            // the num of neighbors
		int m_maximum_evalutions;    // the maximum num of evaluations
		bool m_use_history_nearest = false;
		std::vector<std::list<std::pair<Real, size_t>>> m_dis;
		void setBestPos(int idx_ind, Environment *env, Random *rnd);
		void sortDistance(int idx_ind, Environment *env);
	public:
		SwarmLIP() = default;
		SwarmLIP(size_t size_pop, Environment *env, int max_evals);
		void resize(size_t size_pop, Environment *env, int max_evals);
		void clear() override;
		int evolve(Environment *env, Random *rnd) override;
		void setUseHistoryNearest(bool flag) { m_use_history_nearest = flag; }
	};
}

#endif // !OFEC_LIPS_POP_H

