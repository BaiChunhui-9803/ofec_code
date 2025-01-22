/********* Begin Register Information **********
{
	"name": "kNPSO",
	"identifier": "kNPSO",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#ifndef OFEC_K_NEAREST_PSO_H
#define OFEC_K_NEAREST_PSO_H

#include "../../../../../core/algorithm/algorithm.h"
#include "../../../template/classic/particle_swarm_optimization/swarm.h"
#include "../../../template/classic/particle_swarm_optimization/particle.h"

namespace ofec {
	class kNSwarm : public Swarm<Particle> {
	public:
		void setNumNeighbors(size_t num_nbrs) { m_num_nbrs = num_nbrs; }
		void resize(size_t size_pop, Environment *env) override;
		void setNeighborhood(Random *rnd) override;

	protected:
		size_t m_num_nbrs;
		size_t m_number_variables;
		std::vector<std::vector<Real>> m_swarm_vars;
	};

	class kNPSO : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(kNPSO)
	protected:
		size_t m_swarm_size;
		size_t m_num_nbrs;
		kNSwarm m_swarm;

		void addInputParameters();
		void run_(Environment *env) override;
	};
}

#endif // !OFEC_K_NEAREST_PSO_H
