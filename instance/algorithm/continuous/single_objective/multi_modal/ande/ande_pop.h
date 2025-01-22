#ifndef OFEC_ANDE_POP_H
#define OFEC_ANDE_POP_H

#include "../../../../template/classic/differential_evolution/population.h"
#include "../../../../../../utility/clustering/apc.h"

namespace ofec {
	class PopANDE : public PopulationDE<> {
	protected:
		APC m_apc;
		size_t m_MaxFEs = 1000000;

	public:
		PopANDE(size_t size_pop, Environment *env, Real lambda, size_t Mits, size_t Cits);
		void clustering(Environment *env);
		int evolve(Environment *env, Random *rnd) override;
		const std::vector<std::vector<size_t>>& clusters() { return m_apc.clusters(); }

	private:
		void CPA(const std::vector<size_t> &cluster, Environment *env);
		void TLLS(const std::vector<std::vector<size_t>> &clusters, Environment *env, Random *rnd);
	};
}

#endif // !OFEC_ANDE_POP_H
