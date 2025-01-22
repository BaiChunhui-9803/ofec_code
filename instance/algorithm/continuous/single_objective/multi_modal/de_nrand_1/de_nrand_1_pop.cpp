#include "de_nrand_1_pop.h"

namespace ofec {
	PopDE_nrand_1::PopDE_nrand_1(size_t size_pop, Environment *env) : 
		PopulationDE(size_pop, env) {}

	void PopDE_nrand_1::mutate(int idx, Random *rnd, Environment *env) {
		std::vector<size_t> r;
		select(idx, 2, r, rnd);
		size_t nn = nearestNeighbour(idx, env);
		m_individuals[idx]->mutate(
			m_scaling_factor, 
			m_individuals[nn].get(),
			m_individuals[r[0]].get(), 
			m_individuals[r[1]].get(), 
			env
		);
	}

	size_t PopDE_nrand_1::nearestNeighbour(int idx, Environment *env) {
		size_t nn = idx == 0 ? 1 : 0;
		Real min_dis = m_individuals[idx]->variableDistance(*m_individuals[nn], env);
		for (size_t i = 0; i < m_individuals.size(); ++i) {
			if (i == idx) {
				continue;
			}
			Real dis = m_individuals[idx]->variableDistance(*m_individuals[i], env);
			if (dis < min_dis) {
				nn = i;
				min_dis = dis;
			}
		}
		return nn;
	}
}

