#include "knpso.h"
#include "../../../../../utility/kd-tree/data_adaptor.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void kNSwarm::resize(size_t size_pop, Environment *env) {
		Swarm::resize(size_pop, env);
		m_number_variables = env->problem()->numberVariables();
		m_swarm_vars.assign(size_pop, std::vector<Real>(m_number_variables));
	}

	void kNSwarm::setNeighborhood(Random *rnd) {
		for (size_t i = 0; i < m_individuals.size(); ++i) {
			m_swarm_vars[i] = m_individuals[i]->pbest().variable().vector();
		}
		using MyKdTree = KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<Real>>, Real>;
		MyKdTree mat_index(m_number_variables, m_swarm_vars);
		std::vector<size_t> ret_indexes(m_num_nbrs);
		std::vector<Real> out_dists_sqr(m_num_nbrs);
		nanoflann::KNNResultSet<Real> result_set(m_num_nbrs);
		for (size_t i = 0; i < m_individuals.size(); ++i) {
			result_set.init(&ret_indexes[0], &out_dists_sqr[0]);
			mat_index.index->findNeighbors(result_set, &m_swarm_vars[i][0], nanoflann::SearchParams(10));
			m_link[i].assign(m_individuals.size(), false);
			for (size_t k = 0; k < m_num_nbrs; ++k) {
				m_link[i][ret_indexes[k]] = true;
			}
		}
	}

	void kNPSO::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_swarm_size, 5, 1000, 20));
		m_input_parameters.add("number of neighbors", new RangedSizeT(m_num_nbrs, 2, 5000, 2));
	}

	void kNPSO::run_(Environment *env) {
		m_swarm.setNumNeighbors(m_num_nbrs);
		m_swarm.resize(m_swarm_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_swarm.size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_swarm.at(i).pbest());
		}
#endif
		m_swarm.initialize(env, m_random.get());
		m_swarm.initVelocity(env, m_random.get());
		m_swarm.initVelocityMax(env, m_random.get());
		m_swarm.evaluate(env);
		m_swarm.initPbest(env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
		while (!terminating()) {
			m_swarm.evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
		}
	}
}
