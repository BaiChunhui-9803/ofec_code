#include "ncde.h"
#include "../../../../template/classic/differential_evolution/population.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	void NCDE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 20, 1000, 50));
	}

	void NCDE::initialize_(Environment* env) {
		Algorithm::initialize_(env);
		if (m_maximum_evaluations <= 0) {
			throw Exception("Maximum number of evaluations must be provided.");
		}
	}

	void NCDE::run_(Environment *env) {
		PopulationDE<> pop(m_pop_size, env);
		pop.setParameter(0.1, 0.9);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < pop.size(); ++i) {
			g_multi_pop.pops[0].push_back(&pop[i]);
		}
#endif
		pop.initialize(env, m_random.get());
		pop.evaluate(env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
		while (!terminating()) {
			int rf = kNormalEval;
			for (size_t i = 0; i < m_pop_size; i++) {
				size_t m = 0;
				if (m_pop_size <= 200) {
					m = 5 + 5.0 * (m_maximum_evaluations - m_evaluations) / m_maximum_evaluations;
				}
				else {
					m = 20 + 30.0 * (m_maximum_evaluations - m_evaluations) / m_maximum_evaluations;
				}
				std::vector<std::pair<Real, size_t>> distances;
				for (size_t j = 0; j < m_pop_size; j++) {
					distances.emplace_back(pop[i].variableDistance(pop[j], env), j);
				}
				std::nth_element(distances.begin(), distances.begin() + m, distances.end());
				std::vector<size_t> mnearest;
				for (size_t j = 0; j < m; j++) {
					mnearest.push_back(distances[j].second);
				}
				std::vector<size_t> ridx(3);
				pop.selectInCandidates(3, mnearest, ridx, m_random.get());
				pop[i].mutate(pop.scalingFactor(), &pop[ridx[0]], &pop[ridx[1]], &pop[ridx[2]], env);
				pop.recombine(i, m_random.get(), env);
				size_t nearest_neighbor = 0; 
				Real min_distance = pop[i].trial().variableDistance(pop[nearest_neighbor], env);
				for (size_t j = 1; j < m_pop_size; j++) {
					Real distance = pop[i].trial().variableDistance(pop[j], env);
					if (distance < min_distance) {
						nearest_neighbor = j;
						min_distance = distance;
					}
				}
				rf = pop[i].trial().evaluate(env);
				if (dominate(pop[i].trial(), pop[nearest_neighbor], env->problem()->optimizeMode())) {
					pop[nearest_neighbor].solution() = pop[i].trial();
				}
				if (rf & kTerminate) {
					break;
				}
			}
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
			if (rf & kTerminate) {
				break;
			}
		}
	}
}



