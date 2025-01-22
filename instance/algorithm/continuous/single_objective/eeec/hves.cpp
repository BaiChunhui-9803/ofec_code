#include "hves.h"
#include <numeric>
#include "../multi_modal/hill_vallea/hill_vallea.h"
#include "../../../../../instance/algorithm/continuous/single_objective/global/cmsa_es/cmsa_es_pop.h"
#include "../../../../../core/algorithm/multi_population.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void HVES::addInputParameters() {
		m_input_parameters.add("number of samples", new RangedSizeT(m_num_samples, 5, 1000, 20));
		m_input_parameters.add("number of gradations", new RangedSizeT(m_num_gradations, 1, 15, 2));
	}

	void HVES::run_(Environment *env) {
		std::vector<std::shared_ptr<Solution<>>> sols;
		for (size_t i = 0; i < m_num_samples; i++) {
			auto new_sol = dynamic_cast<Solution<>*>(env->problem()->createSolution());
			new_sol->initialize(env, m_random.get());
			new_sol->evaluate(env);
			sols.emplace_back(new_sol);
		}
		std::vector<size_t> sequence(sols.size());
		std::iota(sequence.begin(), sequence.end(), 0);
		if (env->problem()->optimizeMode(0) == OptimizeMode::kMinimize) {
			std::sort(sequence.begin(), sequence.end(),
				[&sols](size_t i, size_t j) { return sols[i]->objective(0) > sols[j]->objective(0);  }
			);
		}
		else {
			std::sort(sequence.begin(), sequence.end(),
				[&sols](size_t i, size_t j) { return sols[i]->objective(0) < sols[j]->objective(0); }
			);
		}
		std::vector<std::vector<std::shared_ptr<Solution<>>>> clusters;
		while (!sequence.empty()) {
			auto cur_best = sequence.back();
			std::list<size_t> belongings;
			for (size_t k = 0; k < clusters.size(); ++k) {
				if (HillVallEA::test(*clusters[k].front(), *sols[cur_best], m_num_gradations, env)) {
					belongings.push_back(k);
				}
			}
			if (belongings.empty()) {
				clusters.resize(clusters.size() + 1);
				clusters.back().emplace_back(new Solution<>(*sols[cur_best]));
			}
			else {
				int nearest = -1; Real min_distance = std::numeric_limits<Real>::max();
				for (auto k : belongings) {
					auto distance = sols[cur_best]->variableDistance(*clusters[k].front(), env);
					if (distance < min_distance) {
						min_distance = distance;
						nearest = k;
					}
				}
				clusters[nearest].emplace_back(new Solution<>(*sols[cur_best]));
			}
			sequence.pop_back();
		}
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(clusters.size());
		for (size_t k = 0; k < clusters.size(); k++) {
			for (size_t i = 0; i < clusters[k].size(); ++i) {
				g_multi_pop.pops[k].push_back(clusters[k][i].get());
			}
		}
		datumUpdated(env, g_multi_pop);
#endif
		MultiPopulation<PopCMSA_ES> multi_pops;
		Real D = env->problem()->numberVariables();
		for (auto &cluster : clusters) {
			if (cluster.size() > 1) {
				auto pop = std::make_unique<PopCMSA_ES>(3 * sqrt(D), env);
				pop->initializeBySamples(cluster, env);
				multi_pops.append(pop);
			}
		}
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(multi_pops.size());
		for (size_t k = 0; k < multi_pops.size(); k++) {
			for (size_t i = 0; i < multi_pops[k].size(); ++i) {
				g_multi_pop.pops[k].push_back(&multi_pops[k][i]);
			}
		}
#endif
		for (auto &pop : multi_pops) {
			pop->sampleNewPopulation(env, m_random.get());
		}
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
		while (!terminating()) {
			multi_pops.evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
		}
	}
}
