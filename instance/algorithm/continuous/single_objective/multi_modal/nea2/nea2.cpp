#include "nea2.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../utility/clustering/nbc.h"
#include "../../../../../../utility/functional.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	void NEA2::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
		m_input_parameters.add("without restart", new Bool(m_without_restart, false));
	}

	void NEA2::initialize_(Environment *env) {
		Algorithm::initialize_(env);
	}

	void NEA2::run_(Environment *env) {
		while (!terminating()) {
			addSubpops(m_pop_size, env);
			while (!terminating() && m_subpops.isActive()) {
				int rf = m_subpops.evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
				datumUpdated(env, g_multi_pop);
#endif
				if (rf & kTerminate)
					break;
				for (size_t k = 0; k < m_subpops.size(); ++k) {
					if (m_subpops[k].isActive() && stopTolFun(m_subpops[k]))
						m_subpops[k].setActive(false);
				}
			}
			if (m_without_restart) {
				break;
			}
		}
	}

	void NEA2::addSubpops(size_t num_samples, Environment *env) {
		size_t num_vars = env->problem()->numberVariables();
		Population<Solution<>> global_sample(num_samples, env, num_vars);
		global_sample.initialize(env, m_random.get());
		global_sample.evaluate(env);
		NBC nbc(2.0, NBC::kByKDTree, NBC::kByMean, true);
		nbc.setData(global_sample);
		nbc.clustering(env);
		auto &clusters = nbc.clusters();
#ifdef OFEC_DATUM_SEEDS_H
		g_seeds.sols.clear();
		for (size_t k = 0; k < clusters.size(); ++k) {
			SolutionBase *new_seed = nullptr;
			for (size_t i = 0; i < clusters[k].size(); ++i) {
				if (!new_seed || dominate(global_sample[clusters[k][i]], *new_seed, env->problem()->optimizeMode())) {
					new_seed = &global_sample[clusters[k][i]];
				}
			}
			if (new_seed) {
				g_seeds.sols.push_back(new_seed);
			}
		}
		datumUpdated(env, g_seeds);
#endif // OFEC_DATUM_SEEDS_H
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(clusters.size());
		for (size_t k = 0; k < clusters.size(); k++) {
			for (size_t i = 0; i < clusters[k].size(); ++i) {
				g_multi_pop.pops[k].push_back(&global_sample[clusters[k][i]]);
			}
		}
		datumUpdated(env, g_multi_pop);
#endif
		std::vector<size_t> filtered_start_Solutions;
		for (auto& cluster : clusters) {
			auto iter = cluster.begin();
			size_t center = *iter;
			while (++iter != cluster.end()) {
				if (dominate(global_sample[*iter], global_sample[center], env->problem()->optimizeMode())) {
					center = *iter;
				}
			}
			filtered_start_Solutions.push_back(center);
		}

		std::vector<Real> step_size(num_vars);
		std::vector<std::pair<Real, Real>> bnd(num_vars, { std::numeric_limits<Real>::max(),-std::numeric_limits<Real>::max() });
		Solution<> tmp_sol(1, 0, num_vars);
		for (size_t i = 0; i < 1000; ++i) {
			tmp_sol.initialize(env, m_random.get());
			for (size_t j = 0; j < num_vars; ++j) {
				if (tmp_sol.variable()[j] < bnd[j].first) {
					bnd[j].first = tmp_sol.variable()[j];
				}
				if (tmp_sol.variable()[j] > bnd[j].second) {
					bnd[j].second = tmp_sol.variable()[j];
				}
			}
		}
		for (size_t j = 0; j < num_vars; ++j) {
			step_size[j] = 0.05 * (bnd[j].second - bnd[j].first);
		}
		//auto &domain = CAST_CONOP(env->problem())->domain();
		//Real diag = 0;
		//for (size_t j = 0; j < num_vars; ++j)
		//	diag += pow(domain[j].limit.second - domain[j].limit.first, 2);
		//diag = sqrt(diag);
		size_t size_subpop = 4 + int(3 * log(num_vars));
		m_subpops.clear();
		for (size_t start_ind : filtered_start_Solutions) {
			auto pop = std::make_unique<PopCMA_ES>(size_subpop, env);
			//Real step_size = std::max(0.025 * diag, m_random->normal.nextNonStd(0.05 * diag, 0.025 * diag));
			//pop->initializeByNonStd(env, m_random.get(), global_sample[start_ind].variable().vect(), std::vector<Real>(num_vars, step_size));
			pop->initializeByNonStd(env, m_random.get(), global_sample[start_ind].variable().vector(), step_size);
			pop->reproduce(env, m_random.get());
			pop->evaluate(env);
			m_subpops.append(pop);
		}
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(m_subpops.size());
		for (size_t k = 0; k < m_subpops.size(); k++) {
			for (size_t i = 0; i < m_subpops[k].size(); ++i)
				g_multi_pop.pops[k].push_back(&m_subpops[k][i]);
		}
		datumUpdated(env, g_multi_pop);
#endif
	}

	bool NEA2::stopTolFun(PopCMA_ES& subpop) {
		Real min_obj, max_obj;
		min_obj = max_obj = subpop[0].objective(0);
		for (size_t i = 1; i < subpop.size(); i++) {
			if (subpop[i].objective(0) < min_obj)
				min_obj = subpop[i].objective(0);
			if (subpop[i].objective(0) > max_obj)
				max_obj = subpop[i].objective(0);
		}
		if ((max_obj - min_obj) < 1e-4)
			return true;
		else
			return false;
	}
}
