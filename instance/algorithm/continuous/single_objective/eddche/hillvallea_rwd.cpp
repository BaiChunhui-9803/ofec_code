#include "hillvallea_rwd.h"
#include "../../../../../core/environment/environment.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../datum/datum_inclusion.h"

namespace ofec {
	void HillVallEA_RWD::addInputParameters() {

	}

	void HillVallEA_RWD::run_(Environment *env) {
		while (!terminating()) {
			//if (m_maximum_evaluations > 0) {
			//	if (m_maximum_evaluations - m_evaluations < m_N) {
			//		break;
			//	}
			//}
			auto P = uniformlySample(m_N, env, m_random.get());
#ifdef OFEC_DATUM_NUM_SAMPLES_H
			g_num_samples.value = m_N;
			datumUpdated(env, g_num_samples);
#endif // OFEC_DATUM_NUM_SAMPLES_H
			for (auto &e : m_E) { P.push_back(e); }
			auto S = truncationSelection(P, m_tau, env);
			auto K = clustering(S, env);
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(K.size());
			for (size_t k = 0; k < K.size(); k++) {
				for (size_t i = 0; i < K[k].size(); ++i) {
					g_multi_pop.pops[k].push_back(K[k][i].get());
				}
			}
			datumUpdated(env, g_multi_pop);
#endif
			// Second phase - niche optimization
			for (auto &C : K) {
				if (duplicate(*best(C, env), m_E, env)) {
					continue;
				}
				auto x = coreSearch(C, m_N_c, env, m_random.get());
				if (terminating()) {
					break;
				}
				if (distinct(*x, m_E, env)) {
					if (betterThanAll(*x, m_E, m_TOL, env)) {
						m_E.clear();
					}
					m_E.push_back(x);
				}
			}
		}
	}
}
