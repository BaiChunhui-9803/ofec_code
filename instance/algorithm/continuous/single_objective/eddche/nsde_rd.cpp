#include "nsde_rd.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void NSDE_RD::addInputParameters() {
		m_input_parameters.add("maximum stagnant iterations", new RangedSizeT(m_max_stagant_iters, 5, 100, 10));
	}

	void NSDE_RD::initialize_(Environment *env) {
		NSDE::initialize_(env);

	}

	void NSDE_RD::run_(Environment *env) {
		size_t pop_size = m_pop_size;
		initPop(pop_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
		size_t counter_stagnant = 0;
		while (!terminating()) {
			int rf = m_pop->evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
			if (rf == kNormalEval) {
				bool flag_stagnant = true;
				for (size_t i = 0; i < m_pop->size(); ++i) {
					if (m_pop->at(i).isImproved()) {
						flag_stagnant = false;
						break;
					}
				}
				if (flag_stagnant) {
					counter_stagnant++;
				}
				else {
					counter_stagnant = 0;
				}
				if (counter_stagnant > m_max_stagant_iters) {
					pop_size *= 2;
					initPop(pop_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
					datumUpdated(env, g_multi_pop);
#endif
				}
			}
		}
	}
}
