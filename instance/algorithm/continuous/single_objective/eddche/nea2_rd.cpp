#include "nea2_rd.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void NEA2_RD::addInputParameters() {

	}

	void NEA2_RD::initialize_(Environment *env) {
		NEA2::initialize_(env);

	}

	void NEA2_RD::run_(Environment *env) {
		size_t num_samples = m_pop_size;
		while (!terminating()) {
			addSubpops(num_samples, env);
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
			num_samples *= 2;
		}
	}
}
