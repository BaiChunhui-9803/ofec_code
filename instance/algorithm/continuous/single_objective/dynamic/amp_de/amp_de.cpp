#include "amp_de.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void AMP_DE::addInputParamters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
	}

	void AMP_DE::run_(Environment *env) {
		m_multi_pop.reset(new ContAMP<PopulationType>(m_pop_size));
		m_multi_pop->initialize(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
		datum::g_multi_pop.clear();
		datum::g_multi_pop.resize(m_multi_pop->size());
		for (size_t k = 0; k < m_multi_pop->size(); k++) {
			for (size_t i = 0; i < m_multi_pop->at(k).size(); ++i)
				datum::g_multi_pop[k].push_back(&m_multi_pop->at(k)[i]);
		}
#endif // OFEC_DATUM_MULTI_POP_H
		datumUpdated(env);
		while (!terminating()) {
			m_multi_pop->evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			datum::g_multi_pop.clear();
			datum::g_multi_pop.resize(m_multi_pop->size());
			for (size_t k = 0; k < m_multi_pop->size(); k++) {
				for (size_t i = 0; i < m_multi_pop->at(k).size(); ++i)
					datum::g_multi_pop[k].push_back(&m_multi_pop->at(k)[i]);
			}
#endif // OFEC_DATUM_MULTI_POP_H
			datumUpdated(env);
		}
#ifdef OFEC_DATUM_MULTI_POP_H
		datum::g_multi_pop.clear();
#endif // OFEC_DATUM_MULTI_POP_H
	}
}