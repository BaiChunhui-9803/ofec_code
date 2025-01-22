#include "fep.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"


namespace ofec {
	void PopFEP::mutate(Random *rnd, Environment *env) {
		for (size_t i = 0; i < m_individuals.size(); ++i)
			*m_offspring[i] = *m_individuals[i];
		for (size_t i = 0; i < m_individuals.size(); i++) {
			Real N = rnd->normal.next();
			for (size_t j = 0; j < m_offspring[i]->variable().size(); ++j) {
				Real delta_j = rnd->cauchy.next();
				m_offspring[i+ m_individuals.size()]->variable()[j] = m_individuals[i]->variable()[j] + m_individuals[i]->eta()[j] * delta_j;
				Real N_j = rnd->normal.next();
				m_offspring[i + m_individuals.size()]->eta()[j] = m_individuals[i]->eta()[j] * exp(m_tau_prime * N + m_tau * N_j);
			}
		}
	}

	void FEP::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
	}

	void FEP::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		m_tau_prime = 1 / (sqrt(2 / sqrt(m_pop_size)));
		m_tau = 1 / (sqrt(2 * m_pop_size));
		m_q = 10;
	}

	void FEP::run_(Environment *env) {
		m_pop.resize(m_pop_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop.size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop[i]);
		}
#endif
#ifdef OFEC_DATUM_EP_POP_H
		g_ep_pop.value = &m_pop;
#endif // OFEC_DATUM_EP_POP_H
		m_pop.tauPrime() = m_tau_prime;
		m_pop.tau() = m_tau;
		m_pop.q() = m_q;
		m_pop.initialize(env, m_random.get());
		m_pop.evaluate(env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
#ifdef OFEC_DATUM_EP_POP_H
		datumUpdated(env, g_ep_pop);
#endif
		while (!terminating()) {
			m_pop.evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
#ifdef OFEC_DATUM_EP_POP_H
			datumUpdated(env, g_ep_pop);
#endif
		}
	}
}