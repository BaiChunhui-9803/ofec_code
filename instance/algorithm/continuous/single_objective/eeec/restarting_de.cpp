#include "restarting_de.h"
#include "../../../../template/classic/differential_evolution/population.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void RestartingDE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
	}

	void RestartingDE::run_(Environment *env) {
		PopulationDE<> pop(m_pop_size, env);
		pop.initialize(env, m_random.get());
		pop.evaluate(env);
		while (!terminating()) {
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(1);
			for (size_t i = 0; i < pop.size(); ++i) {
				g_multi_pop.pops[0].push_back(&pop[i]);
			}
#endif
			datumUpdated(env);
			pop.evolve(env, m_random.get());
			if (pop.best(env)->objectiveDistance(*pop.worst(env)) < 1e-8) {
#ifdef OFEC_DATUM_MULTI_POP_H
				g_multi_pop.pops.clear();
#endif
				datumUpdated(env);
				pop.initialize(env, m_random.get());
				pop.evaluate(env);
			}
		}
	}
}
