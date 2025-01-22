#include "slpso.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"


namespace ofec {
    void SLPSO::addInputParameters() {
        m_input_parameters.add("population size", new RangedSizeT(m_swarm_size, 5, 1000, 20));
    }

    void SLPSO::initialize_(Environment *env) {
        Algorithm::initialize_(env);
        if (m_maximum_evaluations <= 0) {
            throw Exception("Maximum number of evaluations must be provided.");
        }
        m_maxW = 0.89;
        m_minW = 0.4;
        m_accelerator1 = 1.496;
        m_accelerator2 = 1.496;
    }

    void SLPSO::run_(Environment *env) {
        m_swarm.resize(m_swarm_size, env);
#ifdef OFEC_DATUM_SL_SWARM_H
        g_sl_swarm.value = &m_swarm;
#endif // OFEC_DATUM_SL_SWARM_H
#ifdef OFEC_DATUM_MULTI_POP_H
        g_multi_pop.pops.clear();
        g_multi_pop.pops.resize(1);
        for (size_t i = 0; i < m_swarm.size(); ++i) {
            g_multi_pop.pops[0].push_back(&m_swarm[i].pbest());
        }
#endif // OFEC_DATUM_MULTI_POP_H
#ifdef OFEC_DATUM_OFFSPRING_H
        g_offspring.sols.clear();
        for (size_t i = 0; i < m_swarm.size(); i++) {
            g_offspring.sols.push_back(&m_swarm[i].solution());
        }
#endif // OFEC_DATUM_OFFSPRING_H
        m_swarm.weight() = m_maxW;
        m_swarm.accelerator1() = m_accelerator1;
        m_swarm.accelerator2() = m_accelerator2;
        m_swarm.initialize(env, m_random.get());
        m_swarm.evaluate(env);
        m_swarm.initPbest(env);
        m_swarm.initVelocity(env, m_random.get());
        m_swarm.initVelocityMax(env, m_random.get());
        m_swarm.setParameters(m_random.get());
#ifdef OFEC_DATUM_SL_SWARM_H
        datumUpdated(env, g_sl_swarm);
#endif // OFEC_DATUM_SL_SWARM_H
#ifdef OFEC_DATUM_MULTI_POP_H
        datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
        while (!terminating()) {
            m_swarm.calculateNumLearning(0, this);
            m_swarm.evolve(env, m_random.get());
            m_swarm.weight() = m_maxW - (m_maxW - m_minW) * m_evaluations / m_maximum_evaluations;
#ifdef OFEC_DATUM_SL_SWARM_H
            datumUpdated(env, g_sl_swarm);
#endif // OFEC_DATUM_SL_SWARM_H
#ifdef OFEC_DATUM_OFFSPRING_H
            datumUpdated(env, g_offspring);
#endif // OFEC_DATUM_OFFSPRING_H
#ifdef OFEC_DATUM_MULTI_POP_H
            datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
        }
    }
}