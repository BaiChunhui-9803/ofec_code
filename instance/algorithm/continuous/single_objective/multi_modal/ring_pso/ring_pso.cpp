#include "ring_pso.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	void RingParticle::nextVelocity(const Solution<> *lbest, Real w, Real c1, Real c2, Random *rnd) {
		for (size_t j = 0; j < variable().size(); j++) {
			m_vel[j] = w * (m_vel[j] + 
				c1 * rnd->uniform.next() * (m_pbest.variable()[j] - variable()[j]) + 
				c2 * rnd->uniform.next() * (lbest->variable()[j] - variable()[j]));
		}
	}

	RingSwarm::RingSwarm(size_t size_pop, Topology topology, Environment *env) : 
		Swarm(size_pop, env),
		m_is_neighbor_set(false),
		m_topology(topology) {}

	void RingSwarm::setNeighborhood(Random *rnd) {
		if (m_is_neighbor_set) 
			return;
		for (size_t i = 0; i < m_individuals.size(); ++i) {
			switch (m_topology)
			{
			case Topology::R2:
				m_link[i][i] = true;
				m_link[i][(i + 1) % m_individuals.size()] = true;
				break;
			case Topology::R3:
				m_link[i][(i - 1) % m_individuals.size()] = true;
				m_link[i][i] = true;
				m_link[i][(i + 1) % m_individuals.size()] = true;
				break;
			case Topology::LHC_R2:
				for (size_t j = 0; j < 2 && (i / 2 * 2 + j) < m_individuals.size(); ++j)
					m_link[i][(i / 2 * 2 + j)] = true;
				break;
			case Topology::LHC_R3:
				for (size_t j = 0; j < 3 && (i / 3 * 3 + j) < m_individuals.size(); ++j)
					m_link[i][(i / 3 * 3 + j)] = true;
				break;
			default:
				break;
			}
		}
		m_is_neighbor_set = true;
	}

	void RingPSO::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
		m_input_parameters.add("ring topology", new Enumeration(m_topology,
			{ "R2","R3","LHC-R2","LHC-R3" }, RingSwarm::Topology::R2));
		m_input_parameters.add("weight", new RangedReal(m_weight, 0, 2, 0.7298));
		m_input_parameters.add("accelerator1", new RangedReal(m_accelerator1, 0, 3, 2.05));
		m_input_parameters.add("accelerator2", new RangedReal(m_accelerator2, 0, 3, 2.05));
	}

	void RingPSO::initialize_(Environment *env) {
		Algorithm::initialize_(env);
	}

	void RingPSO::run_(Environment *env) {
		initSwarm(env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
		while (!this->terminating()) {
			m_pop->evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
		}
	}

	void RingPSO::initSwarm(Environment *env) {
		m_pop.reset(new RingSwarm(m_pop_size, m_topology, env));
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop->size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop->at(i).pbest());
		}
#endif
		m_pop->weight() = m_weight;
		m_pop->accelerator1() = m_accelerator1;
		m_pop->accelerator2() = m_accelerator2;
		m_pop->initialize(env, m_random.get());
		m_pop->initVelocity(env, m_random.get());
		m_pop->evaluate(env);
		m_pop->initPbest(env);
		m_pop->setNeighborhood(m_random.get());
	}
}