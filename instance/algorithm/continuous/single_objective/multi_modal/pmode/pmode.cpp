#include "pmode.h"
#include <algorithm>

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	void PMODE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
	}

	void PMODE::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		if (m_maximum_evaluations <= 0) {
			throw Exception("Maximum evaluations must be provided.");
		}
		m_epsilon = 0.1;
		m_freq = 100;
		m_tau = 0.3;
		m_phi = CAST_CONOP(env->problem())->numberVariables() < 5 ? 1e-10 : 1e-8;
		m_tilde_R = std::numeric_limits<Real>::max();
		for (size_t j = 0; j < CAST_CONOP(env->problem())->numberVariables(); ++j) {
			auto &range = CAST_CONOP(env->problem())->range(j);
			if ((range.second - range.first) / 2 < m_tilde_R) {
				m_tilde_R = (range.second - range.first) / 2;
			}
		}	
		m_S.clear();
		m_pop.clear();
	}

	void PMODE::run_(Environment *env) {
		Real currentbest = -std::numeric_limits<Real>::max();
		size_t g = 0;
		std::vector<size_t> flag = { 0 };
		while (!terminating()) {
			updatePenaltyRadius(g, flag, env);
			std::vector<Real> index;
			size_t t = 0;
			m_pop.resize(m_pop_size, env);
			m_pop.initialize(env, m_random.get());
			m_pop.evaluate(env);
			for (size_t i = 0; i < m_pop_size; ++i) {
				m_pop[i].penalizedFitness() = getPenalizedFitness(m_pop[i], env);
			}
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(1);
			for (size_t i = 0; i < m_pop_size; ++i) {
				g_multi_pop.pops[0].push_back(&m_pop[i]);
			}
			datumUpdated(env, g_multi_pop);
#endif
			while (!terminating()) {
				if (t > 15 && abs(index[t - 1] - index[t - 16]) < m_phi) {
					break;
				}
				else {
					m_pop.evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
					g_multi_pop.pops.clear();
					g_multi_pop.pops.resize(1);
					for (size_t i = 0; i < m_pop_size; ++i) {
						g_multi_pop.pops[0].push_back(&m_pop[i]);
					}
					datumUpdated(env, g_multi_pop);
#endif
					Real max_Pf = m_pop[0].penalizedFitness();
					for (size_t i = 1; i < m_pop.size(); ++i) {
						if (max_Pf < m_pop[i].penalizedFitness()) {
							max_Pf = m_pop[i].penalizedFitness();
						}
					}
					index.push_back(max_Pf);
					t++;
				}
			}
			if (index.empty()) {
				break;
			}

			Real max_index = *std::max_element(index.begin(), index.end());
			if (max_index > currentbest) {
				currentbest = max_index;
			}
			g++;

			size_t best = 0;
			if (env->problem()->optimizeMode(0) == OptimizeMode::kMaximize) {
				for (size_t i = 1; i < m_pop.size(); ++i) {
					if (m_pop[i].objective(0) > m_pop[best].objective(0)) {
						best = i;
					}
				}
			}
			else {
				for (size_t i = 1; i < m_pop.size(); ++i) {
					if (m_pop[i].objective(0) < m_pop[best].objective(0)) {
						best = i;
					}
				}
			}
			
			updateEliteSet(m_pop[best], currentbest, env);
			flag.push_back(m_S.size());
		}
	}

	PenaltyIndDE::PenaltyIndDE(size_t num_objs, size_t num_cons, size_t num_vars) :
		IndividualDE(num_objs, num_cons, num_vars) {}

	int PenaltyIndDE::select(Environment *env) {
		int tag = m_pu.evaluate(env);
		Real pu_penalized_fitness = dynamic_cast<PMODE*>(env->algorithm())->getPenalizedFitness(m_pu, env);
		if (pu_penalized_fitness > m_penalized_fitness) {
			variable() = m_pu.variable();
			m_objectives = m_pu.objective();
			m_constraints = m_pu.constraint();
			m_penalized_fitness = pu_penalized_fitness;
			m_improved = true;
		}
		else {
			m_improved = false;
		}
		return tag;
	}
	
	void PMODE::updatePenaltyRadius(size_t g, const std::vector<size_t> &flag, Environment *env) {
		if (g == 0)
			m_R_g = 0.5 * m_tilde_R;
		else if (g > 0 && flag[g] == flag[g - 1]) {
			Real ratio = 1.0 - (Real)m_evaluations / m_maximum_evaluations;
			Real ampl = pow(ratio, 2);
			Real r_g = abs(ampl * sin(m_freq * OFEC_PI * ratio)) + m_tau * ampl;
			m_R_g = r_g * m_tilde_R;
		}
#ifdef OFEC_DATUM_RADIUS_REPULSION_H
		g_radius_repulsion.value = m_R_g;
		datumUpdated(env, g_radius_repulsion);
#endif // OFEC_DATUM_RADIUS_REPULSION_H

	}

	void PMODE::updateEliteSet(const PenaltyIndDE &best_ind, Real currentbest, Environment *env) {
		if (currentbest - best_ind.penalizedFitness() < m_phi) {
			bool add = true;
			for (auto &s : m_S) {
				if (best_ind.variableDistance(s, env) <= m_phi) {
					add = false;
					break;
				}
			}
			if (add) {
				m_S.push_back(best_ind);
				for (size_t i = 0; i < m_pop_size; ++i) {
					m_pop[i].penalizedFitness() = getPenalizedFitness(m_pop[i], env);
				}
			}
		}
	}

	Real PMODE::getPenalizedFitness(const Solution<> &s, Environment *env) const {
		Real penalized_fitness = s.objective(0);
		if (env->problem()->optimizeMode(0) == OptimizeMode::kMinimize) {
			penalized_fitness = -penalized_fitness;
		}
		for (auto &s_ : m_S) {
			Real d_n = s.variableDistance(s_, env);
			Real D_r = 1.0;
			if (d_n <= m_R_g) {
				if (penalized_fitness >= 0) {
					D_r = atan(m_epsilon * d_n);
				}
				else {
					D_r = 1.0 / atan(m_epsilon * d_n);
				}
			}
			penalized_fitness *= D_r;
		}
		return penalized_fitness;
	}
}