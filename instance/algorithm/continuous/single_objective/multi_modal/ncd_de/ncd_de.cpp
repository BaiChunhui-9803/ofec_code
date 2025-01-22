#include "ncd_de.h"
#include <numeric>

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	void NCD_DE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
	}

	void NCD_DE::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		m_scaling_factor = 0.5;
		m_crossover_rate = 0.9;
		m_lambda = 5;
	}

	void NCD_DE::run_(Environment *env) {
		m_pop.resize(m_pop_size, env);
		m_pop.initialize(env, m_random.get());
		m_pop.evaluate(env);
		Real positive = env->problem()->optimizeMode(0) == OptimizeMode::kMaximize ? 1 : -1;
		size_t gen = 0;
		while (!terminating()) {
			for (size_t i = 0; i < m_pop_size; ++i) {
				m_pop[i].setFitness(positive * m_pop[i].objective(0));
			}
			std::shared_ptr<Environment> internal_environment(Environment::create());
			auto NCD = NicheCenterDistinguish::create();
			NCD->setData(m_pop, env->problem());
			internal_environment->setProblem(NCD);
			internal_environment->initializeProblem();
			internal_environment->setAlgorithm(InternalGA::create());
			internal_environment->initializeAlgorithm(m_random);
			internal_environment->runAlgorithm();
			auto &Ch_best = dynamic_cast<InternalGA*>(internal_environment->algorithm())->
				bestChromosome(internal_environment.get());
			std::vector<std::vector<size_t>> niches;
			for (size_t i = 0; i < m_pop_size; ++i) {
				if (Ch_best[i] == 1) {
					niches.push_back({ i });
				}
			}
			if (!niches.empty()) {
				for (size_t i = 0; i < m_pop_size; ++i) {
					if (Ch_best[i] == 0) {
						size_t id_nearest_niche = 0;
						Real min_dis = m_pop[i].variableDistance(m_pop[niches[0].front()], env);
						if (niches.size() >= 1) {
							for (size_t k = 1; k < niches.size(); ++k) {
								Real dis = m_pop[i].variableDistance(m_pop[niches[k].front()], env);
								if (dis < min_dis) {
									min_dis = dis;
									id_nearest_niche = k;
								}
							}
						}
						niches[id_nearest_niche].push_back(i);
					}
				}
			}
#ifdef OFEC_DATUM_SEEDS_H
			g_seeds.sols.clear();
			for (size_t k = 0; k < niches.size(); ++k) {
				SolutionBase *new_seed = nullptr;
				for (size_t i = 0; i < niches[k].size(); ++i) {
					if (!new_seed || dominate(m_pop.at(niches[k][i]), *new_seed, env->problem()->optimizeMode())) {
						new_seed = &m_pop.at(niches[k][i]);
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
			g_multi_pop.pops.resize(niches.size());
			for (size_t k = 0; k < niches.size(); ++k) {
				for (size_t i = 0; i < niches[k].size(); ++i) {
					g_multi_pop.pops[k].push_back(&m_pop[niches[k][i]]);
				}
			}
			datumUpdated(env, g_multi_pop);
#endif
			for (auto &niche : niches) {
				if (niche.size() >= 3) {
					std::vector<size_t> nr;
					for (size_t i : niche) {
						m_pop.selectInCandidates(3, niche, nr, m_random.get());
						m_pop[i].mutate(m_scaling_factor, &m_pop[nr[0]], &m_pop[nr[1]], &m_pop[nr[2]], env);
						m_pop[i].recombine(m_crossover_rate, de::RecombineStrategy::kBinomial, m_random.get(), env);
						m_pop[i].trial().evaluate(env);
						size_t nn = nearestInd(m_pop[i].trial(), env);
						if (dominate(m_pop[i].trial(), m_pop[nn], env->problem()->optimizeMode())) {
							m_pop[nn].solution() = m_pop[i].trial();
						}
					}
				}
				else {
					for (size_t i : niche) {
						size_t nn = nearestInd(m_pop[i], env);
						if (m_random->uniform.next() <= 0.5) {
							Real dis = m_pop[i].variableDistance(m_pop[nn], env);
							for (size_t j = 0; j < env->problem()->numberVariables(); ++j) {
								m_pop[i].trial().variable()[j] = m_pop[i].variable()[j];
								if (m_random->uniform.next() <= 0.5) {
									m_pop[i].trial().variable()[j] += m_random->normal.next() * dis;
								}
							}
						}
						else {
							for (size_t j = 0; j < env->problem()->numberVariables(); ++j) {
								m_pop[i].trial().variable()[j] = m_pop[i].variable()[j]
									+ 0.5 * m_random->normal.next() * (m_pop[nn].variable()[j] - m_pop[i].variable()[j]);
							}
						}
						m_pop[i].trial().evaluate(env);
						if (dominate(m_pop[i].trial(), m_pop[i], env->problem()->optimizeMode())) {
							m_pop[i].solution() = m_pop[i].trial();
						}
					}
				}
			}
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(niches.size());
			for (size_t k = 0; k < niches.size(); ++k) {
				for (size_t i = 0; i < niches[k].size(); ++i) {
					g_multi_pop.pops[k].push_back(&m_pop[niches[k][i]]);
				}
			}
			datumUpdated(env, g_multi_pop);
#endif
			if (gen % m_lambda == 0) {
				std::vector<size_t> r;
				for (size_t i = 0; i < m_pop_size; ++i) {
					m_pop.select(i, 3, r, m_random.get());
					m_pop[i].mutate(m_scaling_factor, &m_pop[r[0]], &m_pop[r[1]], &m_pop[r[2]], env);
					m_pop[i].recombine(m_crossover_rate, de::RecombineStrategy::kBinomial, m_random.get(), env);
					m_pop[i].trial().evaluate(env);
					size_t nn = nearestInd(m_pop[i].trial(), env);
					if (dominate(m_pop[i].trial(), m_pop[nn], env->problem()->optimizeMode())) {
						m_pop[nn].solution() = m_pop[i].trial();
					}
				}
#ifdef OFEC_DATUM_MULTI_POP_H
				g_multi_pop.pops.clear();
				g_multi_pop.pops.resize(1);
				for (size_t i = 0; i < m_pop.size(); ++i) {
					g_multi_pop.pops[0].push_back(&m_pop[i]);
				}
				datumUpdated(env, g_multi_pop);
#endif
			}
			gen++;
		}
	}

	size_t NCD_DE::nearestInd(const Solution<> &s, Environment *env) {
		size_t nn = 0;
		Real min_dis = s.variableDistance(m_pop[0], env);
		for (size_t j = 1; j < m_pop_size; ++j) {
			Real dis = s.variableDistance(m_pop[j], env);
			if (dis < min_dis) {
				min_dis = dis;
				nn = j;
			}
		}
		return nn;
	}

	void NicheCenterDistinguish::setData(const PopulationDE<> &pop, Problem *original_problem) {
		m_pop = &pop;
		m_original_problem = original_problem;
		m_dis.assign(pop.size(), std::vector<Real>(pop.size(), 0));
		for (size_t i = 1; i < pop.size(); i++) {
			for (size_t j = 0; j < i; j++) {
				m_dis[i][j] = m_dis[j][i] = m_original_problem->variableDistance(
					(*m_pop)[i].variableBase(), (*m_pop)[j].variableBase());
			}
		}
	}

	void NicheCenterDistinguish::initializeVariables(VariableBase &x_, Random *rnd) const {
		auto &x = dynamic_cast<VariableVector<int>&>(x_);
		bool all_zero;
		do {
			all_zero = true;
			for (size_t i = 0; i < m_number_variables; ++i) {
				x[i] = rnd->uniform.next() < 0.5 ? 0 : 1;
				if (x[i] == 1) {
					all_zero = false;
				}
			}
		} while (all_zero);
	}

	void NicheCenterDistinguish::initialize_(Environment *env) {
		Problem::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
		m_number_variables = (*m_pop).size();
	}

	void NicheCenterDistinguish::evaluate(const VariableBase &vars, std::vector<Real> &objs, std::vector<Real> &cons) const {
		auto &Ch_i = dynamic_cast<const VariableVector<int>&>(vars);
		std::vector<size_t> niche_centers;
		size_t NP = m_number_variables;
		for (size_t j = 0; j < NP; ++j) {
			if (Ch_i[j] == 1) {
				niche_centers.push_back(j);
			}
		}
		Real FEM = 0;
		for (size_t j : niche_centers) {
			Real sum_p = 0;
			std::vector<Real> ps;
			for (size_t k : niche_centers) {
				if (k == j) {
					continue;
				}
				ps.push_back(exp(-m_dis[j][k]));
				sum_p += ps.back();
			}
			for (Real &p : ps)
				p /= sum_p;
			Real entropy = 0;
			for (Real p : ps) {
				entropy += -p * log(p);
			}
			entropy /= niche_centers.size();
			FEM += (*m_pop)[j].fitness() * entropy;
		}
		FEM /= pow(niche_centers.size(), 2);
		objs[0] = FEM;
		//Real FDM = 0;
		//for (size_t j : niche_centers) {
		//	Real distance = 0;
		//	for (size_t k : niche_centers) {
		//		if (m_dis[j][k] == -1) {
					//m_dis[j][k] = m_dis[k][j] = m_original_problem->variableDistance(
					//	(*m_pop)[j].variableBase(), (*m_pop)[k].variableBase());
		//		}
		//		distance += m_dis[j][k];
		//	}
		//	FDM += (*m_pop)[j].fitness() * distance;
		//}
		//FDM /= pow(niche_centers.size(), 2);
		//objs[0] = FDM;
	}

	void InternalGA::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		m_pop_size = 30;
		m_crossover_rate = 0.9;
		m_mutation_rate = 0.1;
		m_iteration_epochs = 5;
	}

	void InternalGA::run_(Environment *env) {
		m_pop.resize(m_pop_size, env, env->problem()->numberVariables());
		m_pop.initialize(env, m_random.get());
		m_pop.evaluate(env);
		Population<Solution<VariableVector<int>>> pop_and_off(m_pop_size * 2, env, env->problem()->numberVariables());
		std::vector<size_t> p(2);
		std::vector<Solution<VariableVector<int>>> o(2, Solution<VariableVector<int>>(1, 0, env->problem()->numberVariables()));
		for (size_t t = 0; t < m_iteration_epochs; ++t) {
			for (size_t i = 0; i < m_pop_size; ++i) {
				pop_and_off[i] = m_pop[i];
				p[0] = tournamentSelection(m_pop, env);
				do { p[1] = tournamentSelection(m_pop, env); } while (p[1] == p[0]);
				if (m_random->uniform.next() < m_crossover_rate) {
					crossover(m_pop[p[0]], m_pop[p[1]], o[0], o[1], env);
					pop_and_off[i + m_pop_size] = m_random->uniform.next() < 0.5 ? o[0] : o[1];
				}
				else {
					pop_and_off[i + m_pop_size] = dominate(m_pop[p[0]], m_pop[p[1]], env->problem()->optimizeMode()) ?
						m_pop[p[0]] : m_pop[p[1]];
				}
				mutate(pop_and_off[i + m_pop_size]);
				pop_and_off[i + m_pop_size].evaluate(env);
			}
			size_t id_best = 0;
			for (size_t i = 1; i < pop_and_off.size(); ++i) {
				if (dominate(pop_and_off[i], pop_and_off[id_best], env->problem()->optimizeMode())) {
					id_best = i;
				}
			}
			m_pop[0] = pop_and_off[id_best];
			for (size_t i = 1; i < m_pop_size; ++i) {
				m_pop[i] = pop_and_off[tournamentSelection(pop_and_off, env)];
			}
		}
	}

	size_t InternalGA::tournamentSelection(const Population<Solution<VariableVector<int>>> &pop, Environment *env) {
		if (m_rand_seq.size() != pop.size()) {
			m_rand_seq.resize(pop.size());
			std::iota(m_rand_seq.begin(), m_rand_seq.end(), 0);
		}
		m_random->uniform.shuffle(m_rand_seq.begin(), m_rand_seq.end());
		size_t idx_best = m_rand_seq[0];
		if (dominate(pop[m_rand_seq[1]], pop[idx_best], env->problem()->optimizeMode())) {
			idx_best = m_rand_seq[1];
		}
		return idx_best;
	}

	void InternalGA::crossover(const Solution<VariableVector<int>> &p1, const Solution<VariableVector<int>> &p2,
		Solution<VariableVector<int>> &o1, Solution<VariableVector<int>> &o2, Environment *env)
	{
		size_t j_rnd = m_random->uniform.nextNonStd<size_t>(0, env->problem()->numberVariables());
		for (size_t j = 0; j < j_rnd; ++j) {
			o1.variable()[j] = p1.variable()[j];
			o2.variable()[j] = p2.variable()[j];
		}
		for (size_t j = j_rnd; j < env->problem()->numberVariables(); ++j) {
			o1.variable()[j] = p2.variable()[j];
			o2.variable()[j] = p1.variable()[j];
		}
	}

	void InternalGA::mutate(Solution<VariableVector<int>> &s) {
		for (size_t j = 0; j < s.variable().size(); ++j) {
			if (m_random->uniform.next() < m_mutation_rate) {
				s.variable()[j] = 1 - s.variable()[j];
			}
		}
	}

	const VariableVector<int>& InternalGA::bestChromosome(Environment *env) {
		size_t id_best = 0;
		for (size_t i = 1; i < m_pop_size; ++i) {
			if (dominate(m_pop[i], m_pop[id_best], env->problem()->optimizeMode())) {
				id_best = i;
			}
		}
		return m_pop[id_best].variable();
	}
}
