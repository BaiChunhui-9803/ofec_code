#include "espde.h"
#include "../hts/hts_cde.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	void ESPDE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_NP, 5, 1000, 20));
		m_input_parameters.add("speciation frequency", new RangedSizeT(m_t, 5, 1000, 20));
		m_input_parameters.add("maximum size of cluster", new RangedSizeT(m_K, 1, 100, 8));
	}

	void ESPDE::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		if (m_maximum_evaluations <= 0) {
			throw Exception("Maximum number of evaluations must be provided.");
		}
		m_c = 0.1;
		m_positive = env->problem()->optimizeMode(0) == OptimizeMode::kMaximize ? 1 : -1;
		m_f_min = std::numeric_limits<Real>::max();
		m_f_max = -m_f_min;
	}

	void ESPDE::run_(Environment *env) {
		if ((initPop(env) & kNormalEval) == false)
			return;
		size_t G = 1;
		speciation(env);
#ifdef OFEC_DATUM_SEEDS_H
		g_seeds.sols.clear();
		for (size_t k = 0; k < m_S.size(); ++k) {
			SolutionBase *new_seed = nullptr;
			for (size_t i = 0; i < m_S[k].size(); ++i) {
				if (!new_seed || dominate(m_S[k][i], *new_seed, env->problem()->optimizeMode())) {
					new_seed = &m_S[k][i];
				}
			}
			if (new_seed) {
				g_seeds.sols.push_back(new_seed);
			}
		}
		datumUpdated(env, g_seeds);
#endif // OFEC_DATUM_SEEDS_H
		std::array<Real, 2> mu_F = { 0.5, 0.5 }, mu_CR = { 0.5, 0.5 };
		while (!terminating()) {
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(m_S.size());
			for (size_t k = 0; k < m_S.size(); k++) {
				for (size_t i = 0; i < m_S[k].size(); i++) {
					g_multi_pop.pops[k].push_back(&m_S[k][i]);
				}
			}
			datumUpdated(env, g_multi_pop);
#endif
			std::vector<std::vector<Real>> F(m_S.size()), CR(m_S.size());
			std::array<std::list<Real>, 2> A_F, A_CR;
			for (size_t i = 0; i < m_S.size(); ++i) {
				if (!m_S[i].isActive())
					continue;
				State state;
				if (sigmaNegInf(m_S[i], env) > 0.1) {
					state = kExplore;
					generateParam(mu_F[state], mu_CR[state], m_S[i], F[i], CR[i]);
					m_S[i].mutationStrategy() = de::MutateStrategy::kRand1;
				}
				else {
					if (m_positive * m_S[i].best(env)->objective(0) > gamma()) {
						state = kExploit;
						generateParam(mu_F[state], mu_CR[state], m_S[i], F[i], CR[i]);
						m_S[i].mutationStrategy() = de::MutateStrategy::kBest1;
					}
					else {
						m_S[i].setActive(false);
						continue;
					}
				}
				for (size_t j = 0; j < m_S[i].size(); ++j) {
					m_S[i].setParameter(CR[i][j], F[i][j]);
					m_S[i].mutate(j, m_random.get(), env);
					m_S[i].recombine(j, m_random.get(), env);
					m_S[i][j].trial().evaluate(env);
					updateFMinMax(m_S[i][j].trial());
				}
				for (size_t j = 0; j < m_S[i].size(); ++j) {
					int id_nearest = nearest(m_S[i], m_S[i][j].trial(), env);
					if (dominate(m_S[i][j].trial(), m_S[i][id_nearest], env->problem()->optimizeMode())) {
						m_EI.emplace_back(new Solution<>(m_S[i][id_nearest]));
						m_S[i][id_nearest].solution() = m_S[i][j].trial();
						A_F[state].push_back(F[i][j]);
						A_CR[state].push_back(CR[i][j]);
					}
					else {
						m_EI.emplace_back(new Solution<>(m_S[i][j].trial()));
					}
				}
			}
			updateMuParam(A_F[kExplore], A_CR[kExplore], mu_F[kExplore], mu_CR[kExplore]);
			updateMuParam(A_F[kExploit], A_CR[kExploit], mu_F[kExploit], mu_CR[kExploit]);
			if (G % m_t == 0) {
				m_P.clear();
				for (size_t i = 0; i < m_S.size(); ++i) {
					if (m_S[i].isActive()) {
						for (size_t j = 0; j < m_S[i].size(); ++j) {
							m_P.emplace_back(new Solution<>(m_S[i][j]));
						}
					}
				}
				speciation(env);
#ifdef OFEC_DATUM_SEEDS_H
				g_seeds.sols.clear();
				for (size_t k = 0; k < m_S.size(); ++k) {
					SolutionBase *new_seed = nullptr;
					for (size_t i = 0; i < m_S[k].size(); ++i) {
						if (!new_seed || dominate(m_S[k][i], *new_seed, env->problem()->optimizeMode())) {
							new_seed = &m_S[k][i];
						}
					}
					if (new_seed) {
						g_seeds.sols.push_back(new_seed);
					}
				}
				datumUpdated(env, g_seeds);
#endif // OFEC_DATUM_SEEDS_H
			}
			G++;
		}
	}

	int ESPDE::initPop(Environment *env) {
		m_P.clear();
		m_EI.clear();
		std::vector<std::shared_ptr<Solution<>>> samples;
		for (size_t i = 0; i < 2 * m_NP; i++) {
			samples.emplace_back(dynamic_cast<Solution<>*>(env->problem()->createSolution()));
			samples.back()->initialize(env, m_random.get());
			int rf = samples.back()->evaluate(env);
			updateFMinMax(*samples.back());
			if ((rf & kNormalEval) == false) {
				return rf;
			}
		}
		std::sort(samples.begin(), samples.end(),
			[this, env](const std::shared_ptr<Solution<>> &s1, const std::shared_ptr<Solution<>> &s2) {
				return dominate(*s1, *s2, env->problem()->optimizeMode());
			});
		for (size_t i = 0; i < m_NP; ++i) {
			m_P.push_back(samples[i]);
		}
		for (size_t i = m_NP; i < 2 * m_NP; ++i) {
			m_EI.push_back(samples[i]);
		}
		return kNormalEval;
	}

	int ESPDE::speciation(Environment *env) {
		std::set<size_t> id_P;
		for (size_t i = 0; i < m_P.size(); ++i) {
			id_P.insert(i);
		}
		std::vector<std::list<size_t>> S;
		while (!id_P.empty()) {
			auto it_best = id_P.begin();
			for (auto it = id_P.begin(); ++it != id_P.end();) {
				if (dominate(*m_P[*it], *m_P[*it_best], env->problem()->optimizeMode())) {
					it_best = it;
				}
			}
			size_t id_best = *it_best;
			id_P.erase(it_best);
			S.push_back({ id_best });
			std::vector<size_t> K_nearest;
			if (id_P.size() > m_K) {
				std::vector<std::pair<size_t, Real>> seq;
				for (size_t id : id_P) {
					seq.emplace_back(id, m_P[id_best]->variableDistance(*m_P[id], env));
				}
				std::nth_element(seq.begin(), seq.begin() + m_K - 1, seq.end(),
					[](const decltype(seq)::value_type &p1, const decltype(seq)::value_type &p2) {
						return p1.second < p2.second;
					});
				for (size_t i = 0; i < m_K; ++i) {
					K_nearest.push_back(seq[i].first);
				}
			}
			else {
				for (size_t id : id_P) {
					K_nearest.push_back(id);
				}
			}
			for (size_t i : K_nearest) {
				if (test(*m_P[id_best], *m_P[i], m_EI, env)) {
					S.back().push_back(i);
					id_P.erase(i);
				}
			}
		}
		m_S.clear();
		m_S.resize(S.size());
		for (size_t i = 0; i < S.size(); ++i) {
			for (size_t j : S[i]) {
				m_S[i].append(IndividualDE(*m_P[j]));
			}
			while (m_S[i].size() < 5) {
				Real sigma = pow(10, -1 - (10.0 / env->problem()->numberVariables() + 3) * m_evaluations / m_maximum_evaluations);
				m_S[i].append(IndividualDE(env));
				for (size_t j = 0; j < env->problem()->numberVariables(); ++j) {
					m_S[i].back().variable()[j] = m_random->normal.nextNonStd(m_S[i].front().variable()[j], sigma);
				}
				m_S[i].back().validate(env, Validation::kSetToBound);
				int rf = m_S[i].back().evaluate(env);
				updateFMinMax(m_S[i].back());
				if ((rf & kNormalEval) == false) {
					return rf;
				}
			}
		}
		return kNormalEval;
	}

	Real ESPDE::sigmaNegInf(const PopulationDE<> &S_i, Environment *env) const {
		std::vector<Real> mu(env->problem()->numberVariables(), 0);
		for (size_t j = 0; j < env->problem()->numberVariables(); ++j) {
			for (size_t k = 0; k < S_i.size(); ++k) {
				mu[j] += S_i[k].variable()[j];
			}
			mu[j] /= S_i.size();
		}
		std::vector<Real> delta(env->problem()->numberVariables(), 0);
		for (size_t j = 0; j < env->problem()->numberVariables(); ++j) {
			for (size_t k = 0; k < S_i.size(); ++k) {
				delta[j] += pow(S_i[k].variable()[j] - mu[j], 2);
			}
			delta[j] = sqrt(delta[j] / S_i.size());
		}
		return *std::min_element(delta.begin(), delta.end());
	}

	bool ESPDE::test(const Solution<> &x_a, const Solution<> &x_b,
		const std::list<std::shared_ptr<const Solution<>>> &EI, Environment *env
	) {
		std::list<std::shared_ptr<const Solution<>>> omega;
		std::vector<std::pair<Real, Real>> R = HTS_CDE::cubeRegion(x_a, x_b);
		for (auto &ei : EI) {
			if (HTS_CDE::isPointInOpenCubeRegion(*ei, R)) {
				omega.push_back(ei);
			}
		}
		if (omega.empty()) {
			return true;
		}
		else {
			for (auto &x_c : omega) {
				if (dominate(x_a, *x_c, env->problem()->optimizeMode()) && dominate(x_b, *x_c, env->problem()->optimizeMode())) {
					return false;
				}
			}
			return true;
		}
	}

	void ESPDE::updateMuParam(const std::list<Real> &A_F, const std::list<Real> &A_CR,
		Real &mu_F, Real &mu_CR
	) {
		if (A_F.empty()) {
			mu_F = m_random->uniform.next();
			mu_CR = m_random->uniform.next();
		}
		else {
			Real mean_F = 0, pow_mean_F = 0;
			for (Real F : A_F) {
				mean_F += F;
				pow_mean_F += pow(F, 2);
			}
			mean_F = pow_mean_F / mean_F;
			mu_F = (1 - m_c) * mu_F + m_c * mean_F;
			Real mean_CR = 0;
			for (Real CR : A_CR) {
				mean_CR += CR;
			}
			mean_CR /= A_CR.size();
			mu_CR = (1 - m_c) * mu_CR + m_c * mean_CR;
		}
	}

	void ESPDE::updateFMinMax(const Solution<> &s) {
		if (s.objective(0) < m_f_min) {
			m_f_min = s.objective(0);
		}
		if (s.objective(0) > m_f_max) {
			m_f_max = s.objective(0);
		}
	}

	size_t ESPDE::nearest(const PopulationDE<> &S_i, const Solution<> &s, Environment *env) const {
		size_t id_nearest = 0;
		Real dis_to_nearest = s.variableDistance(S_i[id_nearest], env);
		for (size_t k = 1; k < S_i.size(); ++k) {
			Real dis = s.variableDistance(S_i[k], env);
			if (dis < dis_to_nearest) {
				id_nearest = k;
				dis_to_nearest = dis;
			}
		}
		return id_nearest;
	}

	void ESPDE::generateParam(Real mu_F, Real mu_CR, const PopulationDE<> &S_i, 
		std::vector<Real> &F_i, std::vector<Real> &CR_i
	) {
		F_i.resize(S_i.size());
		CR_i.resize(S_i.size());
		for (size_t j = 0; j < S_i.size(); ++j) {
			do {
				F_i[j] = m_random->cauchy.nextNonStd(mu_F, 0.1);
			} while (F_i[j] <= 0);
			if (F_i[j] > 1) F_i[j] = 1;
			CR_i[j] = m_random->normal.nextNonStd(mu_CR, 0.1);
			if (CR_i[j] < 0) CR_i[j] = 0;
			else if (CR_i[j] > 1) CR_i[j] = 1;
		}
	}

	Real ESPDE::gamma() const {
		return (m_f_max - m_f_min) * m_evaluations / m_maximum_evaluations * 0.8 + m_f_min;
	}
}
