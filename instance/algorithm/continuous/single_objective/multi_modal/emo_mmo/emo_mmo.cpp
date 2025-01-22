#include "emo_mmo.h"
#include "../../../../template/selection/multi_objective/nsgaii.h"
#include "../../../../../../utility/functional.h"
#include <cmath>
#include <deque>
#include <fstream>

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	void EMO_MMO::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_size_pop, 5, 1000, 20));
	}

	void EMO_MMO::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		m_ratio_MOFLA_FEs = 0.5;
		m_eta = 0.1;
		m_D.clear();
		m_D_x.clear();
		m_P.clear();
		m_LS_pops.clear();
	}

	void EMO_MMO::run_(Environment *env) {
		m_MOFLA_pop.reset();
		size_t max_evals = m_maximum_evaluations > 0 ? m_maximum_evaluations : 1e6;
		m_t_max = m_ratio_MOFLA_FEs * ceil(max_evals / m_size_pop);
		MOFLA(env);
		binaryCuttingAPD(env);
		localSearch(env);
	}

	void EMO_MMO::MOFLA(Environment *env) {
		m_MOFLA_pop.reset(new PopSBX<>(m_size_pop, env));
		m_MOFLA_pop->setRate(0.9, 0.01);
		m_MOFLA_pop->setEta(50, 20);
		m_MOFLA_pop->initialize(env, m_random.get());
		m_MOFLA_pop->evaluate(env);
		for (size_t i = 0; i < m_MOFLA_pop->size(); ++i) {
			if (m_D_x.count(m_MOFLA_pop->at(i).variable().vector()) == 0) {
				m_D.emplace_back(std::make_unique<Solution<>>(m_MOFLA_pop->at(i)));
				m_D_x.insert(m_MOFLA_pop->at(i).variable().vector());
			}
		}
		for (size_t i = 0; i < m_MOFLA_pop->size(); i++) {
			m_MOFLA_pop->at(i).resizeObjective(2);
		}
		NSGAII nsgaii(2, { env->problem()->optimizeMode(0), OptimizeMode::kMaximize});
		Population<Solution<>> pop_combined(2 * m_size_pop, env, CAST_CONOP(env->problem())->numberVariables());
		for (size_t i = 0; i < pop_combined.size(); i++) {
			pop_combined[i].resizeObjective(2);
		}
		transferedMOP mop(m_t_max, CAST_CONOP(env->problem())->numberVariables(), 2 * m_size_pop);
		size_t t = 0;
		std::vector<size_t> rand_perm(m_MOFLA_pop->size());
		std::iota(rand_perm.begin(), rand_perm.end(), 0);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_D.size(); i++) {
			g_multi_pop.pops[0].push_back(m_D[i].get());
		}
		datumUpdated(env, g_multi_pop);
#endif
		while (!terminating() && t < m_t_max) {
			m_random->uniform.shuffle(rand_perm.begin(), rand_perm.end());
			for (size_t i = 0; i < m_MOFLA_pop->size(); i += 2) {
				m_MOFLA_pop->crossover(rand_perm[i], rand_perm[i + 1], pop_combined[i], pop_combined[i + 1], env, m_random.get());
				m_MOFLA_pop->mutate(pop_combined[i], env, m_random.get());
				m_MOFLA_pop->mutate(pop_combined[i + 1], env, m_random.get());
			}
			for (size_t i = 0; i < m_MOFLA_pop->size(); ++i) {
				pop_combined[i].evaluate(env);
				pop_combined[i + m_MOFLA_pop->size()] = m_MOFLA_pop->at(i);
			}
			mop.updateNormX(pop_combined);
			mop.updateDelta(t);
			mop.evaluate2ndObj(pop_combined);
			nsgaii.survivorSelection(*m_MOFLA_pop, pop_combined);
			for (size_t i = 0; i < m_MOFLA_pop->size(); ++i) {
				if (m_D_x.count(m_MOFLA_pop->at(i).variable().vector()) == 0) {
					m_D.emplace_back(std::make_unique<Solution<>>(m_MOFLA_pop->at(i)));
					m_D.back()->resizeObjective(1);           // only preserve the original objetive value
					m_D_x.insert(m_MOFLA_pop->at(i).variable().vector());
				}
			}
			t++;
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(1);
			for (size_t i = 0; i < m_D.size(); i++) {
				g_multi_pop.pops[0].push_back(m_D[i].get());
			}
			datumUpdated(env, g_multi_pop);
#endif
		}
	}

	void EMO_MMO::binaryCuttingAPD(Environment *env) {
		std::unordered_set<size_t> D_c;
		for (size_t i = 0; i < m_D.size(); i++)
			D_c.insert(i);
		landscapeCutting(D_c, m_eta, env);
		std::vector<std::vector<Real>> manh_dis_mat(m_D.size(), std::vector<Real>(m_D.size()));
		updateManhDisMat(manh_dis_mat, D_c);
		while (!D_c.empty()) {
			APD(manh_dis_mat, D_c, env);
			landscapeCutting(D_c, 0.5, env);
		}
	}
	 
	void EMO_MMO::localSearch(Environment *env) {
		auto &domain = CAST_CONOP(env->problem())->domain();
		size_t size_var = CAST_CONOP(env->problem())->numberVariables();
		std::vector<Real> size_dim(size_var);
		for (size_t j = 0; j < size_var; j++) {
			size_dim[j] = (domain[j].limit.second - domain[j].limit.first) * 0.05;
		}
		std::set<size_t> idxs_best;
		for (size_t k = 0; k < m_P.size(); k++) {
			size_t idx_best = *m_P[k].begin();
			for (size_t idx_ind : m_P[k]) {
				if (dominate(*m_D[idx_ind], *m_D[idx_best], env->problem()->optimizeMode())) {
					idx_best = idx_ind;
				}
			}
			if (idxs_best.count(idx_best) == 1) {
				continue;
			}
			else {
				idxs_best.insert(idx_best);
			}
			auto pop = std::make_unique<PopulationDE<>>(10, env);
			(*pop)[0].solution() = *m_D[idx_best];
			for (size_t j = 0; j < size_var; j++) {
				Real min_x_j = m_D[idx_best]->variable()[j] - size_dim[j] / 2;
				Real max_x_j = m_D[idx_best]->variable()[j] + size_dim[j] / 2;
				for (size_t i = 1; i < 10; i++) {
					(*pop)[i].variable()[j] = m_random->uniform.nextNonStd(min_x_j, max_x_j);
				}
			}
			for (size_t i = 1; i < 10; i++) {
				(*pop)[i].evaluate(env);
			}
			m_LS_pops.append(pop);
		}
#ifdef OFEC_DATUM_SEEDS_H
		g_seeds.sols.clear();
		for (size_t k = 0; k < m_LS_pops.size(); ++k) {
			SolutionBase *new_seed = nullptr;
			for (size_t i = 0; i < m_LS_pops[k].size(); ++i) {
				if (!new_seed || dominate(m_LS_pops[k][i], *new_seed, env->problem()->optimizeMode())) {
					new_seed = &m_LS_pops[k][i];
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
		g_multi_pop.pops.resize(m_LS_pops.size());
		for (size_t k = 0; k < m_LS_pops.size(); k++) {
			for (size_t i = 0; i < m_LS_pops[k].size(); i++) {
				g_multi_pop.pops[k].push_back(&m_LS_pops[k][i]);
			}
		}
		datumUpdated(env, g_multi_pop);
#endif
		while (!terminating()) {
			m_LS_pops.evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
		}
	}

	void EMO_MMO::landscapeCutting(std::unordered_set<size_t> &D_c, Real eta, Environment *env) {
		Real min_obj, max_obj;
		auto iter = D_c.begin();
		min_obj = max_obj = m_D[*iter]->objective(0);
		while (++iter != D_c.end()) {
			if (m_D[*iter]->objective(0) > max_obj) {
				max_obj = m_D[*iter]->objective(0);
			}
			if (m_D[*iter]->objective(0) < min_obj) {
				min_obj = m_D[*iter]->objective(0);
			}
		}
		Real threshold;
		if (env->problem()->optimizeMode(0) == OptimizeMode::kMaximize) {
			threshold = max_obj - eta * (max_obj - min_obj);
			for (auto iter = D_c.begin(); iter != D_c.end();) {
				if (m_D[*iter]->objective(0) > threshold) {
					++iter;
				}
				else {
					iter = D_c.erase(iter);
				}
			}
		}
		else {
			threshold = min_obj + eta * (max_obj - min_obj);
			for (auto iter = D_c.begin(); iter != D_c.end();) {
				if (m_D[*iter]->objective(0) < threshold) {
					++iter;
				}
				else {
					iter = D_c.erase(iter);
				}
			}
		}
		if (D_c.size() > 5000) {
			std::vector<size_t> D_c_cp;
			for (size_t idx_ind : D_c) {
				D_c_cp.push_back(idx_ind);
			}
			m_random->uniform.shuffle(D_c_cp.begin(), D_c_cp.end());
			D_c.clear();
			for (size_t i = 0; i < 5000; ++i) {
				D_c.insert(D_c_cp[i]);
			}
		}
	}

	void EMO_MMO::APD(std::vector<std::vector<Real>> &manh_dis_mat, std::unordered_set<size_t> D_c, Environment *env) {
		int k = 0;
		size_t num_B = D_c.size();
		normcManhDisMat(manh_dis_mat, D_c, env);
		while (num_B > 1) {
			k++;
			auto iter = D_c.begin();
			auto iter_sigma = iter;
			Real min_dis, sigma = minManhDis(manh_dis_mat, D_c, *iter);
			while (++iter != D_c.end()) {
				min_dis = minManhDis(manh_dis_mat, D_c, *iter);
				if (min_dis > sigma) {
					sigma = min_dis;
					iter_sigma = iter;
				}
			}
			std::list<size_t> psi_k;
			std::deque<size_t> neighbors = { *iter_sigma };
			D_c.erase(iter_sigma);
			while (!neighbors.empty()) {
				psi_k.push_back(neighbors.front());
				for (iter = D_c.begin(); iter != D_c.end();) {
					if (manh_dis_mat[neighbors.front()][*iter] == 0) {
						std::cout << "error" << std::endl;
						std::cout << "neighbors.front(): " << neighbors.front() << std::endl;
						std::cout << "*iter: " << *iter << std::endl;
						std::cout << m_D[neighbors.front()]->variable()[0] << ", " << m_D[neighbors.front()]->variable()[1] << std::endl;
						std::cout << m_D[*iter]->variable()[0] << ", " << m_D[*iter]->variable()[1] << std::endl;
					}
					if (manh_dis_mat[neighbors.front()][*iter] <= std::max<Real>(sigma, 1e-3)) {
						neighbors.push_back(*iter);
						iter = D_c.erase(iter);
					}
					else {
						iter++;
					}
				}
				num_B--;
				neighbors.pop_front();
			}
			m_P.push_back(std::move(psi_k));
		}
	}

	void EMO_MMO::updateManhDisMat(std::vector<std::vector<Real>> &manh_dis_mat, const std::unordered_set<size_t> &D_c) {
		for (auto iter1 = D_c.begin(); iter1 != D_c.end(); iter1++) {
			auto iter2 = iter1;
			iter2++;
			for (; iter2 != D_c.end(); iter2++) {
				manh_dis_mat[*iter1][*iter2] = manhattanDistance(m_D[*iter1]->variable().begin(), m_D[*iter1]->variable().end(), m_D[*iter2]->variable().begin());
				manh_dis_mat[*iter2][*iter1] = manh_dis_mat[*iter1][*iter2];
			}
		}
	}

	void EMO_MMO::normcManhDisMat(std::vector<std::vector<Real>> &manh_dis_mat, const std::unordered_set<size_t> &D_c, Environment *env) {
		std::vector<std::vector<Real>> normc_x(D_c.size(), std::vector<Real>(env->problem()->numberVariables()));
		for (size_t j = 0; j < env->problem()->numberVariables(); j++) {
			Real sum = 0;
			for (size_t idx_ind : D_c) {
				sum += pow(m_D[idx_ind]->variable()[j], 2);
			}
			Real beta = sqrt(1 / sum);
			size_t i = 0;
			for (size_t idx_ind : D_c) {
				normc_x[i][j] = beta * m_D[idx_ind]->variable()[j];
				i++;
			}
		}
		size_t i = 0;
		for (auto iter1 = D_c.begin(); iter1 != D_c.end(); iter1++) {
			size_t k = i + 1;
			auto iter2 = iter1;
			iter2++;
			for (; iter2 != D_c.end(); iter2++) {
				manh_dis_mat[*iter1][*iter2] = manhattanDistance(normc_x[i].begin(), normc_x[i].end(), normc_x[k].begin());
				manh_dis_mat[*iter2][*iter1] = manh_dis_mat[*iter1][*iter2];
				k++;
			}
			i++;
		}
	}

	Real EMO_MMO::minManhDis(const std::vector<std::vector<Real>> &manh_dis_mat, const std::unordered_set<size_t> &D_c, size_t i) {
		auto iter = D_c.begin();
		if (*iter == i) iter++;
		Real min_dis = manh_dis_mat[i][*iter];
		while (++iter != D_c.end()) {
			if (*iter == i) continue;
			if (manh_dis_mat[i][*iter] < min_dis) {
				min_dis = manh_dis_mat[i][*iter];
			}
		}
		return min_dis;
	}

	transferedMOP::transferedMOP(size_t t_max, size_t size_var, size_t size_pop) :
		m_t_max(t_max),
		m_delta(0),
		m_N(size_pop),
		m_size_var(size_var),
		m_norm_x(size_pop, std::vector<int>(size_var)),
		m_manh_dis(size_pop, std::vector<int>(size_pop, 0)) {}

	void transferedMOP::evaluate2ndObj(Population<Solution<>> &offspring) {
		for (size_t i = 0; i < m_N; i++) {
			int sum = 0, num = 0;
			for (size_t k = 0; k < m_N; k++) {
				if (m_manh_dis[i][k] < m_delta) {
					sum += m_manh_dis[i][k];
					num++;
				}
			}
			offspring[i].objective(1) = sum / m_delta - num;
		}
	}

	void transferedMOP::updateDelta(int t) {
		updateManhDis();
		int max = minManhDis(0);
		for (size_t i = 1; i < m_N; i++) {
			int min_dis = minManhDis(i);
			if (max < min_dis)
				max = min_dis;
		}
		m_delta = (1 - Real(t - 1) / m_t_max) * max;
	}

	void transferedMOP::updateNormX(const Population<Solution<>> &pop) {
		std::vector<Real> x_max(m_size_var), x_min(m_size_var);
		for (size_t j = 0; j < m_size_var; j++) {
			x_max[j] = x_min[j] = pop[0].variable()[j];
			for (size_t i = 1; i < m_N; i++) {
				if (x_max[j] < pop[i].variable()[j])
					x_max[j] = pop[i].variable()[j];
				if (x_min[j] > pop[j].variable()[j])
					x_min[j] = pop[j].variable()[j];
			}
		}
		for (size_t i = 0; i < m_N; i++) {
			auto &x = pop[i].variable();
			for (size_t j = 0; j < m_size_var; j++) {
				m_norm_x[i][j] = floor((m_N - 1) * (x[j] - x_min[j]) / (x_max[j] - x_min[j])) + 1;
			}
		}
	}

	void transferedMOP::updateManhDis() {
		for (size_t i = 0; i < m_N; i++) {
			for (size_t k = i + 1; k < m_N; k++) {
				m_manh_dis[i][k] = 0;
				for (size_t j = 0; j < m_size_var; j++) {
					m_manh_dis[i][k] += abs(m_norm_x[i][j] - m_norm_x[k][j]);
				}
				m_manh_dis[k][i] = m_manh_dis[i][k];
			}
		}
	}

	int transferedMOP::minManhDis(size_t i) {
		int min_dis = i != 0 ? m_manh_dis[i][0] : m_manh_dis[i][1];
		for (size_t k = 0; k < m_N; k++) {
			if (i == k) continue;
			if (min_dis > m_manh_dis[i][k])
				min_dis = m_manh_dis[i][k];
		}
		return min_dis;
	}
}
