#include "sade_pop.h"

#include <algorithm>

namespace ofec {
	PopSaDE::PopSaDE(size_t size_pop, Environment *env) : 
		PopulationDE(size_pop, env), 
		m_num_strategy(4),
		m_LP(20),
		m_epsilon(0.01),
		mv_F(size_pop),
		mvv_CR(size_pop, std::vector<Real>(m_num_strategy)), 
		m_CRm(m_num_strategy, 0.5), 
		m_probability(m_num_strategy, 1. / m_num_strategy),
		m_strategy_selection(size_pop) 
	{
		for (size_t i = 0; i < m_probability.size(); ++i) {
			if (i > 0) m_probability[i] += m_probability[i - 1];
		}
	}

	int PopSaDE::evolve(Environment *env, Random *rnd) {
		std::vector<size_t> ridx;
		Real K;		
		updateCR(rnd);
		updateF(rnd);
		updateBest(env);
		for (size_t i = 0; i < size(); ++i) {
			Real p = rnd->uniform.next() * m_probability[m_num_strategy - 1];
			m_strategy_selection[i] = lower_bound(m_probability.begin(), m_probability.end(), p) - m_probability.begin();
			switch (m_strategy_selection[i]) {
			case 0:		// DE/rand/1/bin
				select(i, 3, ridx, rnd);
				m_individuals[i]->mutate(
					m_scaling_factor, 
					m_individuals[ridx[0]].get(), 
					m_individuals[ridx[1]].get(), 
					m_individuals[ridx[2]].get(), 
					env
				);
				m_individuals[i]->recombine(mvv_CR[i][m_strategy_selection[i]], de::RecombineStrategy::kBinomial, rnd, env);
				break;
			case 1:		// DE/rand-to-best/2/bin
				select(i, 4, ridx, rnd);
				m_individuals[i]->mutate(
					m_scaling_factor, 
					m_individuals[i].get(),
					m_best_individual,
					m_individuals[i].get(), 
					env, 
					m_individuals[ridx[0]].get(), 
					m_individuals[ridx[1]].get(),
					m_individuals[ridx[2]].get(), 
					m_individuals[ridx[3]].get()
				);
				m_individuals[i]->recombine(mvv_CR[i][m_strategy_selection[i]], de::RecombineStrategy::kBinomial, rnd, env);
				break;
			case 2:		// DE/rand/2/bin
				select(i, 5, ridx, rnd);
				m_individuals[i]->mutate(
					m_scaling_factor, 
					m_individuals[ridx[0]].get(), 
					m_individuals[ridx[1]].get(), 
					m_individuals[ridx[2]].get(), 
					env,
					m_individuals[ridx[3]].get(),
					m_individuals[ridx[4]].get()
				);
				m_individuals[i]->recombine(mvv_CR[i][m_strategy_selection[i]], de::RecombineStrategy::kBinomial, rnd, env);
				break;
			case 3:		// DE/current-to-rand/1
				K = rnd->uniform.next();
				select(i, 3, ridx, rnd);
				m_individuals[i]->mutate(
					K, 
					m_scaling_factor, 
					m_individuals[i].get(),
					m_individuals[ridx[0]].get(),
					m_individuals[i].get(), 
					m_individuals[ridx[1]].get(), 
					m_individuals[ridx[2]].get(), 
					env
				);
				break;
			}
		}
		int tag = kNormalEval;
		for (size_t i = 0; i < size(); ++i) {
			tag = m_individuals[i]->select(env);
			if (!(tag & kNormalEval)) return tag;
		}
		if (tag & kNormalEval) {
			m_iteration++;
		}
		updateMemory();
		return tag;
	}

	void PopSaDE::updateF(Random *rnd) {
		for (size_t i = 0; i < size(); i++) {
			mv_F[i] = rnd->normal.nextNonStd(0.5, 0.3);
		}
	}

	void PopSaDE::updateCR(Random *rnd) {
		if (m_iteration >= m_LP) {
			for (size_t k = 0; k < m_num_strategy; ++k) {
				std::vector<Real> t;
				for (auto it = m_CRsuc.begin(); it != m_CRsuc.end(); ++it) {
					for (auto i : it->at(k)) t.push_back(i);
				}
				if (!t.empty()) {
					std::nth_element(t.begin(), t.begin() + t.size() / 2, t.end());
					m_CRm[k] = t[t.size() / 2];
				}
			}
		}
		for (size_t i = 0; i < size(); i++) {
			for (size_t k = 0; k < m_num_strategy; ++k) {
				do {
					mvv_CR[i][k] = rnd->normal.nextNonStd(m_CRm[k], 0.01);
				} while (mvv_CR[i][k] < 0 || mvv_CR[i][k]>1);
			}
		}
	}

	void PopSaDE::updateMemory() {
		std::vector < std::list<Real>> curmem(m_num_strategy);
		std::vector<int> curSuc(m_num_strategy), curFail(m_num_strategy);
		for (size_t i = 0; i < size(); i++) {
			if (m_individuals[i]->isImproved()) {
				curmem[m_strategy_selection[i]].push_back(mvv_CR[i][m_strategy_selection[i]]);
				curSuc[m_strategy_selection[i]]++;
			}
			else {
				curFail[m_strategy_selection[i]]++;
			}
		}
		m_cnt_success.push_back(move(curSuc));
		m_CRsuc.push_back(move(curmem));
		m_cnt_fail.push_back(move(curFail));
		if (m_iteration >= m_LP) {
			m_cnt_success.pop_front();
			m_CRsuc.pop_front();
			m_cnt_fail.pop_front();
			//update probability for all stategies
			for (size_t k = 0; k < m_num_strategy; ++k) {
				m_probability[k] = 0;
				int fail = 0;
				for (auto& j : m_cnt_success) m_probability[k] += j[k];
				for (auto& j : m_cnt_fail) fail += j[k];
				m_probability[k] = m_probability[k] / (m_probability[k] + fail) + m_epsilon;
				if (k > 0) m_probability[k] += m_probability[k - 1];
			}
		}
	}
}
