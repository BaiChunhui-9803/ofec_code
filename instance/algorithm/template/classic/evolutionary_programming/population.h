/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created by Junchen Wang on Oct. 28, 2018.

#ifndef OFEC_EP_POPULATION_H
#define OFEC_EP_POPULATION_H

#include "../../../../../core/algorithm/population.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../../utility/functional.h"
#include "individual.h"
#include <numeric>

namespace ofec {
	template<typename TInd = IndEP>
	class PopEP : public Population<TInd> {
	public:
		using Population<TInd>::m_individuals;
		PopEP() = default;
		PopEP(size_t size_pop, Environment *env);
		PopEP(const PopEP &rhs);
		PopEP(PopEP &&rhs) noexcept;
		PopEP& operator=(const PopEP &rhs);
		PopEP& operator=(PopEP &&rhs) noexcept;
		void resize(size_t size, Environment *env);
		void initialize(Environment *env, Random *rnd) override;
		int evolve(Environment *env, Random *rnd) override;
		Real& tau() { return m_tau; }
		Real& tauPrime() { return m_tau_prime; }
		size_t& q() { return m_q; }
	protected:
		virtual void mutate(Random *rnd, Environment *env);
		virtual void select(Random *rnd, Environment *env);
		void resizeOffspring(Environment *env);
	protected:
		std::vector<std::unique_ptr<TInd>> m_offspring;
		Real m_tau = 0, m_tau_prime = 0;
		size_t m_q = 0;
	};

	template<typename TInd>
	PopEP<TInd>::PopEP(size_t size_pop, Environment *env) : 
		Population<TInd>(size_pop, CAST_CONOP(env->problem())->numberVariables())
	{
		resizeOffspring(env);
	}

	template<typename TInd>
	PopEP<TInd>::PopEP(const PopEP &rhs) : Population<TInd>(rhs) {
		for (size_t i = 0; i < rhs.m_offspring.size(); ++i) {
			m_offspring.emplace_back(new TInd(*rhs.m_offspring[i]));
		}
	}

	template<typename TInd>
	PopEP<TInd>::PopEP(PopEP &&rhs) noexcept : 
		Population<TInd>(std::move(rhs)),
		m_offspring(std::move(rhs.m_offspring)) {}

	template<typename TInd>
	PopEP<TInd> &PopEP<TInd>::operator=(const PopEP &rhs) {
		if (this == &rhs) return *this;
		Population<TInd>::operator=(rhs);
		m_offspring.resize(rhs.m_offspring.size());
		for (size_t i = 0; i < rhs.m_offspring.size(); ++i) {
			m_offspring[i].reset(new TInd(*rhs.m_offspring[i]));
		}
		return *this;
	}

	template<typename TInd>
	PopEP<TInd> &PopEP<TInd>::operator=(PopEP &&rhs) noexcept {
		if (this == &rhs) return *this;
		Population<TInd>::operator=(std::move(rhs));
		m_offspring = std::move(rhs.m_offspring);
		return *this;
	}

	template<typename TInd>
	void PopEP<TInd>::resize(size_t size, Environment *env) {
		Population<TInd>::resize(size, env, env->problem()->numberVariables());
	}

	template<typename TInd>
	void PopEP<TInd>::initialize(Environment *env, Random *rnd) {
		Population<TInd>::initialize(env, rnd);
		for (auto &ind : m_individuals) {
			ind->initializeEta(env);
		}
	}

	template<typename TInd>
	void PopEP<TInd>::mutate(Random *rnd, Environment *env) {
		for (size_t i = 0; i < m_individuals.size(); ++i)
			*m_offspring[i] = *m_individuals[i];
		for (size_t i = 0; i < m_individuals.size(); i++) {
			Real N = rnd->normal.next();
			for (size_t j = 0; j < m_offspring[i]->variable().size(); ++j) {
				m_offspring[i + m_individuals.size()]->variable()[j] = m_individuals[i]->variable()[j] + m_individuals[i]->eta()[j] * N;
				Real N_j = rnd->normal.next();
				m_offspring[i + m_individuals.size()]->eta()[j] = m_individuals[i]->eta()[j] * exp(m_tau_prime * N + m_tau * N_j);
			}
		}
	}

	template<typename TInd>
	void PopEP<TInd>::select(Random *rnd, Environment *env) {
		std::vector<size_t> wins(m_offspring.size(), 0);
		std::vector<size_t> rand_seq(m_offspring.size());
		std::iota(rand_seq.begin(), rand_seq.end(), 0);
		for (size_t i = 0; i < m_offspring.size(); ++i) {
			rnd->uniform.shuffle(rand_seq.begin(), rand_seq.end());
			for (size_t idx = 0; idx < m_q; idx++) {
				if (!dominate(*m_offspring[rand_seq[idx]], *m_offspring[i], env->problem()->optimizeMode())) {
					wins[i]++;
				}
			}
		}
		std::vector<int> index;
		mergeSort(wins, wins.size(), index, false);
		for (size_t i = 0; i < m_individuals.size(); ++i) {
			*m_individuals[i] = *m_offspring[index[i]];
		}
	}

	template<typename TInd>
	void PopEP<TInd>::resizeOffspring(Environment *env) {
		if (m_offspring.size() > 2 * m_individuals.size()) {
			m_offspring.resize(2 * m_individuals.size());
		}
		else if (m_offspring.size() < 2 * m_individuals.size()) {
			size_t num_vars = env->problem()->numberVariables();
			size_t num_objs = env->problem()->numberObjectives();
			size_t num_cons = env->problem()->numberConstraints();
			for (size_t i = m_offspring.size(); i < 2 * m_individuals.size(); ++i) {
				m_offspring.emplace_back(new TInd(num_objs, num_cons, num_vars));
			}
		}
	}

	template<typename TInd>
	int PopEP<TInd>::evolve(Environment *env, Random *rnd) {
		if (m_offspring.size() != 2 * m_individuals.size()) {
			resizeOffspring(env);
		}
		mutate(rnd, env);
		for (size_t i = m_individuals.size(); i < m_offspring.size(); ++i) {
			if (CAST_CONOP(env->problem())->boundaryViolated(*m_offspring[i])) {
				CAST_CONOP(env->problem())->validateSolution(*m_offspring[i], Validation::kSetToBound, rnd);
			}
			int tag = m_offspring[i]->evaluate(env);
			if (!(tag & kNormalEval)) {
				return tag;
			}
		}
		select(rnd, env);
		++this->m_iteration;
		return kNormalEval;
	}
}

#endif // !OFEC_EP_POPULATION_H
