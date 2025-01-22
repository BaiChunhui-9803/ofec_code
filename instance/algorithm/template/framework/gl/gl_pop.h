/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
* Framework of genetic learning (GL) algorithm
*
*********************************************************************************/

#ifndef OFEC_GL_H
#define OFEC_GL_H

#include "gl_adaptor.h"
#include <deque>
#include <limits>
#include <cmath>
#include <algorithm>
#include "../../../../../core/algorithm/population.h"

namespace ofec {
	template<typename TInd>
	class PopGL : virtual public Population<TInd> {
	public:
		// bsf : best so far of each Solution
		// hb: all historical best solutions of each Solution
		// ci: improved Solutions in the best so far population
		// c: all Solutions in the current population
		enum UpdateScheme { bsf, hb, ci, c };
	protected:
		UpdateScheme m_ms = UpdateScheme::bsf;			//memory scheme

		std::vector<Real> m_fitness;
		std::vector<Real> m_weight;
		std::vector<Real> m_objectives;
		std::vector<bool> m_improved;
		std::vector<bool> m_inds_active;

		Real m_memoryMaxObj, m_memoryMinObj, m_preMemoryMaxObj, m_preMemoryMinObj;
		Real m_wt = 1.e-3;			// weight threshold, Solutions with weight less than m_wt are removed from memory

		std::vector<std::deque<int>> m_exMemory; // explicit memory: the index of all personal best memory 

		std::vector<TInd> m_his;		//historical Solutions
		std::vector<TInd> m_offspring;
		std::unique_ptr<AdaptorGL<TInd>> m_adaptor;

	public:
		template<typename ... Args>
		void resize(size_t size_pop, Environment *env, Args&& ... args);
		void initializeMemory(Environment *env);
		//bool terminating();
		void updateMemory(Environment *env, bool excute = false);
		void updateMemoryHB(const std::vector<int> &index, Environment *env);
		void updateMemoryBSF(const std::vector<int> &index, Environment *env);
		void updateMemoryCI(const std::vector<int> &index, Environment *env);
		void updateMemoryC(Environment *env);
		int numImprove() {
			return m_num_improve;
		}
		void resetFlag() {
			for (int i = 0; i < this->size(); ++i)
				this->m_inds_active[i] = false;
		}
		void setAlpha(Real alpha) { m_alpha = alpha; }
		void setBeta(Real beta) { m_beta = beta; }
		void setGamma(Real gamma) { m_gamma = gamma; }
		void setUpdateScheme(int uspl) { m_ms = (UpdateScheme)uspl; }

		void setImproved(size_t id_ind, bool flag) { m_improved[id_ind] = flag; }
		bool isIndActive(size_t id_ind) const { return m_inds_active[id_ind]; }

		template<typename Adaptor>
		void setAdaptor() {
			m_adaptor.reset(new Adaptor());
		}
	protected:
		int update(Environment *env);
		int m_num_improve = 0;
		Real m_alpha = 0.5, m_beta = 3, m_gamma = 6;
	};

	template<typename TInd>
	template<typename ...Args>
	void PopGL<TInd>::resize(size_t size_pop, Environment *env, Args&& ...args) {
		Population<TInd>::resize(size_pop, env, std::forward<Args>(args)...);
		m_fitness.resize(size_pop);
		m_weight.resize(size_pop);
		m_objectives.resize(size_pop);
		m_exMemory.resize(size_pop);
		m_improved.resize(size_pop);
		m_inds_active.assign(size_pop, true);
	}

	//template<typename TInd>
	//void PopGL<TInd>::initialize_curpop() {
	//	m_offspring.resize(this->size());
	//	for (int i = 0; i < this->size(); i++) {
	//		m_offspring[i].initialize(this->m_pop[i]->id());
	//	}
	//}

	template<typename TInd>
	void PopGL<TInd>::initializeMemory(Environment *env) {
		m_memoryMaxObj = m_memoryMinObj = this->m_individuals[0]->objective(0);
		for (int i = 0; i < this->size(); ++i) {
			Real obj = this->m_individuals[i]->objective(0);
			if (obj > m_memoryMaxObj) m_memoryMaxObj = obj;
			if (obj < m_memoryMinObj) m_memoryMinObj = obj;
		}
		//std::vector<int> indiv(this->m_number_variables);
		Real gap = m_memoryMaxObj - m_memoryMinObj + 1e-5;
		for (int i = 0; i < this->size(); ++i) {
			m_objectives[i] = this->m_individuals[i]->objective(0);
			if (env->problem()->optimizeMode(0) == OptimizeMode::kMinimize)
				m_fitness[i] = (m_memoryMaxObj - m_objectives[i] + 1e-5) / gap;
			else
				m_fitness[i] = (m_objectives[i] - m_memoryMinObj + 1e-5) / gap;
			m_weight[i] = 1. / (1 + exp(-m_fitness[i]));
			m_exMemory[i].push_front(i);
			//for (int j = 0; j < this->m_number_variables; j++)
			//	indiv[j] = this->m_pop[i]->variable()[j];
			//m_hisIndi.push_back(indiv);
			m_his.push_back(*this->m_individuals[i]);
		}
		m_adaptor->updateProbability(env, *this, m_weight);
		m_preMemoryMaxObj = m_memoryMaxObj;
		m_preMemoryMinObj = m_memoryMinObj;
	}

	template<typename TInd>
	void PopGL<TInd>::updateMemory(Environment *env, bool excute) {

		std::vector<int> index;  //the updated Solution in the best so far

		if (env->problem()->optimizeMode(0) == OptimizeMode::kMinimize) {
			if (!excute) {
				m_memoryMinObj = m_preMemoryMinObj;
				for (int i = 0; i < this->size(); ++i) {
					if (this->m_improved[i]) {
						index.push_back(i);
						if (m_memoryMinObj > this->m_individuals[i]->objective(0))
							m_memoryMinObj = this->m_individuals[i]->objective(0);
					}
				}
			}
			else {
				m_memoryMinObj = std::numeric_limits<Real>::max();
				for (int i = 0; i < this->size(); ++i) {
					index.push_back(i);
					if (m_memoryMinObj > this->m_individuals[i]->objective(0))
						m_memoryMinObj = this->m_individuals[i]->objective(0);
				}
			}

		}
		else {
			if (!excute) {
				m_memoryMaxObj = m_preMemoryMaxObj;
				for (int i = 0; i < this->size(); ++i) {
					if (this->m_improved[i]) {
						index.push_back(i);
						if (m_memoryMaxObj < this->m_individuals[i]->objective(0))
							m_memoryMaxObj = this->m_individuals[i]->objective(0);
					}
				}
			}
			else {
				m_memoryMaxObj = -std::numeric_limits<Real>::max();
				for (int i = 0; i < this->size(); ++i) {
					index.push_back(i);
					if (m_memoryMaxObj < this->m_individuals[i]->objective(0))
						m_memoryMaxObj = this->m_individuals[i]->objective(0);
				}
			}
		}

		//if (m_ms != PopGL::c)
		//	if (index.empty())
		//		return;

		switch (m_ms)
		{
		case PopGL::bsf:
			updateMemoryBSF(index, env);
			break;
		case PopGL::hb:
			updateMemoryHB(index, env);
			break;
		case PopGL::ci:
			updateMemoryCI(index, env);
			break;
		case PopGL::c:
			updateMemoryC(env);
			break;
		default:
			break;
		}

	}

	template<typename TInd>
	void PopGL<TInd>::updateMemoryHB(const std::vector<int> &index, Environment *env) {
		if (m_memoryMinObj < m_preMemoryMinObj)
			m_preMemoryMinObj = m_memoryMinObj;
		Real gap = m_preMemoryMaxObj - m_preMemoryMinObj + 1e-5;
		if (env->problem()->optimizeMode(0) == OptimizeMode::kMinimize) {
			if (m_memoryMinObj < m_preMemoryMinObj) {
				for (int i = 0; i < m_weight.size(); ++i) {
					m_fitness[i] = (m_preMemoryMaxObj - m_objectives[i] + 1e-5) / gap;
					m_weight[i] = 1 / (1 + exp(-m_fitness[i]));
				}
			}
		}
		else {
			if (m_memoryMinObj < m_preMemoryMinObj) {
				for (int i = 0; i < m_weight.size(); ++i) {
					m_fitness[i] = (m_preMemoryMaxObj - m_objectives[i] + 1e-5) / gap;
					m_weight[i] = 1 / (1 + exp(-m_fitness[i]));
				}
			}
		}
		std::vector<Real> weight(index.size()), fitness(index.size());
		for (int i = 0; i < weight.size(); ++i) {
			fitness[i] = (m_preMemoryMaxObj - this->m_individuals[index[i]]->objective(0)+1e-5) / gap;
			weight[i] = 1 / (1 + exp(-fitness[i]));
		}
		int z = 0; //改进个体按顺序1，2，...出现
		for (int i = 0; i < m_exMemory.size(); ++i) {
			int exMemorySize = m_exMemory[i].size();
			for (int j = 0; j < exMemorySize; ++j) {
				if (this->m_improved[i]) {
					if (m_memoryMinObj < m_preMemoryMinObj)
						m_weight[m_exMemory[i][j]] = m_weight[m_exMemory[i][j]] / (pow((j + 2), 2));
					else
						m_weight[m_exMemory[i][j]] = m_weight[m_exMemory[i][j]] * pow((j + 1), 2) / pow((j + 2), 2);
					if (j == m_exMemory[i].size() - 1) {
						if (m_weight[m_exMemory[i][j]] / weight[z] < m_wt) {
							int pos = m_exMemory[i][j];
							m_weight[pos] = weight[z];
							m_objectives[pos] = this->m_individuals[i]->objective(0);
							m_fitness[pos] = fitness[z];
							m_exMemory[i].pop_back();
							m_exMemory[i].push_front(pos);
							for (int e1 = 0; e1 < this->m_individuals[i]->variable().size(); e1++)
								m_his[pos].variable()[e1] = this->m_individuals[i]->variable()[e1];
						}
						else {
							m_objectives.push_back(this->m_individuals[i]->objective(0));
							m_fitness.push_back(fitness[z]);
							m_weight.push_back(weight[z]);
							m_exMemory[i].push_front(m_weight.size() - 1);
							m_his.push_back(*this->m_individuals[i]);
						}
						z++;
					}
				}
				else {
					if (m_memoryMinObj < m_preMemoryMinObj)
						m_weight[m_exMemory[i][j]] = m_weight[m_exMemory[i][j]] / pow((j + 1), 2);
				}
			}
		}
		m_adaptor->updateProbability(env, m_his, m_weight);
		std::vector<Real>::iterator maxObj = std::max_element(m_objectives.begin(), m_objectives.end());
		m_preMemoryMaxObj = *maxObj;
	}

	template<typename TInd>
	void PopGL<TInd>::updateMemoryBSF(const std::vector<int> &index, Environment *env) {
		/*if (global::ms_global->m_problem->opt_mode(0)== optimization_mode::Minimization) {
			m_memoryMaxObj = this->m_pop[0]->objective(0);
			for (int i = 1; i < this->size(); ++i) {
				if (m_memoryMaxObj < this->m_pop[i]->objective(0))
					m_memoryMaxObj = this->m_pop[i]->objective(0);
			}
			if ((m_memoryMinObj < m_preMemoryMinObj) || (m_memoryMaxObj != m_preMemoryMaxObj))
			{
				m_preMemoryMinObj = m_memoryMinObj;
				m_preMemoryMaxObj = m_memoryMaxObj;

				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj + 1;
				for (int i = 0; i < this->size(); i++) {
					m_fitness[i] = (m_preMemoryMaxObj - this->m_pop[i]->objective(0) + 1) / gap;
					m_weight[i] = 1 / (1 + exp(-m_fitness[i]));
				}
			}
			else {
				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj + 1;
				for (int i = 0; i < index.size(); ++i)
				{
					m_fitness[index[i]] = (m_preMemoryMaxObj - this->m_pop[index[i]]->objective(0) + 1) / gap;
					m_weight[index[i]] = 1 / (1 + exp(-m_fitness[index[i]]));
				}
			}
		}
		else {
			m_memoryMinObj = this->m_pop[0]->objective(0);
			for (int i = 1; i < this->size(); ++i) {
				if (m_memoryMinObj > this->m_pop[i]->objective(0))
					m_memoryMinObj = this->m_pop[i]->objective(0);
			}

			if (m_memoryMaxObj > m_preMemoryMaxObj || m_memoryMinObj != m_preMemoryMinObj )
			{
				m_preMemoryMinObj = m_memoryMinObj;
				m_preMemoryMaxObj = m_memoryMaxObj;

				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj + 1;
				for (int i = 0; i < this->size(); i++) {
					m_fitness[i] = (this->m_pop[i]->objective(0) - m_preMemoryMinObj + 1) / gap;
					m_weight[i] = 1 / (1 + exp(-m_fitness[i]));
				}
			}
			else {
				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj + 1;
				for (int i = 0; i < index.size(); ++i)
				{
					m_fitness[index[i]] = (this->m_pop[index[i]]->objective(0) - m_preMemoryMinObj + 1) / gap;
					m_weight[index[i]] = 1 / (1 + exp(-m_fitness[index[i]]));
				}
			}
		}
		*/
		/*if (global::ms_global->m_problem->opt_mode(0)== optimization_mode::Minimization) {
			m_memoryMaxObj = this->m_pop[0]->objective(0);
			for (int i = 1; i < this->size(); ++i) {
				if (m_memoryMaxObj < this->m_pop[i]->objective(0))
					m_memoryMaxObj = this->m_pop[i]->objective(0);
			}
			if ((m_memoryMinObj < m_preMemoryMinObj) || (m_memoryMaxObj != m_preMemoryMaxObj))
			{
				m_preMemoryMinObj = m_memoryMinObj;
				m_preMemoryMaxObj = m_memoryMaxObj;

				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj;
				for (int i = 0; i < this->size(); i++) {
					if(gap>0)				m_fitness[i] = (m_preMemoryMaxObj - this->m_pop[i]->objective(0)) / gap;
					else m_fitness[i] = 0;
					m_weight[i] = 1 / (1 + exp(-m_fitness[i]));
				}
			}
			else {
				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj;
				for (int i = 0; i < index.size(); ++i)
				{
					if(gap>0)	m_fitness[index[i]] = (m_preMemoryMaxObj - this->m_pop[index[i]]->objective(0)) / gap;
					else m_fitness[index[i]] = 0;
					m_weight[index[i]] = 1 / (1 + exp(-m_fitness[index[i]]));
				}
			}
		}
		else {
			m_memoryMinObj = this->m_pop[0]->objective(0);
			for (int i = 1; i < this->size(); ++i) {
				if (m_memoryMinObj > this->m_pop[i]->objective(0))
					m_memoryMinObj = this->m_pop[i]->objective(0);
			}

			if (m_memoryMaxObj > m_preMemoryMaxObj || m_memoryMinObj != m_preMemoryMinObj)
			{
				m_preMemoryMinObj = m_memoryMinObj;
				m_preMemoryMaxObj = m_memoryMaxObj;

				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj;
				for (int i = 0; i < this->size(); i++) {
					if(gap>0)
					m_fitness[i] = (this->m_pop[i]->objective(0) - m_preMemoryMinObj) / gap;
					else m_fitness[i] = 0;
					m_weight[i] = 1 / (1 + exp(-m_fitness[i]));
				}
			}
			else {
				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj;
				for (int i = 0; i < index.size(); ++i)
				{
					if(gap>0)				m_fitness[index[i]] = (this->m_pop[index[i]]->objective(0) - m_preMemoryMinObj) / gap;
					else m_fitness[index[i]] = 0;
					m_weight[index[i]] = 1 / (1 + exp(-m_fitness[index[i]]));
				}
			}
		}*/
		if (env->problem()->optimizeMode(0) == OptimizeMode::kMinimize) {
			m_memoryMaxObj = this->m_individuals[0]->objective(0);
			for (int i = 1; i < this->size(); ++i) {
				if (m_memoryMaxObj < this->m_individuals[i]->objective(0))
					m_memoryMaxObj = this->m_individuals[i]->objective(0);
			}
			if ((m_memoryMinObj < m_preMemoryMinObj) || (m_memoryMaxObj != m_preMemoryMaxObj)) {
				m_preMemoryMinObj = m_memoryMinObj;
				m_preMemoryMaxObj = m_memoryMaxObj;

				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj;
				if (gap > 0) {
					for (int i = 0; i < this->size(); i++) {
						m_fitness[i] = (m_preMemoryMaxObj - this->m_individuals[i]->objective(0)) / gap;
						m_weight[i] = pow(m_fitness[i], m_beta) / (1 + exp(m_gamma - 12 * m_fitness[i]));
					}
				}
				else {
					for (int i = 0; i < this->size(); i++) {
						m_fitness[i] = 0;
						m_weight[i] = 1 / (1 + exp(-m_fitness[i]));
					}
				}

			}
			else {
				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj;
				if (gap > 0) {
					for (int i = 0; i < index.size(); ++i) {
						m_fitness[index[i]] = (m_preMemoryMaxObj - this->m_individuals[index[i]]->objective(0)) / gap;
						m_weight[index[i]] = pow(m_fitness[index[i]], m_beta) / (1 + exp(m_gamma - 12 * m_fitness[index[i]]));
					}
				}
				else {
					for (int i = 0; i < index.size(); ++i) {
						m_fitness[index[i]] = 0;
						m_weight[index[i]] = 1 / (1 + exp(-m_fitness[index[i]]));
					}
				}

			}
		}
		else {
			m_memoryMinObj = this->m_individuals[0]->objective(0);
			for (int i = 1; i < this->size(); ++i) {
				if (m_memoryMinObj > this->m_individuals[i]->objective(0))
					m_memoryMinObj = this->m_individuals[i]->objective(0);
			}

			if (m_memoryMaxObj > m_preMemoryMaxObj || m_memoryMinObj != m_preMemoryMinObj) {
				m_preMemoryMinObj = m_memoryMinObj;
				m_preMemoryMaxObj = m_memoryMaxObj;

				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj;
				if (gap > 0) {
					for (int i = 0; i < this->size(); i++) {
						m_fitness[i] = (this->m_individuals[i]->objective(0) - m_preMemoryMinObj) / gap;
						m_weight[i] = pow(m_fitness[i], m_beta) / (1 + exp(m_gamma - 12 * m_fitness[i]));
					}
				}
				else {
					for (int i = 0; i < this->size(); i++) {
						m_fitness[i] = 0;
						m_weight[i] = 1 / (1 + exp(-m_fitness[i]));
					}
				}
			}
			else {
				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj;
				if (gap > 0) {
					for (int i = 0; i < index.size(); ++i) {
						m_fitness[index[i]] = (this->m_individuals[index[i]]->objective(0) - m_preMemoryMinObj) / gap;
						m_weight[index[i]] = pow(m_fitness[index[i]], m_beta) / (1 + exp(m_gamma - 12 * m_fitness[index[i]]));
					}
				}
				else {
					for (int i = 0; i < index.size(); ++i) {
						m_fitness[index[i]] = 0;
						m_weight[index[i]] = 1 / (1 + exp(-m_fitness[index[i]]));
					}
				}
			}
		}
		m_adaptor->updateProbability(env, *this, m_weight);
	}

	template<typename TInd>
	void PopGL<TInd>::updateMemoryCI(const std::vector<int> &index, Environment *env) {
		m_preMemoryMaxObj = this->m_individuals[index[0]]->objective(0);
		m_preMemoryMinObj = this->m_individuals[index[0]]->objective(0);
		for (int i = 1; i < index.size(); ++i) {
			if (m_preMemoryMaxObj < this->m_individuals[index[i]]->objective(0))
				m_preMemoryMaxObj = this->m_individuals[index[i]]->objective(0);
			if (m_preMemoryMinObj > this->m_individuals[index[i]]->objective(0))
				m_preMemoryMinObj = this->m_individuals[index[i]]->objective(0);
		}
		Real gap = m_preMemoryMaxObj - m_preMemoryMinObj + 1e-5;
		for (int i = 0; i < index.size(); ++i) {
			if (env->problem()->optimizeMode(0) == OptimizeMode::kMinimize)
				m_fitness[index[i]] = (m_preMemoryMaxObj - this->m_individuals[index[i]]->objective(0) + 1e-5) / gap;
			else
				m_fitness[index[i]] = (this->m_individuals[index[i]]->objective(0) - m_preMemoryMinObj + 1e-5) / gap;
			m_weight[index[i]] = 1 / (1 + exp(-m_fitness[index[i]]));
		}
		m_adaptor->updateProbability(env, *this, m_weight, &index);
	}

	template<typename TInd>
	void PopGL<TInd>::updateMemoryC(Environment *env) {
		m_preMemoryMaxObj = m_offspring[0].objective(0);
		m_preMemoryMinObj = m_offspring[0].objective(0);
		for (int i = 1; i < this->size(); ++i) {
			if (m_preMemoryMaxObj < m_offspring[i].objective(0))
				m_preMemoryMaxObj = m_offspring[i].objective(0);
			if (m_preMemoryMinObj > m_offspring[i].objective(0))
				m_preMemoryMinObj = m_offspring[i].objective(0);
		}
		Real gap = m_preMemoryMaxObj - m_preMemoryMinObj + 1e-5;
		for (int i = 0; i < this->size(); ++i) {
			if (env->problem()->optimizeMode(0) == OptimizeMode::kMinimize)
				m_fitness[i] = (m_preMemoryMaxObj - m_offspring[i].objective(0) + 1e-5) / gap;
			else
				m_fitness[i] = (m_offspring[i].objective(0) - m_preMemoryMinObj + 1e-5) / gap;
			m_weight[i] = 1 / (1 + exp(-m_fitness[i]));
		}
		m_adaptor->updateProbability(env, m_offspring, m_weight);
	}

	//template<typename TInd>
	//bool PopGL<TInd>::terminating() {
	//	//TermMean * term = dynamic_cast<TermMean*>(this->m_termination.get());
	//	if (this->m_termination->terminating()) {
	//		if (this->m_iteration == 0) return false;
	//		return true;
	//	}
	//	return false;
	//}

	template<typename TInd>
	int PopGL<TInd>::update(Environment *env) {
		int rf = kNormalEval;
		rf = m_adaptor->updateSolution(env, *this, m_offspring, m_num_improve);
		return rf;
	}
}
#endif // !GENETIC_LEARNING_H

