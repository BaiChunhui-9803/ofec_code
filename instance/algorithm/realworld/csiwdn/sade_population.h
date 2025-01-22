/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Li Zhou
* Email: 441837060@qq.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*********************************************************************************/
// updated Jan 19, 2019 by Li Zhou

#ifndef OFEC_SADE_POPULATION
#define OFEC_SADE_POPULATION

#include "csiwdn_individual.h"
#include "../../../../core/algorithm/population.h"
#include <list>
#include <algorithm>
#include "../../../../utility/functional.h"

namespace ofec {

	class SaDEPopulation : public Population<IndCSIWDN>
	{
		enum DEMutationSratgy { rand_1, best_1, target_to_best_1, best_2, rand_2, rand_to_best_1, target_to_rand_1 };
		DEMutationSratgy m_mutation_strategy;
	protected:
		const int m_num_strategy = 4;
		std::vector<std::vector<Real>> m_F;
		std::vector<std::vector<std::vector<Real>>> mvv_CR;

		std::vector<std::list<std::vector<std::list<Real>>>> m_CRsuc;

		std::vector<std::vector<Real>> m_CRm, m_pro_strategy;
		std::vector<std::list<std::vector<int>>> m_cnt_success, m_cnt_fail;
		int m_LP = 5;
		Real m_epsilon = 0.01;
		std::vector<std::vector<int>> m_strategy_selection;
		Violation m_mode = Violation::kBoundary;

		//for ST and duration
		std::vector<std::vector<Real>> m_probability;
		std::vector<std::pair<Real, size_t>> m_ST_data_obj;
		std::vector<std::pair<Real, size_t>> m_duration_data_obj;
		std::list<std::unique_ptr<SolutionType>> m_best;
		
		

	public:
		SaDEPopulation(size_t no, Environment *env, size_t dim);
		SaDEPopulation & operator=(SaDEPopulation & pop) {
			Population::operator=(pop);

			m_F = pop.m_F;
			mvv_CR = pop.mvv_CR;
			m_CRsuc = pop.m_CRsuc;
			m_CRm = pop.m_CRm;
			m_pro_strategy = pop.m_pro_strategy;
			m_cnt_success = pop.m_cnt_success;
			m_cnt_fail = pop.m_cnt_fail;
			m_strategy_selection = pop.m_strategy_selection;
			m_mode = pop.m_mode;
			m_probability = pop.m_probability;
			m_ST_data_obj = pop.m_ST_data_obj;
			m_duration_data_obj = pop.m_duration_data_obj;

			return *this;
		}
		SaDEPopulation & operator=(SaDEPopulation && pop) {
			Population::operator=(std::move(pop));

			m_F = std::move(pop.m_F);
			mvv_CR = std::move(pop.mvv_CR);
			m_CRsuc = std::move(pop.m_CRsuc);
			m_CRm = std::move(pop.m_CRm);
			m_pro_strategy = std::move(pop.m_pro_strategy);
			m_cnt_success = std::move(pop.m_cnt_success);
			m_cnt_fail = std::move(pop.m_cnt_fail);
			m_strategy_selection = std::move(pop.m_strategy_selection);
			m_mode = pop.m_mode;
			m_probability = std::move(pop.m_probability);
			m_ST_data_obj = std::move(pop.m_ST_data_obj);
			m_duration_data_obj = std::move(pop.m_duration_data_obj);

			return *this;
		}
		IndCSIWDN& operator[](size_t i) {
			return *m_individuals[i];
		}

		const IndCSIWDN& operator[](size_t i) const {
			return *m_individuals[i];
		}


		std::list<std::unique_ptr<SolutionType>>& best() {
			return m_best;
		}

		void mutate(Environment *env, Random *rnd, const std::pair<int, int>& source_index);
		void mutate_(Environment *env, Random *rnd, const int idx, const int jdx, Real F, const std::vector<std::vector<Real>> & prob);
		void recombine(Environment *env, Random *rnd, const std::pair<int, int>& source_index);

		void select(Random *rnd, int base, int number, std::vector<int>& result);
		void select(Environment *env, bool is_stable, const std::pair<int, int>& source_index);
		void setMutationStrategy(DEMutationSratgy rS);

		void updateF(Environment *env, Random *rnd, const std::pair<int, int>& source_index);
		void updateCR(Environment *env, Random *rnd, const std::pair<int, int>& source_index);
		void updateProStrategy(Environment *env, const std::pair<int, int>& source_index);
		void updateProbability(Environment *env);

		void fillSolution(VarCSIWDN &indi, Environment *env, const std::pair<int, int>& source_index);
		bool isFeasiblePopulation(Environment *env, const Real tar);


		void updateBest(Environment* env)  {
			for (auto& ptr_ind : m_individuals) {
				updateBest(*ptr_ind, env);
			}
		}


		bool updateBest(const SolutionType& sol, Environment* env) {
			bool is_best = true;
			auto pro = env->problem();
			for (auto iter = m_best.begin(); iter != m_best.end();) {
			//	auto result = dominate(sol,**iter, pro->optimizeMode());

				auto result = objectiveCompare(sol.objective(), (**iter).objective(), pro->optimizeMode());
				if (result == Dominance::kDominated || result == Dominance::kEqual) {
					is_best = false;
					break;
				}
				else if (result == Dominance::kDominant)
					iter = m_best.erase(iter);
				else
					iter++;
			}
			if (is_best) {
				m_best.emplace_back(new SolutionType(sol));
			}
			return is_best;
		}

	};

}

#endif // !OFEC_SADE_POPULATION

