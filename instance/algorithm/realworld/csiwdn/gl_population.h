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

#ifndef OFEC_GL_POPULATION
#define OFEC_GL_POPULATION

#include "csiwdn_individual.h"
#include "../../../../core/algorithm/population.h"
#include "../../../../core/environment/environment.h"

namespace ofec {

	class GLPopulation : public Population<IndCSIWDN>
	{
	protected:
		std::vector<std::vector<std::pair<Real, size_t>>> m_node_data_obj;
		std::vector<std::vector<Real>> m_probability;
		//std::vector<std::vector<int>> m_offspring;
		std::vector<size_t> m_node_cluster;
		std::list<std::unique_ptr<SolutionType>> m_best;

	public:
		GLPopulation(size_t no,  Environment *env, size_t dim);
		GLPopulation & operator=(GLPopulation & pop) {
			Population::operator=(pop);	
			m_probability = pop.m_probability;
			m_node_data_obj = pop.m_node_data_obj;

			return *this;
		}
		GLPopulation & operator=(GLPopulation && pop) {
			Population::operator=(std::move(pop));
			m_probability = std::move(pop.m_probability);
			m_node_data_obj = std::move(pop.m_node_data_obj);
			
			return *this;
		}


		std::list<std::unique_ptr<SolutionType>>& best() {
			return m_best;
		}
		void initialize(Environment *env, Random *rnd) override;
		void evolve(Environment *env, Random *rnd,  bool is_stable, const std::pair<int, int>& source_index);

		void mutate(Environment *env, Random *rnd, const std::pair<int, int>& source_index);
		void updateProbability(Environment *env, const std::pair<int, int>& source_index);
		
		void select(Environment *env,  bool is_stable, const std::pair<int, int>& source_index);
		void fillSolution(VarCSIWDN &indi, Environment *env, const std::pair<int, int>& source_index);
		bool isFeasiblePopulation(Environment *env, const Real tar);



		void updateBest(Environment* env) {
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

#endif // !OFEC_GL_POPULATION

