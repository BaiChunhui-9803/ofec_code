/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
* class Algorithm is an abstract for all algorithms.
*
*********************************************************************************/
#ifndef OFEC_POP_UNCERTIANTY_H
#define OFEC_POP_UNCERTIANTY_H

#include<vector>
#include<memory>
#include "../../../../../core/algorithm/population.h"
#include "../../../../../core/algorithm/Solution.h"
#include "distance_calculator.h"
#include "evaluation_strategy.h"


namespace ofec {
	template <
		template<class> class TDistanceCalculator,
		class TIndi = Solution<>>
	class PopUncertianty: virtual public Population<TIndi> {
	public:
		using SolutionType = typename TIndi;
		//using DistanceCalculatorType = typename TDistanceCalculator;
	protected:
		int this = -1;
		int m_problem.get() = -1;
		int m_random.get() = -1;
		std::shared_ptr<EvaluationStrategyBase<SolutionType>> m_eval_strategy;
		std::unique_ptr<DistanceCalculatorBase<SolutionType>> m_distance_calculator;
		std::vector<SolutionType*> m_survivors;
		//bool m_flag_survivors_update = false;
		std::vector<SolutionType> m_offsprings;
		int m_num_improve = 0;
		int m_iteration = 0;
		int m_memory_iter = 0;

	public:


		PopUncertianty() = default;
		const std::vector<SolutionType*>& getSurvivors() {
			return m_survivors;
		}
		virtual void setSolution(int id, const SolutionType& indi) {
			*m_individuals[id] = indi;
			m_offsprings[id] = indi;
			m_survivors[id] = m_individuals[id].get();
		}
		



		template<typename ... Args>
		void resize(size_t size, Problem *pro, Args&& ...args){
			int origin_size(m_offsprings.size());
			Population::resize(size, pro,std::forward<Args>(args)...);
			updateSurvivors();
			m_offsprings.resize(size);
			if (size > origin_size) {
				for (int idx(origin_size); idx < size; ++idx) {
					m_offsprings[idx] = *m_individuals[idx];
				}
			}
		}


		virtual void setEvalStrategy(
			const std::shared_ptr<EvaluationStrategyBase<SolutionType>>& eval_stra) {
			m_eval_strategy = eval_stra;
		}
		const std::unique_ptr<DistanceCalculatorBase<SolutionType>>&
			getDistanceCalculator()const {
			return m_distance_calculator;
		}
		virtual void initParamenters(Problem *pro, Algorithm *alg) {
			this = alg;
			m_problem.get() = pro;
			m_random.get() = alg.idRandom();
			int id_param = alg.idParam();
			auto& param(GET_PARAM(id_param));
			//Population<TIndi>::assign(std::get<int>(param.at("population size")), m_problem.get());
			m_distance_calculator.reset(new TDistanceCalculator<SolutionType>());
			m_distance_calculator->initialize(pro, alg);
		}
		virtual void initialize(Problem *pro, Random *rnd)override {
			Population<TIndi>::initialize(pro, m_random.get());
			m_num_improve = 0;
			m_memory_iter = 0;
			m_iteration = 0;
		}



		virtual int evolve(Problem *pro, Algorithm *alg, Random *rnd) override {
			int rf = 0;
			++m_iteration;
			//updateMemory();
			updateSurvivors();
			generateOffstrpings();
			rf|= calculateFitness();
			rf|= selectSurvivors();
		//	updateBest(pro);
		//	updateMemory();
			return rf;
		}
		virtual void resizeOffsprings(int size) {
			m_offsprings.resize(size);
			for (int idx(0); idx < size; ++idx) {
				m_offsprings[idx] = *m_individuals.front();
			}
		}

		virtual void updateSurvivors() {
			m_survivors.resize(m_individuals.size());
			for (int idx(0); idx < m_individuals.size(); ++idx) {
				m_survivors[idx] = m_individuals[idx].get();
			}
			//if (!m_flag_survivors_update) {
			//	m_survivors.resize(m_individuals.size());
			//	for (int idx(0); idx < m_individuals.size(); ++idx) {
			//		m_survivors[idx] = m_individuals[idx].get();
			//	}
			//	m_flag_survivors_update = true;
			//}
		}
		virtual void updateMemory() {

			updateSurvivors();
			++m_memory_iter;
			m_distance_calculator->updateMemory(m_survivors);
			m_eval_strategy->updateMemory(m_survivors,m_id);
		}

		virtual int calculateFitness() {
			//updateSurvivors();
			return m_eval_strategy->calFitness(
				m_survivors, m_offsprings, m_id
			);
		}
		virtual void generateOffstrpings() = 0;
		virtual int selectSurvivors() {
			auto& pop(m_individuals);
			m_num_improve = 0;
			for (int idx(0); idx < pop.size(); ++idx) {
				
				if (m_offsprings[idx].fitness() >= pop[idx]->fitness()) {
					bool sameFlag(m_offsprings[idx].same(*pop[idx], m_problem.get()));
					setSolution(idx, m_offsprings[idx]);
					if (!sameFlag) {
						pop[idx]->setImproved(true);
						++m_num_improve;
					}
				}
				else pop[idx]->setImproved(false);
			}
			return 0;
		}
	};




}

#endif !OFEC_POP_UNCERTIANTY_H