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
* Framework of genetic learning (PopGL) algorithm
*
*********************************************************************************/

#ifndef OFEC_POP_MP_SUB_UNCERTIANTY_SEQ_ENSEMBLE_H
#define OFEC_POP_MP_SUB_UNCERTIANTY_SEQ_ENSEMBLE_H


#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

#include "../pop_uncertainty_seq_ensemble.h"
#include <algorithm>


// for test


#include "../../../combination/sequence/SP/sp_operator.h"
#include "../../../template/framework/uncertianty/distance_calculator.h"
#include "../../../template/framework/uncertianty/evaluation_strategy.h"



namespace ofec {
	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		class PopulationMPSubUncertiantySeqEnsemble :
		virtual public PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator> {
		public:
			using OpType = typename TSequenceOperator;
			using SolutionType = typename TSequenceOperator::SolutionType;
			using InterpreterType = typename TSequenceOperator::InterpreterType;
			using PopulationType = typename PopulationMPSubUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>;
		protected:
			int m_min_size = 50;
		public:


			virtual int generateNewSolution() {
				if (size() >= m_min_size) {
					return 0;
				}
				PopUncertianty::initParamenters(m_problem.get(), this);
				SequenceActionEnsembleParametersBase par;
				par.initGL_GlobalNone(m_problem.get(), this);
				m_actions.initialize(par);
				int newlySize = m_min_size- size();
				resizeOffsprings(newlySize);
				int rf=evolve(m_problem.get(), this, m_random.get());
				resetCurState();


				par.initGL_GA_LS(m_problem.get(), this);
				m_actions.initialize(par);
				return rf;
			}

			virtual int selectSurvivors() override {



				if (size() > m_min_size) {
					//return PopulationUncertiantySeqEnsemble::selectSurvivors();
					auto& pop(m_individuals);
					int origin_size = size();
					mergeIndis(m_offsprings, m_min_size);
					return 0;
				}
				else if (size() == m_min_size) {
					return PopulationUncertiantySeqEnsemble::selectSurvivors();
				}
				else {

					
					auto& pop(m_individuals);
					int origin_size = size();
					mergeIndis(m_offsprings, m_min_size);
					m_num_improve = m_min_size - origin_size;
					for (int idx(0); idx < origin_size; ++idx) {
						pop[idx]->setImproved(false);
					}
					for (int idx(origin_size); idx < size(); ++idx) {
						pop[idx]->setImproved(true);
					}
					return 0;
				}
			}

			//void mergePopByFitnessDistance(PopulationType& other) {
			//	std::vector<std::unique_ptr<SolutionType>> inds;

			//	for (int idx(0); idx < m_individuals.size(); ++idx) {
			//		inds.emplace_back(m_individuals[idx].release());
			//		//swap(inds[idx], m_individuals[idx]);
			//	}
			//	int origin_size(m_individuals.size());
			//	for (int idx(0); idx < other.size(); ++idx) {
			//		//swap(inds[idx], other.m_individuals[idx]);
			//		inds.emplace_back(other.m_individuals[idx].release());
			//	}

			//	std::vector<int> order_idxs;
			//	std::vector<int> belong_idx;
			//	sortByFitnessDis(inds, order_idxs, belong_idx, m_problem.get());

			//	for (int idx(0); idx < m_individuals.size(); ++idx) {
			//		m_individuals[idx] = std::move(inds[order_idxs[idx]]);
			//	}
			//}
			void mergePopByFitnessDis(PopulationType& other);
		//	void mergePop(PopulationMPSubUncertiantySeqEnsemble& other);
		//	void mergeIndis(std::unique_ptr<SolutionType>& indis);
			void mergeIndis(std::vector<SolutionType>& indis, int size);
			

			virtual bool judgeSearchInside(PopulationType& other);
			virtual bool judgeSearchInside(const SolutionType& indi);
			virtual bool judgeSearchAreaInside(PopulationType& other);

			virtual int judgeRelationship(PopulationType& other);

			virtual bool judgeRelationship2(PopulationType& other);

			virtual void getMinDis() {
				
			}



		protected:
			void eraseEmptyIndis();
		protected:

			int m_min_iter = 30;

			double m_inner_ratio = 0.1;
			std::vector<double> m_min_dis;
	};


	//template <
	//	template<class> class TDistanceCalculator,
	//	class TSequenceOperator>
	//	inline void PopulationMPSubUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>
	//	::mergePop(PopulationMPSubUncertiantySeqEnsemble& other)
	//{
	//	int originSize(m_individuals.size());
	//	
	//	for (int idx(0); idx < other.size(); ++idx) {
	//		if (judgeInside(other[idx])) {
	//			m_individuals.emplace_back(other.m_individuals[idx].release());
	//		}
	//	}
	//	
	//	deleteIndisFit(originSize);
	//}


	//template <
	//	template<class> class TDistanceCalculator,
	//	class TSequenceOperator>
	//inline void PopulationMPSubUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
	//	mergeIndis(std::unique_ptr<SolutionType>& indi)
	//{
	//	if (judgeSearchInside(*indi)) {
	//		int worst_id(0);
	//		for (int idx(1); idx < m_individuals.size(); ++idx) {
	//			if (m_individuals[worst_id]->fitness() >= m_individuals[idx]->fitness()) {
	//				worst_id = idx;
	//			}
	//		}
	//		setSolution(worst_id, *indi);
	//		//swap(m_individuals[worst_id], indi);
	//	}

	//	resetCurState();
	//}
	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline bool PopulationMPSubUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
	 judgeRelationship2(PopulationType& other) {
		int total_num(0);
		for (auto& it : other.m_individuals) {
			if (m_distance_calculator->disToPop(*it) <= m_distance_calculator->radiusThreadhold()) {
				++total_num;
			}
		}

		bool flag= double(total_num) / double(other.m_individuals.size()) >= m_inner_ratio;
		if (flag) {
			int stop = -1;
		}
		return flag;
	}


	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline int PopulationMPSubUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
	    judgeRelationship(PopulationType& other) {
		int ret = 0; 

		if (this->m_cur_iters < m_min_iter || other.m_cur_iters < m_min_iter) {
			return ret;
		}

		if (m_distance_calculator->disToPop(other.getRepresentative()) <= m_distance_calculator->radiusThreadhold()&& judgeRelationship2(other)) {
			if (getRepresentative().fitness() >= other.getRepresentative().fitness()) {
				ret = 1;
			}
			else {
				ret = -1;
			}
		}
		return ret;
	}


	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline void PopulationMPSubUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
		mergePopByFitnessDis(PopulationType& other) {
		std::vector<std::unique_ptr<SolutionType>> inds;

		for (int idx(0); idx < m_individuals.size(); ++idx) {
			inds.emplace_back(m_individuals[idx].release());
			//swap(inds[idx], m_individuals[idx]);
		}
		int origin_size(m_individuals.size());
		for (int idx(0); idx < other.size(); ++idx) {
			//swap(inds[idx], other.m_individuals[idx]);
			inds.emplace_back(other.m_individuals[idx].release());
		}





		std::vector<int> fitness_order;
		std::vector<int> belong_idx;
		std::vector<int> distance_order;
		sortByFitnessDis(inds, fitness_order, belong_idx, distance_order,  m_problem.get());

		for (int idx(0); idx < m_individuals.size(); ++idx) {
			m_individuals[idx] = std::move(inds[distance_order[idx]]);
		}
		resetCurState();




	}


	// 	void mergeIndis(std::vector< SolutionType>& indis, int size);
	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline void PopulationMPSubUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
		mergeIndis(std::vector<SolutionType>& offspring, int size) {
		
		updateSurvivors();
		m_distance_calculator->updateMemory(m_survivors);

		double before_radius = m_distance_calculator->radiusThreadhold();

		
		std::vector<std::unique_ptr<SolutionType>> inds;
		for (int idx(0); idx < m_individuals.size(); ++idx) {
			inds.emplace_back(m_individuals[idx].release());
			//swap(inds[idx], m_individuals[idx]);
		}
		int origin_size(m_individuals.size());
		for (int idx(0); idx < offspring.size(); ++idx) {
			//swap(inds[idx], other.m_individuals[idx]);
			inds.emplace_back(new SolutionType(offspring[idx]));
		}




		std::vector<int> fitness_order;
		std::vector<int> belong_idx;
		std::vector<int> distance_order;
		sortByFitnessDis(inds, fitness_order, belong_idx, distance_order, m_problem.get());

		double before_avg_dis(0);
		for (int idx(0); idx < inds.size(); ++idx) {
			before_avg_dis += inds[idx]->variableDistance(*inds[belong_idx[idx]],m_problem.get());
		}
		before_avg_dis /= inds.size();

		size = std::min<size_t>(size, inds.size());
		m_individuals.resize(size);
		for (int idx(0); idx < size; ++idx) {
			//setSolution(idx, )
			m_individuals[idx] = std::move(inds[distance_order[idx]]);
		}

		sortByFitnessDis(m_individuals, fitness_order, belong_idx, distance_order, m_problem.get());

		double mavg(0);
		for (int idx(0); idx < m_individuals.size(); ++idx) {
			mavg += m_individuals[idx]->variableDistance(*m_individuals[belong_idx[idx]],m_problem.get());
		}
		mavg /= m_individuals.size();



		resize(size,m_problem.get(),m_problem->numberVariables());


		sortByFitnessDis(m_individuals, fitness_order, belong_idx, distance_order, m_problem.get());

		double after_avg_dis(0);
		for (int idx(0); idx < m_individuals.size(); ++idx) {
			after_avg_dis += m_individuals[idx]->variableDistance(*m_individuals[belong_idx[idx]],m_problem.get());
		}
		after_avg_dis /= m_individuals.size();



		resetCurState();

		updateSurvivors();
		m_distance_calculator->updateMemory(m_survivors);

		double after = m_distance_calculator->radiusThreadhold();
	}




	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
	inline bool PopulationMPSubUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
		judgeSearchInside(PopulationType& other)
	{
		return other.getRepresentative().fitness() < getRepresentative().fitness() && m_distance_calculator->disToPop(other.getRepresentative()) <= m_distance_calculator->radiusThreadhold();
	}


	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline bool PopulationMPSubUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
		judgeSearchAreaInside(PopulationType& other)
	{

		return other.getRepresentative().fitness() > getRepresentative().fitness() && m_distance_calculator->disToPop(other.getRepresentative()) <= m_distance_calculator->radiusThreadhold();

	}



	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline bool PopulationMPSubUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
		judgeSearchInside(const SolutionType& indi)
	{
		return m_distance_calculator->disToPop(indi) <= m_distance_calculator->radiusThreadhold();
	}

	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
	inline void PopulationMPSubUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::eraseEmptyIndis()
	{
		int before_idx(0);
		int cur_idx(0);
		for (; cur_idx < m_individuals.size(); ++cur_idx) {
			if (m_individuals[cur_idx] != nullptr && cur_idx != before_idx) {
				swap(m_individuals[cur_idx], m_individuals[before_idx++]);
			}
		}
		m_individuals.resize(before_idx);
	}



}
#endif

