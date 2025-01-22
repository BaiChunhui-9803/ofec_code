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

#ifndef OFEC_POP_MP_DIV_MAIN_UNCERTIANTY_SEQ_ENSEMBLE_H
#define OFEC_POP_MP_DIV_MAIN_UNCERTIANTY_SEQ_ENSEMBLE_H


#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

#include "../pop_uncertainty_seq_ensemble.h"
#include <algorithm>
#include <memory>


// for test
#include "../../../combination/sequence/SP/sp_operator.h"
#include "../../../template/framework/uncertianty/distance_calculator.h"
#include "../../../template/framework/uncertianty/evaluation_strategy.h"



namespace ofec {
	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		class PopulationMPDivMainUncertiantySeqEnsemble :
		virtual public PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator> {
		public:
			using OpType = typename TSequenceOperator;
			using SolutionType = typename TSequenceOperator::SolutionType;
			using InterpreterType = typename TSequenceOperator::InterpreterType;

		public:
			virtual void insert(Population<SolutionType> &pops);
			virtual void insert(const SolutionType &sol);

			virtual int selectSurvivors() override {

				if (m_flag_save_his) {
					for (int idx(0); idx < m_offsprings.size(); ++idx) {
						m_his.emplace_back(new SolutionType(m_offsprings[idx]));
					}
				}
				mergeByAttraction();
				if (m_random->uniform.next() > m_global_learn_pro) {
					m_type = 1;
				}
				else m_type = 0;

				return 0;
			}

			virtual void setEvolvoType(int type) {
				m_type = type;
			}
			virtual void initialize(Problem *pro, Random *rnd)override {
				PopulationUncertiantySeqEnsemble::initialize(pro, rnd);
				for (int idx(0); idx < m_individuals.size(); ++idx) {
					m_individuals[idx]->setId(idx);
				}
			}

			virtual int evolve(Problem *pro, Algorithm *alg, Random *rnd) override {
				//return 0;

				m_actions.setPro(m_action_pro[m_type]);
				return PopulationUncertiantySeqEnsemble::evolve(pro, alg, rnd);
			}

			void setSaveFlag(int flag) {
				m_flag_save_his = (flag);
				m_his.clear();
				if (m_flag_save_his) {
					m_his.resize(m_individuals.size());
					for (int idx(0); idx < m_individuals.size(); ++idx) {
						m_his[idx].reset(new SolutionType(*m_individuals[idx]));
					}
				}
			}

			bool getSaveFlag() {
				return m_flag_save_his;
			}
			std::vector<std::unique_ptr<SolutionType>>& getHis() {
				return m_his;
			}
			
		protected:
	
			virtual void mergeByAttraction();

		protected:

			int m_type = 0;
			std::vector<std::vector<double>> m_action_pro =
			{ {0.2,0.8,0} ,{1.0,0,0} };

			double m_global_learn_pro = 0.8;

			//std::vector<std::unique_ptr<Solution>> m_his;

			//std::vector<double> m_action_pro_0 = { 0.2,0.8,0 };
			//std::vector<double> m_action_pro_1 = {1.0,0,0};

			double m_inner_radius = 0.2;
			int m_numPop = 50;

			//std::vector<int> m_stagnated_times;

			std::vector<std::unique_ptr<SolutionType>> m_his;
			bool m_flag_save_his = false;
	};

	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline void PopulationMPDivMainUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
		insert(Population<SolutionType> &pop)
	{
		for (size_t i = 0; i < pop.size(); ++i) {
			m_individuals.emplace_back(std::move(pop->atPointer(i)));
		}
		reduceByFitness(m_numPop);
	}
	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline void PopulationMPDivMainUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
		insert(const SolutionType& sol) {
		m_individuals.resize(m_individuals.size() + 1);
		{
			m_individuals.back().reset(new SolutionType(sol));
		}
		reduceByFitness(m_numPop);
	}

	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline void PopulationMPDivMainUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>
		::mergeByAttraction()
	{

		int origin_size(m_individuals.size());
		m_individuals.resize(m_individuals.size() + m_offsprings.size());
		for (int idx(origin_size); idx < m_individuals.size(); ++idx) {
			m_individuals[idx].reset(new SolutionType(m_offsprings[idx- origin_size]));
		}
		reduceByFitness(m_numPop);
	}

}
#endif

