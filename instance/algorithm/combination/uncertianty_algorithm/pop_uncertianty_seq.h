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

#ifndef OFEC_POP_UNCERTIANTY_SEQ_H
#define OFEC_POP_UNCERTIANTY_SEQ_H


#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

#include "../../template/framework/uncertianty/pop_uncertianty.h"
//#include "../ensemble_algo/adaptor/adaptor_seq.h"

#include "../ensemble_algo/ensemble/sequence_actions_ensemble.h"


namespace ofec {
	template <
		template<class> class TDistanceCalculator,
		class TSequenceAction>
		class PopulationUncertiantySeq :
		virtual public PopUncertianty<TDistanceCalculator, typename TSequenceAction::SolutionType> {
		public:
			using SequenceActionType = typename TSequenceAction;
			using OpType = typename TSequenceAction::OpType;
			using SolutionType = typename TSequenceAction::SolutionType;
			using InterpreterType = typename TSequenceAction::InterpreterType;

		public:

			
			virtual void setId(int id)override {
				PopUncertianty::setId(id);
				m_actions.setPopId(id);
			}
			virtual void setEvalStrategy(
				const std::shared_ptr<EvaluationStrategyBase<SolutionType>>& eval_stra)override {
				PopUncertianty::setEvalStrategy(eval_stra);
				m_actions.setEvalStrategy(eval_stra);
			}
			virtual void initParamenters(Problem *pro, Algorithm *alg)override {
				PopUncertianty::initParamenters(pro, alg);
				SequenceActionEnsembleParametersBase par;
				par.initGL_GA_LS(pro, alg);
				m_actions.initialize(par);
			}


			//AdaptorType& getAdaptor() {
			//	return dynamic_cast<AdaptorType&>(*m_adaptor);
			//}
			virtual void updateMemory()override;

			virtual void generateOffstrpings()override;

			//virtual int selectSurvivors() override;

			virtual void initialize(Problem *pro, Random *rnd) {
				PopUncertianty::initialize(pro, m_random.get());
				auto& interpreter(GET_ASeq(this).getInterpreter<InterpreterType>());
				auto& op(GET_ASeq(this).getOp<OpType>());
				for (auto& it : m_individuals) {
					it->reset();
					interpreter.stepFinal(pro, *it);
					interpreter.updateEdges(pro, *it);
				}
				updateSurvivors();
				m_eval_strategy->calFitness(m_survivors, m_id);
				m_offsprings.resize(m_individuals.size());
				for (int idx(0); idx < m_individuals.size(); ++idx) {
					m_offsprings[idx] = *m_individuals[idx];
				}

			}

		protected:
			//	int m_numImp = 0;
			int m_NotConvergeTime = 0;
			int m_NotConvergeTimeThread = 50;

			TSequenceAction m_actions;
			//AdaptorType m_adaptor;
			//std::unique_ptr<AdaptorSeq<OpType>> m_adaptor;
	};

	template <
		template<class> class TDistanceCalculator,
		class TIndi>
		void PopulationUncertiantySeq<TDistanceCalculator, TIndi>::updateMemory() {


		PopUncertianty::updateMemory();
		auto& interpreter(GET_ASeq(this).getInterpreter<InterpreterType>());
		for (auto& it : m_survivors) {
			interpreter.updateEdges(m_problem.get(), *it);
		}
		m_actions.globalUpdate(m_problem.get(), m_random.get(), m_survivors);
	}


	template <
		template<class> class TDistanceCalculator,
		class TIndi>
		void PopulationUncertiantySeq<TDistanceCalculator, TIndi>::generateOffstrpings() {
		m_actions.createSolution(m_problem.get(), this, m_random.get(),
			m_survivors, m_offsprings);
	}
}
#endif

