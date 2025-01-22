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

#ifndef SEQUENCE_ACTION_H
#define SEQUENCE_ACTION_H
#include"../../../../../core/definition.h"
#include"sequence_operator.h"
#include<vector>
#include<memory>

#include"../../framework/uncertianty/evaluation_strategy.h"


namespace ofec {



	class SequenceParametersBase {
	public:
		int m_problem.get() = -1;
		int this = -1;
		 SequenceParametersBase() = default;
		 virtual ~SequenceParametersBase() = default;
	};


	template<typename TSequenceOperator>
	class SequenceAction {
	public:
		using OpType = typename TSequenceOperator;
		using SolutionType = typename TSequenceOperator::SolutionType;
		using InterpreterType = typename TSequenceOperator::InterpreterType;
		using ParType = SequenceParametersBase;
	
	protected:
		bool m_flagStepLearn = true;
		std::shared_ptr<EvaluationStrategyBase<SolutionType>> m_eval_stra=nullptr;
		int m_id_pop = -1;
	public:

		bool flagStepLearn() const{
			return m_flagStepLearn;
		}

		SequenceAction() = default;
		~SequenceAction() = default;
		virtual void initialize(const SequenceParametersBase & par) {
			m_eval_stra = nullptr;
			m_id_pop = -1;
			m_flagStepLearn = false;
		}
		virtual void resize(const std::vector<int>& mat_size, double* initVal = nullptr) {}
		// return -1 means failed
		virtual int learn(
			SolutionType& cur,
			Random *rnd, Problem *pro,
			const InterpreterType& interpreter,
			const OpType& op,
			const std::vector<SolutionType*>& pop,
			int id_indiv
		)const {
			return -1;
		}

		virtual int globalLearn(Problem *pro, Algorithm *alg, Random *rnd,
			std::vector<SolutionType*>& pop,
			std::vector<SolutionType>& offspring) {
			return 0;
		}
		virtual void globalUpdate(Problem *pro, Random *rnd,
			const std::vector<SolutionType*>& pop) {}
		virtual void localUpdate(Problem *pro,
			const SolutionType& cur) {}


		virtual void setPopId(int popId) {
			m_id_pop = popId;
		}
		virtual void setEvalStrategy(const std::shared_ptr<EvaluationStrategyBase<SolutionType>>& eval_stra) {
			m_eval_stra = eval_stra;
		}
		virtual void updateMemory(Problem *pro, Random *rnd,
			const std::vector<SolutionType*>& pop) {
			globalUpdate(pro,rnd, pop);
		}
		virtual int createSolution(Problem *pro, Algorithm *alg, Random *rnd,
			std::vector<SolutionType*>& pop,
			std::vector<SolutionType>& offspring);
		//virtual int updateSolution(Problem *pro, Algorithm *alg, Random *rnd,
		//	std::vector<SolutionType*>& pop,
		//	std::vector<SolutionType>& offspring, int& num_improve);
	};

	template<typename TSequenceOperator>
	inline int SequenceAction<TSequenceOperator>::createSolution(Problem *pro, Algorithm *alg, Random *rnd,
		std::vector<SolutionType*>& pop,
		std::vector<SolutionType>& offspring) {
		
		if (m_flagStepLearn) {
			auto& interpreter(GET_ASeq(alg)->getInterpreter<InterpreterType>());
			auto& op(GET_ASeq(alg)->getOp<TSequenceOperator>());
			for (int offId(0); offId < offspring.size(); ++offId)
			{
				auto& it = offspring[offId];
				it.reset();
				interpreter.stepInit(pro, it);
				while (!interpreter.stepFinish(pro, it)) {
					int next = learn(it, rnd, pro, interpreter, op, pop, offId);
					if (next != -1) {
						interpreter.stepNext(pro, it, next);
						localUpdate(pro, it);
					}
					else {
						interpreter.stepBack(pro, it);
					}
				}
				//	op.learn_from_other(it, rnd, pro, interpreter, m_center, m_radius);
				interpreter.stepFinal(pro, it);
			}
		}

		return globalLearn(pro, alg, rnd,
			pop, offspring);

	}
	//template<typename TSequenceOperator>
	//inline int SequenceAction<TSequenceOperator>::updateSolution(Problem *pro, Algorithm *alg, Random *rnd,
	//	std::vector<SolutionType*>& pop,
	//	std::vector<SolutionType>& offspring, int& num_improve) {

	//	for (int idx(0); idx < pop.size(); ++idx) {
	//		if (offspring[idx].fitness() >= pop[idx]->fitness()) {
	//			*pop[idx] = offspring[idx];
	//			pop[idx]->setImprove(true);
	//			++num_improve;
	//		}
	//	}

	//	return 0;
	//}





}

#endif