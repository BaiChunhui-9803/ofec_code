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

#ifndef SEQUENCE_LOCAL_SEARCH_H
#define SEQUENCE_LOCAL_SEARCH_H
#include"../../../template/combination/sequence/sequence_action.h"
namespace ofec {
	template<typename TSequenceOperator>
	class SequenceLocalSearch :public SequenceAction <TSequenceOperator> {
	public:
		using OpType = TSequenceOperator;
		using SolutionType = typename TSequenceOperator::SolutionType;
		using InterpreterType = typename TSequenceOperator::InterpreterType;
		using ParType = typename  SequenceAction<TSequenceOperator>::ParType;
	protected:
		int m_localSearchTimes = 30;
	public:
		SequenceLocalSearch() = default;
		~SequenceLocalSearch() = default;
		virtual void initialize(const SequenceParametersBase& par) override {
			m_flagStepLearn = false;
		}


		virtual int globalLearn(Problem *pro, Algorithm *alg, Random *rnd,
			std::vector<SolutionType*>& pop,
			std::vector<SolutionType>& offspring) {
			auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
			auto& op(GET_ASeq(alg).getOp<TSequenceOperator>());

			int rf(0);
			if (!offspring.empty()) {
				int randId = rnd->uniform.nextNonStd<int>(0, offspring.size());
				rf |= op.localSearch(offspring[randId], m_localSearchTimes,
					pro, alg, alg.idRandom(),
					m_eval_stra, m_id_pop
				);

				auto& it(offspring[randId]);
				it.reset();
				interpreter.stepFinal(pro, it);

			}
			if (!pop.empty()) {

				int randId = rnd->uniform.nextNonStd<int>(0, pop.size());
				rf |= op.localSearch(*pop[randId], m_localSearchTimes,
					pro, alg, alg.idRandom(),
					m_eval_stra, m_id_pop
				);


				auto& it(*pop[randId]);
				it.reset();
				interpreter.stepFinal(pro, it);

			}
			return rf;
		}


	};

}

#endif