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

#ifndef SEQUENCE_LEARN_OTHER_H
#define SEQUENCE_LEARN_OTHER_H
#include"../../../template/combination/sequence/sequence_action.h"

namespace ofec {
	template<typename TSequenceOperator>
	class SequenceLearnOther :public SequenceAction <TSequenceOperator> {
	public:
		using OpType = TSequenceOperator;
		using SolutionType = typename TSequenceOperator::SolutionType;
		using InterpreterType = typename TSequenceOperator::InterpreterType;
		using ParType = typename  SequenceAction<TSequenceOperator>::ParType;

	public:
		SequenceLearnOther() = default;
		~SequenceLearnOther() = default;
		virtual void initialize(const SequenceParametersBase& par) override {
			m_flagStepLearn = true;
		}
		virtual void globalUpdate(Problem *pro, Random *rnd,
			const std::vector<SolutionType*>& pop)override {
			m_learn_idxs.resize(pop.size());
			for (int idx(0); idx < pop.size(); ++idx) {
				m_learn_idxs[idx] = idx;
			}
			rnd->uniform.shuffle(m_learn_idxs.begin(), m_learn_idxs.end());
			//for(int idx(0);idx<)
		}


		virtual int learn(
			SolutionType& cur,
			Random *rnd, Problem *pro,
			const InterpreterType& interpreter,
			const OpType& op,
			const std::vector<SolutionType*>& pop,
			int id_indiv
		)const override {
			return op.learn_from_other(cur, rnd, pro, interpreter, *pop[m_learn_idxs[id_indiv]]);
		}
	protected:
		std::vector<int> m_learn_idxs;
	};

}

#endif