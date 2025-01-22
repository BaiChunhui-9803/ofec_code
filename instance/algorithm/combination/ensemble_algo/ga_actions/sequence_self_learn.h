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

#ifndef SEQUENCE_SELF_LEARN_H
#define SEQUENCE_SELF_LEARN_H
#include"../../../template/combination/sequence/sequence_action.h"
namespace ofec {
	template<typename TSequenceOperator>
	class SequenceSelfLearn :public SequenceAction <TSequenceOperator> {
	public:
		using OpType = TSequenceOperator;
		using SolutionType = typename TSequenceOperator::SolutionType;
		using InterpreterType = typename TSequenceOperator::InterpreterType;
		using ParType = typename  SequenceAction<TSequenceOperator>::ParType;

	public:
		SequenceSelfLearn() = default;
		~SequenceSelfLearn() = default;

		virtual void initialize(const SequenceParametersBase& par) override {
			m_flagStepLearn = true;
		}

		virtual int learn(
			SolutionType& cur,
			Problem *pro, Random *rnd,
			const InterpreterType& interpreter,
			const OpType& op,
			const std::vector<SolutionType*>& pop,
			int id_indiv
		)const override {
			/*
				* 			virtual int learn_from_other(
					const SolutionType& cur, Random *rnd, Problem *pro,
					const InterpreterType& interpreter,
					const SolutionType& other
				)const override;
				*/
			return op.learn_from_other(cur, rnd, pro, interpreter, *pop[id_indiv]);
		}

	};

}

#endif