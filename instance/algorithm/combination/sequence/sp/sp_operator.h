/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Yiya Diao
* Email: diaoyiyacug@163.com
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
* Operator adaptor for combinatorial problem that can be classified into the sequence problem
*
*********************************************************************************/




#ifndef OPERATOR_SELECTION_PROBLEM_H
#define OPERATOR_SELECTION_PROBLEM_H


#include"../../../template/combination/sequence/sequence_operator.h"
#include"SP_interpreter.h"
namespace ofec {
		class SelectionProblemOperator: public SequenceOperator<SP_Interpreter> {
		public:
			enum class EFitnessCalType { Objective, HLS };
			enum class EMemoryType  {Decrease,Change};
			
		public:
			using SolutionType = typename  SequenceOperator::SolutionType;
			using InterpreterType = typename SequenceOperator::InterpreterType;
		public:

			
			SelectionProblemOperator() = default;
            ~SelectionProblemOperator() = default;


			virtual void initialize(Problem *pro, Algorithm *alg)override;
			virtual int learn_from_local(
				const SolutionType& cur, Random *rnd, Problem *pro,
				const InterpreterType& interpreter
			)const override;
			virtual int learn_from_global(
				const SolutionType& cur, Random *rnd, Problem *pro,
				const InterpreterType& interpreter,
				const std::function<Real(const SolutionType& cur, int to)>& proFun
			)const override;
			virtual int learn_from_other(
				const SolutionType& cur, Random *rnd, Problem *pro,
				const InterpreterType& interpreter,
				const SolutionType& other
			)const override;



			virtual bool learn_from_other(
				SolutionType& cur, Random *rnd, Problem *pro,
				const InterpreterType& interpreter,
				const SolutionType& other,
				Real radius
			)const;
			

			virtual void localSearch(SolutionType& cur,
				Random *rnd, Problem *pro, Algorithm *alg, int totalEvals, int curType)override;
			//virtual void localSearch(SolutionType& cur,
			//	int totalEvals, int curType,
			//	const std::function<
			//	void(Problem *pro, Algorithm *alg, Random *rnd,
			//		SolutionType& cur)>& fitnessCal,
			//	Problem *pro, Algorithm *alg, Random *rnd
			//) override;


			virtual void localSearch(SolutionType& cur,
				int totalEvals, int curType,
				const std::function<void(SolutionType&)>& fitnessCal,
				Problem *pro, Algorithm *alg, Random *rnd
			) override;



			virtual int localSearch(SolutionType& cur,
				int total_iters,
				Problem *pro,Algorithm *alg, Random *rnd,
				std::shared_ptr<EvaluationStrategyBase<SolutionType>>& eval_stra,
				int pop_idx = 0
			) override;

			virtual int evaluate(Problem *pro, Algorithm *alg, Random *rnd, SolutionType& curSol, bool effective = true)const override;
		protected:
			int m_sample_num = 100;
			EFitnessCalType m_fitness_type = EFitnessCalType::Objective;
			EMemoryType m_memory_type = EMemoryType::Change;
		};

}

#endif // !OPERATOR  _INTERPRETER_SEQ_H
