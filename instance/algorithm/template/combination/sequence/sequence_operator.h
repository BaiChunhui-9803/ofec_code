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




#ifndef OPERATOR_SEQUENCE_H
#define OPERATOR_SEQUENCE_H

#include<functional>
#include<vector>
#include<memory>
#include "../../../../../core/problem/encoding.h"
#include "../../../../../core/definition.h"
#include "../../../../../utility/functional.h"
#include"../../framework/uncertianty/evaluation_strategy.h"



namespace ofec {


		class SequenceOperatorBase {
		public:
			SequenceOperatorBase() = default;
			~SequenceOperatorBase() = default;
			virtual void initialize(Problem *pro, Algorithm *alg) {}
		};

		template<typename TInterpreter>
		class SequenceOperator : public SequenceOperatorBase {
		public:
			using SolutionType = typename TInterpreter::SolutionType;
			using InterpreterType = typename TInterpreter;

		public:
			SequenceOperator() = default;
			~SequenceOperator() = default;
			virtual int learn_from_local(
				const SolutionType& cur, Random *rnd, Problem *pro,
				const InterpreterType& interpreter
			)const = 0;
			virtual int learn_from_global(
				const SolutionType& cur, Random *rnd, Problem *pro,
				const InterpreterType& interpreter,
				const std::function<Real(const SolutionType& cur, int to)>& proFun
			)const = 0;
			virtual int learn_from_other(
				const SolutionType& cur, Random *rnd, Problem *pro,
				const InterpreterType& interpreter,
				const SolutionType& other
			)const = 0;

			virtual bool learn_from_other(
				SolutionType& cur, Random *rnd, Problem *pro,
				const InterpreterType& interpreter,
				const SolutionType& other,
				Real radius
			)const = 0;

			virtual void localSearch(SolutionType& cur,
				Random *rnd, Problem *pro, Algorithm *alg, int totalEvals, int curType) {}

			virtual void localSearch(SolutionType& cur,
				int totalEvals, int curType,
				const std::function<void (SolutionType&)>& fitnessCal,
				 Problem *pro, Algorithm *alg, Random *rnd
				) {}

			virtual int localSearch(SolutionType& cur,
				int total_iters,
				Problem *pro, Algorithm *alg, Random *rnd,
				std::shared_ptr<EvaluationStrategyBase<SolutionType>>& eval_stra,
				int pop_idx = 0
			) {
				return eval_stra->calFitness(cur, pop_idx);
			}

			virtual int evaluate(Problem *pro, Algorithm *alg, Random *rnd, SolutionType& curSol, bool effective = true)const = 0;

			virtual bool better(Problem *pro,const SolutionType& a, const SolutionType& b) {
				return a.dominate(b, pro);
			}
			
		//	virtual int getIdx(Problem *pro,std::vector<std::unique_ptr<SolutionType>>& curIdx,bool betterFlag=true);
		//	virtual int getIdx(Problem *pro, std::vector<SolutionType>& curIdx, bool betterFlag = true);

		};

		//template<typename TInterpreter>
		//int SequenceOperator<TInterpreter>::getIdx(Problem *pro, std::vector<std::unique_ptr<SolutionType>>& cur,bool betterFlag){
		//	if (betterFlag) {
		//		return calIdx(cur, [&](const std::unique_ptr<SolutionType>& a, const  std::unique_ptr<SolutionType>& b) {
		//			return better(pro, *a, *b);
		//		});
		//	}
		//	else {
		//		return calIdx(cur, [&](const std::unique_ptr<SolutionType>& a, const  std::unique_ptr<SolutionType>& b) {
		//			return !better(pro, *a, *b);
		//		});
		//	}
		//}

		//template<typename TInterpreter>
		//int SequenceOperator<TInterpreter>::getIdx(Problem *pro, std::vector<SolutionType>& cur, bool betterFlag) {
		//	if (betterFlag) {
		//		return calIdx<SolutionType>(cur, [&](const SolutionType& a, const  SolutionType& b) {
		//			return better(pro, a, b);
		//		});
		//	}
		//	else {
		//		return calIdx(cur, [&](const SolutionType& a, const SolutionType& b) {
		//			return !better(pro, a, b);
		//		});
		//	}
		//}




	//template<typename T>
	//int getWorstIdx(const std::vector<T>& indis) {
	//	if (T.size() == 0) return -1;
	//	int worstIdx(0);
	//	for (int curIdx(0); curIdx < T.size(); ++curIdx) {
	//		if (T[worstIdx].dominate(T[curIdx], m_problem.get())) {
	//			worstIdx = curIdx;
	//		}
	//	}
	//	return worstIdx;
	//}


	//template<typename T>
	//int getWorstIdx(Problem *pro,const std::vector<T>& indis) {
	//	if (T.size() == 0) return -1;
	//	int worstIdx(0);
	//	for (int curIdx(0); curIdx < T.size(); ++curIdx) {
	//		if (T[worstIdx].dominate(T[curIdx], m_problem.get())) {
	//			worstIdx = curIdx;
	//		}
	//	}
	//	return worstIdx;
	//}


	//template<typename T>
	//int getWorstPtrIdx(Problem *pro,const std::vector<std::unique_ptr<T>>& indis) {
	//	if (T.size() == 0) return -1;
	//	int worstIdx(0);
	//	for (int curIdx(0); curIdx < T.size(); ++curIdx) {
	//		if (T[worstIdx]->dominate(*T[curIdx], m_problem.get())) {
	//			worstIdx = curIdx;
	//		}
	//	}
	//	return worstIdx;
	//}

}

#endif // !OPERATOR_INTERPRETER_SEQ_H

