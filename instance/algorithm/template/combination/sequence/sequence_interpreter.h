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




#ifndef OPERATOR_INTERPRETER_SEQ_H
#define OPERATOR_INTERPRETER_SEQ_H

#include"../../../../../core/definition.h"
#include"../../../../../core/problem/encoding.h"
#include<vector>
#include<iostream>

// problem interpreter and operator for sequence problem
namespace ofec {


		class SequenceInterpreterBase {

		protected:

			virtual void initializeByProblem(Problem *pro) = 0;
		public:
			SequenceInterpreterBase() = default;
			~SequenceInterpreterBase() = default;

			// for online evaluation
			virtual void remapDimemsion(Problem *pro, Algorithm *alg) {}

			virtual void initialize(Problem *pro, Algorithm *alg) {
				initializeByProblem(pro);
			}
			virtual void updateByProblem(Problem *pro) {}
			const std::vector<int>& getMatrixSize() const {
				return m_matrix_size;
			}

			// for analyse
			virtual bool optimalEdgeGiven(Problem *pro)const {
				return false;
			}
			virtual void calOptEdges(Problem *pro)const {};
		protected:
			std::vector<int> m_matrix_size;
		};


		template<typename TSolution>
		class SequenceInterpreter : public SequenceInterpreterBase {
		public:
			using SolutionType = TSolution;
		protected:
			
		public:
			SequenceInterpreter() = default;
			~SequenceInterpreter() = default;



			virtual Real heursticInfo(Problem *pro,int from, int to, int obj_idx = 0)const = 0;
			virtual Real heursticInfo(Problem *pro,const SolutionType& pop, int to, int obj_idx = 0)const = 0;

			virtual SolutionType heuristicSol(Problem *pro,Algorithm *alg, Random *rnd)const = 0;
			virtual void getNearestNeighbors(Problem *pro, int len, int from_idx, std::vector<int>& neighbors)const = 0;
			virtual int curPositionIdx(Problem *pro,const SolutionType& sol)const = 0;
			virtual void stepInit(Problem *pro, SolutionType& sol)const {
				sol.evaluateTimes() = 0;
			}
			virtual void calNextFeasible(Problem *pro,SolutionType& sol)const = 0;
			virtual bool stepFeasible(Problem *pro,SolutionType& sol) const = 0;
			virtual void fillterFeasible(Problem *pro, SolutionType& sol, std::vector<int>& feasible)const = 0;
			virtual void stepBack(Problem *pro,SolutionType& sol)const = 0;
			virtual bool stepNext(Problem *pro,SolutionType& sol, int next)const = 0;
			virtual void updateCurEdges(Problem *pro,SolutionType& sol)const = 0;
			virtual void updateEdges(Problem *pro, SolutionType& sol, bool all=true)const = 0;
			virtual bool stepFinish(Problem *pro,const SolutionType& sol) const = 0;
			virtual void stepFinal(Problem *pro,SolutionType& sol)const = 0;
		};

}

#endif // !OPERATOR_INTERPRETER_SEQ_H
