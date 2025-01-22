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
* Operator adaptor for TSP
*
*********************************************************************************/




#ifndef OPERATOR_INTERPRETER_SP_H
#define OPERATOR_INTERPRETER_SP_H

#include"../../../template/combination/sequence/sequence_interpreter.h"
#include"../../../template/combination/sequence/sequence_solution.h"
#include"../../../../../core/algorithm/Solution.h"

namespace ofec {
	class SP_Solution : public Solution<VariableVector<int>>, public SequenceSolution {
	public:
		template<typename ... Args>
		SP_Solution(size_t num_obj, size_t num_con, Args&& ... args) :
			Solution(num_obj, num_con, std::forward<Args>(args)...), SequenceSolution(){}
		SP_Solution() = default;
		SP_Solution(const SolutionBase& sol) :Solution(sol), SequenceSolution() {}

		SP_Solution& operator=(const SP_Solution& rhs)  {
			Solution<VariableVector<int>>::operator=(dynamic_cast<const Solution<VariableVector<int>>&>(rhs));
			SequenceSolution::operator=(dynamic_cast<const SequenceSolution&>(rhs));
			return *this;
		}
		void reset() {
			Solution::reset();
			SequenceSolution::reset();
		}
	};


	class SP_Interpreter : public SequenceInterpreter<SP_Solution> {
	private:

	public:
		using SolutionType = typename SP_Solution;
	public:
		SP_Interpreter() = default;
		~SP_Interpreter() = default;

		virtual void initializeByProblem(Problem *pro) override;
		virtual Real heursticInfo(Problem *pro, int from, int to, int obj_idx = 0)const override;
		virtual Real heursticInfo(Problem *pro, const SolutionType& pop, int to, int obj_idx = 0)const override;
		virtual SolutionType heuristicSol(Problem *pro, Algorithm *alg, Random *rnd)const override;
		virtual int  curPositionIdx(Problem *pro, const SolutionType& sol)const override;
		virtual void stepInit(Problem *pro, SolutionType& sol)const override;
		virtual void calNextFeasible(Problem *pro, SolutionType& sol)const override;
		virtual bool stepFeasible(Problem *pro, SolutionType& sol) const override;
		virtual void stepBack(Problem *pro, SolutionType& sol)const override;
		virtual bool stepNext(Problem *pro, SolutionType& sol, int next)const override;
		virtual void updateCurEdges(Problem *pro, SolutionType& sol)const override;
		virtual void updateEdges(Problem *pro, SolutionType& sol, bool all = true)const  override;
		virtual bool stepFinish(Problem *pro, const SolutionType& sol)const override;
		virtual void stepFinal(Problem *pro, SolutionType& sol)const override;
		virtual void getNearestNeighbors(Problem *pro, int len, int from_idx, std::vector<int>& neighbors)const override;
		virtual void fillterFeasible(Problem *pro, SolutionType& sol, std::vector<int>& feasible)const override;

		
		// to do
		//virtual void generateSamplesSolutions(Problem *pro, Random *rnd, const SolutionBase& curSol, std::vector<std::unique_ptr<SolutionBase>>& sols)const override;

	protected:
		std::vector<std::vector<int>> m_cur2idx;
		std::vector<std::pair<int, int>> m_idx2cur;
		int m_pos_size = 0;
		int m_dim_size = 0;
		int m_start_state_idx = 0;
		int m_size = 0;
		Real m_eps = 1e-5;

	};



}

#endif // !OPERATOR_INTERPRETER_SEQ_H
