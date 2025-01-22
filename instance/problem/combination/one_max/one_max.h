/********* Begin Register Information **********
{
	"name": "OneMax",
	"identifier": "OneMax",
	"problem tags": [ "OneMax", "ComOP", "SOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com  Or cugxiayong@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 18 Apr 2016
// Modified: 14 Mar 2018 By Junchen Wang (wangjunchen@cug.edu.cn)
// Modified: 31 May 2024 By Yiya Diao

#ifndef ONE_MAX_H
#define ONE_MAX_H

#include "../../../../core/problem/problem.h"
#include "../../../../core/problem/optima.h"
#include "../../../../core/problem/solution.h"

namespace ofec {
#define CAST_ONEMAX(pro) dynamic_cast<OneMax*>(pro)

	class OneMax : public ProblemVariableVector<int>{
		OFEC_CONCRETE_INSTANCE(OneMax)
	
	protected:
		bool m_if_valid_check = true;

		
	public:
		virtual void evaluate(const VariableBase& vars, std::vector<Real>& objs, std::vector<Real>& cons)override;
	

		bool judgeAccessibility(const SolutionBase& s1, const SolutionBase& s2) {
			return int(variableDistance(s1.variableBase(), s2.variableBase())) != m_number_variables;
		}

		virtual void initializeVariables(VariableBase& x, Random* rnd) const override;
		virtual bool same(const VariableBase& x1, const VariableBase& x2) const override;
		virtual Real variableDistance(const VariableBase& x1, const VariableBase& x2) const override;

		bool same(const SolutionBase &s1, const SolutionBase &s2) const ;

		

		bool isValid(const SolutionBase &s);
		void initSolutionNearBy(const SolutionBase& centerSol, SolutionBase& sol, int neighbork, ofec::Random* rnd);
		void initSolVecNearBy(const SolutionBase& centerSol, SolutionBase& sol, std::vector<bool>& vec);
		void filterSameSols(std::vector<const SolutionBase*>& sols, std::vector<bool>& uniqueFlag);

	protected:
		void addInputParameters();
		virtual void initialize_(Environment* env)override;
		virtual void updateOptima(Environment* env)override;
	};

}

#endif // !ONE_MAX_H

