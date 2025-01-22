/********* Begin Register Information **********
{
	"name": "MKP",
	"identifier": "MultiDimensionalKnapsack",
	"problem tags": [ "MKP", "ComOP", "SOP" ]
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

#ifndef OFEC_MKP_H
#define OFEC_MKP_H

#include "../../../../core/problem/problem.h"
#include "../../../../core/problem/optima.h"
#include "../../../../core/problem/domain.h"
#include "../../../../core/problem/solution.h"

namespace ofec {
#define GET_MKP(pro) dynamic_cast<MultiDimensionalKnapsack&>(pro)

	class MultiDimensionalKnapsack : public Problem {
	public:
		using SolutionType = Solution<VariableVector<int>>;
		using variableType = VariableVector<int>;
	protected:
		size_t m_number_variables;
		std::vector<Real> mv_p;
		std::vector<std::vector<Real>> mvv_r;
		std::vector<Real> mv_b;
		std::string m_file_name;
		int m_m;
		Real m_maxP;
		bool m_if_valid_check = true;
		Real m_opt_obj;

	public:
		SolutionBase* createSolution() override {
			return new SolutionType(numberObjectives(), numberConstraints(), numberVariables());
		}
		SolutionBase* createSolution(const SolutionBase& s) override {
			return new SolutionType(dynamic_cast<const SolutionType&>(s));
		}
		SolutionBase* createSolution(const VariableBase& s) override {
			auto curSol = new SolutionType(numberObjectives(), numberConstraints(), numberVariables());
			curSol->variable() = dynamic_cast<const VariableVector<int>&>(s);
			return curSol;
		}
		bool judgeAccessibility(const SolutionBase& s1, const SolutionBase& s2) override {
			return int(variableDistance(s1, s2)) != m_number_variables;
		}
		virtual void initializeSolution(SolutionBase& s, Random *rnd) const override;
		bool same(const SolutionBase &s1, const SolutionBase &s2) const override;
		Real variableDistance(const SolutionBase &s1, const SolutionBase &s2) const override;
		Real variableDistance(const VariableBase &x1, const VariableBase &x2) const override;	
		virtual Real normalizedVariableDistance(const SolutionBase& s1, const SolutionBase& s2) const override {
			return variableDistance(s1, s2) / m_number_variables;
		}
		virtual Real normalizedVariableDistance(const VariableBase& s1, const VariableBase& s2) const override {
			return variableDistance(s1, s2) / m_number_variables;
		}
		size_t numberVariables() const override { return m_number_variables; }

		bool isValid(const SolutionBase &s) const;
		Real getConstraintValue(const SolutionBase &s, std::vector<Real> &val) const;
		int numInvalidConstraints(SolutionBase &s) const;
		
	protected:
		void initialize_() override;
		void updateOptima() override;
		void evaluate_(SolutionBase &s) override;
		void readProblem();
	};
}

#endif // !OFEC_MKP_H

