/********* Begin Register Information **********
{
	"name": "QAP",
	"identifier": "QuadraticAssignment",
	"problem tags": [ "QAP", "ComOP", "SOP" ]
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

#ifndef OFEC_QAP_H
#define OFEC_QAP_H

#include "../../../../core/problem/problem.h"
#include "../../../../core/problem/optima.h"
#include "../../../../core/problem/domain.h"
#include "../../../../core/problem/solution.h"

namespace ofec {
#define GET_QAP(pro) dynamic_cast<QuadraticAssignment&>(pro)

	class QuadraticAssignment : public Problem {
	public:
		using solutionType = Solution<VariableVector<int>>;
		using variableType = VariableVector<int>;
	protected:
		size_t m_number_variables;
		std::vector<std::vector<Real> > mvv_flow;
		std::vector<std::vector<Real> > mvv_distance;
		std::string m_file_name;
		Domain<int> m_domain;
		bool m_if_valid_check = true;
	public:
		SolutionBase* createSolution() override {
			return new solutionType(numberObjectives(), numberConstraints(), numberVariables());
		}
		SolutionBase* createSolution(const SolutionBase &s) override {
			return new solutionType(dynamic_cast<const solutionType &>(s));
		}
		SolutionBase* createSolution(const VariableBase &s) override {
			auto curSol = new solutionType(numberObjectives(), numberConstraints(), numberVariables());
			curSol->variable() = dynamic_cast<const VariableVector<int>&>(s);
			return curSol;
		}
		bool judgeAccessibility(const SolutionBase &s1, const SolutionBase &s2) override {
			return int(variableDistance(s1, s2)) != m_number_variables;
		}
		void initializeSolution(SolutionBase &s, Random *rnd) const override;
		bool same(const SolutionBase &s1, const SolutionBase &s2) const override;
		Real variableDistance(const SolutionBase &s1, const SolutionBase &s2) const override;
		Real variableDistance(const VariableBase &x1, const VariableBase &x2) const override;
		Real normalizedVariableDistance(const SolutionBase &s1, const SolutionBase &s2) const override {
			return variableDistance(s1, s2) / m_number_variables;
		}
		Real normalizedVariableDistance(const VariableBase &s1, const VariableBase &s2) const override {
			return variableDistance(s1, s2) / m_number_variables;
		}
		void updateCandidates(const SolutionBase &sol, std::list<std::unique_ptr<SolutionBase>> &candidates) const override;
		size_t numOptimaFound(const std::list<std::unique_ptr<SolutionBase>> &candidates) const override;
		int updateEvaluationTag(SolutionBase &s, Algorithm *alg) override;
		size_t numberVariables() const override { return m_number_variables; }

		bool isValid(const SolutionBase &s) const;

	protected:
		void initialize_() override;
		void updateOptima() override;
		virtual void readProblem();    //read source data from file
		virtual void readOptima();	   //read optima data from file
		void evaluate_(SolutionBase &s) override;
	};
}

#endif // !OFEC_QAP_H

