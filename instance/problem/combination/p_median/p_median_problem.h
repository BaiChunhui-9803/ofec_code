/********* Begin Register Information **********
{
	"name": "PMP",
	"identifier": "P_MedianProblem",
	"problem tags": [ "PMP", "ComOP", "SOP" ]
}
*********** End Register Information **********/

#ifndef P_MEDIA_PROBLEM_H
#define P_MEDIA_PROBLEM_H

#include "../../../../core/problem/problem.h"
#include "../../../../core/problem/optima.h"
#include "../../../../core/problem/domain.h"
#include "../../../../core/problem/solution.h"


namespace ofec {
#define GET_PMEDIA_PRO(pro) dynamic_cast<P_MedianProblem&>(pro)

	class P_MedianProblem : public Problem {
	public:

	public:
		using solutionType = Solution<VariableVector<bool>>;
		using variableType = VariableVector<bool>;
	protected:

		int m_p = 10;

	public:

		virtual SolutionBase* createSolution()override {
			return new solutionType(numberObjectives(), numberConstraints(), numberVariables());
		};

		virtual SolutionBase* createSolution(const SolutionBase& s)override {
			return new solutionType(dynamic_cast<const solutionType&>(s));
		}
		virtual SolutionBase* createSolution(const VariableBase& s)override {
			auto curSol = new solutionType(numberObjectives(), numberConstraints(), numberVariables());
			curSol->variable() = dynamic_cast<const VariableVector<bool>&>(s);
			return curSol;
		}


		virtual bool judgeAccessibility(const SolutionBase& s1, const SolutionBase& s2)override {
			return int(variableDistance(s1, s2)) != m_p;
		}

		virtual void initializeSolution(SolutionBase& s, Random *rnd) const override {}
		bool same(const SolutionBase& s1, const SolutionBase& s2) const override {
			return false;
		}
		Real variableDistance(const SolutionBase& s1, const SolutionBase& s2) const override {

			return 0;

		}
		Real variableDistance(const VariableBase& x1, const VariableBase& x2) const override {
			return 0;
		}
		//inline virtual Real maximumVariableDistance()const override {
		//	return m_num_cities;
		//}
		virtual Real normalizedVariableDistance(const SolutionBase& s1, const SolutionBase& s2) const override {
			return variableDistance(s1, s2) / m_p;
		}
		virtual Real normalizedVariableDistance(const VariableBase& s1, const VariableBase& s2) const override {
			return variableDistance(s1, s2) / m_p;
		}


	protected:

		virtual void updateOptima() override {}
		virtual void initialize_() override {}
		virtual void evaluate_(SolutionBase& s) override {}
	};
}

#endif