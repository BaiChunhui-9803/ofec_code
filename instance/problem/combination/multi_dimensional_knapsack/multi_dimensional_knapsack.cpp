#include "multi_dimensional_knapsack.h"
#include "../../../../core/global.h"
#include "../../../../core/problem/solution.h"
#include <fstream>
#include <numeric>

namespace ofec {
	void MultiDimensionalKnapsack::initialize_() {
		Problem::initialize_();
		auto& v = *m_param;;
		m_file_name = v.get<std::string>("dataFile1");
		readProblem();
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
	}

	void MultiDimensionalKnapsack::updateOptima() {
		m_optima.reset(new Optima<VariableVector<int>>());
		if (m_opt_obj != 0) {
			m_optima->appendObj(m_opt_obj);
			m_optima->setObjectiveGiven(true);
		}
	}

	void MultiDimensionalKnapsack::evaluate_(SolutionBase &s) {
		VariableVector<int> &x = dynamic_cast<Solution<VariableVector<int>>&>(s).variable();
		std::vector<Real> &obj = dynamic_cast<Solution<VariableVector<int>>&>(s).objective();
		int m = numInvalidConstraints(s);
		for (size_t i = 0; i < m_number_objectives; ++i)
			obj[i] = 0;
		for (size_t n = 0; n < m_number_objectives; ++n) {
			for (size_t i = 0; i < m_number_variables; ++i)
				obj[n] += mv_p[i] * x[i];
			obj[n] -= m * m_maxP;
		}
	}

	bool MultiDimensionalKnapsack::isValid(const SolutionBase &s) const {
		if (!m_if_valid_check) return true;
		const VariableVector<int> &x = dynamic_cast<const Solution<VariableVector<int>> &>(s).variable();
		for (int i = 0; i < m_number_variables; i++) {
			if (x[i] != 0 && x[i] != 1)
				return false;
		}
		return true;
	}

	void MultiDimensionalKnapsack::readProblem() {
		size_t i;
		std::string Line;
		std::ostringstream oss;
		std::ifstream infile;
		oss << static_cast<std::string>(g_working_dir);
		oss << "/instance/problem/combination/multi_dimensional_knapsack/" << m_file_name << ".mkp";
		infile.open(oss.str().c_str());
		if (!infile) {
			throw MyExcept("read Multidimensional Knapsack data error");
		}
		infile >> Line;
		m_number_variables = atoi(Line.c_str());
		infile >> Line;
		m_m = atoi(Line.c_str());
		infile >> Line;
		m_opt_obj = atof(Line.c_str());
		mv_p.resize(m_number_variables);
		mv_b.resize(m_m);
		mvv_r.resize(m_m);
		for (i = 0; i < m_m; i++) {
			mvv_r[i].resize(m_number_variables);
		}
		for (i = 0; i < m_number_variables; i++)
		{
			infile >> mv_p[i];
			if (i == 0) {
				m_maxP = mv_p[i];
			}
			else if (m_maxP < mv_p[i]) {
				m_maxP = mv_p[i];
			}
		}
		for (i = 0; i < m_m; i++) {
			for (int j = 0; j < m_number_variables; j++) {
				infile >> mvv_r[i][j];
			}
		}
		for (i = 0; i < m_m; i++) {
			infile >> mv_b[i];
		}
		infile.close();
		infile.clear();
	}

	int MultiDimensionalKnapsack::numInvalidConstraints(SolutionBase &s) const {
		const VariableVector<int> &x = dynamic_cast<const Solution<VariableVector<int>> &>(s).variable();
		int n = 0;
		Real sum;
		for (int i = 0; i < m_m; i++) {
			sum = 0;
			for (int j = 0; j < m_number_variables; j++) {
				sum += mvv_r[i][j] * x[j];
			}
			if (sum > mv_b[i]) n++;
		}
		return n;
	}

	Real MultiDimensionalKnapsack::getConstraintValue(const SolutionBase &s, std::vector<Real> &val) const {
		const VariableVector<int> &x = dynamic_cast<const Solution<VariableVector<int>> &>(s).variable();
		val.resize(m_m);
		Real sum = 0;
		for (int i = 0; i < m_m; ++i) {
			val[i] = 0;
			for (int j = 0; j < m_number_variables; ++j)
				val[i] += mvv_r[i][j] * x[j];
			val[i] = sum - mv_b[i] < 0 ? 0 : sum - mv_b[i];
		}
		return std::accumulate(val.begin(), val.end(), 0.0);
	}
	void MultiDimensionalKnapsack::initializeSolution(SolutionBase& s, Random *rnd) const {
		VariableVector<int>& x = dynamic_cast<Solution<VariableVector<int>> &>(s).variable();
		for (int i = 0; i < m_number_variables; ++i)
			x[i] = rnd->uniform.nextNonStd(0, 2);
		if (!isValid(s))
			throw MyExcept("error in @multi_dimensional_knapsack::initialize_solution() in multi_dimensional_knapsack.cpp");
	}




	bool MultiDimensionalKnapsack::same(const SolutionBase &s1, const SolutionBase &s2) const {
		const VariableVector<int> &x1 = dynamic_cast<const Solution<VariableVector<int>> &>(s1).variable();
		const VariableVector<int> &x2 = dynamic_cast<const Solution<VariableVector<int>> &>(s2).variable();
		for (int i = 0; i < m_number_variables; i++)
			if (x1[i] != x2[i])
				return false;
		return true;
	}

	Real MultiDimensionalKnapsack::variableDistance(const SolutionBase &s1, const SolutionBase &s2) const {
		const VariableVector<int> &x1 = dynamic_cast<const Solution<VariableVector<int>> &>(s1).variable();
		const VariableVector<int> &x2 = dynamic_cast<const Solution<VariableVector<int>> &>(s2).variable();
		Real dis = 0;
		for (int i = 0; i < m_number_variables; i++)
			if (x1[i] != x2[i])
				dis++;
		return dis;
	}

	Real MultiDimensionalKnapsack::variableDistance(const VariableBase &x1, const VariableBase &x2) const {
		const VariableVector<int> &x1_ = dynamic_cast<const VariableVector<int>&>(x1);
		const VariableVector<int> &x2_ = dynamic_cast<const VariableVector<int>&>(x2);
		Real dis = 0;
		for (int i = 0; i < m_number_variables; i++)
			if (x1_[i] != x2_[i])
				dis++;
		return dis;
	}
}