#include "mccormick.h"

namespace ofec {
	void McCormick::addInputParameters() {

	}

	void McCormick::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(2);
		setDomain(-3, 4);
	}

	void McCormick::updateOptima(Environment *env) {
		Continuous::updateOptima(env);
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		temp_sol.variable().vector() = { -0.54719, -1.54719 };
		evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
		optima()->appendSolution(temp_sol);
	}

	void McCormick::evaluateObjective(Real *vars, std::vector<Real> &objs) const {
		Real x = vars[0], y = vars[1];
		objs[0] = sin(x + y) + pow(x - y, 2) - 1.5 * x + 2.5 * y + 1;
	}
}
