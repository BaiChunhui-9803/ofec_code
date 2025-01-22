#include "beale.h"

namespace ofec {
	void Beale::addInputParameters() {

	}

	void Beale::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(2);
		setDomain(-5, 5);
	}

	void Beale::updateOptima(Environment *env) {
		Continuous::updateOptima(env);
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		temp_sol.variable().vector() = { 3,0.5 };
		evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
		optima()->appendSolution(temp_sol);
	}

	void Beale::evaluateObjective(Real *vars, std::vector<Real> &objs) const {
		Real x = vars[0], y = vars[1];
		objs[0] = log(pow(1.5 - x + x * y, 2) + pow(2.25 - x + x * y * y, 2) + pow(2.625 - x + x * y * y * y, 2) + 1);
	}
}
