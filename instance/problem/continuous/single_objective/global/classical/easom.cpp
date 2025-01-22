#include "easom.h"

namespace ofec {
	void Easom::addInputParameters() {

	}

	void Easom::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(2);
		setDomain(-1, 7);
	}

	void Easom::updateOptima(Environment *env) {
		Continuous::updateOptima(env);
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		temp_sol.variable().vector() = { OFEC_PI, OFEC_PI };
		evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
		optima()->appendSolution(temp_sol);
	}

	void Easom::evaluateObjective(Real *vars, std::vector<Real> &objs) const {
		Real x = vars[0], y = vars[1];
		objs[0] = -cos(x) * cos(y) * exp(-(pow(x - OFEC_PI, 2) + pow(y - OFEC_PI, 2)));
	}
}
