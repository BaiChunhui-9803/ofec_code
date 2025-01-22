#include "goldstein_price.h"

namespace ofec {
	void GoldsteinPrice::addInputParameters() {

	}

	void GoldsteinPrice::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(2);
		m_domain.setRange(-2, 2, 0);
		m_domain.setRange(-3, 1, 1);
	}

	void GoldsteinPrice::updateOptima(Environment *env) {
		Continuous::updateOptima(env);
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		temp_sol.variable().vector() = { 0, -1 };
		evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
		optima()->appendSolution(temp_sol);
	}

	void GoldsteinPrice::evaluateObjective(Real *vars, std::vector<Real> &objs) const {
		Real x = vars[0], y = vars[1];
		double part1 = 1 + pow(x + y + 1, 2) * (19 - 14 * x + 3 * pow(x, 2) - 14 * y + 6 * x * y + 3 * pow(y, 2));
		double part2 = 30 + pow(2 * x - 3 * y, 2) * (18 - 32 * x + 12 * pow(x, 2) + 48 * y - 36 * x * y + 27 * pow(y, 2));
		objs[0] = log(part1 * part2);
	}
}
