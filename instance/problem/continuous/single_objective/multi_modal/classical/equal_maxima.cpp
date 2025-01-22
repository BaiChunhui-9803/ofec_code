#include "equal_maxima.h"

namespace ofec {
	void EqualMaxima::addInputParameters() {
		m_input_parameters.add("objective accuracy", new EnumeratedReal(
			m_objective_accuracy, { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 }, 1e-4));
	}

	void EqualMaxima::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		resizeVariable(1);
		setDomain(0, 1.);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
	}

	void EqualMaxima::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		std::vector<Real> opt_x = { 0.1,0.3,0.5,0.7,0.9 };
		for (Real x : opt_x) {
			temp_sol.variable()[0] = x;
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
		m_variable_niche_radius = 1.e-2;
	}

	void EqualMaxima::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real s;
		s = pow(sin(5 * OFEC_PI*x[0]), 6);
		obj[0] = s;  // note
	}
}