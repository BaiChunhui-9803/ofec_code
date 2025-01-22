#include "uneven_de_maxima.h"

namespace ofec {
	void UnevenDeMaxima::addInputParameters() {
		m_input_parameters.add("objective accuracy", new EnumeratedReal(
			m_objective_accuracy, { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 }, 1e-4));
	}

	void UnevenDeMaxima::initialize_(Environment *env) { 
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
		resizeVariable(1);
		setDomain(0, 1);
	}

	void UnevenDeMaxima::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		std::vector<Real> opt_x = { 0.079699779582100 };
		for (Real x : opt_x) {
			temp_sol.variable()[0] = x;
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
		m_variable_niche_radius = 0.01;
	}

	void UnevenDeMaxima::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real tmp1 = -2 * log(2) * ((x[0] - 0.08) / 0.854) * ((x[0] - 0.08) / 0.854);
		Real tmp2 = sin(5 * OFEC_PI * (pow(x[0], 3.0 / 4.0) - 0.05));
		Real s = exp(tmp1) * pow(tmp2, 6);
		obj[0] = s;  // note
	}
}