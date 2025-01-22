#include "michalewicz.h"

namespace ofec {
	void Michalewicz::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
		resizeVariable(2);
		setDomain(0, OFEC_PI);
		m_m = 20;
	}

	void Michalewicz::updateOptima(Environment *env) {
		m_objective_accuracy = 1.e-3;
		 m_variable_niche_radius = 0.2;
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		//there is only one global optima 
		std::vector<std::vector<Real>> var_data = { { 2.20291f, 1.5708f } };
		for (auto &i : var_data) {
			temp_sol.variable()[0] = i[0];
			temp_sol.variable()[1] = i[1];
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
	}

	void Michalewicz::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real s = 0;
		for (int i = 0; i < m_number_variables; ++i) {
			s += sin(x[i])*pow(sin((i + 1)*x[i] * x[i] / OFEC_PI), m_m);
		}
		obj[0] = s;
	}
	
}