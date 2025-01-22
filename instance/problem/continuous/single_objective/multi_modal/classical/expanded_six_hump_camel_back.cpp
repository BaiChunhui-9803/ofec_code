#include "expanded_six_hump_camel_back.h"

namespace ofec {
	void ExpandedSixHumpCamelBack::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void ExpandedSixHumpCamelBack::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		if (m_number_variables % 2) {
			throw Exception("Number of variables of the expanded six-hump camel back function must be an even.");
		}
		resizeVariable(m_number_variables);
		for (size_t j = 0; j < m_number_variables / 2; j++) {
			m_domain.setRange(-2, 2, j * 2);
			m_domain.setRange(-0.5, 2, j * 2 + 1);
		}
	}

	void ExpandedSixHumpCamelBack::updateOptima(Environment *env) {
		// 2^(Dim/2) gopt 
		m_optima.reset(new Optima<>());
		size_t num_opts = (int)pow(2, m_number_variables / 2);
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		std::vector<int> binary_code(m_number_variables);
		std::vector<std::vector<Real>> binary_value = { { 0.0,0.0 }, {-0.179684,1.425312} };
		for (size_t i = 0; i < num_opts; i++) {
			int num = i;
			size_t j = 0;
			while (num != 0) {
				binary_code[j++] = num % 2;
				num /= 2;
			}
			for (j = 0; j < m_number_variables / 2; j++) {
				temp_sol.variable()[j * 2] = binary_value[binary_code[j]][0];
				temp_sol.variable()[j * 2 + 1] = binary_value[binary_code[j]][1];
			}
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
		m_variable_niche_radius = 0.01;
		m_objective_accuracy = 1.e-4;
	}

	void ExpandedSixHumpCamelBack::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		size_t i;
		obj[0] = 0.0;
		for (i = 0; i<m_number_variables - 1; i += 2) {
			x[i] += 0.089842;
			x[i + 1] -= 0.712656;
			obj[0] += ((4.0 - 2.1*pow(x[i], 2.0) + pow(x[i], 4.0) / 3.0)*pow(x[i], 2.0) + x[i] * x[i + 1] + ((-4.0 + 4.0*pow(x[i + 1], 2.0))*pow(x[i + 1], 2.0)))*4.0;
		}
		obj[0] += 4.126514*m_number_variables / 2.0;
	}

}