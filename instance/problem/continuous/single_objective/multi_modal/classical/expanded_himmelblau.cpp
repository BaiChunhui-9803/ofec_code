#include "expanded_himmelblau.h"

namespace ofec {
	void ExpandedHimmelblau::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 2, 8, 2));
	}

	void ExpandedHimmelblau::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		if (m_number_variables % 2) {
			throw Exception("Number of variables of the expanded Himmelblau function must be an even.");
		}
		resizeVariable(m_number_variables);
		setDomain(-10., 4.);
	}

	void ExpandedHimmelblau::updateOptima(Environment *env) {
		// 4^(Dim/2) gopt 
		m_optima.reset(new Optima<>());
		size_t num_opts = (int)pow(4, m_number_variables / 2);
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		std::vector<int> quaternary_code(m_number_variables);
		std::vector<std::vector<Real>> quaternary_value = {
			{ 0.0,0.0 }, 
			{0.584428,-3.848126}, 
			{-6.779310,-5.283186}, 
			{-5.805118,1.131312} 
		};
		for (size_t i = 0; i < num_opts; i++) {
			int num = i;
			size_t j = 0;
			while (num != 0) {
				quaternary_code[j++] = num % 4;
				num /= 4;
			}
			for (j = 0; j < m_number_variables / 2; j++) {
				temp_sol.variable()[j * 2] = quaternary_value[quaternary_code[j]][0];
				temp_sol.variable()[j * 2 + 1] = quaternary_value[quaternary_code[j]][1];
			}
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
		m_variable_niche_radius = 0.01;
		m_objective_accuracy = 1.e-4;
	}

	void ExpandedHimmelblau::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		size_t i;
		obj[0] = 0.0;
		for (i = 0; i < m_number_variables - 1; i += 2) {
			x[i] += 3.0;
			x[i + 1] += 2.0;
			obj[0] += pow((x[i] * x[i] + x[i + 1] - 11.0), 2.0) + pow((x[i] + x[i + 1] * x[i + 1] - 7.0), 2.0);
		}
	}

}