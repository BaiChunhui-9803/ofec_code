#include "expanded_equal_maxima.h"

namespace ofec {
	void ExpandedEqualMaxima::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void ExpandedEqualMaxima::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-2, 2);
	}

	void ExpandedEqualMaxima::updateOptima(Environment *env) {
		m_variable_niche_radius = 0.01;
		m_objective_accuracy = 1.e-4;
		// 5^Dim gopt 
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		std::vector<int> quinary_code(m_number_variables);
		std::vector<Real> quinary_value = { 0.0,0.2,0.4,0.6,0.8 };
		size_t num_opts = (int)pow(5, m_number_variables);
		for (size_t i = 0; i < num_opts; i++) {
			int num = i;
			size_t j = 0;
			while (num != 0) {
				quinary_code[j++] = num % 5;
				num /= 5;
			}
			for (j = 0; j < m_number_variables; j++) {
				temp_sol.variable()[j] = quinary_value[quinary_code[j]];
			}
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
	}

	void ExpandedEqualMaxima::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		size_t i;
		obj[0] = 0.0;
		for (i = 0; i < m_number_variables; ++i) {
			x[i] += 0.1;
			if ((x[i] <= 1.0) & (x[i] >= 0.0))
				obj[0] += -pow((sin(5 * OFEC_PI * x[i])), 6.0);
			else
				obj[0] += pow(x[i], 2.0);
		}
		obj[0] += 1.0 * m_number_variables;
	}

}