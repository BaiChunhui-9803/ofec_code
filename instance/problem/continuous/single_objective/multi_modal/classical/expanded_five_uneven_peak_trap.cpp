#include "expanded_five_uneven_peak_trap.h"

namespace ofec {
	void ExpandedFiveUnevenPeakTrap::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void ExpandedFiveUnevenPeakTrap::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-40, 40);
	}

	void ExpandedFiveUnevenPeakTrap::updateOptima(Environment *env) {
		m_variable_niche_radius = 0.01;
		m_objective_accuracy = 1.e-4;
		// 2^Dim gopt 
		m_optima.reset(new Optima<>());
		size_t num_opts = (int)pow(2, m_number_variables);
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		std::vector<int> binary_code(m_number_variables);
		std::vector<Real> binary_value = { 0.0,30.0 };
		for (size_t i = 0; i < num_opts; i++) {
			int num = i;
			size_t j = 0;
			while (num != 0) {
				binary_code[j++] = num % 2;
				num /= 2;
			}
			for (j = 0; j < m_number_variables; j++) {
				temp_sol.variable()[j] = binary_value[binary_code[j]];
			}
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
	}

	void ExpandedFiveUnevenPeakTrap::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		size_t i;
		obj[0] = 0.0;
		for (i = 0; i < m_number_variables; ++i) {
			if (x[i] < 0)
				obj[0] += -200.0 + pow(x[i], 2.0);
			else if (x[i] < 2.5)
				obj[0] += -80.0 * (2.5 - x[i]);
			else if (x[i] < 5.0)
				obj[0] += -64.0 * (x[i] - 2.5);
			else if (x[i] < 7.5)
				obj[0] += -160.0 + pow(x[i], 2.0);
			else if (x[i] < 12.5)
				obj[0] += -28.0 * (x[i] - 7.5);
			else if (x[i] < 17.5)
				obj[0] += -28.0 * (17.5 - x[i]);
			else if (x[i] < 22.5)
				obj[0] += -32.0 * (x[i] - 17.5);
			else if (x[i] < 27.5)
				obj[0] += -32.0 * (27.5 - x[i]);
			else if (x[i] <= 30.0)
				obj[0] += -80.0 * (x[i] - 27.5);
			else
				obj[0] += -200.0 + pow(x[i] - 30.0, 2.0);
		}
		obj[0] += 200.0 * m_number_variables;
	}

}