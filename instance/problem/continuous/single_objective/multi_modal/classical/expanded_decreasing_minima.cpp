#include "expanded_decreasing_minima.h"

namespace ofec {
	void ExpandedDecreasingMinima::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void ExpandedDecreasingMinima::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-40, 40);
		setOriginalGlobalOpt();
	}

	void ExpandedDecreasingMinima::updateOptima(Environment *env) {
		m_variable_niche_radius = 0.01;
		m_objective_accuracy = 1.e-4;
		m_optima.reset(new Optima<>(m_original_optima));
	}

	void ExpandedDecreasingMinima::evaluateOriginalObj(Real *x, std::vector<Real> &obj) {
		size_t i;
		obj[0] = 0.0;
		for (i = 0; i < m_number_variables; ++i) {
			x[i] += 0.1;
			if ((x[i] <= 1.0) & (x[i] >= 0.0)) {
				obj[0] += -exp(-2.0 * log(2.0) * pow((x[i] - 0.1) / 0.8, 2.0)) * pow(sin(5.0 * OFEC_PI * x[i]), 6.0);
			}
			else {
				obj[0] += pow(x[i], 2.0);
			}
		}
		obj[0] += 1.0 * m_number_variables;
	}
}