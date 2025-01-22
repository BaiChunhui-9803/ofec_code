#include "expanded_two_peak_trap.h"

namespace ofec {
	void ExpandedTwoPeakTrap::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void ExpandedTwoPeakTrap::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-40., 40.);
		setOriginalGlobalOpt();
	}

	void ExpandedTwoPeakTrap::updateOptima(Environment *env) {
		m_variable_niche_radius = 0.01;
		m_objective_accuracy = 1.e-4;
		m_optima.reset(new Optima<>(m_original_optima));
	}

	void ExpandedTwoPeakTrap::evaluateOriginalObj(Real *x, std::vector<Real> &obj) {
		size_t i;
		obj[0] = 0.0;
		for (i = 0; i < m_number_variables; ++i) {
			x[i] += 20.0;
			if ((x[i] < 15.0) & (x[i] >= 0.0))
				obj[0] += -(160.0 / 15.0) * (15.0 - x[i]);
			else if ((x[i] <= 20.0) & (x[i] >= 15.0))
				obj[0] += -40.0 * (x[i] - 15.0);
			else if (x[i] < 0.0)
				obj[0] += -160.0 + pow(x[i], 2.0);
			else
				obj[0] += -200.0 + pow(x[i] - 20.0, 2.0);
		}
		obj[0] += 200.0 * m_number_variables;
	}

}