#include "different_powers.h"

namespace ofec {
	void DifferentPowers::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void DifferentPowers::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-100., 100.);
		setOriginalGlobalOpt();
	}

	void DifferentPowers::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1e-8;
	}

	void DifferentPowers::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		obj[0] = 0.0;
		for (size_t i = 0; i<m_number_variables; i++) {
			obj[0] += pow(fabs(x[i]), 2.0 + 4.0*i / (m_number_variables - 1.0));
		}
		obj[0] = pow(obj[0], 0.5);
	}

}