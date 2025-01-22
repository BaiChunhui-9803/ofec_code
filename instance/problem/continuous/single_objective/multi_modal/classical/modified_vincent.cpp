#include "modified_vincent.h"

namespace ofec {
	void ModifiedVincent::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void ModifiedVincent::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-10, 10);
	}

	void ModifiedVincent::updateOptima(Environment *env) {
		// 6^Dim gopt 
		m_optima.reset(new Optima<>());
		size_t num = (int)pow(6, m_number_variables);
		for (size_t i = 0; i < num; i++) {
			m_optima->appendObjective(0.0);
		}
		m_variable_niche_radius = 0.01;
		m_objective_accuracy = 1.e-4;
	}

	void ModifiedVincent::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		size_t i;
		obj[0] = 0.0;
		for (i = 0; i < m_number_variables; i++) {
			x[i] += 4.1112;
			if ((x[i] >= 0.25) & (x[i] <= 10.0))
			{
				obj[0] += -sin(10.0 * log(x[i]));
			}
			else if (x[i] < 0.25)
			{
				obj[0] += pow(0.25 - x[i], 2.0) - sin(10.0 * log(2.5));
			}
			else
			{
				obj[0] += pow(x[i] - 10, 2.0) - sin(10.0 * log(10.0));
			}
		}
		obj[0] /= m_number_variables;
		obj[0] += 1.0;
	}

}