#include "f8_sr_modified_vincent.h"

namespace ofec::cec2015 {
	void F8_SR_modified_vincent::initialize_(Environment *env) {
		Function::initialize_(env);
		setDomain(20, 100);
		m_bias = 800;
		m_scale = 5;
		loadOptima("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		loadTranslation("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		loadRotation("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		// 6^Dim gopt
		evaluateOptima();
	}

	void F8_SR_modified_vincent::evaluateObjective(Real *x, std::vector<Real> &obj) {
		translate(x);
		scale(x);
		rotate(x);
		obj[0] = 0.0;
		for (size_t i = 0; i < m_number_variables; i++) {
			x[i] += 4.1112;//shift to orgin
			if ((x[i] >= 0.25) & (x[i] <= 10.0))
				obj[0] += -sin(10.0 * log(x[i]));
			else if (x[i] < 0.25)
				obj[0] += pow(0.25 - x[i], 2.0) - sin(10.0 * log(2.5));
			else
				obj[0] += pow(x[i] - 10, 2.0) - sin(10.0 * log(10.0));
		}
		obj[0] /= m_number_variables;
		obj[0] += 1 + m_bias;
	}
}