#include "f3_sr_expanded_equal_minima.h"

namespace ofec::cec2015 {
	void F3_SR_expanded_equal_minima::initialize_(Environment *env) {
		Function::initialize_(env);
		m_bias = 300;
		m_scale = 20;
		loadOptima("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		loadTranslation("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		loadRotation("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		// 5^Dim gopt 
		evaluateOptima();
	}

	void F3_SR_expanded_equal_minima::evaluateObjective(Real *x, std::vector<Real> &obj) {
		translate(x);
		scale(x);
		rotate(x);
		obj[0] = 0.0;
		for (size_t i = 0; i < m_number_variables; ++i) {
			x[i] += 0.1;
			if ((x[i] <= 1.0) & (x[i] >= 0.0))
				obj[0] += -pow((sin(5 * OFEC_PI * x[i])), 6.0);
			else
				obj[0] += pow(x[i], 2.0);
		}
		obj[0] += 1.0 * m_number_variables;
		obj[0] += m_bias;
	}
}