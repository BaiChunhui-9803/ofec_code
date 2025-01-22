#include "f5_sr_expanded_uneven_minima.h"

namespace ofec::cec2015 {
	void F5_SR_expanded_uneven_minima::initialize_(Environment *env) {
		Function::initialize_(env);
		m_bias = 500;
		m_scale = 20;
		loadOptima("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		loadTranslation("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		loadRotation("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		// 5^Dim gopt		 
		evaluateOptima();
	}

	void F5_SR_expanded_uneven_minima::evaluateObjective(Real *x, std::vector<Real> &obj) {
		translate(x);
		scale(x);
		rotate(x);
		obj[0] = 0.0;
		for (size_t i = 0; i < m_number_variables; i++) {
			x[i] += 0.079699392688696;//pow(0.15,4.0/3.0);//shift to orgin
			if ((x[i] <= 1.0) & (x[i] >= 0.0))
				obj[0] -= pow(sin(5.0 * OFEC_PI * (pow(x[i], 0.75) - 0.05)), 6.0);
			else
				obj[0] += pow(x[i], 2.0);
		}
		obj[0] += 1.0 * m_number_variables;
		obj[0] += m_bias;
	}
}