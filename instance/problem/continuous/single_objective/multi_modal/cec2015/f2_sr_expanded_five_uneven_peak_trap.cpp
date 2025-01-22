#include "f2_sr_expanded_five_uneven_peak_trap.h"

namespace ofec::cec2015 {
	void F2_SR_expanded_five_uneven_peak_trap::initialize_(Environment *env) {
		Function::initialize_(env);
		m_bias = 200;
		loadOptima("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		loadTranslation("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		loadRotation("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		// 2^Dim gopt and 5^Dim - 2^Dim lopt 
		evaluateOptima();
	}

	void F2_SR_expanded_five_uneven_peak_trap::evaluateObjective(Real *x, std::vector<Real> &obj) {
		translate(x);
		rotate(x);
		obj[0] = 0.0;
		for (size_t i = 0; i < m_number_variables; ++i) {
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
		obj[0] += m_bias;
	}
}