#include "f1_sr_expanded_two_peak_trap.h"

namespace ofec::cec2015 {
	void F1_SR_expanded_two_peak_trap::initialize_(Environment *env) {
		Function::initialize_(env);
		m_bias = 100;
		//E:\Diao_Yiya\code\OFEC\instance\problem\continuous\single_objective\multi_modal\cec2015\data
		loadOptima("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		loadTranslation("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		loadRotation("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
		// 5  Dim : 1 gopt and 15 lopt
		// 10 Dim : 1 gopt and 55 lopt
		// 20 Dim : 1 gopt and 210 lopt
		evaluateOptima();
	}

	void F1_SR_expanded_two_peak_trap::evaluateObjective(Real *x, std::vector<Real> &obj) {
		translate(x);
		rotate(x);
		obj[0] = 0.0;
		for (size_t i = 0; i < m_number_variables; ++i) {
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
		obj[0] += m_bias;
	}
}