#include "f8_modified_rastrigin.h"

namespace ofec::cec2013 {
	void ModifiedRastrigin::initialize_(Environment *env) {
		ofec::ModifiedRastrigin::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
	}

	void ModifiedRastrigin::evaluateObjective(Real* x, std::vector<Real> &obj) const {
		ofec::ModifiedRastrigin::evaluateObjective(x, obj);
		Real s = 0;
		for (int i = 0; i < m_number_variables; ++i) {
			s += 10 + 9 * cos(2 * OFEC_PI * m_k[i] * x[i]);
		}
		obj[0] = -obj[0];
	}
}