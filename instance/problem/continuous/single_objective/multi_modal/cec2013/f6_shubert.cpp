#include  "f6_shubert.h"

namespace ofec::cec2013 {
	void Shubert::initialize_(Environment *env) {
		ofec::Shubert::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
	}

	void Shubert::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		ofec::Shubert::evaluateObjective(x, obj);
		obj[0] = -obj[0];
	}
}