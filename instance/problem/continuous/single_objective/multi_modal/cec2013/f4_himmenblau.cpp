#include "f4_himmenblau.h"

namespace ofec::cec2013 {
	void Himmenblau::initialize_(Environment *env) {
		ofec::Himmenblau::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
	}

	void Himmenblau::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		ofec::Himmenblau::evaluateObjective(x, obj);
		obj[0] = 200 - obj[0];
	}
}