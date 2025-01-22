#include "f5_six_hump_camel_back.h"
namespace ofec::cec2013 {
	void SixHumpCamelBack::initialize_(Environment *env) {
		ofec::SixHumpCamelBack::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
	}

	void SixHumpCamelBack::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		ofec::SixHumpCamelBack::evaluateObjective(x, obj);
		obj[0] = -obj[0];
	}
}