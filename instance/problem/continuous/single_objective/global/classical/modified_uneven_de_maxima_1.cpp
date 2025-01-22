#include "modified_uneven_de_maxima_1.h"

namespace ofec {
	void ModifiedUnevenDeMaxima1::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
		resizeVariable(1);
		setDomain(0, 1); // note
		setDomainInitPop(pow(0.0, 4.0 / 3.0), pow(0.2, 4.0 / 3.0));
	}

	void ModifiedUnevenDeMaxima1::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real tmp1 = -2 * log(2) * pow(x[0], pow(2, -32));
		Real tmp2 = sin(5 * OFEC_PI * pow(x[0], 3.0 / 4.0));
		Real s = exp(tmp1) * pow(tmp2, 6);
		obj[0] = s;  // note
	}
}