#include "messenger.h"
#include "gtoptoolbox/trajobjfuns.h"

namespace ofec::trajectory {
	void Messenger::initialize_(Environment* env) {
		Continuous::initialize_(env);
		resizeVariable(18);
		m_domain.setRange(1000, 4000, 0);
		m_domain.setRange(1, 5, 1);
		m_domain.setRange(0, 1, 2);
		m_domain.setRange(0, 1, 3);
		m_domain.setRange(200, 400, 4);
		for (size_t i = 5; i < 8; ++i)
			m_domain.setRange(30, 400, i);
		for (size_t i = 8; i < 12; ++i)
			m_domain.setRange(0.01, 0.99, i);
		for (size_t i = 12; i < 15; ++i)
			m_domain.setRange(1.1, 6, i);
		for (size_t i = 15; i < 18; ++i)
			m_domain.setRange(-OFEC_PI, OFEC_PI, i);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
	}

	void Messenger::evaluateObjective(Real *x, std::vector<Real> &obj) {
		std::vector<double> vx(x, x + m_number_variables);
		obj[0] = gtop_toolbox::messenger(vx);
	}
}