#include "rosetta.h"
#include "gtoptoolbox/trajobjfuns.h"

namespace ofec::trajectory {
	void Rosetta::initialize_() {
		Continuous::initialize_();
		resizeVariable(22);
		m_domain.setRange(1460, 1825, 0);
		m_domain.setRange(3, 5, 1);
		m_domain.setRange(0, 1, 2);
		m_domain.setRange(0, 1, 3);
		m_domain.setRange(300, 500, 4);
		m_domain.setRange(150, 800, 5);
		m_domain.setRange(150, 800, 6);
		m_domain.setRange(300, 800, 7);
		m_domain.setRange(700, 1850, 8);
		for (size_t i = 9; i < 14; ++i)
			m_domain.setRange(0.01, 0.9, i);
		for (size_t i = 14; i < 18; ++i)
			m_domain.setRange(1.05, 9, i);
		for (size_t i = 18; i < 22; ++i)
			m_domain.setRange(-OFEC_PI, OFEC_PI, i);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
	}

	void Rosetta::evaluateObjective(Real *x, std::vector<Real> &obj) {
		std::vector<double> vx(x, x + m_number_variables);
		obj[0] = gtop_toolbox::rosetta(vx);
	}
}