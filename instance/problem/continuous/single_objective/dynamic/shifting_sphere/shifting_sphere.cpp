#include "shifting_sphere.h"

namespace ofec {
	void ShiftingSphere::change(Random* rnd) {
		std::vector<Real> rand_trans(m_number_variables);
		for (size_t j = 0; j < m_number_variables; j++) {
			Real mean = (m_domain[j].limit.first + m_domain[j].limit.second) / 2;
			Real std = (m_domain[j].limit.second - m_domain[j].limit.first) / 3;
			rand_trans[j] = rnd->normal.nextNonStd(mean, std);
		}
		setTranslation(rand_trans.data(), rnd);
		setGlobalOpt(rand_trans.data());
	}
}