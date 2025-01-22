#include "dtlz2.h"

namespace ofec {
	void DTLZ2::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real g = 0;
		for (size_t i = m_number_objectives - 1; i < m_number_variables; i += 1)
			g += pow((x[i] - 0.5), 2);
		for (size_t m = 0; m < m_number_objectives; m += 1) {
			Real product = (1 + g);
			size_t i = 0;
			for (; i + m <= m_number_objectives - 2; i += 1)
				product *= cos(x[i] * OFEC_PI / 2);
			if (m > 0)
				product *= sin(x[i] * OFEC_PI / 2);
			obj[m] = product;
		}
	}
}