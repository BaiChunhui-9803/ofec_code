#include "dtlz1.h"

namespace ofec {
	void DTLZ1::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real g = 0;
		for (size_t i = m_number_objectives - 1; i < m_number_variables; ++i)
			g += (x[i] - 0.5)*(x[i] - 0.5) - cos(20 * OFEC_PI*(x[i] - 0.5));
		g = (m_number_variables + 1 - m_number_objectives + g) * 100;
		for (size_t m = 0; m < m_number_objectives;++m) {
			Real product = 0.5*(1 + g);
			size_t i = 0;
			for (; m_number_objectives >= 2 + m && i <= m_number_objectives - 2 - m; i += 1)
				product *= x[i];
			if (m > 0)
				product *= (1 - x[i]);
			obj[m] = product;
		}
	}
}