#include "dtlz7.h"
 
namespace ofec {
	void DTLZ7::evaluateObjective(Real *x, std::vector<Real> &obj) {
		int i = 0;
		int j = 0;
		int n = m_number_variables;
		int k = n - m_number_objectives + 1;
		Real g = 0;
		Real h = 0;
		for (i = n - k + 1; i <= n; i++)
			g += x[i - 1];
		g = 1 + 9 * g / k;
		for (i = 1; i <= m_number_objectives - 1; i++)
			obj[i - 1] = x[i - 1];
		for (j = 1; j <= m_number_objectives - 1; j++)
			h += x[j - 1] / (1 + g) * (1 + sin(3 * OFEC_PI * x[j - 1]));
		h = m_number_objectives - h;
		obj[m_number_objectives - 1] = (1 + g) * h;
	}
}