#include "uf02.h"

namespace ofec {
	void UF02::evaluateObjective(Real *x, std::vector<Real> &obj) {
		int count1, count2;
		Real sum1, sum2, yj;
		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		for (int j = 2; j <= m_number_variables; j++) {
			if (j % 2 == 0)	{
				yj = x[j - 1] - 0.3*x[0] * (x[0] * cos(24.0*OFEC_PI*x[0] + 4.0*j*OFEC_PI / m_number_variables) + 2.0)*sin(6.0*OFEC_PI*x[0] + j * OFEC_PI / m_number_variables);
				sum2 += yj * yj;
				count2++;
			}
			else
			{
				yj = x[j - 1] - 0.3*x[0] * (x[0] * cos(24.0*OFEC_PI*x[0] + 4.0*j*OFEC_PI / m_number_variables) + 2.0)*cos(6.0*OFEC_PI*x[0] + j * OFEC_PI / m_number_variables);
				sum1 += yj * yj;
				count1++;
			}
		}
		if (count1 == 0)
			obj[0] = x[0];
		else
			obj[0] = x[0] + 2.0 * sum1 / (Real)count1;

		obj[1] = 1.0 - std::sqrt(x[0]) + 2.0 * sum2 / (Real)count2;
	}
}