#include "uf06.h"

namespace ofec {
	void UF06::evaluateObjective(Real *x, std::vector<Real> &obj) {
		int count1, count2;
		Real sum1, sum2, prod1, prod2, yj, hj, pj, N, E;
		N = 2.0; E = 0.1;
		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		prod1 = prod2 = 1.0;
		for (int j = 2; j <= m_number_variables; j++) {
			yj = x[j - 1] - sin(6.0*OFEC_PI*x[0] + j * OFEC_PI / m_number_variables);
			pj = cos(20.0*yj*OFEC_PI / std::sqrt(j + 0.0));
			if (j % 2 == 0) {
				sum2 += yj * yj;
				prod2 *= pj;
				count2++;
			}
			else {
				sum1 += yj * yj;
				prod1 *= pj;
				count1++;
			}
		}
		hj = 2.0*(0.5 / N + E)*sin(2.0*N*OFEC_PI*x[0]);
		if (hj < 0.0) hj = 0.0;
		if (count1 == 0)
			obj[0] = x[0] + hj;
		else
			obj[0] = x[0] + hj + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (Real)count1;

		obj[1] = 1.0 - x[0] + hj + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (Real)count2;
	}
}