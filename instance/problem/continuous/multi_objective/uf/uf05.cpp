#include "uf05.h"

namespace ofec {
	void UF05::evaluateObjective(Real *x, std::vector<Real> &obj) {
		int count1, count2;
		Real sum1, sum2, yj, hj, N, E;
		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		N = 10.0; E = 0.1;
		for (int j = 2; j <= m_number_variables; j++) {
			yj = x[j - 1] - sin(6.0*OFEC_PI*x[0] + j * OFEC_PI / m_number_variables);
			hj = 2.0*yj*yj - cos(4.0*OFEC_PI*yj) + 1.0;
			if (j % 2 == 0) {
				sum2 += hj;
				count2++;
			}
			else {
				sum1 += hj;
				count1++;
			}
		}
		hj = (0.5 / N + E)*fabs(sin(2.0*N*OFEC_PI*x[0]));
		if (count1 == 0)
			obj[0] = x[0] + hj;
		else
			obj[0] = x[0] + hj + 2.0*sum1 / (Real)count1;

		obj[1] = 1.0 - x[0] + hj + 2.0*sum2 / (Real)count2;
	}
}