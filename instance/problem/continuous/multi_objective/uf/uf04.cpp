#include "uf04.h"

namespace ofec {
	void UF04::evaluateObjective(Real *x, std::vector<Real> &obj) {
		int count1, count2;
		Real sum1, sum2, yj, hj;
		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		for (int j = 2; j <= m_number_variables; j++) {
			yj = x[j - 1] - sin(6.0*OFEC_PI*x[0] + j * OFEC_PI / m_number_variables);
			hj = fabs(yj) / (1.0 + exp(2.0*fabs(yj)));
			if (j % 2 == 0) {
				sum2 += hj;
				count2++;
			}
			else {
				sum1 += hj;
				count1++;
			}
		}
		if (count1 == 0)
			obj[0] = x[0];
		else
			obj[0] = x[0] + 2.0 * sum1 / (Real)count1;

		obj[1] = 1.0 - x[0] * x[0] + 2.0*sum2 / (Real)count2;
	}
}