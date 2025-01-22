#include "uf03.h"

namespace ofec {
	void UF03::evaluateObjective(Real *x, std::vector<Real> &obj) {
		int count1, count2;
		Real sum1, sum2, prod1, prod2, yj, pj;
		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		prod1 = prod2 = 1.0;
		for (int j = 2; j <= m_number_variables; j++) {
			if (m_number_variables == 2)
				yj = x[j - 1] - pow(x[0], 0.5 * (1.0 + 3.0 * (j - 1.0) / (m_number_variables - 1.0)));
			else
				yj = x[j - 1] - pow(x[0], 0.5*(1.0 + 3.0*(j - 2.0) / (m_number_variables - 2.0)));
			pj = cos(20.0*yj*OFEC_PI / std::sqrt(j + 0.0));
			if (j % 2 == 0){
				sum2 += yj * yj;
				prod2 *= pj;
				count2++;
			}
			else{
				sum1 += yj * yj;
				prod1 *= pj;
				count1++;
			}
		}
		if (m_number_variables == 2) {
			obj[0] = x[0];
			obj[1] = 1.0 - std::sqrt(x[0]) + 2.0 * (4.0 * sum2 - 2.0 * prod2 + 2.0) / (Real)count2;
		}
		else {
			obj[0] = x[0] + 2.0 * (4.0*sum1 - 2.0*prod1 + 2.0) / (Real)count1;
			obj[1] = 1.0 - std::sqrt(x[0]) + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (Real)count2;
		}
	}
}