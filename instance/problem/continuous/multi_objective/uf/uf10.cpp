#include "uf10.h"

namespace ofec {
	void UF10::evaluateObjective(Real *x, std::vector<Real> &obj) {
		int count1, count2, count3;
		Real sum1, sum2, sum3, yj, hj;
		sum1 = sum2 = sum3 = 0.0;
		count1 = count2 = count3 = 0;
		for (int j = 2; j <= m_number_variables; j++) {
			yj = x[j - 1] - 2.0*x[1] * sin(2.0*OFEC_PI*x[0] + j * OFEC_PI / m_number_variables);
			hj = 4.0*yj*yj - cos(8.0*OFEC_PI*yj) + 1.0;
			if (j % 3 == 1)	{
				sum1 += hj;
				count1++;
			}
			else if (j % 3 == 2) {
				sum2 += hj;
				count2++;
			}
			else {
				sum3 += hj;
				count3++;
			}
		}
		if (count1 == 0)
			obj[0] = cos(0.5*OFEC_PI*x[0])*cos(0.5*OFEC_PI*x[1]);
		else
			obj[0] = cos(0.5*OFEC_PI*x[0])*cos(0.5*OFEC_PI*x[1]) + 2.0*sum1 / (Real)count1;

		if (count2 == 0)
			obj[1] = cos(0.5*OFEC_PI*x[0])*sin(0.5*OFEC_PI*x[1]);
		else
			obj[1] = cos(0.5*OFEC_PI*x[0])*sin(0.5*OFEC_PI*x[1]) + 2.0*sum2 / (Real)count2;

		if (count3 == 0)
			obj[2] = sin(0.5*OFEC_PI*x[0]);
		else
			obj[2] = sin(0.5*OFEC_PI*x[0]) + 2.0*sum3 / (Real)count3;
	}
}