#include "wfg6.h"

namespace ofec {
	void WFG6::t1(std::vector<Real> &y) {	//as t1 from WFG1
		for (int i = m_k; i < y.size(); i++)
			y[i] = sLinear(y[i], 0.35);
	}

	void WFG6::t2(std::vector<Real> &y) {
		int n = y.size();
		std::vector<Real> t;
		for (int i = 1; i <= m_number_objectives - 1; i++) {
			const int head = (i - 1)*m_k / (m_number_objectives - 1);
			const int tail = i * m_k / (m_number_objectives - 1);
			const std::vector<Real>& y_sub = subvector(y, head, tail);
			t.push_back(rNonsep(y_sub, m_k / (m_number_objectives - 1)));
		}
		const std::vector<Real>& y_sub = subvector(y, m_k, n);
		t.push_back(rNonsep(y_sub, n - m_k));
		y = t;
	}

	void WFG6::shape(std::vector<Real> &y) {
		int M = y.size();
		std::vector<Real> x = calculateX(y);
		for (int m = 1; m <= M; m++) 
			m_h[m - 1] = concave(x, m);
		y = x;
	}
}