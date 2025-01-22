#include "wfg4.h"

namespace ofec {
	void WFG4::t1(std::vector<Real> &y) {
		for (int i = 0; i < y.size(); i++)
			y[i] = sMulti(y[i], 30, 10, 0.35);
	}

	void WFG4::t2(std::vector<Real> &y) {
		int n = y.size();
		std::vector<Real> w(n, 1.0);
		std::vector<Real> t;
		for (int i = 1; i <= m_number_objectives - 1; i++) {
			const int head = (i - 1)*m_k / (m_number_objectives - 1);
			const int tail = i * m_k / (m_number_objectives - 1);
			const std::vector<Real>& y_sub = subvector(y, head, tail);
			const std::vector<Real>& w_sub = subvector(w, head, tail);
			t.push_back(rSum(y_sub, w_sub));
		}
		const std::vector<Real>& y_sub = subvector(y, m_k, n);
		const std::vector<Real>& w_sub = subvector(w, m_k, n);
		t.push_back(rSum(y_sub, w_sub));
		y = t; 
	}
	
	void WFG4::shape(std::vector<Real> &y) {
		int M = y.size();
		std::vector<Real> x = calculateX(y);
		for (int m = 1; m <= M; m++)
			m_h[m - 1] = concave(x, m);
		y = x;
	}
}