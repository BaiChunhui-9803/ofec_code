#include "wfg1.h"

namespace ofec {
	void WFG1::t1(std::vector<Real> &y) {		
		for (int i = m_k; i < y.size(); i++)
			y[i] = sLinear(y[i], 0.35);
	}

	void WFG1::t2(std::vector<Real> &y) {
		for (int i = m_k; i < y.size(); i++)
			y[i] = bFlat(y[i], 0.8, 0.75, 0.85);
	}

	void WFG1::t3(std::vector<Real> &y) {
		for (int i = m_k; i < y.size(); i++)
			y[i] = bPoly(y[i], 0.02);
	}

	void WFG1::t4(std::vector<Real> &y) {
		int n = y.size();
		std::vector<Real> w;
		for (int i = 1; i <= n; i++) {
			w.push_back(2.0*i);
		}
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

	void WFG1::shape(std::vector<Real> &y) {
		int M = y.size();		
		std::vector<Real> x = calculateX(y);
		for (int m = 1; m <= M - 1; m++)
			m_h[m - 1] = convex(x, m);
		m_h[M - 1] = mixed(x, 5, 1.0);
		y = x;
	}
}