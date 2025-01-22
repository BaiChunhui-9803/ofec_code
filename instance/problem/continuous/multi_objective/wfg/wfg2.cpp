#include "wfg2.h"

namespace ofec {
	void WFG2::t1(std::vector<Real> &y) {	//as t1 from WFG1
		for (int i = m_k; i < y.size(); i++) 
			y[i] = sLinear(y[i], 0.35);
	}

	void WFG2::t2(std::vector<Real> &y) {
		int n = y.size();
		int l = n - m_k;
		std::vector<Real> t(y.begin(), y.begin() + m_k);
		for (int i = m_k + 1; i <= m_k + l / 2; i++) {
			const int head = m_k + 2 * (i - m_k) - 2;
			const int tail = m_k + 2 * (i - m_k);
			t.push_back(rNonsep(subvector(y, head, tail), 2));
		}
		y = t;
	}

	void WFG2::t3(std::vector<Real> &y) {
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

	void WFG2::shape(std::vector<Real> &y) {
		int M = y.size();		
		std::vector<Real> x = calculateX(y);
		for (int m = 1; m <= M - 1; m++)
			m_h[m - 1] = convex(x, m);
		m_h[M - 1] = disc(x, 5, 1.0, 1.0);		
		y = x;
	}
}