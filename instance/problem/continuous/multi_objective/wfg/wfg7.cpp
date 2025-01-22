#include "wfg7.h"

namespace ofec {
	void WFG7::t1(std::vector<Real> &y) {
		int n = y.size();
		const std::vector<Real> w(n, 1.0);
		std::vector<Real> t;
		for (int i = 0; i < m_k; i++)	{
			const std::vector<Real>& y_sub = subvector(y, i + 1, n);
			const std::vector<Real>& w_sub = subvector(w, i + 1, n);
			const double u = rSum(y_sub, w_sub);
			t.push_back(bParam(y[i], u, 0.98 / 49.98, 0.02, 50));
		}
		for (int i = m_k; i < n; i++)
			t.push_back(y[i]);
		y = t;
	}

	void WFG7::t2(std::vector<Real> &y) {	//as t1 from WFG1
		for (int i = m_k; i < y.size(); i++)
			y[i] = sLinear(y[i], 0.35);
	}

	void WFG7::t3(std::vector<Real> &y) {	//as t2 from WFG4
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
		y= t;
	}
	
	void WFG7::shape(std::vector<Real> &y) {
		int M = y.size();
		std::vector<Real> x = calculateX(y);
		for (int m = 1; m <= M; m++) 
			m_h[m - 1] = concave(x, m);
		y = x;
	}
}