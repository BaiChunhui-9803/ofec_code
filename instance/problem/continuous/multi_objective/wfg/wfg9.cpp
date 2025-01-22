#include "wfg9.h"

namespace ofec {
	void WFG9::t1(std::vector<Real> &y) {
		int n = y.size();		
		const std::vector<Real> w(n, 1.0);
		std::vector<Real> t;
		for (int i = 0; i < n - 1; i++)	{
			const std::vector<Real>& y_sub = subvector(y, i + 1, n);
			const std::vector<Real>& w_sub = subvector(w, i + 1, n);
			const Real u = rSum(y_sub, w_sub);
			t.push_back(bParam(y[i], u, 0.98 / 49.98, 0.02, 50));
		}
		t.push_back(y.back());
		y = t;
	}

	void WFG9::t2(std::vector<Real> &y) {
		for (int i = m_k; i < y.size(); i++)
			y[i] = sMulti(y[i], 30, 95, 0.35);
	}

	void WFG9::t3(std::vector<Real> &y) {	//as t2 from WFG6
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
	
	void WFG9::shape(std::vector<Real> &y) {
		int M = y.size();
		std::vector<Real> x = calculateX(y);
		for (int m = 1; m <= M; m++)
			m_h[m - 1] = concave(x, m);
		y = x;
	}
}