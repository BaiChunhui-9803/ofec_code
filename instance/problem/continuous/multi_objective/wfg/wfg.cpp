#include "wfg.h"
#include <fstream>
#include <cmath>
#include <algorithm>
#include "../../../../../core/global.h"

namespace ofec {
	void WFG::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;;
		resizeObjective(v.get<int>("number of objectives"));
		for (size_t i = 0; i < m_number_objectives; ++i)
			m_optimize_mode[i] = OptimizeMode::kMinimize;

		int k_factor = v.has("k factor") ? v.get<int>("k factor") : 1;
		m_k = k_factor * (m_number_objectives - 1);

		resizeVariable(v.get<int>("number of variables"));
		if (m_name == "MOP_WFG2" || m_name == "MOP_WFG3") {
			int l_factor = v.has("l factor") ? v.get<int>("l factor") : 2;
			if (m_number_variables != m_k * 2 * l_factor)
				throw MyExcept("Number of variable must be equal to (number of objectives - 1) * (k factor) * 2 * (l factor)");
			if (m_number_variables <= m_number_objectives)
				throw MyExcept("The number of variables should be no greater than the number of objectives.");
		}
		else {
			if (m_number_variables <= m_k)
				throw MyExcept("Number of variable must be no less than (number of objectives - 1) * (k factor)");
		}
		std::vector<std::pair<Real, Real>> r;
		for (int i = 0; i < m_number_variables; ++i)
			r.push_back(std::make_pair(0., 2. * (i + 1.)));
		setDomain(r);

		m_D = 1;
		if (m_name == "MOP_WFG3")
			m_A.assign(m_number_objectives - 1, 1);
		else
			m_A.assign(m_number_objectives - 1, 0);
		m_A[0] = 1;
		m_S.resize(m_number_objectives);
		for (int i = 0; i < m_number_objectives; ++i)
			m_S[i] = 2 * (i + 1.);
		m_h.resize(m_number_objectives);

		loadParetoFront();
	}

	void WFG::loadParetoFront() {
		std::ifstream infile;
		std::stringstream os;
		os << g_working_dir << "/instance/problem/continuous/multi_objective/wfg/data/" << m_name << "." << m_number_objectives << "D.pf";
		infile.open(os.str());
		if (!infile)
			throw MyExcept("open PF file of WFG problem is fail");
		std::string str;
		size_t line = 0;
		while (getline(infile, str))
			++line;
		m_optima->resizeObjectiveSet(line);
		infile.close();
		infile.clear();
		infile.open(os.str());
		for (size_t i = 0; i < line; i++) {
			std::vector<Real> temp_obj(m_number_objectives);
			for (size_t j = 0; j < m_number_objectives; j++)
				infile >> temp_obj[j];
			m_optima->setObjective(temp_obj, i);
		}
		m_optima->setObjectiveGiven(true);
		infile.close();
	}

	void WFG::evaluateObjective(Real *x, std::vector<Real> &obj) {
		std::vector<Real> y = normalise(x);	//normalization
		t1(y);
		t2(y);
		t3(y);
		t4(y);
		shape(y); //calculate x, and set m_h
		for (int i = 0; i < m_number_objectives; i++) {	//calculate objective value
			obj[i] = m_D * y[m_number_objectives - 1] + m_S[i] * m_h[i];
		}
	}

	std::vector<Real> WFG::calculateX(const std::vector<Real> &t_p) {
		std::vector<Real> x(t_p.size());
		int i;
		for (i = 0; i < (t_p.size() - 1); i++) {
			Real tmp1 = std::max<Real>(t_p.back(), m_A[i]);
			x[i] = tmp1 * (t_p[i] - 0.5) + 0.5;
		}
		x[i] = t_p.back();
		return x;
	}

	std::vector<Real> WFG::normalise(Real *z) {
		std::vector<Real> y(m_number_variables);
		for (int i = 0; i < m_number_variables; i++) {
			Real z_max = 2 * (i + 1.);
			y[i] = z[i] / z_max;
		}
		return y;
	}

	std::vector<Real> WFG::subvector(const std::vector<Real> &v, int head, int tail) {
		std::vector<Real> result;
		for (int i = head; i < tail; i++) {
			result.push_back(v[i]);
		}
		return result;
	}

	/*Shape Function*/
	Real WFG::linear(const std::vector<Real> &x, int m) {
		int M = static_cast<int>(x.size());
		Real result = 1.0;
		for (int i = 1; i <= M - m; i++) {
			result *= x[i - 1];
		}
		if (m != 1) {
			result *= 1 - x[M - m];
		}
		return correctTo01(result);
	}

	Real WFG::convex(const std::vector<Real> &x, int m) {
		int M = x.size();
		Real result = 1.0;
		for (int i = 1; i <= M - m; i++) {
			result *= 1.0 - cos(x[i - 1] * OFEC_PI / 2.0);
		}
		if (m != 1) {
			result *= 1.0 - sin(x[M - m] * OFEC_PI / 2.0);
		}
		return correctTo01(result);
	}

	Real WFG::concave(const std::vector<Real> &x, int m) {
		int M = static_cast<int>(x.size());
		Real result = 1.0;
		for (int i = 1; i <= M - m; i++) {
			result *= sin(x[i - 1] * OFEC_PI / 2.0);
		}
		if (m != 1) {
			result *= cos(x[M - m] * OFEC_PI / 2.0);
		}
		return correctTo01(result);
	}

	Real WFG::mixed(const std::vector<Real> &x, int A, Real alpha) {
		Real tmp = 2.0 * A * OFEC_PI;
		return correctTo01(pow(1.0 - x[0] - cos(tmp * x[0] + OFEC_PI / 2.0) / tmp, alpha));
	}

	Real WFG::disc(const std::vector<Real> &x, int A, Real alpha, Real beta) {
		Real tmp1 = A * pow(x[0], beta) * OFEC_PI;
		return correctTo01(1.0 - pow(x[0], alpha) * pow(cos(tmp1), 2.0));
	}

	Real WFG::correctTo01(const Real &a, Real epsilon) {
		Real min = 0.0;
		Real max = 1.0;
		Real min_epsilon = min - epsilon;
		Real max_epsilon = max + epsilon;
		if (a <= min && a >= min_epsilon) {
			return min;
		}
		else if (a >= max && a <= max_epsilon) {
			return max;
		}
		else {
			return a;
		}
	}

	/*TransFunctions*/
	Real WFG::bPoly(Real y, Real alpha) {
		return correctTo01(pow(y, alpha));
	}

	Real WFG::bFlat(Real y, Real A, Real B, Real C) {
		Real tmp1 = std::min<Real>(0.0, floor(y - B)) * A * (B - y) / B;
		Real tmp2 = std::min<Real>(0.0, floor(C - y)) * (1.0 - A) * (y - C) / (1.0 - C);
		return correctTo01(A + tmp1 - tmp2);
	}

	Real WFG::bParam(Real y, Real u, Real A, Real B, Real C) {
		Real v = A - (1.0 - 2.0 * u) * fabs(floor(0.5 - u) + A);
		return correctTo01(pow(y, B + (C - B) * v));
	}

	Real WFG::sLinear(Real y, Real A) {
		return correctTo01(fabs(y - A) / fabs(floor(A - y) + A));
	}

	Real WFG::sDecept(Real y, Real A, Real B, Real C) {
		Real tmp1 = floor(y - A + B) * (1.0 - C + (A - B) / B) / (A - B);
		Real tmp2 = floor(A + B - y) * (1.0 - C + (1.0 - A - B) / B) / (1.0 - A - B);
		return correctTo01(1.0 + (fabs(y - A) - B) * (tmp1 + tmp2 + 1.0 / B));
	}

	Real WFG::sMulti(Real y, int A, Real B, Real C) {
		Real tmp1 = fabs(y - C) / (2.0 * (floor(C - y) + C));
		Real tmp2 = (4.0 * A + 2.0) * OFEC_PI * (0.5 - tmp1);
		return correctTo01((1.0 + cos(tmp2) + 4.0 * B * pow(tmp1, 2.0)) / (B + 2.0));
	}

	Real WFG::rSum(const std::vector<Real> &y, const std::vector<Real> &w) {
		Real numerator = 0.0;
		Real denominator = 0.0;
		for (int i = 0; i < static_cast<int>(y.size()); i++) {
			numerator += w[i] * y[i];
			denominator += w[i];
		}
		return correctTo01(numerator / denominator);
	}

	Real WFG::rNonsep(const std::vector<Real> &y, int A) {
		int y_len = static_cast<int>(y.size());
		Real numerator = 0.0;
		for (int j = 0; j < y_len; j++) {
			numerator += y[j];
			for (int k = 0; k <= A - 2; k++) {
				numerator += fabs(y[j] - y[(j + k + 1) % y_len]);
			}
		}
		Real tmp = ceil(A / 2.0);
		Real denominator = y_len * tmp * (1.0 + 2.0 * A - 2.0 * tmp) / A;
		return correctTo01(numerator / denominator);
	}

}