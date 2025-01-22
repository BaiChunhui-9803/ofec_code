#include "mmea_f7.h"

namespace ofec {
	void MMEA_F7::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real y = (x[0] + x[1] + x[2]) / 2.;
		Real sum = 0.;
		Real g;
		for (size_t j = 3; j < m_number_variables; ++j) {
			Real h;
			if (j % 2) {
				h = 2 * x[j] - std::sin(0.5 * OFEC_PI * y) * std::cos(2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables) - 1;
			}
			else {
				h = 2 * x[j] - std::cos(0.5 * OFEC_PI * y) * std::sin(1. / 3 * (2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables)) - 1;
			}
			sum = sum + h * h;
		}
		if (m_number_variables == 3) {
			g = 1 + 5. / (m_number_variables - 2.) * sum;
		}
		else {
			g = 1 + 5. / (m_number_variables - 3.) * sum;
		}

		obj[0] = g * std::cos(0.25 * OFEC_PI * (x[0] + x[1])) * std::sin(0.5 * OFEC_PI * x[2]);
		obj[1] = g * std::cos(0.25 * OFEC_PI * (x[0] + x[1])) * std::cos(0.5 * OFEC_PI * x[2]);
		obj[2] = g * std::sin(0.25 * OFEC_PI * (x[0] + x[1]));
	}

	void MMEA_F7::updateOptima() {
		m_optima.reset(new Optima<>);
		sampleParetoSols(100);
		//#ifdef OFEC_DEMO
		//		sampleParetoSets(100);
		//#endif //
	}

	void MMEA_F7::sampleParetoSols(size_t sample_num) {
		std::vector<VariableVector<Real>> all_vars;
		std::vector<std::vector<Real>> all_objs;
		size_t div_num = 25;
		for (size_t i = 0; i < div_num; i++) {
			for (size_t j = 0; j < div_num; j++) {
				for (size_t k = 0; k < div_num; k++) {
					VariableVector<Real> temp_var(m_number_variables);
					std::vector<Real> temp_obj(m_number_objectives);
					for (size_t p = 0; p < m_number_variables; p++) {
						if (p == 0) {
							temp_var[p] = (Real)i / (div_num - 1.);
						}
						else if (p == 1) {
							temp_var[p] = (Real)j / (div_num - 1.);
						}
						else if (p == 2) {
							temp_var[p] = (Real)k / (div_num - 1.);
						}
						else {
							Real y = (temp_var[0] + temp_var[1] + temp_var[2]) / 3.;
							if (p % 2) {
								temp_var[p] = 0.5 + 0.5 * std::sin(0.5 * OFEC_PI * y) * std::cos(2 * OFEC_PI * y + (p + 1) * OFEC_PI / m_number_variables);
							}
							else {
								temp_var[p] = 0.5 + 0.5 * std::cos(0.5 * OFEC_PI * y) * std::sin(1. / 3 * (2 * OFEC_PI * y + (p + 1) * OFEC_PI / m_number_variables));
							}
						}
						all_vars.emplace_back(temp_var);
						evaluateObjective(all_vars[i].data(), temp_obj);
						all_objs.emplace_back(temp_obj);
					}
				}
			}
		}
		for (size_t i = 0; i < all_vars.size(); ++i) {
			Solution<> sol(m_number_objectives, m_number_constraints);
			sol.variable() = all_vars[i];
			sol.objective() = all_objs[i];
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(sol);
		}
	}
}