#include "mmea_f4.h"

namespace ofec {
	void MMEA_F4::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real y = (x[0] + x[1]) / 2.;
		Real sum = 0.;
		Real g;
		for (size_t j = 2; j < m_number_variables; ++j) {
			Real h;
			if (j % 2) {
				h = 2 * x[j] - y * std::cos(2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables) - 1;
			}
			else {
				h = 2 * x[j] - y * std::sin(2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables) - 1;
			}
			sum = sum + h * h;
		}
		if (m_number_variables == 2) {
			g = 1 + 5. / (m_number_variables - 1.) * sum;
		}
		else {
			g = 1 + 5. / (m_number_variables - 2.) * sum;
		}

		obj[0] = y;
		obj[1] = g - std::pow(obj[0], 2);
	}

	void MMEA_F4::updateOptima() {
		m_optima.reset(new Optima<>);
		sampleParetoSols(100);
		//#ifdef OFEC_DEMO
		//		sampleParetoSets(100);
		//#endif //
	}

	void MMEA_F4::sampleParetoSols(size_t sample_num) {
		size_t div_num = sample_num;
		std::vector<VariableVector<Real>> all_vars;
		std::vector<std::vector<Real>> all_objs;
		for (size_t i = 0; i < div_num; i++) {
			for (size_t k = 0; k < div_num; ++k) {
				VariableVector<Real> temp_var(m_number_variables);
				std::vector<Real> temp_obj(m_number_objectives);
				for (size_t j = 0; j < m_number_variables; j++) {
					if (j == 0)
						temp_var[j] = (Real)i / (div_num - 1.);
					else if (j == 1) {
						temp_var[j] = (Real)k / (div_num - 1.);
					}
					else {
						Real y = (temp_var[0] + temp_var[1]) / 2.;
						if (j % 2) {
							temp_var[j] = 0.5 + 0.5 * y * std::cos(2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables);
						}
						else {
							temp_var[j] = 0.5 + 0.5 * y * std::sin(2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables);
						}
					}
				}
				all_vars.emplace_back(temp_var);
				evaluateObjective(all_vars[i].data(), temp_obj);
				all_objs.emplace_back(temp_obj);
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