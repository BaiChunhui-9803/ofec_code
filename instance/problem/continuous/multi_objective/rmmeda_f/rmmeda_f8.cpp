#include "rmmeda_f8.h"

namespace ofec {
	void RMMEDA_F8::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real g = 1.;
		Real sum = 0.;
		for (size_t i = 2; i < m_number_variables; ++i) {
			Real temp = std::pow(x[i],2) - x[0];
			sum += (temp * temp);
		}
		g = sum;

		obj[0] = (1. + g) * std::cos(OFEC_PI / 2 * x[0]) * std::cos(OFEC_PI / 2 * x[1]);
		obj[1] = (1. + g) * std::cos(OFEC_PI / 2 * x[0]) * std::sin(OFEC_PI / 2 * x[1]);
		obj[2] = (1. + g) * std::sin(OFEC_PI / 2 * x[0]);
	}

	void RMMEDA_F8::updateOptima() {
		m_optima.reset(new Optima<>);
		sampleParetoSols(m_num_reference_points);
		//#ifdef OFEC_DEMO
		//		sampleParetoSets(10000);
		//#endif //
	}

	void RMMEDA_F8::sampleParetoSols(size_t sample_num) {
		size_t num = sample_num;
		std::vector<VariableVector<Real>> all_vars;
		std::vector<std::vector<Real>> all_objs;
		
		size_t div_num = 100;
		for (size_t i = 0; i < div_num; i++) {
			for (size_t j = 0; j < div_num; j++) {
				VariableVector<Real> temp_var(m_number_variables);
				std::vector<Real> temp_obj(m_number_objectives);
				for (size_t k = 0; k < m_number_variables; k++) {
					if (k == 0) {
						temp_var[k] = (Real)i / (div_num - 1.);
					}
					else if (k == 1) {
						temp_var[k] = (Real)j / (div_num - 1.);
					}
					else {
						temp_var[k] = std::pow(temp_var[0], 0.5);
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