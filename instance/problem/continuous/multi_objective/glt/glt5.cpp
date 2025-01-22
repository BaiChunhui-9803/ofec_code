#include "glt5.h"

namespace ofec {
	void GLT5::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real yj, g = 0.;
		for (size_t j = 2; j < m_number_variables; j++) {
			yj = x[j] - sin(2 * OFEC_PI*x[0] + (1+j) * OFEC_PI / m_number_variables);
			g += yj * yj;
		}

		obj[0] = (1 + g)*(1 - cos((x[0] / 2)*OFEC_PI))*(1 - cos((x[1] / 2)*OFEC_PI));
		obj[1] = (1 + g)*(1 - cos((x[0] / 2)*OFEC_PI))*(1 - sin((x[1] / 2)*OFEC_PI));
		obj[2] = (1 + g) * (1 - sin((x[0] / 2) * OFEC_PI));
		
		/*obj[0] = (1 + g) * (1 - cos((x[0] / 2) * OFEC_PI) * cos((x[1] / 2) * OFEC_PI));
		obj[1] = (1 + g) * (1 - cos((x[0] / 2) * OFEC_PI) * sin((x[1] / 2) * OFEC_PI));
		obj[2] = (1 + g)*(1 - sin((x[0] / 2)*OFEC_PI));*/
	}

	void GLT5::updateOptima() {
		m_optima.reset(new Optima<>);
		sampleParetoSols(m_num_reference_points);
		//#ifdef OFEC_DEMO
		//		sampleParetoSets(m_num_reference_points);
		//#endif //
	}

	void GLT5::sampleParetoSols(size_t sample_num) {
		size_t div_num = 80;
		std::vector<VariableVector<Real>> all_vars;
		std::vector<std::vector<Real>> all_objs;
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
						temp_var[k] = std::sin(2 * OFEC_PI * temp_var[0] + (1 + k) * OFEC_PI / m_number_variables);
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