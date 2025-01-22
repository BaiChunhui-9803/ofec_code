#include "zdt4.h"

namespace ofec {
	void ZDT4::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real g = 0;
		for (size_t n = 1; n < m_number_variables; n++) {
			//g = g + (pow(x[n], 2) - 10 * cos(4 * OFEC_PI * x[n]));
			g = g + (pow(x[n], 2) - 10 * cos(4 * OFEC_PI * x[n]));
		}
		g = 1 + 10 * (m_number_variables - 1) + g;
		obj[0] = x[0];
		obj[1] = g*(1 - std::sqrt(x[0] / g));
	}

	void ZDT4::updateOptima() {
		m_optima.reset(new Optima<>);
		//loadParetoFront(10000);
		sampleParetoSols(m_num_reference_points);
		//#ifdef OFEC_DEMO
		//		sampleParetoSets(10000);
		//#endif //
	}

	void ZDT4::sampleParetoSols(size_t sample_num) {
		std::vector<VariableVector<Real>> all_vars;
		std::vector<std::vector<Real>> all_objs;
		for (size_t i = 0; i < sample_num; i++) {
			VariableVector<Real> temp_var(m_number_variables);
			std::vector<Real> temp_obj(m_number_objectives);
			for (size_t j = 0; j < m_number_variables; j++) {
				if (j == 0) {
					temp_var[j] = (Real)i / (sample_num - 1.);
				}
				else {
					temp_var[j] = 0.;
				}
			}
			Solution<> sol(m_number_objectives, m_number_constraints);
			sol.variable() = temp_var;
			evaluateObjective(temp_var.data(), temp_obj);
			sol.objective() = temp_obj;
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(sol);
		}
	}
}
