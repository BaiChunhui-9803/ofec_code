#include "zdt2.h"

namespace ofec {
	void ZDT2::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real g = 0;
		for (size_t n = 1; n < m_number_variables; n++) {
			g = g + x[n];
		}
		g = 1 + 9 * g / (m_number_variables - 1);
		obj[0] = x[0];
		obj[1] = g*(1 - pow(obj[0] / g, 2));
	}

	void ZDT2::updateOptima() {
		m_optima.reset(new Optima<>);
		//loadParetoFront(10000);
		sampleParetoSols(m_num_reference_points);
		//#ifdef OFEC_DEMO
		//		sampleParetoSets(10000);
		//#endif //
	}

	void ZDT2::sampleParetoSols(size_t sample_num) {
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
