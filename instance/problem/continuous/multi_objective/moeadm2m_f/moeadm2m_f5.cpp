#include "moeadm2m_f5.h"

namespace ofec {
	void MOEADM2M_F5::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real sum = 0.;
		for (size_t i = 1; i < m_number_variables; ++i) {
			Real t = x[i] - std::sin(0.5 * OFEC_PI * x[0]);
			sum = sum + std::pow(std::fabs(t), 0.6) - 0.9 * std::pow(t, 2);
		}
		Real g = 2 * std::fabs(std::cos(OFEC_PI * x[0])) * sum;

		obj[0] = (1 + g) * x[0];
		obj[1] = (1 + g) * (1 - std::pow(x[0], 0.5));
	}

	void MOEADM2M_F5::updateOptima() {
		m_optima.reset(new Optima<>);
		sampleParetoSols(m_num_reference_points);
		//#ifdef OFEC_DEMO
		//		sampleParetoSets(10000);
		//#endif //
	}

	void MOEADM2M_F5::sampleParetoSols(size_t sample_num) {
		size_t num = sample_num;
		std::vector<VariableVector<Real>> all_vars;
		std::vector<std::vector<Real>> all_objs;
		for (size_t i = 0; i < num; i++) {
			VariableVector<Real> temp_var(m_number_variables);
			std::vector<Real> temp_obj(m_number_objectives);
			for (size_t j = 0; j < m_number_variables; j++) {
				if (j == 0)
					temp_var[j] = (Real)i / (num - 1.);
				else {
					temp_var[j] = std::sin(0.5 * OFEC_PI * temp_var[0]);
				}
			}
			all_vars.emplace_back(temp_var);
			evaluateObjective(all_vars[i].data(), temp_obj);
			all_objs.emplace_back(temp_obj);
		}
		for (size_t i = 0; i < all_vars.size(); ++i) {
			Solution<> sol(m_number_objectives, m_number_constraints);
			sol.variable() = all_vars[i];
			sol.objective() = all_objs[i];
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(sol);
		}
	}
}