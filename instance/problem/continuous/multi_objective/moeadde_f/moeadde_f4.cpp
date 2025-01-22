#include "moeadde_f4.h"

namespace ofec {
	void MOEADDE_F4::evaluateObjective(Real* x, std::vector<Real>& obj) {
		std::vector<size_t> obj1_dim_inx;
		std::vector<size_t> obj2_dim_inx;
		for (size_t i = 1; i < m_number_variables; i++) {
			if (i % 2 == 0) {
				obj1_dim_inx.push_back(i);
			}
			else {
				obj2_dim_inx.push_back(i);
			}
		}
		Real beta1 = 0.;
		Real beta2 = 0.;
		for (size_t i = 0; i < obj1_dim_inx.size(); ++i) {
			beta1 += (std::pow(x[obj1_dim_inx[i]] - 0.8 * x[0] * std::cos((6 * OFEC_PI * x[0] + (1. + obj1_dim_inx[i]) * OFEC_PI / m_number_variables)/3.), 2));
		}
		for (size_t i = 0; i < obj2_dim_inx.size(); ++i) {
			beta2 += (std::pow(x[obj2_dim_inx[i]] - 0.8 * x[0] * std::sin(6 * OFEC_PI * x[0] + (1. + obj2_dim_inx[i]) * OFEC_PI / m_number_variables), 2));
		}
		if (m_number_variables < 3) {
			obj[0] = x[0];
		}
		else {
			obj[0] = x[0] + 2. / obj1_dim_inx.size() * beta1;
		}
		obj[1] = 1 - std::sqrt(x[0]) + 2. / obj2_dim_inx.size() * beta2;
	}

	void MOEADDE_F4::updateOptima() {
		m_optima.reset(new Optima<>);
		sampleParetoSols(m_num_reference_points);
//#ifdef OFEC_DEMO
//		sampleParetoSols(m_num_reference_points);
//#endif // 
	}

	void MOEADDE_F4::sampleParetoSols(size_t sample_num) {
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
					if (j % 2 == 1)
						temp_var[j] = 0.8 * temp_var[0] * std::sin(6 * OFEC_PI * temp_var[0] + (j + 1) * OFEC_PI / m_number_variables);
					else
						temp_var[j] = 0.8 * temp_var[0] * std::cos((6 * OFEC_PI * temp_var[0] + (j + 1) * OFEC_PI / m_number_variables) / 3);
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