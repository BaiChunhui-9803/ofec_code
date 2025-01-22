#include "moeadde_f6.h"

namespace ofec {
	void MOEADDE_F6::evaluateObjective(Real* x, std::vector<Real>& obj) {
		std::vector<size_t> obj1_dim_inx;
		std::vector<size_t> obj2_dim_inx;
		std::vector<size_t> obj3_dim_inx;
		for (size_t i = 2; i < m_number_variables; i++) {
			if (i % 3 == 0) {
				obj1_dim_inx.push_back(i);
			}
			else if(i % 3 == 1) {
				obj2_dim_inx.push_back(i);
			}
			else {
				obj3_dim_inx.push_back(i);
			}
		}
		Real beta1 = 0.;
		Real beta2 = 0.;
		Real beta3 = 0.;
		for (size_t i = 0; i < obj1_dim_inx.size(); ++i) {
			beta1 += (std::pow(x[obj1_dim_inx[i]] - 2. * x[1] * std::sin(2 * OFEC_PI * x[0] + (1. + obj1_dim_inx[i]) * OFEC_PI / m_number_variables), 2));
		}
		for (size_t i = 0; i < obj2_dim_inx.size(); ++i) {
			beta2 += (std::pow(x[obj2_dim_inx[i]] - 2. * x[1] * std::sin(2 * OFEC_PI * x[0] + (1. + obj2_dim_inx[i]) * OFEC_PI / m_number_variables), 2));
		}
		for (size_t i = 0; i < obj3_dim_inx.size(); ++i) {
			beta3 += (std::pow(x[obj3_dim_inx[i]] - 2. * x[1] * std::sin(2 * OFEC_PI * x[0] + (1. + obj3_dim_inx[i]) * OFEC_PI / m_number_variables), 2));
		}
		if (m_number_variables < 3) {
			obj[0] = std::cos(0.5*OFEC_PI*x[0]) * std::cos(0.5 * OFEC_PI * x[1]);
			obj[1] = std::cos(0.5 * OFEC_PI * x[0]) * std::sin(0.5 * OFEC_PI * x[1]);
			obj[2] = std::sin(0.5 * OFEC_PI * x[0]);
		}
		else if(m_number_variables < 4){
			obj[0] = std::cos(0.5 * OFEC_PI * x[0]) * std::cos(0.5 * OFEC_PI * x[1]);
			obj[1] = std::cos(0.5 * OFEC_PI * x[0]) * std::sin(0.5 * OFEC_PI * x[1]);
			obj[2] = std::sin(0.5 * OFEC_PI * x[0]) + 2. / obj3_dim_inx.size() * beta3;
		}
		else if(m_number_variables < 5){
			obj[0] = std::cos(0.5 * OFEC_PI * x[0]) * std::cos(0.5 * OFEC_PI * x[1])+ 2. / obj1_dim_inx.size() * beta1;
			obj[1] = std::cos(0.5 * OFEC_PI * x[0]) * std::sin(0.5 * OFEC_PI * x[1]);
			obj[2] = std::sin(0.5 * OFEC_PI * x[0]) + 2. / obj3_dim_inx.size() * beta3;
		}
		else {
			obj[0] = std::cos(0.5 * OFEC_PI * x[0]) * std::cos(0.5 * OFEC_PI * x[1]) + 2. / obj1_dim_inx.size() * beta1;
			obj[1] = std::cos(0.5 * OFEC_PI * x[0]) * std::sin(0.5 * OFEC_PI * x[1])+ 2. / obj2_dim_inx.size() * beta2;
			obj[2] = std::sin(0.5 * OFEC_PI * x[0]) + 2. / obj3_dim_inx.size() * beta3;
		}
	}

	void MOEADDE_F6::updateOptima() {
		m_optima.reset(new Optima<>);
		sampleParetoSols(m_num_reference_points);
//#ifdef OFEC_DEMO
//		sampleParetoSols(m_num_reference_points);
//#endif // 
	}

	void MOEADDE_F6::sampleParetoSols(size_t sample_num) {
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
						temp_var[k] = 2 * temp_var[1] * std::sin(2 * OFEC_PI * temp_var[0] + (1 + k) * OFEC_PI / m_number_variables);
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