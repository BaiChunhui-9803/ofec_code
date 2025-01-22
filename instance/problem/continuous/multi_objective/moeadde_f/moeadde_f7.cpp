#include "moeadde_f7.h"

namespace ofec {
	void MOEADDE_F7::evaluateObjective(Real* x, std::vector<Real>& obj) {
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
			Real temp = x[obj1_dim_inx[i]] - std::pow(x[0], 0.5 * (1. + 3. * (obj1_dim_inx[i] - 1) / (m_number_variables - 2)));
			beta1 += (1. + 4 * temp * temp - std::cos(8 * temp * OFEC_PI));
		}
		for (size_t i = 0; i < obj2_dim_inx.size(); ++i) {
			Real temp = 1.;
			if (m_number_variables == 2) {
				temp = x[obj2_dim_inx[i]] - std::pow(x[0], 0.5 * (1. + 3. * obj2_dim_inx[i] / (m_number_variables - 1)));
			}
			else {
				temp = x[obj2_dim_inx[i]] - std::pow(x[0], 0.5 * (1. + 3. * (obj2_dim_inx[i] - 1) / (m_number_variables - 2)));
			}
			
			beta2 += (1. + 4 * temp * temp - std::cos(8 * temp * OFEC_PI));
		}
		if (m_number_variables < 3) {
			obj[0] = x[0];
		}
		else {
			obj[0] = x[0] + 2. / obj1_dim_inx.size() * beta1;
		}
		obj[1] = 1 - std::sqrt(x[0]) + 2. / obj2_dim_inx.size() * beta2;
	}

	void MOEADDE_F7::updateOptima() {
		m_optima.reset(new Optima<>);
		sampleParetoSols(m_num_reference_points);
//#ifdef OFEC_DEMO
//		sampleParetoSols(m_num_reference_points);
//#endif // 
	}

	void MOEADDE_F7::sampleParetoSols(size_t sample_num) {
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
					if (m_number_variables == 2)
						temp_var[j] = std::pow(temp_var[0], 0.5 * (1. + 3 * j / (m_number_variables - 1)));
					else
						temp_var[j] = std::pow(temp_var[0], 0.5 * (1. + 3 * (j - 1) / (m_number_variables - 2)));
				}
			}
			all_vars.emplace_back(temp_var);
			evaluateObjective(all_vars[i].data(), temp_obj);
			all_objs.emplace_back(temp_obj);
		}
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < all_objs.size(); ++i) {
			objs.emplace_back(&all_objs[i]);
		}
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, m_optimize_mode);
		for (size_t i = 0; i < rank.size(); ++i) {
			if (rank[i] == 0) {
				Solution<> sol(m_number_objectives, m_number_constraints);
				sol.variable() = all_vars[i];
				sol.objective() = all_objs[i];
				dynamic_cast<Optima<>&>(*m_optima).appendSolution(sol);
			}
		}
	}
}