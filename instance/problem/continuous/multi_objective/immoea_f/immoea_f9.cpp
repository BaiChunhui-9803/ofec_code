#include "immoea_f9.h"

namespace ofec {
	void IMMOEA_F9::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real g = 1.;
		Real sum1 = 0.;
		Real sum2 = 1.;
		for (size_t i = 1; i < m_number_variables; ++i) {
			Real temp1 = std::pow(std::pow(x[i], 1. / (1. + m_beta * (i + 1.) / m_number_variables)) - x[0], 2)/4000;
			sum1 += temp1;
			Real temp2 = std::cos((std::pow(x[i], 1. / (1. + m_beta * (i + 1.) / m_number_variables)) - x[0])/(std::sqrt(i)));
			sum2 *= temp2;
		}
		g = 2 + sum1-sum2;

		obj[0] = x[0];
		obj[1] = g * (1 - std::pow(obj[0] / g, 0.5));
	}

	void IMMOEA_F9::updateOptima() {
		m_optima.reset(new Optima<>);
		sampleParetoSols(m_num_reference_points);
		//#ifdef OFEC_DEMO
		//		sampleParetoSets(10000);
		//#endif //
	}

	void IMMOEA_F9::sampleParetoSols(size_t sample_num) {
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
					temp_var[j] = std::pow(temp_var[0], 1. + m_beta * (j + 1.) / m_number_variables);
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