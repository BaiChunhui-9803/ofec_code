#include "glt1.h"

namespace ofec {
	void GLT1::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real yj, g = 0;
		for (size_t j = 2; j <= m_number_variables; j++) {
			yj = x[j - 1] - sin(2 * OFEC_PI*x[0] + j * OFEC_PI / m_number_variables);
			g += yj * yj;
		}	
		
		obj[0] = (1 + g)*x[0];
		obj[1] = (1 + g)*(2 - x[0] - sign(cos(2 * OFEC_PI*x[0])));
	}

	void GLT1::updateOptima() {
		m_optima.reset(new Optima<>);
		sampleParetoSols(m_num_reference_points);
//#ifdef OFEC_DEMO
//		sampleParetoSets(m_num_reference_points);
//#endif //
	}

	void GLT1::sampleParetoSols(size_t sample_num) {
		size_t num = sample_num;
		std::vector<VariableVector<Real>> all_vars;
		std::vector<std::vector<Real>> all_objs;
		for (size_t i = 0; i < num; i++) {
			VariableVector<Real> temp_var(m_number_variables);
			for (size_t j = 0; j < m_number_variables; j++) {
				if (j == 0)
					temp_var[j] = (Real)i / (num - 1.);
				else
					temp_var[j] = std::sin(2 * OFEC_PI * temp_var[0] + (1 + j) * OFEC_PI / m_number_variables);
			}
			all_vars.emplace_back(temp_var);
		}
		//first, nd-sort, then add into optima
		for (size_t i = 0; i < all_vars.size(); ++i) {
			std::vector<Real> temp_obj(m_number_objectives);
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