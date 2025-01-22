#include "zdt3.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../utility/nondominated_sorting/filter_sort.h"

namespace ofec {
	void ZDT3::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real g = 0;
		for (size_t n = 1; n < m_number_variables; n++) {
			g = g + x[n];
		}
		g = 1 + 9 * g / (m_number_variables - 1);
		obj[0] = x[0];
		obj[1] =1+ g*(1 - std::sqrt(x[0] / g) - x[0] * sin(10 * x[0] * OFEC_PI) / g);
	}

	void ZDT3::updateOptima() {
		m_optima.reset(new Optima<>);
		//loadParetoFront(10000);
		sampleParetoSols(m_num_reference_points);
		//#ifdef OFEC_DEMO
		//		sampleParetoSets(10000);
		//#endif //
	}

	void ZDT3::sampleParetoSols(size_t sample_num) {
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
			evaluateObjective(temp_var.data(), temp_obj);
			all_vars.emplace_back(temp_var);
			all_objs.emplace_back(temp_obj);
		}
		//first nd-sort, then add into optima
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
