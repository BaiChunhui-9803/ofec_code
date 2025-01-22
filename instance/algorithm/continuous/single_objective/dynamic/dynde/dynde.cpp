#include "dynde.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void DynDE::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;;
		m_r_excl = v.get<Real>("exlRadius size");
		m_pop_size = v.get<int>("population size");
		m_subpop_size = v.get<int>("subpopulation size");
		m_multi_pop.clear();
	}

	void DynDE::run_() {
		m_multi_pop.resize(m_pop_size / m_subpop_size, m_subpop_size, m_problem.get());
		for (size_t i = 0; i < m_multi_pop.size(); ++i) {
			m_multi_pop[i].initialize(m_problem.get(), m_random.get());
			m_multi_pop[i].evaluate(m_problem.get(), this);
		}
#ifdef OFEC_DEMO
		updateBuffer();
#endif
		int tag = kNormalEval;
		while (!terminating()) {
			for (size_t i = 0; i < m_multi_pop.size(); ++i) {
				if (!m_multi_pop[i].getFlag())
					tag = m_multi_pop[i].evolve(m_problem.get(), this, m_random.get());
				//handle_evaluation_tag(tag);
			}
			if (tag == kTerminate)
				break;
			//exclusion
			exclusion_check();
			for (size_t i = 0; i < m_multi_pop.size(); ++i) {
				if (m_multi_pop[i].getFlag()) {
					m_multi_pop[i].initialize(m_problem.get(), m_random.get());
					m_multi_pop[i].setFlag(false);
				}
			}
#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
	}

#ifdef OFEC_DEMO
	void DynDE::updateBuffer() {
		m_solution.resize(m_multi_pop.size());
		for (size_t i = 0; i < m_multi_pop.size(); i++) {
			for (size_t j = 0; j < m_multi_pop[i].size(); ++j)
				m_solution[i].emplace_back(&m_multi_pop[i][j]);
		}
		ofec_demo::g_buffer->appendAlgBuffer(this);
	}
#endif

	void DynDE::exclusion_check() {
		for (size_t i = 0; i < m_multi_pop.size(); ++i) {
			for (size_t j = i + 1; j < m_multi_pop.size(); ++j) {
				if (m_multi_pop[i].getFlag() == false && m_multi_pop[j].getFlag() == false 
					&& m_multi_pop[i].best().front()->variableDistance(*m_multi_pop[j].best().front(), m_problem.get()) < m_r_excl)
				{
					if (m_multi_pop[i].best().front()->dominate(*m_multi_pop[j].best().front(), m_problem.get())) {
						m_multi_pop[j].setFlag(true);
					}
					else {
						m_multi_pop[i].setFlag(true);
					}
				}
			}
		}
	}

	void DynDE::record() {
		//if (CONTINUOUS_CAST->has_tag(problem_tag::DOP)) {
		//	// ******* Dynamic Optimization ***********
		//}
		//else if (CONTINUOUS_CAST->has_tag(problem_tag::MMOP)) {
		//	// ******* Multi-Modal Optimization *******
		//	size_t evals = CONTINUOUS_CAST->evaluations();
		//	size_t num_opt_found = CONTINUOUS_CAST->num_optima_found();
		//	size_t num_opt_known = CONTINUOUS_CAST->optima().number_objective();
		//	Real peak_ratio = (Real)num_opt_found / (Real)num_opt_known;
		//	Real success_rate = CONTINUOUS_CAST->solved() ? 1 : 0;
		//	measure::get_measure()->record(global::ms_global.get(), evals, peak_ratio, success_rate);
		//}
		//else if (CONTINUOUS_CAST->has_tag(problem_tag::GOP)) {
		//	// ******* Global Optimization ************
		//	size_t evals = CONTINUOUS_CAST->evaluations();
		//	Real err = std::fabs(problem::get_sofar_best<Solution<>>(0)->objective(0) - CONTINUOUS_CAST->optima().objective(0).at(0));
		//	measure::get_measure()->record(global::ms_global.get(), evals, err);
		//}
	}
}
