#include "amp_pso.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void AMP_PSO::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;;
		m_pop_size = v.get<int>("population size");
		m_weight = v.has("weight") ? v.get<Real>("weight") : 0.7298;
		m_accelerator1 = v.has("accelerator1") ? v.get<Real>("accelerator1") : 1.496;
		m_accelerator2 = v.has("accelerator2") ? v.get<Real>("accelerator2") : 1.496;
		m_multi_pop.reset();
		m_minmax_objective_monitored = true;
	}

	void AMP_PSO::record() {}

	void AMP_PSO::initMultiPop() {
		m_multi_pop.reset(new ContAMP<PopulationType>(m_pop_size));
		m_multi_pop->initialize(m_problem.get(), this, m_random.get());
		for (size_t k = 0; k < m_multi_pop->size(); k++) {
			m_multi_pop->at(k).weight() = m_weight;
			m_multi_pop->at(k).accelerator1() = m_accelerator1;
			m_multi_pop->at(k).accelerator2() = m_accelerator2;
			m_multi_pop->at(k).initPbest(m_problem.get());
		}
#ifdef OFEC_DEMO
		updateBuffer();
#endif
	}

	void AMP_PSO::run_() {
		initMultiPop();
		while (!terminating()) {
			m_multi_pop->evolve(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
	}

#ifdef OFEC_DEMO
	void AMP_PSO::updateBuffer() {
		m_solution.clear();
		m_solution.resize(m_multi_pop->size());
		for (size_t k = 0; k < m_multi_pop->size(); k++) {
			for (size_t i = 0; i < m_multi_pop->at(k).size(); ++i)
				m_solution[k].push_back(&m_multi_pop->at(k)[i]);
		}
		ofec_demo::g_buffer->appendAlgBuffer(this);
	}

	std::vector<bool> AMP_PSO::getPopHiberState() const {
		std::vector<bool> hiber_state;
		return hiber_state;
	}
#endif
}