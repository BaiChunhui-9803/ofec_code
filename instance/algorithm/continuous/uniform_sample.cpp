#include "uniform_sample.h"

#ifdef OFEC_PLAYBACK
#include <core/global.h>
#endif

namespace ofec {
	void UniformSample::run_() {
#ifdef OFEC_PLAYBACK
		m_solution.resize(1);
#endif
		size_t num_vars = m_problem->numberVariables();
		size_t number_objectives = m_problem->numberObjectives();
		size_t num_cons = m_problem->numberConstraints();
		while (!terminating()) {
			for (size_t i = 0; i < 10; ++i) {
				m_samples.emplace_back(number_objectives, num_cons, num_vars);
				m_samples.back().initialize(m_problem.get(), m_random.get());
				m_samples.back().evaluate(m_problem.get(), this);
#ifdef OFEC_PLAYBACK
				m_solution[0].push_back(&m_samples.back());
#endif
			}
#ifdef OFEC_PLAYBACK
			ofec_playback::g_buffer->appendAlgBuffer(this);
#endif
		}
	}

	void UniformSample::record() {

	}
}