#include "individual.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../../core/environment/environment.h"

namespace ofec {
	IndEP::IndEP(size_t number_objectives, size_t num_cons, size_t num_vars) :
		Solution<>(number_objectives, num_cons, num_vars), 
		m_eta(num_vars) {}

	void IndEP::initializeEta(Environment *env) {
		for (size_t i = 0; i < m_eta.size(); ++i) {
			m_eta[i] = (CAST_CONOP(env->problem())->range(i).second - CAST_CONOP(env->problem())->range(i).first) / 10;
		}
	}
}