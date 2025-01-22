#include "particle07.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../core/environment/environment.h"


namespace ofec {
	void Particle07::initVelocity(Environment *env, Random *rnd) {
		for (size_t j = 0; j < variable().size(); j++) {
			auto &range = CAST_CONOP(env->problem())->range(j);
			m_vel[j] = (rnd->uniform.nextNonStd(range.first, range.second) - variable()[j]) / 2;
		}
	}

	void Particle07::nextVelocity(const Solution<> *lbest, Real w, Real c1, Real c2, Random *rnd) {
		if (lbest != &m_pbest) {
			for (size_t j = 0; j < variable().size(); j++)
				m_vel[j] = w * m_vel[j]
				+ c1 * rnd->uniform.next() * (m_pbest.variable()[j] - variable()[j])
				+ c2 * rnd->uniform.next() * (lbest->variable()[j] - variable()[j]);
		}
		else {
			for (size_t j = 0; j < variable().size(); j++)
				m_vel[j] = w * m_vel[j] + c1 * rnd->uniform.next() * (m_pbest.variable()[j] - variable()[j]);
		}
	}
}