#include "gl_cont_pop.h"
#include "gl_adaptor_cont.h"
#include "../../../../../../core/problem/continuous/continuous.h"

namespace ofec {
	void PopContGL::resize(size_t size_pop, Environment *env) {
		PopGL<Solution<>>::resize(size_pop, env, CAST_CONOP(env->problem())->numberVariables());
	}

	void PopContGL::initialize(Environment *env, Random *rnd) {
		m_adaptor.reset(new AdaptorContGL(m_alpha, CAST_CONOP(env->problem())->numberVariables(), m_individuals.size()));
		Population<Solution<>>::initialize(env, rnd);
		initializeCurpop();
	}

	void PopContGL::initializeCurpop() {
		for (int i = 0; i < this->size(); i++) {
			m_offspring.push_back(*this->m_individuals[i]);
		}
	}

	int PopContGL::evolve(Environment *env, Random *rnd) {
		dynamic_cast<AdaptorContGL *>(m_adaptor.get())->updateStep(env, *this);
		m_adaptor->createSolution(env, rnd, *this, m_offspring);
		int rf = update(env);
		//updateBest(env);
		updateMemory(env);
		m_iteration++;
		return rf;
	}
}
