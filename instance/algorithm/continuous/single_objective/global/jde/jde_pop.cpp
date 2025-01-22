#include "jde_pop.h"

namespace ofec {
	PopJDE::PopJDE(size_t size_pop, Environment *env) : 
		PopulationDE(size_pop, env), 
		mv_F(size_pop), 
		mv_CR(size_pop),
		mv_tF(size_pop),
		mv_tCR(size_pop),
		m_t1(0.1),
		m_t2(0.1),
		m_Fl(0.1),
		m_Fu(0.9) {}

	void PopJDE::initialize(Environment *env, Random *rnd) {
		PopulationDE::initialize(env, rnd);
		for (size_t i = 0; i < m_individuals.size(); i++) {
			mv_F[i] = 0.5;
			mv_CR[i] = 0.9;
		}
	}

	int PopJDE::evolve(Environment *env, Random *rnd) {
		updateParams(rnd);
		int tag = kNormalEval;
		std::vector<size_t> ridx(3);
		for (size_t i = 0; i < size(); ++i) {
			select(i, 3, ridx, rnd);
			m_individuals[i]->mutate(mv_tF[i], m_individuals[ridx[0]].get(), m_individuals[ridx[1]].get(), m_individuals[ridx[2]].get(), env);
			m_individuals[i]->recombine(mv_tCR[i], m_recombine_strategy, rnd, env);
			tag = m_individuals[i]->select(env);
			if (!(tag & kNormalEval)) return tag;
			if (m_individuals[i]->isImproved()) {
				mv_F[i] = mv_tF[i];
				mv_CR[i] = mv_tCR[i];
			}
		}
		if (tag & kNormalEval) {
			++m_iteration;
		}
		return tag;
	}

	void PopJDE::updateParams(Random *rnd) {
		std::vector<Real> rand(4);
		for (size_t i = 0; i < size(); ++i) {
			for (size_t j = 0; j < 4; ++j) {
				rand[j] = rnd->uniform.next();
			}
			mv_tF[i] = rand[1] < m_t1 ? m_Fl + rand[0] * m_Fu : mv_F[i];
			mv_tCR[i] = rand[3] < m_t2 ? rand[2] : mv_CR[i];
		}
	}
}


