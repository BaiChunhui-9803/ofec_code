#ifndef OFEC_JDE_POP_H
#define OFEC_JDE_POP_H

#include "../../../../template/classic/differential_evolution/population.h"

namespace ofec {
	class PopJDE : public PopulationDE<> {
	public:
		PopJDE(size_t size_pop, Environment *env);
		void initialize(Environment *env, Random *rnd) override;
		int evolve(Environment *env, Random *rnd) override;
		Real F(size_t i) const { return mv_F[i]; }
		Real CR(size_t i) const { return mv_CR[i]; }
	protected:
		virtual void updateParams(Random *rnd);
	protected:
		Real m_t1 = 0, m_t2 = 0, m_Fl = 0, m_Fu = 0;
		std::vector<Real> mv_F, mv_CR, mv_tF, mv_tCR;
	};
}


#endif // OFEC_JDE_POP_H
