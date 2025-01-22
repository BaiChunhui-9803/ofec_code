#ifndef OFEC_GL_CONT_POP_H
#define OFEC_GL_CONT_POP_H

#include "../../../../template/framework/gl/gl_pop.h"
#include "../../../../../../core/problem/solution.h"

namespace ofec {
	class PopContGL : public PopGL<Solution<>> {
	public:
		void resize(size_t size_pop, Environment *env);
		void initialize(Environment *env, Random *rnd) override;
		int evolve(Environment *env, Random *rnd) override;
	protected:
		void initializeCurpop();
	};
}


#endif // !OFEC_GL_CONT_POP_H

