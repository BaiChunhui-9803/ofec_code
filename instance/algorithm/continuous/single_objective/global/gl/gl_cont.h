/********* Begin Register Information **********
{
	"name": "GLC",
	"identifier": "ContGL",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#ifndef OFEC_GL_CONT_H
#define OFEC_GL_CONT_H

#include "gl_cont_pop.h"
#include "../../../../../../core/problem/solution.h"

namespace ofec {
	class ContGL : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(ContGL)
	public:
		using UpdateScheme = PopContGL::UpdateScheme;
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

	protected:
		size_t m_pop_size;
		UpdateScheme m_update_scheme;
		Real m_alpha, m_beta, m_gamma;
		PopContGL m_pop;
	};
}


#endif // !OFEC_GL_CONT_H

