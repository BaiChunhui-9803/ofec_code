#ifndef OFEC_NBN_EAX_POP_V3_H
#define OFEC_NBN_EAX_POP_V3_H

#include "../../../../../core/algorithm/population.h"
#include "../../../../../core/algorithm/individual.h"

#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"
#include "../environment.h"
#include "../../nbn_alg_com/gl_calculator.h"
#include "../eax_nbn_2/nbn_eax_pop_v2.h"

namespace ofec {

	


	class PopNBN_EAX_V3 : public PopNBN_EAX_V2 {




	public:
		PopNBN_EAX_V3() :PopNBN_EAX_V2(){}
		PopNBN_EAX_V3(size_t size_pop, Problem* pro) : PopNBN_EAX_V2(size_pop, pro){};

		virtual ~PopNBN_EAX_V3() = default;
		int evolve(Problem* pro, Algorithm* alg, Random* rnd) override;
		virtual void initialize(Problem* pro, Random* rnd) override;
		//	int evaluate(Problem* pro, Algorithm* alg) override;
	protected:

	public:

	};
}

#endif // !OFEC_PopGL_NBN_COM_ALG_H