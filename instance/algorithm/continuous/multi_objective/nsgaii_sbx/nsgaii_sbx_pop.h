#ifndef OFEC_NSGAII_SBX_POP_H
#define OFEC_NSGAII_SBX_POP_H

#include "../../template/multi_objective/nsgaii/nsgaii.h"
#include "../../template/classic/ga/sbx_pop.h"

namespace ofec {
	class PopNSGAII : public PopSBX<>, public NSGAII {
	public:
		PopNSGAII(size_t size_pop, Problem *pro);
		int evolve(Problem *pro, Algorithm *alg, Random *rnd) override;
	protected:
		Population<Solution<>> m_pop_combined;  // combination of parent and children
	};
}

#endif // !OFEC_NSGAII_SBX_POP_H

