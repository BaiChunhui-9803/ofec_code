#ifndef OFEC_DE_NRAND_1_POP_H
#define OFEC_DE_NRAND_1_POP_H

#include "../../../../template/classic/differential_evolution/population.h"

namespace ofec {
	class PopDE_nrand_1 : public PopulationDE<> {
	public:
		PopDE_nrand_1() = default;
		PopDE_nrand_1(const PopDE_nrand_1 &rhs) = default;
		PopDE_nrand_1(PopDE_nrand_1 &&rhs) noexcept = default;
		PopDE_nrand_1(size_t size_pop, Environment *env);
		void mutate(int idx, Random *rnd, Environment *env) override;
	protected:
		size_t nearestNeighbour(int idx, Environment *env);
	};
}
#endif // OFEC_DE_NRAND_1_POP_H
