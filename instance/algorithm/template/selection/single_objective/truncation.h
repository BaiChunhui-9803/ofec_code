#ifndef OFEC_SELECTION_TRUNCATION_H
#define OFEC_SELECTION_TRUNCATION_H

#include "../../../../../core/problem/problem.h"
#include <algorithm>
#include <numeric>

namespace ofec::selection {
	template<typename TPopulation, typename TPopAndOff>
	void truncation(TPopulation &pop, const TPopAndOff &pop_and_off, Environment *env) {
		std::vector<size_t> order(pop_and_off.size());
		std::iota(order.begin(), order.end(), 0);
		if (env->problem()->optimizeMode(0) == OptimizeMode::kMinimize) {
			std::sort(order.begin(), order.end(),
				[&pop_and_off](size_t i, size_t j) {
					return pop_and_off[i].objective(0) < pop_and_off[j].objective(0);
				}
			);
		}
		else {
			std::sort(order.begin(), order.end(),
				[&pop_and_off](size_t i, size_t j) {
					return pop_and_off[i].objective(0) > pop_and_off[j].objective(0);
				}
			);
		}
		int curr_order = 0, num_filled = 0;
		while (num_filled < pop.size()) {
			pop[num_filled++] = pop_and_off[order[curr_order++]];
			if (num_filled == pop.size()) {
				break;
			}
		}
	}
}

#endif // !OFEC_SELECTION_TRUNCATION_H
