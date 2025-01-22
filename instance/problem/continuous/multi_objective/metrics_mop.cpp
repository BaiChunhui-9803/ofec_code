#include "metrics_mop.h"

namespace ofec {
	void MetricsMOP::updateCandidates(const SolutionBase& sol, std::list<std::unique_ptr<SolutionBase>> &candidates) const {
		if (candidates.empty())
			candidates.emplace_back(createSolution(sol));
		else if (sol.dominate(*candidates.front(), this))
			candidates.front().reset(createSolution(sol));
		/*bool dominated = false;
		Dominance result;
		for (auto it = candidates.begin(); it != candidates.end();) {
			result = sol.compare(**it, m_optimize_mode);
			sol.compare();
			if (result == Dominance::kDominated) {
				dominated = true;
				break;
			}
			if (result == Dominance::kDominant) {
				it = candidates.erase(it);
			}
			else
				it++;
		}
		if (!dominated)
			candidates.emplace_back(createSolution(sol));*/
	}
}