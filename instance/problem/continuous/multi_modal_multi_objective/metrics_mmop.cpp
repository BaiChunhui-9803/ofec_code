#include "metrics_mmop.h"

namespace ofec {
	void MetricsMMOP::updateCandidates(const SolutionBase& sol, std::list<std::unique_ptr<SolutionBase>>& candidates) const {
		bool dominated = false;
		Dominance result;
		for (auto it = candidates.begin(); it != candidates.end();) {
			result = sol.compare(**it, m_optimize_mode);
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
			candidates.emplace_back(createSolution(sol));
	}
}