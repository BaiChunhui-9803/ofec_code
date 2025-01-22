#ifndef OFEC_METRICS_DMOP_H
#define OFEC_METRICS_DMOP_H

#include "../../../../../core/problem/continuous/continuous.h"

namespace ofec {
	class MetricsDMOP : virtual public Continuous {
	public:
		void updateCandidates(const SolutionBase& sol, std::list<std::unique_ptr<SolutionBase>>& candidates) const override;
		size_t numOptimaFound(const std::list<std::unique_ptr<SolutionBase>>& candidates) const override;
	};
}
#endif // !OFEC_METRICS_DMOP_H
