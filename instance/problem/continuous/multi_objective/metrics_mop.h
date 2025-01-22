#ifndef OFEC_METRICS_MOP_H
#define OFEC_METRICS_MOP_H

#include "../../../../core/problem/problem.h"

namespace ofec {
	class MetricsMOP : virtual public Problem {
	public:
		void updateCandidates(const SolutionBase& sol, std::list<std::unique_ptr<SolutionBase>> &candidates) const override;
	};
}
#endif // !OFEC_METRICS_MOP_H