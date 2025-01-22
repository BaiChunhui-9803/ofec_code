#ifndef OFEC_METRICS_MMOP_H
#define OFEC_METRICS_MMOP_H

#include "../../../../core/problem/continuous/continuous.h"

namespace ofec {
	class MetricsMMOP : virtual public Problem {
	public:
		void updateCandidates(const SolutionBase& sol, std::list<std::unique_ptr<SolutionBase>>& candidates) const override;
	};
}
#endif // !OFEC_METRICS_MMOP_H