#ifndef OFEC_COP_BASE_H
#define OFEC_COP_BASE_H

#include "../../../../../core/problem/continuous/function.h"
#include <functional>

namespace ofec {
	class CopBase : public Function {
	public:
		void updateCandidates(const SolBase& sol, std::list<std::unique_ptr<SolBase>>& candidates) const override;
		size_t numOptimaFound(const std::list<std::unique_ptr<SolBase>>& candidates) const override;
		bool constraintViolated(const SolBase& s) const override;
	};
}

#endif // !OFEC_METRICS_GOP_H
