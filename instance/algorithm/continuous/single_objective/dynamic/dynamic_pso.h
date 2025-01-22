#ifndef OFEC_DYNAMIC_PSO
#define OFEC_DYNAMIC_PSO

#include "../../../../../core/problem/solution.h"

namespace ofec {
	class DynamicPSO {
	public:
		virtual std::vector<bool> getPopHiberState() const = 0;
		const std::vector<Solution<>> &getGbestLocation() const { return m_gbestslocation; }

	protected:
		std::vector<Solution<>> m_gbestslocation;
	};
}

#endif // !OFEC_DYNAMIC_PSO
