#ifndef OFEC_HILL_CLIMBING_H
#define OFEC_HILL_CLIMBING_H

#include "../../../../../core/algorithm/algorithm.h"
#include "../../../../../core/problem/solution.h"

/* Hill Climbing */

namespace ofec {
	class BaseHC : virtual public Algorithm {
		OFEC_ABSTRACT_INSTANCE(BaseHC)
	protected:
		void addInputParameters() {}
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

		virtual void pickNeighbour(Environment *env) = 0;
		virtual void replaceCurrent() = 0;

	protected:
		std::unique_ptr<SolutionBase> m_cur, m_neighbour;
	};

	template <typename TVariable>
	class HC : virtual public BaseHC {
		OFEC_ABSTRACT_INSTANCE(HC)
	protected:
		using TypeVar = TVariable;
		void addInputParameters() {}
		virtual void replaceCurrent() override {
			dynamic_cast<Solution<TVariable>&>(*m_cur) = dynamic_cast<Solution<TVariable>&>(*m_neighbour);
		}
	};
}

#endif // !OFEC_HILL_CLIMBING_H

