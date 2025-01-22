#ifndef OFEC_SIMULATED_ANNEALING_H
#define OFEC_SIMULATED_ANNEALING_H

#include "../../../../../core/algorithm/algorithm.h"
#include "../../../../../core/problem/solution.h"

/* Simulated Annealing */

namespace ofec {
	class BaseSA : virtual public Algorithm {
		OFEC_ABSTRACT_INSTANCE(BaseSA)
	protected:
		void addInputParameters() {}
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

		virtual void pickNeighbour(Environment *env) = 0;
		virtual Real temperature(Real r);
		virtual Real criterion(Real t, Environment *env);
		virtual void replaceCurrent() = 0;

	protected:
		std::unique_ptr<SolutionBase> m_cur, m_neighbour;

	public:
		void record();
	};

	template <typename TVariable>
	class SA : public BaseSA {
		OFEC_ABSTRACT_INSTANCE(SA)
	protected:
		using TypeVar = TVariable;
		void addInputParameters() {}
		virtual void replaceCurrent() override {
			dynamic_cast<Solution<TVariable>&>(*m_cur) = dynamic_cast<Solution<TVariable>&>(*m_neighbour);
		}
	};
}

#endif // !OFEC_SIMULATED_ANNEALING_H

