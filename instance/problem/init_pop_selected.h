#ifndef OFEC_INIT_POP_SELECTED_H
#define OFEC_INIT_POP_SELECTED_H

#include "../../core/problem/problem.h"

namespace ofec {
	class InitPopSelectedBase : virtual public Problem {
	public:
		void initializeVariables(VariableBase &x, Random *rnd) const override;
		void addSolInitPop(const SolutionBase &s);
		void clearSolsInitPop();
	protected:
		std::list<std::unique_ptr<const SolutionBase>> m_sols_init_pop; // solutions for intializing population
	};

	template<typename TPro>
	class InitPopSelected : public TPro, public InitPopSelectedBase {
		OFEC_CONCRETE_INSTANCE(InitPopSelected<TPro>)
	protected:
		void addInputParameters() {}

	public:
		void initializeVariables(VariableBase &x, Random *rnd) const override {
			if (m_sols_init_pop.empty()) {
				TPro::initializeVariables(x, rnd);
			}
			else {
				InitPopSelectedBase::initializeVariables(x, rnd);
			}
		}
	};
}

#endif // !OFEC_INIT_POP_SELECTED_H