#ifndef ALGORITHM_UNCERTIANTY_H
#define ALGORITHM_UNCERTIANTY_H

#include "../../../../../core/algorithm/algorithm.h"
#include "pop_uncertianty.h"
#include "evaluation_strategy.h"

namespace ofec {
//#define CAST_AlgUncrty(alg) dynamic_cast<AlgorighmUncertianty*>(alg)


	template <template<class> class TEvaluationStrategy,
		class TIndi = Solution<>>
	class AlgorighmUncertianty:  virtual public Algorithm {
	public:
		using SolutionType = typename TIndi;
		//using EvaluationStrategyType = typename TEvaluationStrategy;
	protected:
		//PopUncertianty<SolutionType> m_pop;
		std::shared_ptr<EvaluationStrategyBase<SolutionType>>
			m_eval_strategy;
		virtual void clearPop(int idx = 0) {
			m_eval_strategy->clearMemory(idx);
		}
		virtual void setPop(int popIdx = 0) {
			getPop(popIdx).setId(popIdx);
			getPop(popIdx).setEvalStrategy(m_eval_strategy);
		}

		//virtual void evolvePop(int popIdx = 0) {
		//	getPop(popIdx).generateOffstrpings();
		//	calculateFitness(popIdx);               
		//	getPop(popIdx).selectSurvivors(m_eval_strategy);
		//	getPop(popIdx).updateInfo();
		//}




		virtual void initialize_()override {
			Algorithm::initialize_();
			m_eval_strategy.reset(new EvaluationStrategyType<SolutionType>());
			m_eval_strategy->initialize(m_problem.get(), this);
		}

	public:
		//virtual const std::vector<TIndi>* 
		//	getIndis(int idx = 0, bool flag_parents=true) {
		//	return nullptr;
		//}
		////virtual const std::vector<std::shared_ptr<TIndi>>* 
		////	getIndisSP(int idx, bool flag_parents=true) {
		////	return nullptr;
		////}
		//virtual const std::vector<std::unique_ptr<TIndi>>* 
		//	getIndisUP(int idx, bool flag_parents=true) {
		//	return nullptr;
		//}
	};
}

#endif // ! ALGORITHM_UNCERTIANTY_H
