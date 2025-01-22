#ifndef OFEC_HGHE_H
#define OFEC_HGHE_H

#include "../../../../../core/algorithm/algorithm.h"
#include <memory>
#include "../../../../../core/problem/solution.h"
#include <list>

namespace ofec {
	template <typename TVariable>
	class HGHE : virtual public Algorithm {
		OFEC_ABSTRACT_INSTANCE(HGHE)
	protected:
		enum class TaskEval { kExplore, kExploit };
		using SolutionType = Solution<TVariable>;

		std::list<std::unique_ptr<const SolutionType>> m_his_sols;
		std::list<const SolutionType*> m_his_sols_explore, m_his_sols_exploit;

		virtual const SolutionType* archiveSolution(const SolutionType &sol, TaskEval task, Environment *env) {
			m_his_sols.emplace_back(new SolutionType(sol));
			if (task == TaskEval::kExplore)
				m_his_sols_explore.push_back(m_his_sols.back().get());
			else if (task == TaskEval::kExploit)
				m_his_sols_exploit.push_back(m_his_sols.back().get());
			return m_his_sols.back().get();
		}

		void addInputParameters() {}

		void initialize_(Environment *env) override {
			Algorithm::initialize_(env);
			m_his_sols.clear();
			m_his_sols_explore.clear();
			m_his_sols_exploit.clear();
		}

	public:
		const std::list<std::unique_ptr<const SolutionType>>& hisSols() const { return m_his_sols; }
		const std::list<const SolutionType*>& hisSolsExplore() const { return m_his_sols_explore; }
		const std::list<const SolutionType*>& hisSolsExploit() const { return m_his_sols_exploit; }
	};
}

#endif // !OFEC_HGHE_H
