#include "metrics_dynamic.h"
#include "../../../../problem/continuous/single_objective/dynamic/uncertianty_continuous.h"

namespace ofec {
	void MetricsDynamicConOEA::updateCandidates(const SolutionBase& sol){
		Algorithm::updateCandidates(sol);
		if (m_problem->hasTag(ProblemTag::kDOP)) {
			calculateOfflineError(m_candidates.back().get());
			calculateBBCError(m_candidates.back().get());
		}
	}

	void MetricsDynamicConOEA::initialize_(){
		Algorithm::initialize_();
		m_offline_error = 0.;
		m_best_before_change_error = 0.;
		m_num_envirs = 0;
	}

	
	void MetricsDynamicConOEA::calculateOfflineError(const SolutionBase * best_sol) {
		if (CAST_CONOP(m_problem.get())->numberObjectives() == 1) {
			m_offline_error = ((m_problem->optimaBase()->objective(0)[0] - best_sol->objective()[0]) + (m_evaluations - 1) * m_offline_error) / m_evaluations;
		}
	}

	void MetricsDynamicConOEA::calculateBBCError(const SolutionBase* best_sol){
		if (CAST_CONOP(m_problem.get())->numberObjectives() == 1){
			if (m_evaluations % GET_DOP(m_problem.get())->getFrequency() == 0){
				m_num_envirs++;
				m_best_before_change_error = ((m_problem.get()->optimaBase()->objective(0)[0] - best_sol->objective()[0]) + (m_num_envirs - 1) * m_best_before_change_error) / m_num_envirs;
			}
		}
	}
}