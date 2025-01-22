#include"evaluation_dynamic_environment.h"

#ifdef OFEC_PLAYBACK
#include <playback/global.h>
#include <statistics/global.h>
#endif // OFEC_PLAYBACK

#include "dynamic_problem.h"


namespace ofec {
	void EvaluationDynamicEnvironment::addInputParameters() {
		m_input_parameters.add("changeFre", new RangedInt(m_frequency, 500, 100000, 1000));
	}

	int EvaluationDynamicEnvironment::evaluate_(SolutionBase &sol) {
		m_algorithm->increaseEvaluations();
		sol.timeEvaluate() = m_algorithm->evaluations();
		int tag = updateEvaluationTag(sol, m_algorithm.get());
		m_problem->evaluate(sol.variableBase(), sol.objective(), sol.constraint());
		handleEvaluateTag(tag);
		return tag;
	}

	int EvaluationDynamicEnvironment::updateEvaluationTag(SolutionBase &s, Algorithm *alg) {
		int rf = kNormalEval;
		if (alg!=nullptr) {
			if ((alg->evaluations()) % (m_frequency) == 0) {
				rf |= kChangeNextEval;
			}
			else if (alg->evaluations() != 1 &&
				(alg->evaluations() - 1) % m_frequency == 0) {
				rf |= kChangeCurEval;
				if (m_flag_objective_memory_change) {
					rf |= kChangeObjectiveMemory;
				}
				if (m_flag_variable_memory_change) {
					rf |= kChangeVariableMemory;
				}
			}
		}
		return rf;
	}

	void EvaluationDynamicEnvironment::handleEvaluateTag(int tag) {
		if (tag & kChangeNextEval) {
			if (auto dp = dynamic_cast<DynamicProblem *>(m_problem.get())) {
				dp->change(m_random.get());
			}

#ifdef OFEC_PLAYBACK
			ofec_playback::g_buffer_manager->handleAlgorithmDatumUpdated(this);
#endif // OFEC_PLAYBACK
#ifdef OFEC_STATISTICS
		//	ofec_statistics::g_record_task->handleAlgorithmDatumUpdated(this);
#endif // OFEC_STATISTICS

		}
	}
}
