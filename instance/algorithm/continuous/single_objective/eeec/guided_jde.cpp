#include "guided_jde.h"
#include "../../../../problem/continuous/free_peaks/free_peaks.h"

namespace ofec {
	void GuidedJDE::addInputParameters() {}

	void GuidedJDE::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		m_pop.reset(new PopGuidedJDE(m_pop_size, env));
	}

	PopGuidedJDE::PopGuidedJDE(size_t size_pop, Environment* env) :
		PopJDE(size_pop, env) {}

	void PopGuidedJDE::initialize(Environment *env, Random *rnd) {
		PopulationDE::initialize(env, rnd);
		auto kd_tree = CAST_FPs(env->problem())->subspaceTree().tree.get();
		for (size_t i = 0; i < size(); ++i) {
			std::string peak_name = kd_tree->getRegionName(
				m_individuals[i]->variable().vector()
			);
			if (peak_name == "a1") {
				do {
					mv_F[i] = rnd->normal.nextNonStd(0.5, 1);
				} while (mv_F[i] < 0 || mv_F[i] > 1);
				do {
					mv_CR[i] = rnd->normal.nextNonStd(0.9, 1.8);
				} while (mv_CR[i] < 0 || mv_CR[i] > 1);
			}
			else {
				mv_F[i] = 0.5;
				mv_CR[i] = 0.9;
			}
		}
	}
}
