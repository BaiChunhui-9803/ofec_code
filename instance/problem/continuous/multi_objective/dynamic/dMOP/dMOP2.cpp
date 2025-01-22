#include "dMOP2.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void dMOP2::initialize_() {//all variables are in [0,1], recommend the number of x is 10
		DMOPs::initialize_();
		generateAdLoadPF();
	}

	//f1=x1, f2=g*h=1-f1^Ht
	void dMOP2::generateAdLoadPF() {
		int num = 2000;
		Real t = get_time();
		Real Ht = 1.25 + 0.75*std::sin(0.5*OFEC_PI*t);
		Matrix o_point(num, 2);
		for (int i = 0; i < num; i++) {
			o_point[i][0] = 0 + i * 1. / (num - 1);
			o_point[i][1] = 1 - std::pow(o_point[i][0], Ht);
			m_optima->appendObj(o_point[i].vect());
		}
	}

	int dMOP2::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
		auto eval = alg->evaluations();
		set_effective_eval(eval);
		if (alg != nullptr && eval % get_change_fre() == 0 && eval > 0) {
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
			return kChangeNextEval;
		}
		else if (alg != nullptr && (eval + 1) % get_change_fre() == 0 && eval > 0) {
			//m_optima.reset(new Optima<>());
			set_updated_state(true);
			return kNormalEval;
		}
		else {
			return Continuous::updateEvaluationTag(s, alg);
		}
	}

	void dMOP2::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(false);
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		
		Real Ht = 0.75*std::sin(0.5*OFEC_PI*t) + 1.25;
		Real Gt = std::sin(0.5*OFEC_PI*t);
		Real gt = 1;
		for (size_t i = 1; i < m_number_variables; ++i) {
			gt += std::pow(x[i] - Gt, 2);
		}
		for (size_t m = 0; m < m_number_objectives; ++m) {
			if (m < 1)
				obj[m] = x[0];
			else
				obj[m] = gt * (1 - std::pow(obj[0] / gt, Ht));
		}
	}
}