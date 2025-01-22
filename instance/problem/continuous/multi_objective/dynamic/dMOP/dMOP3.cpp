#include "dMOP3.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void dMOP3::initialize_() {//all variables are in [0,1], recommend the number of x is 10
		DMOPs::initialize_();
		generateAdLoadPF();
	}

	//f1=x1, f2=g*h=1-f1^(1/2)
	void dMOP3::generateAdLoadPF() {
		int num = 2000;
		Matrix o_point(num, 2);
		for (int i = 0; i < num; i++) {
			o_point[i][0] = 0 + i * 1. / (num - 1);
			o_point[i][1] = 1 - std::pow(o_point[i][0], (Real)0.5);
			m_optima->appendObj(o_point[i].vect());
		}
	}

	int dMOP3::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void dMOP3::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {//防止不计数评价重复更新问题和重复采样PF
			size_t temp = std::floor(m_number_variables * m_random->uniform.next());//choose one dimension as the first objective randomly
			set_r(temp);
			set_updated_state(false);
			generateAdLoadPF();
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}

		Real Ht = 0.75 * std::sin(0.5 * OFEC_PI * t) + 1.25;
		Real Gt = std::sin(0.5 * OFEC_PI * t);
		Real gt = 1;
		for (size_t i = 0; i <m_number_variables; ++i) {
			if (i != get_r())
				gt += pow(x[i] - Gt, 2);
		}
		obj[0] = x[get_r()];
		obj[1] = gt * (1 - sqrt(obj[0] / gt));
	}
}