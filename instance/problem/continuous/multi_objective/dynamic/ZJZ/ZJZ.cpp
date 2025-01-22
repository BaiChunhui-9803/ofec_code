#include "ZJZ.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void ZJZ::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the first variable is in [0,1], others are in [-1,2]
			if (i < 1)
				m_domain.setRange(0., 1., i);
			else
				m_domain.setRange(-1., 2., i);
		}
		generateAdLoadPF();
	}

	//f1=x1,f2=g*h=1-f1^Ht, recommend n=10;
	void ZJZ::generateAdLoadPF() {
		int num = 2000;
		Real t = get_time();
		Real temp = 1.5 + std::sin(0.5*OFEC_PI*t);
		Matrix o_point(num, 2);
		for (int i = 0; i < num; i++) {
			o_point[i][0] = 0 + i * 1. / (num-1);
			o_point[i][1] = 1 - std::pow(o_point[i][0], temp);
			m_optima->appendObj(o_point[i].vect());
		}
	}

	int ZJZ::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	//f1=x1,f2=g*h=1-(f1/g)^Ht;
	void ZJZ::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real t = get_time();
		Real gt = 1;
		Real Gt = sin(0.5 * OFEC_PI * t);
		Real Ht = 1.5 + Gt;
		for (size_t i = 1; i < m_number_variables; ++i)
		{
			gt += std::pow(x[i] + Gt - std::pow(x[0], Ht), 2);
		}
		if (time_changed() && t != 0. && get_updated_state()) {
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(false);
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}		
		for (size_t m = 0; m < m_number_objectives; ++m)
		{
			if (m < 1)
				obj[m] = x[m];
			else
				obj[m] = 1 - std::pow(obj[0] / gt, Ht);
		}
	}
}