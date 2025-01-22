#include "CEC2018_2.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void CEC2018_2::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the variable is in (0,1)
			m_domain.setRange(0., 1., i);
		}
		generateAdLoadPF();
	}

	void CEC2018_2::generateAdLoadPF() {
		int num = 2000;
		Matrix o_point(num, 2);
		for (int i = 0; i < num; i++) {
			o_point[i][0] = 0 + i * 1. / (num - 1);
			o_point[i][1] = 1 - pow(o_point[i][0], (Real)0.5);
			m_optima->appendObj(o_point[i].vect());
		}
		m_optima->setVariableGiven(true);
	}

	 int CEC2018_2::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
		auto eval = alg->evaluations();
		set_effective_eval(eval);
		if (alg != nullptr && eval % get_change_fre() == 0 && eval > 0) {
			size_t temp = (size_t)std::ceil(m_number_variables * m_random->uniform.next());//choose one dimension as the first objective randomly
			set_r(temp);
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

	void CEC2018_2::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			set_updated_state(false);
			generateAdLoadPF();
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}

		Real Gt = fabs(sin(0.5*OFEC_PI*t));
		set_r(1 + floor(Gt*(m_number_variables-1)));
		Real gt = 1;
		for (size_t i = 0; i < m_number_variables; ++i) {
			if (i != get_r() - 1)
				gt += pow(x[i] - Gt, 2);
		}
		obj[0] = x[get_r() - 1];
		obj[1] = gt * (1 - sqrt(obj[0] / gt));
	}
}