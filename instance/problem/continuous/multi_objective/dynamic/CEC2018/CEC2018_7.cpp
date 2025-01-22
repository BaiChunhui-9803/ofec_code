#include "CEC2018_7.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void CEC2018_7::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the variable is in (0,1)
			if (i < 1)
				m_domain.setRange(1., 4., i);
			else
				m_domain.setRange(0., 1., i);
		}
		generateAdLoadPF();
	}

	void CEC2018_7::generateAdLoadPF() {
		int num = 2000;
		Real t = get_time();
		std::vector<std::vector<Real>*> point(num);
		Matrix o_point(num, 2);
		for (int i = 0; i < num; i++) {
			o_point[i][0] = (1. + t) / 4+ (1+t-(1.+t)/4) * i / (num - 1);
			o_point[i][1] = 1./o_point[i][0];
			point[i] = &o_point[i].vect();
			m_optima->appendObj(o_point[i].vect());
		}
		m_optima->setVariableGiven(true);//given optima
	}

	 int CEC2018_7::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void CEC2018_7::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			set_updated_state(false);
			generateAdLoadPF();
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		Real gt = 1; //g function
		Real alpha_t=5*cos(0.5*OFEC_PI*t);//change factor
		for (size_t i = 1; i < m_number_variables; ++i)
		{
			gt += pow(x[i] - 1./(1+exp(alpha_t*(x[0]-2.5))), 2);
		}

		obj[0] = gt*(1+t)/x[0];
		obj[1] = gt * x[0]/(1+t);
	}
}