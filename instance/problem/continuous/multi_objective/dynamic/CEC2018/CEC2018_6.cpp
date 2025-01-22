#include "CEC2018_6.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void CEC2018_6::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the variable is in (0,1)
			if (i < 1)
				m_domain.setRange(0., 1., i);
			else
				m_domain.setRange(-1., 1., i);
		}
		generateAdLoadPF();
	}

	void CEC2018_6::generateAdLoadPF() {
		int num = 2000;
		Real t = get_time();
		Real Gt = sin(0.5*OFEC_PI*t);//change factor
		Real At = 0.1;
		Real Wt = 3;
		Real alpha_t = 0.2 + 2.8* fabs(Gt);
		std::vector<std::vector<Real>*> point(num);
		Matrix o_point(num, 2);
		for (int i = 0; i < num; i++) {
			Real x = i * 1. / (num - 1);
			o_point[i][0] = pow((x + At * sin(Wt*OFEC_PI * x)), alpha_t);
			o_point[i][1] = pow((1 - x + At * sin(Wt*OFEC_PI * x)), alpha_t);
			point[i] = &o_point[i].vect();
		}

		std::vector<std::vector<Real>*> &temp = get_nondominated_set(point, OptimizeMode::kMinimize);
		for (int i = 0; i < temp.size(); i++) {
			m_optima->appendObj(*temp[i]);
		}
		m_optima->setVariableGiven(true);//given optima
	}

	 int CEC2018_6::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void CEC2018_6::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			set_updated_state(false);
			generateAdLoadPF();
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		Real gt = 1; //g function
		Real Gt = sin(0.5*OFEC_PI*t);//change factor
		Real At = 0.1;
		Real Wt = 3;
		Real alpha_t = 0.2 + 2.8* fabs(Gt);
		for (size_t i = 1; i < m_number_variables; ++i) {
			gt += fabs(Gt)*pow(x[i] - Gt, 2) - 10 * cos(2 * OFEC_PI*(x[i] - Gt)) + 10;
		}

		obj[0] = gt*pow(x[0] + At * sin(Wt*OFEC_PI*x[0]), alpha_t);
		obj[1] = gt*pow(1 - x[0] + At * sin(Wt*OFEC_PI*x[0]), alpha_t);
	}
}