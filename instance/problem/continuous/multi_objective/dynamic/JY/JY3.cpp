#include "JY3.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void JY3::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the first variable is in (0,1), others are in (-1,1)
			if (i < 1)
				m_domain.setRange(0., 1., i);
			else
				m_domain.setRange(-1., 1., i);
		}
		generateAdLoadPF();
	}

	void JY3::generateAdLoadPF() {
		int num = 2000;
		Real t = get_time();
		Real At = 0.05;
		Real Wt = floor(6 * sin(0.5*OFEC_PI*(t - 1)));
		Real alpha_t = floor(100 * pow(sin(0.5*OFEC_PI*t), 2));
		std::vector<std::vector<Real>*> point(num);
		Matrix o_point(num, 2);
		for (int i = 0; i < num; i++) {
			Real y = fabs(i * 1. / (num - 1)*sin((2 * alpha_t + 0.5)*OFEC_PI*i * 1. / (num - 1)));
			o_point[i][0] = y + At * sin(Wt*OFEC_PI*y);
			o_point[i][1] = 1 - y + At * sin(Wt*OFEC_PI*y);
			point[i] = &o_point[i].vect();
		}

		std::vector<std::vector<Real>*> &temp = get_nondominated_set(point, OptimizeMode::kMinimize);
		for (int i = 0; i < temp.size(); i++) {
			m_optima->appendObj(*temp[i]);
		}
		m_optima->setVariableGiven(true);//given optima
	}

	int JY3::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void JY3::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(false);
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		
		Real gt = 0; //g function
		Real At = 0.05;
		Real alpha_t = floor(100 * std::pow(sin(0.5*OFEC_PI*t),2));
		Real Wt = floor(6 * sin(0.5*OFEC_PI*(t - 1)));
		Real y1 = fabs(x[0]*sin((2*alpha_t+0.5)*OFEC_PI*x[0]));
		for (size_t i = 1; i < m_number_variables; ++i){
			if (i < 2)
				gt += pow(pow(x[i], 2) - y1, 2);
			else
				gt += pow(pow(x[i],2)-x[i-1],2);
		}
		obj[0] = (1 + gt) * (y1 + At * sin(Wt * OFEC_PI * y1));
		obj[1] = (1 + gt) * (1 - y1 + At * sin(Wt * OFEC_PI * y1));
	}
}