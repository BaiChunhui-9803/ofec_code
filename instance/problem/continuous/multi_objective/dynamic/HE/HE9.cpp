#include "HE9.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void HE9::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//all in [0,1]
			m_domain.setRange(0., 1., i);
		}
		generateAdLoadPF();
	}

	//recommend n=20, f1 is in [0,1]
	void HE9::generateAdLoadPF() {//the method to get Real PF is over sampling in 2 dimension
		int num = 2000;
		Real t = get_time();
		Real temp = 0.75*sin(0.5*OFEC_PI*t) + 1.25;
		Matrix o_point(num, 2);
		for (int i = 0; i < num; i++) {
			o_point[i][0] = 0 + i * 1. / (num - 1);
			o_point[i][1] = (2 - sqrt(o_point[i][0]))*(1 - pow(o_point[i][0] / 2, temp));
			m_optima->appendObj(o_point[i].vect());
		}
		m_optima->setVariableGiven(true);//given optima
	}

	int HE9::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void HE9::evaluateObjective(Real *x, std::vector<Real> &obj) {//recommend the number of variables is 20
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(false);
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		
		Real Ht = 0.75*std::sin(0.5 * OFEC_PI*t) + 1.25;

		Real tmp1 = 0;
		if (m_number_variables > 2) {
			Real temp1 = 0, temp2 = 1;
			for (size_t i = 2; i < m_number_variables; i += 2) {
				Real y = x[i] - pow(x[0], 0.5 * (1. + 3. * (i - 1) / (m_number_variables - 2)));
				temp1 = temp1 + y * y;
				temp2 *= cos(20 * y * OFEC_PI / sqrt(i + 1.));
			}
			tmp1 = tmp1 + 4 * temp1 - 2 * temp2 + 2.;
		}
		obj[0] = x[0] + 2. / (ceil(m_number_variables / 2.) - 1) + tmp1;

		Real gt = 0;
		if (m_number_variables > 2) {
			Real temp1 = 0, temp2 = 1;
			for (size_t i = 1; i < m_number_variables; i += 2) {
				Real y = x[i] - pow(x[0], 0.5 * (1. + 3. * (i - 1) / (m_number_variables - 2.)));
				temp1 += y * y;
				temp2 *= cos(20 * y * OFEC_PI / sqrt(i + 1.));
			}
			gt = gt + 4 * temp1 - 2 * temp2 + 2.;
		}
		gt = 2 - sqrt(x[0]) + 2. / floor(m_number_variables / 2.) * gt;
		obj[1] = gt * (1 - pow(obj[0] / gt, Ht));
	}
}