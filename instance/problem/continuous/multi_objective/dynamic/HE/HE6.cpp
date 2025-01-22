#include "HE6.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void HE6::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {
			if (i < 1)
				m_domain.setRange(0., 1., i);
			else
				m_domain.setRange(-1., 1., i);
		}
		generateAdLoadPF();
	}

	//recommend n=20;f1 is in [0,1]
	void HE6::generateAdLoadPF() {//the method to get Real PF is over sampling in 2 dimension
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

	int HE6::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void HE6::evaluateObjective(Real* x, std::vector<Real>& obj) {//recommend the number of variables is 20
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

		for (size_t m = 0; m < m_number_objectives; ++m)
		{
			if (m < 1) {
				if (m_number_variables > 2) {
					obj[m] = 0;
					for (size_t i = 2; i < m_number_variables; i += 2) {
						obj[m] += pow(x[i] - 0.8 * x[0] * std::cos((6 * OFEC_PI * x[0] + (i + 1) * OFEC_PI / m_number_variables) / 3), 2);
					}
					obj[m] = x[0] + 2. / (ceil(m_number_variables / 2.) - 1) * obj[m];
				}
				obj[m] = x[0];
			}
			else {
				Real gt = 0;
				for (size_t i = 1; i < m_number_variables; i += 2) {
					gt += pow(x[i] - 0.8 * cos(6 * OFEC_PI * x[0] + (i + 1) * OFEC_PI / m_number_variables), 2);
				}
				gt = 2 - sqrt(x[0]) + 2. / floor(m_number_variables / 2.) * gt;
				obj[m] = gt * (1 - pow(obj[0] / gt, Ht));
			}
		}
	}
}