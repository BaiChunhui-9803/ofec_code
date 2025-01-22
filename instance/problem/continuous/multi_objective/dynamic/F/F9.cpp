#include "F9.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec::DMOP {
	void F9::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//all variables are in [0,5], recommend n=20;f1=[0,1]
			m_domain.setRange(0., 5., i);
		}
		generateAdLoadPF();
	}

	//f1=s^H, f2=(1-s)^H, f2=(1-f1^(1/H))^H
	void F9::generateAdLoadPF() {
		int num = 2000;
		Real t = get_time();
		Real temp1 = 1.25 + 0.75*std::sin(OFEC_PI*t);
		Matrix o_point(num, 2);
		for (int i = 0; i < num; i++) {
			Real s = 0 + i * 1. / (num - 1);
			o_point[i][0] = pow(s, temp1);
			o_point[i][1] = pow(1 - s, temp1);
			m_optima->appendObj(o_point[i].vect());
		}
		m_optima->setVariableGiven(true);//given optima
	}

	int F9::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	//f1=x1, f2=g*h
	void F9::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(false);
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		
		Real a = 2 + 2 * std::cos(OFEC_PI*(t - std::floor(t)));
		Real b = 2 + 2 * std::sin(OFEC_PI*2*(t - std::floor(t)));
		Real Ht = 1.25 + 0.75*std::sin(OFEC_PI*t);

		for (size_t m = 0; m < m_number_objectives; ++m)
		{
			if (m < 1) {
				if (m_number_variables < 2) {
					obj[m] = pow(fabs(x[0] - a), Ht);
				}
				else {
					obj[m] = pow(fabs(x[0] - a), Ht);
					for (size_t i = 2; i < m_number_variables; i += 2)
						obj[m] += pow(x[i] - b - 1 + pow(fabs(x[0] - a), Ht + (i + 1.) / m_number_variables), 2);
				}
			}
			else {
				obj[m] = pow(fabs(x[0] - a - 1), Ht);
				for (size_t i = 1; i < m_number_variables; i += 2)
					obj[m] += pow(x[i] - b - 1 + pow(fabs(x[0] - a), Ht + (i + 1.) / m_number_variables), 2);
			}
		}
	}
}