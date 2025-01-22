#include "UDF2.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void UDF2::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the first one is in [0,1], others are in [-2,2]
			if (i < 1)
				m_domain.setRange(0., 1., i);
			else
				m_domain.setRange(-2., 2., i);
		}
		generateAdLoadPF();
	}

	//recommend n=20;
	void UDF2::generateAdLoadPF() {
		int num = 2000;
		Real t = get_time();
		Real Gt = sin(0.5*OFEC_PI*t);
		Real temp = std::fabs(Gt);
		Matrix o_point(num, 2);
		for (int i = 0; i < num; i++) {
			o_point[i][0] = temp + i * 1. / (num - 1);
			o_point[i][1] = 1 - o_point[i][0] + 2 * temp;
			m_optima->appendObj(o_point[i].vect());
		}
		m_optima->setVariableGiven(true);//given optima
	}

	int UDF2::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void UDF2::evaluateObjective(Real *x, std::vector<Real> &obj) {//recommend the number of variables is 20
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(false);
		}
		
		Real Gt = std::sin(0.5 * OFEC_PI*t);

		if (m_number_variables < 3) {
			obj[0] = x[0] + fabs(Gt);
			obj[1] = 1 - x[0] + fabs(Gt);
		}
		else {
			obj[0] = 0;
			for (size_t i = 2; i < m_number_variables; i += 2) {
				obj[0] += pow(x[i] - pow(x[0], 0.5 * (2 + Gt + 3. * (i - 1) / (m_number_variables - 2))) - Gt, 2);
			}
			obj[0] = x[0] + 2. / (ceil(m_number_variables / 2.) - 1) * obj[0] + fabs(Gt);

			obj[1] = 0;
			for (size_t i = 1; i < m_number_variables; i += 2) {
				obj[1] += pow(x[i] - pow(x[0], 0.5 * (2 + Gt + 3. * (i - 1) / (m_number_variables - 2))) - Gt, 2);
			}
			obj[1] = 1 - x[0] + 2. / (floor(m_number_variables / 2.)) * obj[1] + fabs(Gt);
		}
	}
}