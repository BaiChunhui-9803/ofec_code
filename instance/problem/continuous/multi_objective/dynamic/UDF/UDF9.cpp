#include "UDF9.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void UDF9::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the first one is in [0,1], others are in [-2,2]
			if (i < 1)
				m_domain.setRange(0., 1., i);
			else
				m_domain.setRange(-1., 2., i);
		}
		generateAdLoadPF();
	}

	//recommend n=20;
	void UDF9::generateAdLoadPF() {
		int num = 2000;
		Real Gt3 = sin(0.5*OFEC_PI*m_count[2]);;
		Real Ht4 = 0.5 + std::fabs(std::sin(0.5*OFEC_PI*m_count[3]));
		Real Ht5 = 0.5 + std::fabs(std::sin(0.5*OFEC_PI*m_count[4]));
		Matrix o_point(num, 2);
		for (int i = 0; i < num; i++) {
			o_point[i][0] = std::fabs(Gt3) + i * 1. / (num - 1);
			o_point[i][1] = 1 + std::fabs(Gt3) - Ht4 * std::pow(o_point[i][0], Ht5);
			m_optima->appendObj(o_point[i].vect());
		}
		m_optima->setVariableGiven(true);//given optima
	}

	int UDF9::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
		auto eval = alg->evaluations();
		set_effective_eval(eval);
		if (alg != nullptr && eval % get_change_fre() == 0 && eval > 0) {
			Real temp = m_random->uniform.next();
			if (temp >= 0 && temp <= 0.2)
				m_count[0]++;
			else if (temp > 0.2 && temp <= 0.4)
				m_count[1]++;
			else if (temp > 0.4 && temp <= 0.6)
				m_count[2]++;
			else if (temp > 0.6 && temp <= 0.8)
				m_count[3]++;
			else
				m_count[4]++;
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
			return kChangeNextEval;
		}
		else if (alg != nullptr && (eval + 1) % get_change_fre() == 0 && eval > 0) {
			set_updated_state(true);
			return kNormalEval;
		}
		else {
			return Continuous::updateEvaluationTag(s, alg);
		}
	}

	void UDF9::evaluateObjective(Real *x, std::vector<Real> &obj) {//recommend the number of variables is 20
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			set_updated_state(false);
			generateAdLoadPF();
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}

		if (m_number_variables < 3) {
			obj[0] = x[0] + fabs(sin(0.5 * OFEC_PI * m_count[2]));
			obj[1] = 1 - (0.5 + fabs(sin(0.5 * OFEC_PI * m_count[3]))) * pow(x[0], 0.5 + fabs(sin(0.5 * OFEC_PI * m_count[4]))) + fabs(sin(0.5 * OFEC_PI * m_count[2]));
		}
		else {
			obj[0] = 0;
			for (size_t i = 2; i < m_number_variables; i += 2) {
				obj[0] += pow(x[i] - pow(x[0], 0.5 * (2 + 3. * (i - 1) / (m_number_variables - 2)\
					+ sin(0.5 * OFEC_PI * m_count[1]))) - sin(0.5 * OFEC_PI * m_count[2]), 2);
			}
			obj[0] = x[0] + fabs(sin(0.5 * OFEC_PI * m_count[2])) + 2. / (ceil(m_number_variables / 2.) - 1) * obj[0];

			obj[1] = 0;
			for (size_t i = 1; i < m_number_variables; i += 2) {
				obj[1] += pow(x[i] - pow(x[0], 0.5 * (2 + 3. * (i - 1) / (m_number_variables - 2)\
					+ sin(0.5 * OFEC_PI * m_count[1]))) - sin(0.5 * OFEC_PI * m_count[2]), 2);
			}
			obj[1] = 1 - (0.5 + fabs(sin(0.5 * OFEC_PI * m_count[3]))) * pow(x[0], 0.5 + \
				fabs(sin(0.5 * OFEC_PI * m_count[4]))) + fabs(sin(0.5 * OFEC_PI * m_count[2])) + 2. / (floor(m_number_variables / 2.)) * obj[1];
		}
	}
}