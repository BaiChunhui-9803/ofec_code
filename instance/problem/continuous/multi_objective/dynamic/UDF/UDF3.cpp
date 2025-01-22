#include "UDF3.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void UDF3::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the first one is in [0,1], others are in [-1,1]
			if (i < 1)
				m_domain.setRange(0., 1., i);
			else
				m_domain.setRange(-1., 1., i);
		}
		generateAdLoadPF();
	}

	//recommend n=20;
	void UDF3::generateAdLoadPF() {
		int num = 2000;
		Real t = get_time();
		Real Gt = sin(0.5 * OFEC_PI * t);
		size_t N = 10;
		int m_num = num / N;
		Matrix o_point(num, 2);
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < m_num; j++) {
				o_point[i][0] = (2. * i + 1.) / 2 / N + fabs(Gt) + i * 1. / 2 / N / (m_num - 1);
				o_point[i][1] = 1 - o_point[i][0];
				m_optima->appendObj(o_point[i].vect());
			}
		}
		m_optima->setVariableGiven(true);//given optima
	}

	int UDF3::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void UDF3::evaluateObjective(Real *x, std::vector<Real> &obj) {//recommend the number of variables is 20
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(false);
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		
		Real Gt = std::sin(0.5 * OFEC_PI*t);
		Real alpha = 0.1;
		size_t N = 10;
		Real temp = std::max(0., (1. / 2 / N + alpha) * (sin(2 * N * OFEC_PI * x[0]) - 2 * N * fabs(Gt)));

		for (size_t m = 0; m < m_number_objectives; ++m) {
			if (m < 1) {
				if (m_number_variables < 3) {
					obj[m] = x[0] + temp;
				}
				else {
					Real temp1 = 0, temp2 = 1;
					for (size_t i = 2; i < m_number_variables; i += 2) {
						temp1 += 2 * pow(x[i] - sin(6 * OFEC_PI * x[0] + (i + 1) * OFEC_PI / m_number_variables), 2);
						temp2 *= cos(20 * OFEC_PI * (x[i] - sin(6 * OFEC_PI * x[0] + (i + 1) * OFEC_PI / m_number_variables)) / sqrt(i + 1));
					}
					obj[m] = x[0] + 2. / (ceil(m_number_variables / 2.) - 1) * pow(4 * temp1 - 2 * temp2 + 2, 2) + temp;
				}

			}
			else {
				Real temp1 = 0, temp2 = 1;
				for (size_t i = 1; i < m_number_variables; i += 2) {
					temp1 += 2 * pow(x[i] - sin(6 * OFEC_PI * x[0] + (i + 1) * OFEC_PI / m_number_variables), 2);
					temp2 *= cos(20 * OFEC_PI * (x[i] - sin(6 * OFEC_PI * x[0] + (i + 1) * OFEC_PI / m_number_variables)) / sqrt(i + 1));
				}
				obj[m] = 1 - x[0] + 2. / (floor(m_number_variables / 2.)) * pow(4 * temp1 - 2 * temp2 + 2, 2) + temp;
			}
		}
	}
}