#include "UDF6.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void UDF6::initialize_() {
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
	void UDF6::generateAdLoadPF() {
		int num = 2000;
		Real t = get_time();
		Real Gt = sin(0.5 * OFEC_PI * t);
		Real Mt = 0.5 + fabs(sin(0.5 * OFEC_PI * t));
		Real alpha = 0.1;
		size_t N = 10;
		Real tmp;
		Matrix o_point(num, 2);//
		Real x;
		for (int i = 0; i < num; i++) {
			x = i * 1. / (num - 1);
			tmp = (1. / 2 / N + alpha) * (sin(2 * N * OFEC_PI * x) - 2 * N * fabs(Gt));
			o_point[i][0] = x + tmp;
			o_point[i][1] = 1 - Mt * o_point[i][0] + tmp;
			m_optima->appendObj(o_point[i].vect());
		}
		m_optima->setVariableGiven(true);//given optima
	}

	int UDF6::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void UDF6::evaluateObjective(Real *x, std::vector<Real> &obj) {//recommend the number of variables is 20
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(false);
		}
		
		Real Gt = std::sin(0.5 * OFEC_PI*t);
		Real Mt = 0.5 + std::fabs(std::sin(0.5 * OFEC_PI*t));//Angular shift in PF
		Real alpha = 0.1;
		size_t N = 10;
		Real temp = (1. / 2 / N + alpha)*std::fabs(std::sin(2 * N*OFEC_PI*x[0]) - 2 * N*std::fabs(Gt));
		for (size_t m = 0; m < m_number_objectives; ++m){
			if (m < 1) {
				if (m_number_variables < 3) {
					obj[m] = x[0] + temp;
				}
				else {
					obj[m] = 0;
					for (size_t i = 2; i < m_number_variables; i += 2) {
						obj[m] += pow(2 * pow(x[i] - sin(6 * OFEC_PI * x[0] + (i + 1) * OFEC_PI / m_number_variables), 2)\
							- cos(4 * OFEC_PI * (x[i] - sin(6 * OFEC_PI * x[0] + (i + 1) * OFEC_PI / m_number_variables))) + 1, 2);
					}
					obj[m] = x[0] + 2. / (ceil(m_number_variables / 2.) - 1) * obj[m] + temp;
				}
			}
			else {
				Real temp1 = 0, temp2 = 1;
				for (size_t i = 1; i < m_number_variables; i += 2) {
					obj[m] += pow(2 * pow(x[i] - sin(6 * OFEC_PI * x[0] + (i + 1) * OFEC_PI / m_number_variables), 2)\
						- cos(4 * OFEC_PI * (x[i] - sin(6 * OFEC_PI * x[0] + (i + 1) * OFEC_PI / m_number_variables))) + 1, 2);
				}
				obj[m] = 1 - Mt * x[0] + 2. / (floor(m_number_variables / 2.)) * obj[m] + temp;
			}
		}
	}
}