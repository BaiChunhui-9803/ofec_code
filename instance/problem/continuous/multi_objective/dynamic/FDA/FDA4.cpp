#include "FDA4.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void FDA4::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {
			m_domain.setRange(0., 1., i);
		}
		generateAdLoadPF();
	}

	//f1^2+f2^2+f3^2=1, recommend n=20
	void FDA4::generateAdLoadPF() {
		int n = 100;
		int num = n * (1 + n) / 2;
		Matrix o_point(num, 3);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < i + 1; j++) {
				o_point[i*(1 + i) / 2 + j][0] = 1 - (Real)i / (n - 1);
				o_point[i*(1 + i) / 2 + j][1] = (Real)i*j / (n - 1);
				o_point[i*(1 + i) / 2 + j][2] = 1 - o_point[i*(1 + i) / 2 + j][0] - o_point[i*(1 + i) / 2 + j][1];
			}
		}
		//map points on x+y+z=1 into x^2+y^2+z^2=1
		for (int i = 0; i < num; i++) {
			Real temp = sqrt(pow(o_point[i][0], 2) + pow(o_point[i][1], 2) + pow(o_point[i][2], 2));
			o_point[i][0] = o_point[i][0] / temp;
			o_point[i][1] = o_point[i][1] / temp;
			o_point[i][2] = o_point[i][2] / temp;
			m_optima->appendObj(o_point[i].vect());
		}
	}

	int FDA4::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
		size_t eval = alg->evaluations();
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

	//f1^2+f2^2+...+fn^2=1, f = [0, 11],
	void FDA4::evaluateObjective(Real *x, std::vector<Real> &obj) {//recommend n=M+9, x2 has 10 variables
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(false);
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		Real Gt = std::fabs(std::sin(0.5*OFEC_PI*t));
		Real gt = 0;
		for (size_t i = m_number_variables - m_number_objectives+1; i < m_number_variables; ++i) {
			gt += std::pow(x[i] - Gt, 2);
		}
		for (size_t m = 0; m < m_number_objectives; ++m) {
			obj[m] = 1 + gt;
			if (m < m_number_objectives - 1)
			{
				size_t M = 0;
				for (size_t j=0; j < m_number_objectives - m - 1; ++j, ++M) {
					obj[m] *= cos(0.5 * x[M] * OFEC_PI);
				}
				if (m > 0) {
					obj[m] = obj[m] * sin(0.5*x[M + 1] * OFEC_PI);
				}
			}
			else
				obj[m] = obj[m] * std::sin(0.5*x[0] * OFEC_PI);
		}
	}
}