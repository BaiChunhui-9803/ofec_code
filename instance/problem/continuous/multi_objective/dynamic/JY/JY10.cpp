#include "JY10.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void JY10::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the first variable is in (0,1), others are in (-1,1)
			if (i < 1)
				m_domain.setRange(0., 1., i);
			else
				m_domain.setRange(-1., 1., i);
		}
		generateAdLoadPF();
	}

	void JY10::generateAdLoadPF() {
		int num = 2000;
		Real t = get_time();
		Real At = 0.05;
		Real Wt = 6;
		Real Gt = fabs(sin(0.5*OFEC_PI*t));//change factor
		Real alpha_t = 1 + m_sigma * Gt;
		Real beta_t = alpha_t;
		std::vector<std::vector<Real>*> point(num);
		Matrix o_point(num, 2);
		if (m_sigma == 0 || m_sigma == 1) {
			for (int i = 0; i < num; i++) {
				Real x = i * 1. / (num - 1);
				o_point[i][0] = pow((x + At * sin(Wt * OFEC_PI * x)), alpha_t);
				o_point[i][1] = pow((1 - x + At * sin(Wt * OFEC_PI * x)), beta_t);
				point[i] = &o_point[i].vect();
			}
		}
		else {
			for (int i = 0; i < num; i++) {
				Real x = i * 1. / (num - 1);
				o_point[i][0] = (1 + (m_number_variables - 1) * pow((1 - Gt), 2)) * pow((x + At * sin(Wt * OFEC_PI * x)), alpha_t);
				o_point[i][1] = (1 + (m_number_variables - 1) * pow((1 - Gt), 2)) * pow((1 - x + At * sin(Wt * OFEC_PI * x)), beta_t);
				point[i] = &o_point[i].vect();
			}
		}

		std::vector<std::vector<Real>*>& temp = get_nondominated_set(point, OptimizeMode::kMinimize);
		for (int i = 0; i < temp.size(); i++) {
			m_optima->appendObj(*temp[i]);
		}
		m_optima->setVariableGiven(true);//given optima
	}

	int JY10::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void JY10::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real t = get_time();
		size_t rho_t = 5;//the time steps a type of problem lasts
		size_t changefre = get_change_fre();
		size_t effective_eval = get_effective_eval();
		if (effective_eval / changefre / rho_t == 0) {
			size_t temp = (size_t)std::ceil(3 * m_random->uniform.next());//choose one dimension as the first objective randomly
			set_r(temp);
		}
		m_sigma = (effective_eval / changefre / rho_t + get_r()) % 3;
		if (time_changed() && t != 0. && get_updated_state()) {
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(false);
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		Real gt = 0; //g function
		Real Gt = fabs(std::sin(0.5*OFEC_PI*t));//change factor
		Real At = 0.05;
		Real Wt = 6;
		Real alpha_t = 1 + get_sigma() * Gt;
		Real beta_t = alpha_t;
		for (size_t i = 1; i < m_number_variables; ++i) {
			gt += pow(x[i] + get_sigma() - Gt, 2);
		}
		obj[0] = (1 + gt) * pow(x[0] + At * sin(Wt * OFEC_PI * x[0]), alpha_t);
		obj[1] = (1 + gt) * pow(1 - x[0] + At * sin(Wt * OFEC_PI * x[0]), beta_t);
	}
}