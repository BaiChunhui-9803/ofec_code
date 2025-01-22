#include "CEC2018_13.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void CEC2018_13::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the variable is in (0,1)
			if (i < 2)
				m_domain.setRange(0., 1., i);
			else
				m_domain.setRange(-1., 1., i);
		}
		generateAdLoadPF();
	}

	void CEC2018_13::generateAdLoadPF() {
		int n = 100;
		int num = n * n;
		Matrix o_point(num, 3);
		std::vector<std::vector<Real>*> point(num);
		Real t = get_time();
		Real Gt = sin(0.5*OFEC_PI*t);
		Real Pt = floor(6 *Gt);
		Real x1, x2;
		for (size_t i = 0; i < n; i++) {
			x1 = i * 1. / (n - 1);
			for (size_t j = 0; j < n; j++) {
				x2 = j * 1. / (n - 1);
				o_point[i*n + j][0] = cos(0.5*OFEC_PI*x1)*cos(0.5*OFEC_PI*x1);
				o_point[i*n + j][1] = cos(0.5*OFEC_PI*x2)*cos(0.5*OFEC_PI*x2);
				o_point[i*n + j][2] = pow(sin(0.5*OFEC_PI*x1), 2) + sin(0.5*OFEC_PI*x1)*pow(cos(Pt*OFEC_PI*x1), 2)\
					+pow(sin(0.5*OFEC_PI*x2), 2) + sin(0.5*OFEC_PI*x2)*pow(cos(Pt*OFEC_PI*x2), 2);
				point[i*n + j] = &o_point[i*n + j].vect();
				//m_optima->append(o_point[i*n + j].vect());
			}
		}
		std::vector<std::vector<Real>*> &temp = get_nondominated_set(point, OptimizeMode::kMinimize);
		for (int i = 0; i < temp.size(); i++) {
			m_optima->appendObj(*temp[i]);
		}
		m_optima->setVariableGiven(true);//given optima
	}

	 int CEC2018_13::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void CEC2018_13::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			set_updated_state(false);
			generateAdLoadPF();
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		Real gt = 1; //g function
		Real Gt = sin(0.5*OFEC_PI*t);//change factor
		Real Pt = floor(6*Gt);
		for (size_t i = 2; i < m_number_variables; ++i) {
			gt += pow(x[i] - Gt, 2);
		}

		obj[0] = gt*cos(0.5*OFEC_PI*x[0])*cos(0.5*OFEC_PI*x[0]);
		obj[1] = gt * cos(0.5*OFEC_PI*x[1])*cos(0.5*OFEC_PI*x[1]);
		obj[2] = gt * (pow(sin(0.5*OFEC_PI*x[0]),2)+ sin(0.5*OFEC_PI*x[0])*pow(cos(Pt*OFEC_PI*x[0]),2)\
			+pow(sin(0.5*OFEC_PI*x[1]), 2) + sin(0.5*OFEC_PI*x[1])*pow(cos(Pt*OFEC_PI*x[1]), 2));
	}
}