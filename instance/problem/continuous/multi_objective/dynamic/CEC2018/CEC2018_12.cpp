#include "CEC2018_12.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void CEC2018_12::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the variable is in (0,1)
			if (i < 2)
				m_domain.setRange(0., 1., i);
			else
				m_domain.setRange(-1., 1., i);
		}
		generateAdLoadPF();
	}

	void CEC2018_12::generateAdLoadPF() {
		int n = 100;
		int num = n * n;
		Matrix o_point(num, 3);
		std::vector<std::vector<Real>*> point(num);
		Real t = get_time();
		int Kt = (int)floor(10 * sin(OFEC_PI*t));//change factor
		Real r = 1 - Kt % 2;
		Real x1, x2;
		Real gt;
		for (size_t i = 0; i < n; i++) {
			x1 = i * 1. / (n - 1);
			for (size_t j = 0; j < n; j++) {
				x2 = j * 1. / (n - 1);
				if (((int)fabs((int)floor(Kt*(2 * x1 - r))) % 2)*((int)fabs((int)floor(Kt*(2 * x2 - r))) % 2) == 0) {
					gt = 1 + fabs(sin(floor(Kt*(2 * x1 - r))*OFEC_PI / 2)*sin(floor(Kt*(2 * x2 - r))*OFEC_PI / 2));
					o_point[i*n + j][0] = gt * cos(0.5*OFEC_PI*x1)*cos(0.5*OFEC_PI*x2);
					o_point[i*n + j][1] = gt * cos(0.5*OFEC_PI*x1)*sin(0.5*OFEC_PI*x2);
					o_point[i*n + j][2] = gt * sin(0.5*OFEC_PI*x1);
					//point[i*n + j] = &o_point[i*n + j].vect();
					m_optima->appendObj(o_point[i*n + j].vect());
				}
			}
		}
		/*std::vector<std::vector<Real>*> &temp = get_nondominated_set(point, optimization_mode::Minimization);
		for (int i = 0; i < temp.size(); i++) {
			m_optima->append(*temp[i]);
		}*/
		m_optima->setVariableGiven(true);
	}

	 int CEC2018_12::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void CEC2018_12::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			set_updated_state(false);
			generateAdLoadPF();
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		Real gt = 1; //g function
		int Kt = (int)floor(10*sin(OFEC_PI*t));//change factor
		Real r = 1 - Kt%2;
		for (size_t i =2; i < m_number_variables; ++i){
			gt += pow(x[i] - sin(t*x[0]), 2);
		}
		gt += fabs(sin(floor(Kt*(2*x[0]-r))*OFEC_PI/2)*sin(floor(Kt*(2 * x[1] - r))*OFEC_PI / 2));
		obj[0] = gt*cos(0.5*OFEC_PI*x[0])*cos(0.5*OFEC_PI*x[1]);
		obj[1] = gt * cos(0.5*OFEC_PI*x[0])*sin(0.5*OFEC_PI*x[1]);
		obj[2] = gt*sin(0.5*OFEC_PI*x[0]);
	}
}