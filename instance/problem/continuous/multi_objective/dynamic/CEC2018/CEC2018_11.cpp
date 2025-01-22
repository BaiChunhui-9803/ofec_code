#include "CEC2018_11.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void CEC2018_11::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the variable is in (0,1)
			m_domain.setRange(0., 1., i);
		}
		generateAdLoadPF();
	}

	void CEC2018_11::generateAdLoadPF() {//PF还没设计好
		int n = 100;
		int num = n * n;
		Matrix o_point(num, 3);
		//std::vector<std::vector<Real>*> point(num);
		Real t = get_time();
		Real Gt = fabs(sin(0.5*OFEC_PI*t));
		Real y1,y2;
		for (size_t i = 0; i < n; i++) {
			y1 = OFEC_PI*Gt/6+(OFEC_PI/2- OFEC_PI*Gt/3)*(i * 1. / (n - 1));
			for (size_t j = 0; j < n; j++) {
				y2 = OFEC_PI * Gt / 6 + (OFEC_PI / 2 - OFEC_PI * Gt / 3)*(j * 1. / (n - 1));
				o_point[i*n + j][0] = (1+Gt)*sin(y1);
				o_point[i*n + j][1] = (1 + Gt)*sin(y2)*cos(y1);
				o_point[i*n + j][2] = (1 + Gt)*cos(y2)*cos(y1);
				//point[i*n + j] = &o_point[i*n + j].vect();
				m_optima->appendObj(o_point[i*n + j].vect());
			}
		}
		/*std::vector<std::vector<Real>*> &temp = get_nondominated_set(point, optimization_mode::Minimization);
		for (int i = 0; i < temp.size(); i++) {
			m_optima->append(*temp[i]);
		}*/
		m_optima->setVariableGiven(true);//given optima
	}

	 int CEC2018_11::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void CEC2018_11::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			set_updated_state(false);
			generateAdLoadPF();
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		
		Real Gt = fabs(sin(0.5*OFEC_PI*t));//change factor
		Real gt = 1+Gt; //g function
		for (size_t i = 2; i < m_number_variables; ++i)
		{
			gt += pow(x[i] - 0.5*Gt*x[0], 2);
		}

		Real y1 = OFEC_PI/6*Gt+(OFEC_PI/2- OFEC_PI/3*Gt)*x[0];
		Real y2 = OFEC_PI / 6 * Gt + (OFEC_PI / 2 - OFEC_PI / 3 * Gt)*x[1];
		obj[0] = gt*sin(y1);
		obj[1] = gt * sin(y2)*cos(y1);
		obj[2] = gt * cos(y2)*cos(y1);
	}
}