#include "UDF7.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void UDF7::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the first two is in [0,1], others are in [-2,2]
			if (i < 2)
				m_domain.setRange(0., 1., i);
			else
				m_domain.setRange(-2., 2., i);
		}
		generateAdLoadPF();
	}

	//(f1-Gt)^2+(f2-Gt)^2+(f3-Gt)^2=Rt^2, recommend n=20;
	void UDF7::generateAdLoadPF() {
		int n = 100;
		int num = n * n;
		Real t = get_time();
		Real Gt = sin(0.5 * OFEC_PI * t);
		Real Rt = 1 + Gt;
		Matrix o_point(num, 3);
		//std::vector<std::vector<real>*> point(num);
		Real x1, x2;
		for (int i = 0; i < n; i++) {
			x1 = i * 1. / (n - 1);
			for (int j = 0; j < n; j++) {
				x2 = j * 1. / (n - 1);
				o_point[i * n + j][0] = Rt * cos(x1 * OFEC_PI / 2) * cos(x2 * OFEC_PI / 2) + Gt;
				o_point[i * n + j][1] = Rt * cos(x1 * OFEC_PI / 2) * sin(x2 * OFEC_PI / 2) + Gt;
				o_point[i * n + j][2] = Rt * sin(x1 * OFEC_PI / 2) + Gt;
				m_optima->appendObj(o_point[i * n + j].vect());
			}
		}
		/*std::vector<std::vector<real>*> &temp = get_nondominated_set(point, optimization_mode::Minimization);
		for (int i = 0; i < temp.size(); i++) {
			m_optima->append(*temp[i]);
		}*/
		m_optima->setVariableGiven(true);//given optima
	}

	int UDF7::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void UDF7::evaluateObjective(Real *x, std::vector<Real> &obj) {//recommend the number of variables is 20
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(false);
		}
		
		Real Gt = std::sin(0.5 * OFEC_PI*t);//Vertical or Horizontal shift in PF or PS
		Real Rt = 1 + std::fabs(Gt);//Radius variation in three dimensional PF

		if (m_number_variables < 4) {
			obj[0] = Gt + Rt * cos(0.5 * OFEC_PI * x[0]) * cos(0.5 * OFEC_PI * x[1]);
		}
		else {
			obj[0] = 0;
			for (size_t i = 3; i < m_number_variables; i += 3) {
				obj[0] += pow(x[i] - 2 * x[1] * sin(2 * OFEC_PI * x[0] + (i + 1) * OFEC_PI / m_number_variables), 2);
			}
			obj[0] = Rt * cos(0.5 * OFEC_PI * x[0]) * cos(0.5 * OFEC_PI * x[1]) + 2. / (ceil(m_number_variables / 3.) - 1) * obj[0] + Gt;
		}

		if (m_number_variables < 5) {
			obj[1] = Gt + Rt * cos(0.5 * OFEC_PI * x[0]) * sin(0.5 * OFEC_PI * x[1]);
		}
		else {
			obj[1] = 0;
			for (size_t i = 4; i < m_number_variables; i += 3) {
				obj[1] += pow(x[i] - 2 * x[1] * sin(2 * OFEC_PI * x[0] + (i + 1) * OFEC_PI / m_number_variables), 2);
			}
			obj[1] = Rt * cos(0.5 * OFEC_PI * x[0]) * sin(0.5 * OFEC_PI * x[1]) + 2. / (ceil((m_number_variables + 1) / 3.) - 1) * obj[1] + Gt;
		}

		if (m_number_variables < 3) {
			obj[2] = Gt + Rt * sin(0.5 * OFEC_PI * x[0]);
		}
		else {
			obj[2] = 0;
			for (size_t i = 2; i < m_number_variables; i += 3) {
				obj[2] += pow(x[i] - 2 * x[1] * sin(2 * OFEC_PI * x[0] + (i + 1) * OFEC_PI / m_number_variables), 2);
			}
			obj[2] = Gt * sin(0.5 * OFEC_PI * x[0]) + 2. / (ceil((m_number_variables - 1) / 3.) - 1) * obj[2] + Gt;
		}
	}
}