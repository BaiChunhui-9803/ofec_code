#include "CEC2018_4.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void CEC2018_4::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the variable is in (0,1)
			m_domain.setRange(-2., 2., i);
		}
		generateAdLoadPF();
	}

	void CEC2018_4::generateAdLoadPF() {
		int num = 2000;
		Real t = get_time();
		Real a = sin(0.5*OFEC_PI*t);
		Real b = 1 + fabs(cos(0.5*OFEC_PI*t));
		Real Ht = 1.5 + a;
		std::vector<std::vector<Real>*> point(num);
		Matrix o_point(num, 2);
		Real x;
		for (int i = 0; i < num; i++) {
			x = a + b * i / (num-1);
			o_point[i][0] = pow(fabs(x-a),Ht);
			o_point[i][1] = pow(fabs(x - a-b), Ht);
			point[i] = &o_point[i].vect();
		}

		std::vector<std::vector<Real>*> &temp = get_nondominated_set(point, OptimizeMode::kMinimize);
		for (int i = 0; i < temp.size(); i++) {
			m_optima->appendObj(*temp[i]);
		}
		m_optima->setVariableGiven(true);//given optima
	}

	 int CEC2018_4::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void CEC2018_4::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			set_updated_state(false);
			generateAdLoadPF();
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		Real a = sin(0.5*OFEC_PI*t);
		Real b = 1 + fabs(cos(0.5*OFEC_PI*t));
		Real c = std::max(fabs(a),a+b);
		Real Ht = 1.5 + a;
		Real gt = 1;
		for (size_t i = 1; i < m_number_variables; ++i)
		{
			gt += pow(x[i] -a*x[0]*x[0]/(i+1)/c/c, 2);
		}

		obj[0] = gt * pow(fabs(x[0]-a),Ht);
		obj[1] = gt * pow(fabs(x[0]-a-b), Ht);
	}
}