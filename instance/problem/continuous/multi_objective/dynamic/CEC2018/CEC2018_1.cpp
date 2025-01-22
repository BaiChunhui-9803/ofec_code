#include "CEC2018_1.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void CEC2018_1::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the variable is in (0,1)
			m_domain.setRange(0., 1., i);
		}
		generateAdLoadPF();
	}

	void CEC2018_1::generateAdLoadPF() {
		int num = 2000;
		Real t = get_time();
		Real Ht = 1.25 + 0.75*sin(0.5*OFEC_PI*t);
		std::vector<std::vector<Real>*> point(num);
		Matrix o_point(num, 2);
		for (int i = 0; i < num; i++) {
			o_point[i][0] = i * 1. / (num - 1);
			o_point[i][1] = 1 - pow(o_point[i][0],Ht);
			point[i] = &o_point[i].vect();
		}

		std::vector<std::vector<Real>*> &temp = get_nondominated_set(point, OptimizeMode::kMinimize);
		for (int i = 0; i < temp.size(); i++) {
			m_optima->appendObj(*temp[i]);
		}
		m_optima->setVariableGiven(true);//given optima
	}

	int CEC2018_1::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	void CEC2018_1::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			set_updated_state(false);
			generateAdLoadPF();
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		Real gt = 1; //g function
		Real Gt = fabs(sin(0.5*OFEC_PI*t));//change factor
		Real Ht = 1.25 + 0.75*sin(0.5*OFEC_PI*t);
		for (size_t i = 1; i < m_number_variables; ++i)
		{
			gt += pow(x[i] - Gt, 2);
		}
		obj[0] = x[0];
		obj[1] = gt * (1-pow(obj[0]/gt,Ht));
	}
}