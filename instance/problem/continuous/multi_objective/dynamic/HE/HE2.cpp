#include "HE2.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void HE2::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//all variable is in (0,1)
			m_domain.setRange(0., 1., i);
		}
		generateAdLoadPF();
	}

	//f1 = x1; f2 = g * h; f2 = 1 - (f1)^Ht-f1^Ht*sin(10¦Ðtf1); f1 = [0, 1]
	void HE2::generateAdLoadPF() {
		int num = 2000;
		Real t = get_time();
		Real temp1 = 0.75*std::sin(0.5*OFEC_PI*t) + 1.25;
		std::vector<std::vector<Real>*> point(num);
		Matrix o_point(num, 2);
		for (int i = 0; i < num; i++) {
			o_point[i][0] = 0 + i * 1. / (num - 1);
			o_point[i][1] = 1 - pow(sqrt(o_point[i][0]), temp1) - pow(o_point[1][0], temp1)*sin(10 * OFEC_PI*t*o_point[1][0]);
			point[i] = &o_point[i].vect();
		}

		std::vector<std::vector<Real>*> &temp = get_nondominated_set(point, OptimizeMode::kMinimize);
		for (int i = 0; i < temp.size(); i++) {
			m_optima->appendObj(*temp[i]);
		}
		m_optima->setVariableGiven(true);//given optima
	}

	int HE2::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	//f1 = x1; f2 = g * h; f2 = 1 - (f1)^Ht-f1^Ht*sin(10¦Ðtf1)
	void HE2::evaluateObjective(Real* x, std::vector<Real>& obj) {//recommend the number of variables is 20
		Real t = get_time();
		if (time_changed() && t != 0. && get_updated_state()) {
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(false);
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		Real gt = 0; //g function
		Real Ht = 0.75*std::sin(0.5*OFEC_PI*t) + 1.25;
		for (size_t i = 1; i < m_number_variables; ++i)
		{
			gt += x[i];
		}
		gt = 1 + 9 * gt / (m_number_variables - 1);
		obj[0] = x[0];
		obj[1] = gt * (1 - pow(obj[0] / gt, Ht / 2) - pow(obj[0] / gt, Ht) * sin(10 * OFEC_PI * obj[0]));
	}
}