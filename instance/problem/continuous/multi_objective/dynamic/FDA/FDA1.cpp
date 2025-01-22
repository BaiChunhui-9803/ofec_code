#include "FDA1.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void FDA1::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the first variable is in [0,1], others are in [-1,1]
			if (i < 1)
				m_domain.setRange(0., 1., i);
			else
				m_domain.setRange(-1., 1., i);
		}
		generateAdLoadPF();
		//m_initialized = true;
	}

	//f1 = x1, f2 = g * h, f2 = 1 - sqrt(f1), f1 = [0, 1]
	void FDA1::generateAdLoadPF(){
		int num = 2000;
		Matrix o_point(num,2);
		for (int i = 0; i < num; i++) {
			o_point[i][0] = 0+i*1./(num-1);
			o_point[i][1] = 1 - std::sqrt(o_point[i][0]);
			m_optima->appendObj(o_point[i].vect());
		}
	}

	 int FDA1::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	//f1=x1, f2=g*h, f2=1-sqrt(f1), f1=[0,1]
	void FDA1::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real t = get_time();
		Real gt = 1; //g function
		Real Gt = std::sin(0.5 * OFEC_PI * t);//change factor
		for (size_t i = 1; i < m_number_variables; ++i)
		{
			gt += std::pow(x[i] - Gt, 2);
		}
		if (time_changed() && t != 0. && get_updated_state()) {//avoid re-sample
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(false);
//#ifdef OFEC_DEMO
//			ofec_demo::g_buffer->appendProBuffer(this);
//#endif
		}
		obj[0] = x[0];
		obj[1] = gt * (1 - std::sqrt(obj[0] / gt));
		//std::cout << obj[0] << "  " << obj[1] << std::endl;
	}
}