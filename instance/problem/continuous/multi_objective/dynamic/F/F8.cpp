#include "F8.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec::DMOP {
	void F8::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the first two is in[0,1]��others are in [-1,2], recommend n=20
			if (i < 2)
				m_domain.setRange(0., 1., i);
			else
				m_domain.setRange(-1., 2.,i);
		}
		generateAdLoadPF();
	}

	//f1^2+f2^2+f3^2=1
	void F8::generateAdLoadPF() {
		int n = 100;//the number of sample points in each dimension
		int num = n*n;//total number of sample points in Real PF
		Matrix o_point(num, 3);
		std::vector<std::vector<Real>*> point(num);
		Real x1, x2;
		for (int i = 0; i < n; i++) {
			x1 = i * 1. / (n - 1);
			for (int j = 0; j < n; j++) {
				x2 = j * 1. / (n - 1);
				o_point[i * n + j][0] = cos(x1 * OFEC_PI / 2) * cos(x2 * OFEC_PI / 2);
				o_point[i * n + j][1] = cos(x1 * OFEC_PI / 2) * sin(x2 * OFEC_PI / 2);
				o_point[i * n + j][2] = sin(x1 * OFEC_PI / 2);
				point[i * n + j] = &o_point[i * n + j].vect();
				m_optima->appendObj(o_point[i].vect());
			}
		}
		//std::vector<std::vector<Real>*>& temp = get_nondominated_set(point, OptimizeMode::kMinimize);
		m_optima->setVariableGiven(true);//given optima
	}

	int F8::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	//f1=x1, f2=g*h
	void F8::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real t = get_time();
		Real Gt = std::sin(0.5*OFEC_PI*t);
		Real Ht = 1.25 + 0.75*std::sin(OFEC_PI*t);
		Real g = 0;
		for (size_t i = 2; i < m_number_variables; i++) {
			g += std::pow(x[i] - std::pow((x[0] + x[1]) / 2, Ht) - Gt, 2);
		}

		obj[0] = (1 + g)*std::cos(0.5*OFEC_PI*x[1])*std::cos(0.5*OFEC_PI*x[0]);
		obj[1] = (1 + g)*std::cos(0.5*OFEC_PI*x[1])*std::sin(0.5*OFEC_PI*x[0]);
		obj[2] = (1 + g)*std::sin(0.5*OFEC_PI*x[1]);
	}
}