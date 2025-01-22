#include "DTF.h"
#include "../../../../../../core/algorithm/algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void DTF::initialize_() {
		DMOPs::initialize_();
		for (size_t i = 0; i < m_number_variables; ++i) {//the first one is in [0,1], others are in [-1,1]
			if (i < 1)
				m_domain.setRange(0., 1., i);
			else
				m_domain.setRange(-1., 1., i);
		}
		generateAdLoadPF();
	}

	void DTF::generateAdLoadPF() {
		int num = 2000;
		Real t = get_time();
		Real temp1 = 0.2 + 4.8*t*t;
		size_t s = 3;
		Real temp2 = t * s;
		std::vector<std::vector<Real>*> point(num);
		Matrix o_point(num, 2);
		for (int i = 0; i < num; i++) {
			o_point[i][0] = 0 + i * 1. / (num - 1);
			o_point[i][1] = 2 - std::pow(o_point[i][0], temp1) - o_point[i][0] * std::pow(std::fabs(std::sin(temp2*OFEC_PI*o_point[i][0])), temp1);
			point[i] = &o_point[i].vect();
		}
		
		std::vector<std::vector<Real>*> &temp = get_nondominated_set(point,OptimizeMode::kMinimize);
		for (int i = 0; i < temp.size(); i++) {
			m_optima->appendObj(*temp[i]);
		}
		//std::vector<int> m_rank(num);
		//std::vector<optimization_mode> opt_mode({optimization_mode::Minimization,optimization_mode::Minimization });
		////nd_sort::deductive_sort(point, m_rank, opt_mode);
		//int lin=nd_sort::fast_sort(point,m_rank,opt_mode);//∑µªÿ≈≈–Ú≤„ ˝
		/*std::ofstream out("../../../../../../DTF.txt",std::ios::app);
		if (out) {
			out << "obj1" << " " << "obj2" << std::endl;
			for (int i = 0; i < temp.size(); i++) {
				out << (*point[i])[0] << " " << (*point[i])[1] << std::endl;
			}
			out.close();
		}*/
	}

	int DTF::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
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

	//f1=x1, f2=g*h, recommend the number of decision variable is 20
	void DTF::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real t = get_time();
		Real gt = 1;
		Real alpha_t = 0.2 + 4.8 * t * t;
		Real beta_t = std::pow(10, 2 * std::sin(0.5 * OFEC_PI * t));
		Real gamma_t = std::sin(0.5 * OFEC_PI * t);
		size_t s = 3;
		Real pha_t = s * t;//the number of sections of the PF
		Real omega_t = 0.4 * pha_t;
		for (size_t i = 1; i < m_number_variables; ++i) {
			gt += std::pow(x[i] - gamma_t, 2) - std::cos(omega_t * OFEC_PI * (x[i] - gamma_t)) + 1;
		}
		if (time_changed() && t != 0. && get_updated_state()) {
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(false);
#ifdef OFEC_DEMO
			ofec_demo::g_buffer->appendProBuffer(this);
#endif
		}
		for (size_t m = 0; m < m_number_objectives; ++m){
			if (m < 1)
				obj[m] = std::pow(x[0],beta_t);
			else
				obj[m] = gt * (2 - std::pow(obj[0] / gt,alpha_t)- obj[0] / gt*std::pow(std::fabs(std::sin(pha_t*OFEC_PI*obj[0])),alpha_t));
		}
	}
}