#include "T3.h"

namespace ofec {
	T3::T3(param_map & v) : T3(v.at("problem name"), v.at("number of variables")) {

	}
	T3::T3(const std::string & name, size_t size_var) : problem(name, size_var, 2), DMOPs(name, size_var, 2) {

	}

	void T3::initialize() {//all variables are in（0,1）,need revise
		setInitialDomain(0., 1.);
		setDomain(0., 1.);
		
		//set_duration_gen(global::ms_arg.find("durationGeneration")->second);
		set_change_severity(global::ms_arg.find("changeSeverity")->second);
		set_change_fre(global::ms_arg.find("changeFre")->second);

		generateAdLoadPF();
		m_initialized = true;
	}

	void T3::generateAdLoadPF() {

	}

	int T3::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real t = get_time();
		if (time_changed() && t != 0. && (!get_updated_state())) {//防止不计数评价重复更新问题和重复采样PF
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(true);
		}
		else if (m_evaluations % (size_t)global::ms_arg.find("changeFre")->second != 0) {
			set_updated_state(false);
		}
		Real Gt = floor(m_number_objectives*fabs(sin(0.5*OFEC_PI*t)));
		Real r_t = floor(m_number_variables*fabs(pow(cos(2 * t), 3)));
		Real g = 0;
		for (size_t i = 0; i < m_number_variables; i++) {
			g += pow(x[i] - 0.5, 2);
		}

		for (size_t m = 0; m < m_number_objectives; ++m) {
			obj[m] = 1 + g;
			if (m < m_number_objectives - 1)
			{
				size_t j = 0;
				for (size_t M=0; M <m_number_variables - m - 1;++j, ++M)
					obj[m] *= cos(0.5*x[M] * OFEC_PI);
				if (m > 0)
					obj[m] = obj[m] * sin(0.5*x[j + 1] * OFEC_PI);
			}
			else
				obj[m] = obj[m] * sin(0.5*x[0] * OFEC_PI);
		}
		return kNormalEval;
	}
}