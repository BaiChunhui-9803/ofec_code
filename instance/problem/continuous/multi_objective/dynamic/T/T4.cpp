#include "T4.h"

namespace ofec {
	T4::T4(param_map & v) : T4(v.at("problem name"), v.at("number of variables")) {

	}
	T4::T4(const std::string & name, size_t size_var) : problem(name, size_var, 2), DMOPs(name, size_var, 2) {

	}

	void T4::initialize() {//all variables are in（0,1）,need revise
		setInitialDomain(0., 1.);
		setDomain(0., 1.);
		
		//set_duration_gen(global::ms_arg.find("durationGeneration")->second);
		set_change_severity(global::ms_arg.find("changeSeverity")->second);
		set_change_fre(global::ms_arg.find("changeFre")->second);

		generateAdLoadPF();
		m_initialized = true;
	}

	void T4::generateAdLoadPF() {

	}

	int T4::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real t = get_time();
		if (time_changed() && t != 0. && (!get_updated_state())) {//防止不计数评价重复更新问题和重复采样PF
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(true);
		}
		else if (m_evaluations % (size_t)global::ms_arg.find("changeFre")->second != 0) {
			set_updated_state(false);
		}
		Real r_t = 0;
		for (size_t i = 0; i < m_number_variables; i++) {
			r_t+= x[i];
		}
		r_t = r_t / m_number_variables;

		for (size_t m = 0; m < m_number_objectives; ++m) {
			if (m < 1){
				obj[m] = 0;
				for (size_t j=0; j < m_number_variables; ++j)
					obj[m] += pow(x[j],2)-10*cos(2*x[j] * OFEC_PI)+10;
			}
			else {
				obj[m] = pow(x[0] - r_t, 2);
				for (size_t i = 1; i < m_number_variables; i++)
					obj[m] += pow(pow(x[i], 2) - x[i - 1], 2);
			}
		}
		return kNormalEval;
	}
}