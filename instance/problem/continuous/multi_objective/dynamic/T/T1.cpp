#include "T1.h"

namespace ofec {
	namespace DMOP {
		T1::T1(param_map& v) : T1(v.at("problem name"), v.at("number of variables")) { //param_numDim = 20 is suggested

		}
		T1::T1(const std::string& name, size_t size_var) : problem(name, size_var, 2), DMOPs(name, size_var, 2) {

		}

		void T1::initialize() {//all variables are in（0,1）,need revise
			setInitialDomain(0., 1.);
			setDomain(0., 1.);

			//set_duration_gen(global::ms_arg.find("durationGeneration")->second);
			set_change_severity(global::ms_arg.find("changeSeverity")->second);
			set_change_fre(global::ms_arg.find("changeFre")->second);

			generateAdLoadPF();
			m_initialized = true;
		}

		void T1::generateAdLoadPF() {

		}

		//the number of variables are dynamic
		int T1::evaluateObjective(Real* x, std::vector<Real>& obj) {
			Real t = get_time();
			if (time_changed() && t != 0. && (!get_updated_state())) {//防止不计数评价重复更新问题和重复采样PF
				m_optima.reset(new Optima<>());
				generateAdLoadPF();
				set_updated_state(true);
			}
			else if (m_evaluations % (size_t)global::ms_arg.find("changeFre")->second != 0) {
				set_updated_state(false);
			}
			Real d1_t = floor(m_number_variables * fabs(sin(t)));
			Real d2_t = floor(m_number_variables * fabs(pow(cos(2 * t), 3)));

			for (size_t m = 0; m < m_number_objectives; ++m)
			{
				if (m < 1)
				{
					obj[m] = 0;
					for (size_t i = 0; i < d1_t; i++) {
						obj[m] += pow(x[i], 2) - 10 * cos(2 * OFEC_PI * x[i]) + 10;
					}
				}
				else
				{
					obj[m] = pow(x[0] - 1, 2);
					for (size_t i = 1; i < d2_t; i++) {
						obj[m] += pow(pow(x[i], 2) - x[i - 1], 2);
					}
				}
			}
			return kNormalEval;
		}
	}
}