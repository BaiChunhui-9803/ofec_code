#include "DMOPs.h"

namespace ofec {
	void DMOPs::initialize_() {
		auto& v = *m_param;
		Continuous::initialize_();
		resizeObjective(v.get<int>("number of objectives"));
		resizeVariable(v.get<int>("number of variables"));
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		m_optimize_mode[1] = OptimizeMode::kMinimize;
		set_change_fre(v.get<int>("change frequency"));
		set_change_severity(v.get<Real>("change severity"));
		//setInitialDomain(0., 1.);
		setDomain(0., 1.);
	}

	Real DMOPs::get_time() {
		size_t eval = get_change_fre();
		Real severity = get_change_severity();
		size_t effective_eval = get_effective_eval();
		return (effective_eval / eval) / severity;
	}

	bool DMOPs::time_changed() {
		size_t eval = get_change_fre();
		size_t effective_eval = get_effective_eval();
		if (effective_eval % eval == 0)
			return true;
		else
			return false;
	}
}