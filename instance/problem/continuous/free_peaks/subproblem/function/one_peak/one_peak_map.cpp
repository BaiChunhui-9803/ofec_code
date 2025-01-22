#include "one_peak_map.h"



namespace ofec::free_peaks {
	OnePeakBBOB_Rug::OnePeakBBOB_Rug(Problem *pro, const std::string& subspace_name, const ParameterMap& param) :
		OnePeakBase(pro, subspace_name, param) {

		m_condition_number = 100.;
		m_F0 = 0.;
		m_aK.resize(12);
		m_bK.resize(12);
		loadRotation(1. / sqrt(m_condition_number));
		for (size_t i = 0; i < 12; ++i) /// number of summands, 20 in CEC2005, 10/12 saves 30% of time
		{
			m_aK[i] = pow(0.5, (Real)i);
			m_bK[i] = pow(3., (Real)i);
			m_F0 += m_aK[i] * cos(2 * OFEC_PI * m_bK[i] * 0.5);
		}
	}

	Real OnePeakBBOB_Rug::evaluate_(Real dummy, size_t var_size) {




		return m_height - dummy;
	}
}