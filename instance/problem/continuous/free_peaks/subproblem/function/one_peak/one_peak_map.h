#ifndef OFEC_FREE_PEAKS_ONE_PEAK_BBOB_RUG_H
#define OFEC_FREE_PEAKS_ONE_PEAK_BBOB_RUG_H

#include "one_peak_base.h"
#include "../../../../../../../instance/problem/continuous/single_objective/global/bbob/bbob.h"
namespace ofec::free_peaks {
	class OnePeakBBOB_Rug : public OnePeakBase, public BBOB {

	protected:
		std::unique_ptr<OnePeakBase> m_peak_fun;
	public:
		OnePeakBBOB_Rug(Problem *pro, const std::string& subspace_name, const ParameterMap& param);
		Real evaluate_(Real dummy, size_t var_size) override;
	};
}

#endif // !OFEC_FREE_PEAKS_ONE_PEAK_S1_H
