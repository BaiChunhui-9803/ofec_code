#ifndef OFEC_FREE_PEAKS_X_ASSYMETRIX_H
#define OFEC_FREE_PEAKS_X_ASSYMETRIX_H

#include "transform_base.h"
//#include "map_x_base.h"

namespace ofec::free_peaks {
	class MapXAssymetrix : public TransformBase {
	private:
		Real m_belta = 0.1;
	public:
		MapXAssymetrix(Problem *pro, const std::string& subspace_name, const ParameterMap& param, Random *rnd);
		void transfer(std::vector<Real>& x, const std::vector<Real>& var) override;
		virtual void bindData() override;
	};
}

#endif // !OFEC_FREE_PEAKS_LINEAR_MAP_H
