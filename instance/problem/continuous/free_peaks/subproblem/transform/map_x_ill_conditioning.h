#ifndef OFEC_FREE_PEAKS_X_ILL_CONDITIONING_H
#define OFEC_FREE_PEAKS_X_ILL_CONDITIONING_H

#include "transform_base.h"
#include "../../../../../../utility/linear_algebra/matrix.h"

namespace ofec::free_peaks {
	class MapXIllConditioning : public TransformBase {
	private:
		Matrix m_rotation_matrix;
		int m_idx = 0;
		Real m_condition = 1.0;
	public:
		MapXIllConditioning(Problem *pro, const std::string& subspace_name, const ParameterMap& param, Random *rnd);
		void transfer(std::vector<Real>& x, const std::vector<Real>& var) override;
		void bindData() override;
	};
}

#endif // !OFEC_FREE_PEAKS_LINEAR_MAP_H
