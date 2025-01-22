#ifndef OFEC_FREE_PEAKS_WEIERSTRASS_RUG_H
#define OFEC_FREE_PEAKS_WEIERSTRASS_RUG_H

#include "transform_base.h"

#include "../../../../../../instance/problem/continuous/single_objective/global/bbob/bbob.h"


namespace ofec::free_peaks {
	class WeierstrassRugenessMap : public TransformBase, public BBOB {
	private:
		std::vector<std::pair<Real, Real>> m_to, m_from;
		double m_alpha = 0.1;
	public:
		WeierstrassRugenessMap(Problem *pro, const std::string& subspace_name, const ParameterMap& param);
		void transfer(std::vector<Real>& obj, const std::vector<Real>& var) override;
		void bindData() override;
	};
}

#endif // !OFEC_FREE_PEAKS_LINEAR_MAP_H
