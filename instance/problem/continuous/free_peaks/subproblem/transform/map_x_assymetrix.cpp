#include "map_x_assymetrix.h"
#include "../../free_peaks.h"

namespace ofec::free_peaks {
	MapXAssymetrix::MapXAssymetrix(Problem *pro, const std::string &subspace_name, const ParameterMap &param, Random *rnd) :
		TransformBase(pro, subspace_name, param, rnd)
	{
		auto generation = param.get<std::string>("generation");
		if (generation == "random") {
			Real max = param.get<Real>("max_belta");
			Real min = param.get<Real>("min_belta");
			m_belta = rnd->uniform.nextNonStd<Real>(min, max);
		}
		else if (generation == "specified") {
			if (param.has("belta")) {
				m_belta = param.get<Real>("belta");
			}
		}
	}

	void MapXAssymetrix::transfer(std::vector<Real>& x_, const std::vector<Real>& var) {
		if (x_.size()==1) return;

		for (int i = 0; i < x_.size(); ++i) {
			if (x_[i] > 0) {
				x_[i] = pow(x_[i], 1. + m_belta * (i / (x_.size() - 1)) * log(x_[i] + 1));
			}
		}
	}

	void MapXAssymetrix::bindData() {

	}
}