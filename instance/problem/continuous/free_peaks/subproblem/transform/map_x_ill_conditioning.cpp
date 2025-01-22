#include "map_x_ill_conditioning.h"
#include "../../free_peaks.h"

namespace ofec::free_peaks {
	MapXIllConditioning::MapXIllConditioning(Problem *pro, const std::string &subspace_name, const ParameterMap &param, Random *rnd) :
		TransformBase(pro, subspace_name, param, rnd)
	{
		auto generation = param.get<std::string>("generation");
		if (generation == "random") {
			Real max = param.get<Real>("max_condition");
			Real min = param.get<Real>("min_condition");
			m_condition = rnd->uniform.nextNonStd<Real>(min, max);
		}
		else if (generation == "specified") {
			m_condition = param.get<Real>("condition");
		}
		int numVar = pro->numberVariables();
		m_rotation_matrix.resize(numVar, numVar);
		m_rotation_matrix.randomize(&rnd->uniform);
		m_rotation_matrix.generateRotationClassical(&rnd->normal, m_condition);
	}

	void MapXIllConditioning::transfer(std::vector<Real> &x_, const std::vector<Real> &var) {
		std::vector<double> x(x_);
		for (int i = 0; i < x.size(); i++) {
			x[i] = 0;

			for (int j = 0; j < x.size(); j++) {
				x[i] += m_rotation_matrix[j][i] * x_[j];
			}
		}
		swap(x, x_);
	}

	void MapXIllConditioning::bindData() {

	}
}