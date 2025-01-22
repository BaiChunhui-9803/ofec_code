#include "map_x_rotation.h"
#include "../../free_peaks.h"

namespace ofec::free_peaks {
	MapXRotation::MapXRotation(Problem *pro, const std::string &subspace_name, const ParameterMap &param, Random *rnd) :
		TransformBase(pro, subspace_name, param, rnd) 
	{
		auto generation = param.get<std::string>("generation");
		if (generation == "random") {
			int numVar = pro->numberVariables();
			m_rotation_matrix.resize(numVar, numVar);

			auto &uniform = pro->random()->uniform;
			auto &normal = pro->random()->normal;

			m_rotation_matrix.randomize(&uniform);
			m_rotation_matrix.generateRotationClassical(&normal, 1.0);
		}
	}

	void MapXRotation::transfer(std::vector<Real> &x_, const std::vector<Real> &var) {
		std::vector<double> x(x_);
		for (int i = 0; i < x.size(); i++) {
			x[i] = 0;

			for (int j = 0; j < x.size(); j++) {
				x[i] += m_rotation_matrix[j][i] * x_[j];
			}
		}
		swap(x, x_);
	}

	void MapXRotation::bindData() {

	}
}