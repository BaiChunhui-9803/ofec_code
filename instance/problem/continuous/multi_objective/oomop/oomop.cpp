#include"oomop.h"

namespace ofec {
	void OOMOP::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;
		resizeObjective(v.get<int>("number of objectives"));
		for (size_t i = 0; i < numberObjectives(); ++i) {
			m_optimize_mode[i] = OptimizeMode::kMinimize;
		}
		//resizeVariable(v.get<int>("number of variables"));//recomend n=10
	}

	std::vector<std::pair<size_t, bool>> OOMOP::attachPS(std::vector<Real>& sol) {
		std::vector<std::pair<size_t, bool>> attached_flag;
		std::pair<size_t, bool> flag;
		for (size_t i = 0; i < m_ps_shapes.size(); ++i) {
			flag.first = i;
			flag.second = m_ps_shapes[i]->if_in_range(sol);
			attached_flag.emplace_back(flag);
		}
		return attached_flag;
	}

	bool OOMOP::checkBoundary(std::vector<std::pair<Real, Real>>& dim_bounds, std::vector<Real>& point) {
		for (size_t i = 0; i < dim_bounds.size(); ++i) {
			if (point[i]<dim_bounds[i].first || point[i] > dim_bounds[i].second) {
				return false;
			}
		}
		return true;
	}

	void OOMOP::sample2HighDim(std::vector<std::vector<Real>>& low_dim_samples, Real offset, size_t sample_num, size_t dim_index) {
		Real step = offset / sample_num;
		size_t pre_num = low_dim_samples.size();
		for (size_t i = 0; i < pre_num; ++i) {
			auto temp_point0 = low_dim_samples[i];
			auto temp_point = low_dim_samples[i];
			for (size_t j = 1; j < sample_num + 1; ++j) {
				temp_point[dim_index] = temp_point0[dim_index] + j * step;
				low_dim_samples.push_back(temp_point);
			}
			temp_point[dim_index] = temp_point0[dim_index] + offset;
			low_dim_samples.push_back(temp_point);
		}
	}
}