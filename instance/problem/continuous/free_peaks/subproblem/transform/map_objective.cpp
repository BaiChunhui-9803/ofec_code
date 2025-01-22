#include "map_objective.h"
#include "../../free_peaks.h"

namespace ofec::free_peaks {
	MapObjective::MapObjective(Problem *pro, const std::string &subspace_name, const ParameterMap &param, Random *rnd) :
		TransformBase(pro, subspace_name, param, rnd) {}
	
	void MapObjective::transfer(std::vector<Real> &obj, const std::vector<Real>& var) {
		if (obj[0] < m_from[0].first*10) {
			return;
		}
		for (size_t j = 0; j < m_pro->numberObjectives(); ++j) {
			if (m_from[j].second != m_from[j].first) {
				obj[j] = m_to[j].first + (m_to[j].second - m_to[j].first)
					* (obj[j] - m_from[j].first) / (m_from[j].second - m_from[j].first);
			}
		}
	}

	void MapObjective::bindData() {
		auto &subpro = CAST_FPs(m_pro)->subspaceTree().name_box_subproblem.at(m_subspace_name).second;
		m_to = m_from = subpro->function()->objRanges();
		for (size_t j = 0; j < m_to.size(); ++j)
			m_to[j].first = 0;
	}
}