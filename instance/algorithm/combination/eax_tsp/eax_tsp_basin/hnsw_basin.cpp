#include "hnsw_basin.h"

void ofec::HNSWbasin::initialize(const nbn::HnswModel& model, const std::vector<int>& belongBasinId,
	const std::vector<std::vector<int>>& basinIds
	)
{
	m_model.copy(model);
	m_belongBasinId = belongBasinId;
	m_basinIds = basinIds;
}

int ofec::HNSWbasin::calculateBasinId(SolutionBase& sol)const
{
	using namespace nbn;
	std::vector<Info> neighbors;
	m_model.searchBySolution(sol, neighbors, m_numNeighbors);
	auto& sols = m_model.solutions();

	int belongId = -1;

	Info neighborInfo(-1, std::numeric_limits<Real>::max());
	for (auto& it : neighbors) {
		if (sols[it.nodeId()]->fitness() > sol.fitness()) {
			if (neighborInfo.value() > it.value()) {
				neighborInfo = it;
			}
		}
	}

	if (neighborInfo.nodeId() == -1) {

		for (auto& it : neighbors) {
			if (neighborInfo.value() > it.value()) {
				neighborInfo = it;
			}
		}
	}

	belongId = m_belongBasinId[neighborInfo.nodeId()];
	return belongId;
	//return belongId == curBasinId;
}
