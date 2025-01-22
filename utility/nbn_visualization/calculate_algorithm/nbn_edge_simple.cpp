#include "nbn_edge_simple.h"
#include "../function/custom_function.h"
#include "../../core/problem/solution.h"

#include "../../instance/problem/combination/travelling_salesman/travelling_salesman.h"

void ofec::NBN_EdgeSimpleDivison::updateSol(int solId) {
	using namespace ofec;
	//m_sols[solId] = sols[solId];
	//	m_eval_fun(*m_sols[solId], m_id_pro);

	ofec::Real pos = (m_problem->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
	m_sols[solId]->setFitness(pos * m_sols[solId]->objective(0));

	TravellingSalesman::transferEdgeSol(*m_sols[solId], m_sols_edges[solId]);
	//dynamic_cast<TravellingSalesman&>(m_problem).transferEdgeSol(*m_sols[solId], m_sols_edges[solId]);
	//	GET_TSP(m_id_pro)->transferEdgeSol(*m_sols[solId], m_sols_edges[solId]);
	m_fitness[solId] = m_sols[solId]->fitness();
	//m_insertId[solId] = solId;
}


void ofec::NBN_EdgeSimpleDivison::updateSols() {
	int num_task = std::thread::hardware_concurrency();
	int num_samples = m_sols.size();
	std::vector<std::thread> thrds;
	std::vector<int> tasks;
	UTILITY::assignThreads(num_samples, num_task, tasks);
	std::pair<int, int> from_to;
	for (size_t i = 0; i < num_task; ++i) {
		from_to.first = tasks[i];
		from_to.second = tasks[i + 1];
		thrds.push_back(std::thread(
			&NBN_EdgeSimpleDivison::updateSolsFromTo, this,
			tasks[i], tasks[i + 1]));
	}
	for (auto& thrd : thrds)
		thrd.join();
}

void ofec::NBN_EdgeSimpleDivison::udpateNBNmultithread() {
	int num_task = m_numTaskUpdateNetwork;

	std::vector<std::vector<Node>> networkBetter(num_task);
	std::vector<std::vector<Node>> networkBetterEqual(num_task);
	{

		std::vector<std::thread> thrds;
		for (size_t i = 0; i < num_task; ++i) {
			thrds.push_back(std::thread(
				&NBN_EdgeSimpleDivison::updateNetWorkLoop, this,
				m_numLoop,std::ref(networkBetter[i]),std::ref(networkBetterEqual[i])));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}
	{
		int num_task = std::thread::hardware_concurrency();
		int num_samples = m_sols.size();
		std::vector<std::thread> thrds;
		std::vector<int> tasks;
		UTILITY::assignThreads(num_samples, num_task, tasks);
		std::pair<int, int> from_to;
		for (size_t i = 0; i < num_task; ++i) {
			from_to.first = tasks[i];
			from_to.second = tasks[i + 1];
			thrds.push_back(std::thread(
				&NBN_EdgeSimpleDivison::mergeNBNmultithread, this,
				tasks[i], tasks[i + 1], std::cref(networkBetter),std::cref(networkBetterEqual)));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}


}

void ofec::NBN_EdgeSimpleDivison::mergeNBNmultithread(
	int from,int to,
	const std::vector<std::vector<Node>>& networkBetter,
const std::vector<std::vector<Node>>& networkBetterEqual
) {
	for (int idx(from); idx < to; ++idx) {
		for (auto& curNetwork : networkBetter) {
			updateSolNetworkInfo(idx, curNetwork[idx].m_belong, curNetwork[idx].m_dis2parent,
				m_solInfoBetter);
		}
		for (auto& curNetwork : networkBetterEqual) {
			updateSolNetworkInfo(idx, curNetwork[idx].m_belong, curNetwork[idx].m_dis2parent,
				m_solInfoBetterEqual);
		}
	}
}


void ofec::NBN_EdgeSimpleDivison::updateNetWorkLoop(int numLoop, std::vector<Node>& networkBetter, std::vector<Node>& networkBetterEqual) {
	
	initNetwork(networkBetter);
	initNetwork(networkBetterEqual);
	while (numLoop--) {
		updateNetwork(networkBetter, networkBetterEqual);
	}
}


void ofec::NBN_EdgeSimpleDivison::updateNetwork(std::vector<Node>& networkBetter, std::vector<Node>& networkBetterEqual) {
	using namespace ofec;
	//initNetwork(networkBetter);
	//initNetwork(networkBetterEqual);
	std::vector<int> sortedEdgeId(m_dim, 0);
	for (int idx(0); idx < sortedEdgeId.size(); ++idx) {
		sortedEdgeId[idx] = idx;
	}
	m_random->uniform.shuffle(sortedEdgeId.begin(), sortedEdgeId.end());

	std::vector<int> solDirection(m_sols.size(), 0);
	for (int idx(0); idx < solDirection.size(); ++idx) {
		solDirection[idx] = m_random->uniform.nextNonStd<int>(0, 2);
	}

	//int sortedEdgeIdFrom(0);

	std::vector<int> divisons(m_sols.size());
	for (int idx(0); idx < divisons.size(); ++idx) {
		divisons[idx] = idx;
	}


	updateDivision(divisons,
		sortedEdgeId,
		0,
		solDirection,
		networkBetter, networkBetterEqual);

}



int ofec::NBN_EdgeSimpleDivison::updateDivision(const std::vector<int> & totalSol,
	std::vector<int>& sortedEdgeId,
	int sortedEdgeIdFrom, 
	const std::vector<int>& solDirection,
	std::vector<Node>& networkBetter, std::vector<Node>& networkBetterEqual
	) {
	std::vector<int> edgeId2DivisonId(m_dim + 1, -1);
	//std::map<int, int> edgeId2DivisonId;
	std::vector<int> divisionId2EdgeId;
	std::vector<std::vector<int>> divisons;
	std::vector<int> representativeId;

	for (auto& idx:totalSol) {
		auto& curEdgeId = m_sols_edges[idx][sortedEdgeId[sortedEdgeIdFrom]][solDirection[idx]];
		int divisionId(-1);
		if (edgeId2DivisonId[curEdgeId] == -1) {
			divisionId = edgeId2DivisonId[curEdgeId] = divisionId2EdgeId.size();
			divisionId2EdgeId.push_back(curEdgeId);
			divisons.emplace_back(std::vector<int>());
			representativeId.push_back(-1);
		}
		else divisionId = edgeId2DivisonId[curEdgeId];
		divisons[divisionId].push_back(idx);
	}

	for (int divisionId(0); divisionId<divisons.size(); ++divisionId) {
		if (divisons[divisionId].size() == 1) {
			representativeId[divisionId] = divisons[divisionId].front();
		}
		else if (divisons[divisionId].size() == 2) {


			int solIdx = divisons[divisionId].front();
			int solIdy = divisons[divisionId].back();
			//double curDis = ofec::TravellingSalesman::tspNorDis(m_sols_edges[solIdx], m_sols_edges[solIdy]);
			representativeId[divisionId] = updateRelation(solIdx, solIdy
				, networkBetter, networkBetterEqual);
		}
		else {
			//m_random->uniform.shuffle(sortedEdgeId.begin()+ sortedEdgeIdFrom + 1, sortedEdgeId.end());
			representativeId[divisionId] = updateDivision(divisons[divisionId],
				sortedEdgeId,
				sortedEdgeIdFrom + 1,
				solDirection,
				networkBetter, networkBetterEqual);
		}
	}


	
	m_random->uniform.shuffle(representativeId.begin(), representativeId.end());
	int bestId(representativeId.front());
	for (int idx(1); idx < representativeId.size(); ++idx) {
		bestId= updateRelation(bestId, representativeId[idx]
			, networkBetter, networkBetterEqual);
	}
	return bestId;
	//return representativeId[divisionId];
}

int ofec::NBN_EdgeSimpleDivison::updateRelation(int idx, int idy
	, std::vector<Node>& networkBetter, std::vector<Node>& networkBetterEqual) {
	using namespace ofec;
	double curDis = ofec::TravellingSalesman::tspNorDis(m_sols_edges[idx], m_sols_edges[idy]);
	int returnId(-1);
	bool updateBetterX = false, updateBetterY = false;
	bool updateBetterEqualX = false, updateBetterEqualY = false;
	if (m_fitness[idx] == m_fitness[idy]) {
		updateBetterEqualX = true;
		updateBetterEqualY = true;
		if (m_random->uniform.next() < 0.5) {
			returnId = idx;
		}
		else returnId = idy;
	}
	else if (m_fitness[idx] < m_fitness[idy]) {
		updateBetterX = true;
		updateBetterEqualX = true;
		returnId = idy;
	}
	else {
		updateBetterY = true;
		updateBetterEqualY = true;
		returnId = idx;
	}

	if (updateBetterX) {
		updateSolNetworkInfo(idx, idy, curDis, networkBetter);
	}
	if (updateBetterEqualX) {
		updateSolNetworkInfo(idx, idy, curDis, networkBetterEqual);
	}

	if (updateBetterY) {
		updateSolNetworkInfo(idy, idx, curDis, networkBetter);
	}
	if (updateBetterEqualY) {
		updateSolNetworkInfo(idy, idx, curDis, networkBetterEqual);
	}

	return returnId;
}

int ofec::NBN_EdgeSimpleDivison::updateSolNetworkInfo(int idx, int belong,
	double curDis,
	std::vector<Node>& network) {
	if (network[idx].m_dis2parent > curDis) {
		network[idx].m_dis2parent = curDis;
		network[idx].m_belong = belong;
	}
	else if (network[idx].m_dis2parent == curDis && m_random->uniform.next() < 0.5) {
		network[idx].m_belong = belong;
	}
	return idx;
}

