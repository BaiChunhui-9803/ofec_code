#include "nbn_edge_division.h"
#include "../function/custom_function.h"
#include "../../core/problem/solution.h"

#include "../../instance/problem/combination/travelling_salesman/travelling_salesman.h"


void ofec::NBN_EdgeDivision::initialize_(bool flag_grid_sample )
{
	using namespace ofec;
	int dim = m_pro->numVariables();
	int powDim = 0;
	double divSolNum(m_maxSample);
	m_divide_times = 0;
	m_dim = dim;

	while (divSolNum > m_divide_threadhold ) {
		divSolNum /= double(dim);
		++m_divide_times;
	}


	m_edge_seq.resize(m_dim);
	for (int idx(0); idx < m_dim; ++idx) {
		m_edge_seq[idx] = idx;
	}

	updateNetwork();
	generateSols();
	updateDivision();
	

	//for (int loop(0); loop < m_num_division; ++loop) {
	//	updateDivision();
	//}

	//int stop = -1;
}

void ofec::NBN_EdgeDivision::solToIdx(int solId, int &nodeId) {
	auto& cursol(m_sols_edges[solId]);
	nodeId = 0;
	int edgeId(m_sol_direction[solId]);
	for (int idx(0); idx < m_divide_times; ++idx) {
		auto& curEdge(m_edge_seq[idx]);
		//edgeId = m_random->uniform.nextNonStd<int>(0, 2);
		nodeId = nodeId* m_dim + cursol[curEdge][edgeId];
	}

}


void ofec::NBN_EdgeDivision::idxToVec(int idx, std::vector<int>& cur) {

	for (int idDim(cur.size()-1); idDim >= 0; --idDim) {
		cur[idDim] = idx % m_dim;
		idx /=  m_dim;
	}

}
void ofec::NBN_EdgeDivision::vecToIdx(const std::vector<int>& cur, int& idx) {
	//for(int idDim)
	idx = 0;
	for (int idDim(0); idDim < cur.size(); ++idDim) {
		idx *= m_dim;
		idx += cur[idDim];
	}

}

void ofec::NBN_EdgeDivision::Cnm(std::vector<unsigned>& selectedIdx, unsigned cur,
	int curSize,
	int m, int from) {
	if (m == 0) {
		selectedIdx.push_back(cur);
		return;
	}
	if (from + m > curSize)return;
	for (int idx(from); idx < curSize; ++idx) {
		cur ^= (1 << idx);
		//cur[idx] = !cur[idx];
		Cnm(selectedIdx, cur, curSize ,m - 1, idx + 1);
		cur ^= (1 << idx);
		//cur[idx] = !cur[idx];
	}
}
void ofec::NBN_EdgeDivision::Cnm(std::vector<std::vector<bool>>& selectedIdx,
	std::vector<bool>& cur, int m,int from) {
	if (m == 0) {
		selectedIdx.push_back(cur);
		return;
	}
	if (from + m > cur.size())return;
	for (int idx(from); idx < cur.size(); ++idx) {
		cur[idx] = !cur[idx];
		Cnm(selectedIdx, cur, m - 1, idx + 1);
		cur[idx] = !cur[idx];
	}
	
}

void ofec::NBN_EdgeDivision::clearNetworkSol() {
	for (auto& it : m_network) {
		for (auto& it2 : it) {
			it2.m_best_solId = -1;
			it2.m_indis_ids.clear();
		}
	}
}

void ofec::NBN_EdgeDivision::updateNetwork()
{
	//m_random->uniform.shuffle(m_edge_seq.begin(), m_edge_seq.end());
	//std::vector<std::vector<NetworkNode>> network(m_divide_times + 2);
	auto& network(m_network);
	network.resize(m_divide_times + 1);

	for (int idx(0); idx < network.size(); ++idx) {
		for (int idy(0); idy < network[idx].size(); ++idy) {
			network[idx][idy].m_network_id = { idx,idy };
		}
	}
	

//	std::vector<int> selected_;
	std::vector<std::vector<unsigned>> selectedV(m_divide_times + 1);
	std::vector<std::vector<int>>  fromId(m_divide_times + 1);
	std::vector<int> sumNum(m_divide_times + 1);
	
	unsigned curVal((1<<m_divide_times) -1);
	//for (int idx(0); idx < m_divide_times; ++idx) {
	//	curVal |= (1 << idx);
	//}
	
	selectedV[m_divide_times].resize(1);
	selectedV[m_divide_times].front() = curVal;
	fromId[m_divide_times].resize(1, 0);
	sumNum[m_divide_times] = pow(m_dim, m_divide_times);
	//network[m_divide_times].resize(sumNum[m_divide_times]);

	for (int dividedId(m_divide_times-1); dividedId >= 1; --dividedId) {
		Cnm(selectedV[dividedId], curVal, m_divide_times, m_divide_times - dividedId, 0);
		sumNum[dividedId]= std::pow(m_dim, dividedId);
		fromId[dividedId].resize(selectedV[dividedId].size(),0);
		for (int idx(1); idx < fromId[dividedId].size(); ++idx) {
			fromId[dividedId][idx] = fromId[dividedId][idx-1] + sumNum[dividedId];
		}
		//
	}
	network[0].resize(1);
	for (int idx(1); idx <= m_divide_times; ++idx) {
		network[idx].resize(fromId[idx].size() * sumNum[idx]);
	}

	// divied  = m_divide_times;


	std::vector<int> vec(m_divide_times,0);
	for (int idx(0); idx < network[m_divide_times].size(); ++idx) {
		idxToVec(idx, vec);
		network[m_divide_times][idx].m_repre_edge = vec;
	}

	for (int idDiv(m_divide_times); idDiv > 1; --idDiv) {
		for (int idNode(0); idNode < network[idDiv].size(); ++idNode) {
			int curfromId = idNode / sumNum[idDiv];
			auto& curNode(network[idDiv][idNode]);
			for (int headFromId(0); headFromId < selectedV[idDiv - 1].size(); ++ headFromId) {
				if ((selectedV[idDiv][curfromId] & selectedV[idDiv - 1][headFromId]) == selectedV[idDiv - 1][headFromId]) {
					int index(0);
					auto& curV(selectedV[idDiv - 1][headFromId]);
					auto& curEdge(m_network[idDiv][idNode].m_repre_edge);
					for (int idx(0); idx < m_divide_times; ++idx) {
						if (curV & (1 << idx)) {
							index = index * m_dim + curEdge[idx];
						}
					}
					auto& headNode(network[idDiv - 1][fromId[idDiv - 1][headFromId] + index]);
					if (headNode.m_repre_edge.empty()) headNode.m_repre_edge = curEdge;
					curNode.m_parents.push_back(&headNode);
				}
			}
		}
	}

	
	for (auto& curNode : network[1]) {
		curNode.m_parents.push_back(&network[0].front());
	}

	
	
}

void ofec::NBN_EdgeDivision::generateSols()
{
	resize(m_maxSample);

	if (!m_flag_multiThread || m_pro->hasTag(ProTag::kCSIWDN)) {
		for (int idx(0); idx < m_maxSample; ++idx) {
			generateSol(idx, m_random.get());
		}
	}
	else {
		int num_task = std::thread::hardware_concurrency();
		int num_samples = m_sols.size();
		std::vector<std::thread> thrds;
		std::vector<int> tasks;
		UTILITY::assignThreads(num_samples, num_task, tasks);
		std::pair<int, int> from_to;
		for (size_t i = 0; i < num_task; ++i) {
			from_to.first = tasks[i];
			from_to.second = tasks[i + 1];
			double randomSeed(m_random->uniform.nextNonStd<double>(0.01, 1.0));
			thrds.push_back(std::thread(
				&NBN_EdgeDivision::generateSolsThreadTask, this,
				tasks[i], tasks[i + 1], randomSeed));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}
}
void ofec::NBN_EdgeDivision::generateSolsThreadTask(int from, int to, double seed) {
	Random rand(0.5);
	for (int idx(from); idx < to; ++idx) {
		generateSol(idx, &rand);
	}
	
}
void ofec::NBN_EdgeDivision::generateSol(int solId, Random *rnd)
{/*
 		std::vector<std::vector<std::array<int, 2>>> m_sols_edges;
		std::vector<std::shared_ptr<SolBase>> m_sols;
		std::vector<double> m_fitness;
		std::vector<int> m_belong;
		std::vector<double> m_dis2parent;
			std::vector<bool> m_flagOpt;
		*/
	using namespace ofec;
	m_sols[solId].reset(m_pro->createSolution());
	m_pro->initSolution(*m_sols[solId], rnd);
	m_eval_fun(*m_sols[solId], m_pro.get());
//	GET_TSP(m_id_pro)->transferEdgeSol()
	dynamic_cast<TravellingSalesman&>(*m_pro).transferEdgeSol(*m_sols[solId], m_sols_edges[solId]);
//	GET_TSP(m_id_pro)->transferEdgeSol(*m_sols[solId], m_sols_edges[solId]);
	m_fitness[solId] = m_sols[solId]->fitness();
	m_belong[solId] = solId;
	m_dis2parent[solId] = std::numeric_limits<double>::max();
	m_flagOpt[solId] = false;
	
	m_sol_direction[solId] = 0;
}



void ofec::NBN_EdgeDivision::divideByEdge( NetworkNode& node, int divDim,bool divFlag)
{
	if (node.m_indis_ids.empty())return;
	if (node.m_indis_ids.size() >= m_divide_threadhold&&divFlag) {
		std::vector<NetworkNode> treeDiv(m_dim);
		int edgeIdx = m_edge_seq[divDim];
		for (auto solId : node.m_indis_ids) {
			int edgeDi = m_sol_direction[solId];
			treeDiv[m_sols_edges[solId][edgeIdx][edgeDi]].m_indis_ids.push_back(solId);
		}
		for (auto& node : treeDiv) {
			divideByEdge(node, divDim+1, true);
		}
		node.m_indis_ids.clear();
		for (auto& it : treeDiv) {
			if(it.m_best_solId!=-1)
			node.m_indis_ids.push_back(it.m_best_solId);
		}
	}
	{
		int curdiv = std::max<int>(1,node.m_indis_ids.size() / m_edge_division);
		m_random->uniform.shuffle(node.m_indis_ids.begin(), node.m_indis_ids.end());
		std::vector<std::vector<int>> edgeDivision(curdiv);
		std::vector<int>  edgeDivisionBest;
		for (int nodeId(0); nodeId < node.m_indis_ids.size(); ++nodeId) {
			edgeDivision[nodeId % curdiv].push_back(node.m_indis_ids[nodeId]);
		}

		for (int divId(0); divId < curdiv; ++divId) {
			auto& space(edgeDivision[divId]);
			//auto curBest(edgeDivisionBest[divId]);
			auto curBest = udpateNeighbor(space);
			if (curBest != -1) edgeDivisionBest.push_back(curBest);
		}

		node.m_best_solId = udpateNeighbor(edgeDivisionBest);
	}
}

void ofec::NBN_EdgeDivision::updateDivision() {
	clearNetworkSol();
	m_random->uniform.shuffle(m_edge_seq.begin(), m_edge_seq.end());
	for (int solId(0); solId < m_sols.size(); ++solId) {
		m_sol_direction[solId] = m_random->uniform.nextNonStd<int>(0, 2);
	}
	int nodeId(0);
	auto& networkLevel(m_network[m_divide_times]);
	for (int solId(0); solId < m_sols.size(); ++solId) {
		solToIdx(solId, nodeId);
		networkLevel[nodeId].m_indis_ids.push_back(solId);
	}

	bool divFlag(true);
	for (int idDiv(m_divide_times); idDiv >= 0; --idDiv) {
		auto& curNetwork(m_network[idDiv]);
		for (auto& node : curNetwork) {
			divideByEdge(node, idDiv + 1, divFlag);
			if (node.m_best_solId != -1) {
				for (auto& parentNode : node.m_parents) {
					parentNode->m_indis_ids.push_back(node.m_best_solId);
				}
			}
		}
		divFlag = false;
	}

	for (auto& curNetwork : m_network) {
		for (auto& node : curNetwork) {
			if (node.m_best_solId == -1) continue;
			int parentId = m_belong[node.m_best_solId];
			for (auto& curSolId : node.m_indis_ids) {
				if (m_belong[curSolId] == curSolId) {
					if (m_fitness[curSolId] < m_fitness[parentId]) {
						double curDis = solDis(curSolId, parentId);
						m_belong[curSolId] = parentId;
						m_dis2parent[curSolId] = curDis;
					}
				}
			}
		}
	}


//	for(auto& node: d)
	

}
int ofec::NBN_EdgeDivision::udpateNeighbor(std::vector<int>& indis)
{
	if (indis.empty())return -1;

	std::sort(indis.begin(), indis.end(), [&](int a,int b) {
		return m_fitness[a] > m_fitness[b];
	});

	for (int idx(0); idx < indis.size(); ++idx) {
		int solIdx(indis[idx]);
		for (int idy(0); idy < idx; ++idy) {
			int solIdy(indis[idy]);
			if (m_fitness[solIdy] > m_fitness[solIdx]) {
				double curDis = solDis(solIdx, solIdy);
				if (m_dis2parent[solIdx] > curDis) {
					m_belong[solIdx] = solIdy;
					m_dis2parent[solIdx] = curDis;
				}
				else if (m_dis2parent[solIdx] == curDis && m_random->uniform.next() < 0.5) {
					m_belong[solIdx] = solIdy;					
				}
			}
		}
	}

	return indis.front();
}
double ofec::NBN_EdgeDivision::solDis(int solIdx, int solIdy) {
	auto& edgex(m_sols_edges[solIdx]);
	auto& edgey(m_sols_edges[solIdy]);
	double dis(0);
	for (int idv(0); idv < m_dim; ++idv) {
		for (auto& to1 : edgex[idv]) {
			for (auto& to2 : edgey[idv]) {
				if (to1 == to2) {
					dis += 1;
					break;
				}
			}
		}
	}
	return dis / double(m_dim * edgex.front().size());
}


size_t ofec::NBN_EdgeDivision::size() const
{
	return m_sols.size();
//return size_t();
}

void ofec::NBN_EdgeDivision::addSol(const SolBase& new_sol, int& belong_id, bool flag_opt , int popIter , int popSolId , int algId )
{

	int solId(m_sols.size());
	resize(m_sols.size() + 1);
	m_sols[solId].reset(m_pro->createSolution(new_sol));
	dynamic_cast<TravellingSalesman&>(*m_pro).transferEdgeSol(*m_sols[solId], m_sols_edges[solId]);

	m_fitness[solId] = m_sols[solId]->fitness();
	m_belong[solId] = solId;
	m_dis2parent[solId] = std::numeric_limits<double>::max();
	m_flagOpt[solId] = false;
	m_sol_direction[solId] = 0;

	belong_id = -1;
}

void ofec::NBN_EdgeDivision::getSharedNBN(
	std::vector<std::shared_ptr<SolBase>>& sols, 
	std::vector<double>& fitness, 
	std::vector<int>& belong, 
	std::vector<double>& dis2parent, 
	std::vector<bool>& flagOpt) const
{
	sols = m_sols;
	fitness = m_fitness;
	belong = m_belong;
	dis2parent = m_dis2parent;
	flagOpt = m_flagOpt;
}

void ofec::NBN_EdgeDivision::getSharedNBN(std::vector<std::shared_ptr<SolBase>>& sols, std::vector<double>& fitness, std::vector<int>& belong, std::vector<double>& dis2parent, std::vector<int> popIters, std::vector<int> popSolIds, std::vector<int> algIds) const
{

	sols = m_sols;
	fitness = m_fitness;
	belong = m_belong;
	dis2parent = m_dis2parent;
//	flagOpt = m_flagOpt;
}
