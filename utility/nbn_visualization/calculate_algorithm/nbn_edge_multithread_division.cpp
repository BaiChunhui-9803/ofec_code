#include "nbn_edge_multithread_division.h"
#include "../../function/custom_function.h"
#include "../../../core/problem/solution.h"
#include "../../../instance/problem/combination/travelling_salesman/travelling_salesman.h"



void display(const std::vector<int>& ids) {
	for (auto& it : ids) {
		std::cout << it << "\t";
	}
	std::cout << std::endl;
}

void ofec::NBN_EdgeMultiThreadDivision::initialize_(bool flag_grid_sample)
{


	//m_edge_seq.resize(m_dim);
	//for (int idx(0); idx < m_dim; ++idx) {
	//	m_edge_seq[idx] = idx;
	//}
	if (m_sols.empty())return;
	generateSols(m_maxSample);
	updateNetwork();

	//updateDivision();
	//for (int loop(0); loop < m_num_division; ++loop) {
	//	updateDivision();
	//}
}

void ofec::NBN_EdgeMultiThreadDivision::solToIdx(int solId, 
	const std::vector<int>& edge_seq,
	const std::vector<int>& sol_direction,
	int& nodeId
) const{
	auto& cursol(m_sols_edges[solId]);
	nodeId = 0;
	int edgeId(sol_direction[solId]);
	for (int idx(0); idx < m_divide_times; ++idx) {
		auto& curEdge(edge_seq[idx]);
		//edgeId = m_random->uniform.nextNonStd<int>(0, 2);
		nodeId = nodeId * m_dim + cursol[curEdge][edgeId];
	}

}


void ofec::NBN_EdgeMultiThreadDivision::idxToVec(int idx, std::vector<int>& cur)const {

	for (int idDim(cur.size() - 1); idDim >= 0; --idDim) {
		cur[idDim] = idx % m_dim;
		idx /= m_dim;
	}

}
void ofec::NBN_EdgeMultiThreadDivision::vecToIdx(const std::vector<int>& cur, int& idx) const{
	idx = 0;
	for (int idDim(0); idDim < cur.size(); ++idDim) {
		idx *= m_dim;
		idx += cur[idDim];
	}
}

void ofec::NBN_EdgeMultiThreadDivision::Cnm(std::vector<unsigned>& selectedIdx, unsigned cur,
	int curSize,int m, int from) const{
	if (m == 0) {
		selectedIdx.push_back(cur);
		return;
	}
	if (from + m > curSize)return;
	for (int idx(from); idx < curSize; ++idx) {
		cur ^= (1 << idx);
		//cur[idx] = !cur[idx];
		Cnm(selectedIdx, cur, curSize, m - 1, idx + 1);
		cur ^= (1 << idx);
		//cur[idx] = !cur[idx];
	}
}
void ofec::NBN_EdgeMultiThreadDivision::Cnm(std::vector<std::vector<bool>>& selectedIdx,
	std::vector<bool>& cur, int m, int from) const{
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


void ofec::NBN_EdgeMultiThreadDivision::clearNetworkSols(std::vector<std::vector<NetworkNodeSolBatch>>& network_sols)const {
	for (auto& it : network_sols) {
		for (auto& it2 : it) {
			it2.m_best_solId = -1;
			it2.m_indis_ids.clear();
		}
	}
}


void ofec::NBN_EdgeMultiThreadDivision::updateNetwork()
{

	using namespace ofec;
	auto pro = m_environment->problem();
	


	int dim = pro->numberVariables();
	int powDim = 0;
	double divSolNum(m_sols.size());
	m_divide_times = 0;
	m_dim = dim;
	while (divSolNum > m_divide_threadhold) {
		divSolNum /= double(dim);
		++m_divide_times;
	}


	//m_random->uniform.shuffle(m_edge_seq.begin(), m_edge_seq.end());
	//std::vector<std::vector<NetworkNode>> network(m_divide_times + 2);


	//	std::vector<int> selected_;
	std::vector<std::vector<unsigned>> selectedV(m_divide_times + 1);
	std::vector<std::vector<int>>  fromId(m_divide_times + 1);
	std::vector<int> sumNum(m_divide_times + 1);

	unsigned curVal((1 << m_divide_times) - 1);
	//for (int idx(0); idx < m_divide_times; ++idx) {
	//	curVal |= (1 << idx);
	//}

	selectedV[m_divide_times].resize(1);
	selectedV[m_divide_times].front() = curVal;
	fromId[m_divide_times].resize(1, 0);
	sumNum[m_divide_times] = pow(m_dim, m_divide_times);
	//network[m_divide_times].resize(sumNum[m_divide_times]);

	for (int dividedId(m_divide_times - 1); dividedId >= 1; --dividedId) {
		Cnm(selectedV[dividedId], curVal, m_divide_times, m_divide_times - dividedId, 0);
		sumNum[dividedId] = std::pow(m_dim, dividedId);
		fromId[dividedId].resize(selectedV[dividedId].size(), 0);
		for (int idx(1); idx < fromId[dividedId].size(); ++idx) {
			fromId[dividedId][idx] = fromId[dividedId][idx - 1] + sumNum[dividedId];
		}
		//
	}
	auto& network(m_network);
	network.resize(m_divide_times + 1);
	network[0].resize(1);
	for (int idx(1); idx <= m_divide_times; ++idx) {
		network[idx].resize(fromId[idx].size() * sumNum[idx]);
	}
	for (int idx(0); idx < network.size(); ++idx) {
		for (int idy(0); idy < network[idx].size(); ++idy) {
			network[idx][idy].m_network_id = { idx,idy };
		}
	}
	// divied  = m_divide_times;


	std::vector<int> vec(m_divide_times, 0);
	for (int idx(0); idx < network[m_divide_times].size(); ++idx) {
		idxToVec(idx, vec);
		network[m_divide_times][idx].m_repre_edge = vec;
	}

	for (int idDiv(m_divide_times); idDiv > 1; --idDiv) {
		for (int idNode(0); idNode < network[idDiv].size(); ++idNode) {
			int curfromId = idNode / sumNum[idDiv];
			auto& curNode(network[idDiv][idNode]);
			for (int headFromId(0); headFromId < selectedV[idDiv - 1].size(); ++headFromId) {
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

void ofec::NBN_EdgeMultiThreadDivision::generateSols(int numSample)
{
	auto pro = m_environment->problem();

	int startedId = m_sols.size();
	resize(numSample+startedId);

	if (!m_flag_multiThread) {
		for (int idx(0); idx < m_maxSample; ++idx) {
			generateSol(idx+startedId, m_random.get());
		}
	}
	else {

		std::cout << "generate solutions by multithread" << std::endl;
		int num_task = m_numThread;
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
				&NBN_EdgeMultiThreadDivision::generateSolsThreadTask, this,
				startedId,
				tasks[i], tasks[i + 1], randomSeed));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}
}
void ofec::NBN_EdgeMultiThreadDivision::generateSolsThreadTask(int startedId, int from, int to, double seed) {
	Random rand(0.5);

	for (int idx(from); idx < to; ++idx) {
		generateSol(startedId+idx, &rand);
	}
	
}
void ofec::NBN_EdgeMultiThreadDivision::generateSol(int solId, Random *rnd)
{
	using namespace ofec;
	auto pro = m_environment->problem();
	m_sols[solId].reset(pro->createSolution());
	pro->initializeVariables(m_sols[solId]->variableBase(), rnd);
	m_eval_fun(*m_sols[solId], m_environment.get());
	//	GET_TSP(m_id_pro)->transferEdgeSol()
	
	TravellingSalesman::transferEdgeSol(*m_sols[solId], m_sols_edges[solId]);
	//dynamic_cast<TravellingSalesman&>(m_pro).transferEdgeSol(*m_sols[solId], m_sols_edges[solId]);
	//	GET_TSP(m_id_pro)->transferEdgeSol(*m_sols[solId], m_sols_edges[solId]);
	m_fitness[solId] = m_sols[solId]->fitness();
	m_belong[solId] = solId;
	m_dis2parent[solId] = 1.0;
	m_flagOpt[solId] = false;

	//m_sol_direction[solId] = 0;
}


void ofec::NBN_EdgeMultiThreadDivision::updateSol(int solId, 
	const std::vector<std::shared_ptr<SolutionBase>>& sols) {
	using namespace ofec;
	m_sols[solId] = sols[solId];
	m_eval_fun(*m_sols[solId], m_environment.get());

	//ofec::Real pos = (m_pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
	//m_sols[solId]->setFitness(pos * m_sols[solId]->objective(0));

	
	
	TravellingSalesman::transferEdgeSol(*m_sols[solId], m_sols_edges[solId]);
	//dynamic_cast<TravellingSalesman&>(m_pro).transferEdgeSol(*m_sols[solId], m_sols_edges[solId]);
	//	GET_TSP(m_id_pro)->transferEdgeSol(*m_sols[solId], m_sols_edges[solId]);
	m_fitness[solId] = m_sols[solId]->fitness();
	m_belong[solId] = solId;
	m_dis2parent[solId] = std::numeric_limits<double>::max();
	m_flagOpt[solId] = false;
}

void ofec::NBN_EdgeMultiThreadDivision::updateSolsThreadTask(
	int startedId, 
	const std::vector<std::shared_ptr<SolutionBase>>& sols, int from, int to) {
	for (int idx(from); idx < to; ++idx) {
		updateSol(idx+startedId, sols);
	}
}

void ofec::NBN_EdgeMultiThreadDivision::updateSols(const std::vector<std::shared_ptr<SolutionBase>>& sols) {
	int num_task = m_numThread;
	int startedId = m_sols.size();
	resize(sols.size()+startedId);
	int num_samples = m_sols.size();
	std::vector<std::thread> thrds;
	std::vector<int> tasks;
	UTILITY::assignThreads(num_samples, num_task, tasks);
	std::pair<int, int> from_to;
	for (size_t i = 0; i < num_task; ++i) {
		from_to.first = tasks[i];
		from_to.second = tasks[i + 1];
		thrds.push_back(std::thread(
			&NBN_EdgeMultiThreadDivision::updateSolsThreadTask, this,
			startedId,
			std::cref(sols),
			tasks[i], tasks[i + 1]));
	}
	for (auto& thrd : thrds)
		thrd.join();


	m_sortedIds.resize(m_sols.size());
	for (int idx(0); idx < m_sortedIds.size(); ++idx) {
		m_sortedIds[idx] = idx;
	}
	std::sort(m_sortedIds.begin(), m_sortedIds.end(), [&](int a,int  b) {
		return m_fitness[a] < m_fitness[b];
		});

	//display(m_sortedIds);
	
}


void ofec::NBN_EdgeMultiThreadDivision::calHashThreadTask(int from, int to,
	std::vector<unsigned>& hashVal, const utility::HashRandomValue& hash_table) {
	for (int idx(from); idx < to; ++idx) {
		hashVal[idx] = 0;

		for (int idEdge(0); idEdge < m_dim; ++idEdge) {
			hashVal[idx] ^= (hash_table.m_hash_randNum[idEdge] * hash_table.m_hash_randNum[m_sols_edges[idx][idEdge][1]]);
				//			Hash ^= Rand[t1->Id] * Rand[t2->Id];
		}
	}
}



void ofec::NBN_EdgeMultiThreadDivision::setCenterSol(int solId) {
	m_centerSolId = solId;
	TravellingSalesman::transferEdgeSol(*m_sols[solId], m_centerSolEdge);
	
}


void ofec::NBN_EdgeMultiThreadDivision::setCenterSol(const std::shared_ptr<ofec::SolutionBase>& centerSol) {
	//m_centerSolId = centerSol;
	TravellingSalesman::transferEdgeSol(*centerSol, m_centerSolEdge);

}




void ofec::NBN_EdgeMultiThreadDivision::filterSameSolutions(int from,int to) {
	utility::HashRandomValue hash_table;
	hash_table.initialize(m_dim, m_random.get(), 0, std::sqrt(std::numeric_limits<unsigned>::max()));
	
	
	
	std::vector<unsigned> hashVal(m_sols.size(),0);
	//cal hash_vals by multithread

	{
		int num_task = m_numThread;
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
				&NBN_EdgeMultiThreadDivision::calHashThreadTask, this,
				tasks[i], tasks[i + 1], std::ref(hashVal),std::cref(hash_table)));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}


	//std::vector<int> sortedId(m_sols.size());
	//std::vector<int> solIds(m_sols.size());
	std::vector<bool> flagSame(m_sols.size(), false);
	//for (int idx(0); idx < m_sols.size(); ++idx) {
	//	sortedId[idx] = idx;
	//}


	// filter the same solutions

	std::map<std::pair<int, unsigned>, int> mapper;
	std::vector<std::vector<int>> bucket_ids(m_sols.size());
	
	int bucket_size(0);
	std::pair<int, unsigned> cur_id;
	int bucket_id(0);
	for (int idx(0); idx < m_sols.size(); ++idx) {
		cur_id.first = m_fitness[idx];
		cur_id.second = hashVal[idx];
		if (mapper.find(cur_id) != mapper.end()) {
			bucket_id = mapper[cur_id];
		}
		else {
			mapper[cur_id]  = bucket_id = bucket_size++;
		}
		auto& neighbor = bucket_ids[bucket_id];
		for (auto& it : neighbor) {
			if (ofec::TravellingSalesman::tspDis(m_sols_edges[idx], m_sols_edges[it]) == 0) {
				flagSame[idx] = true;
				break;
			}
		}
		if (!flagSame[idx]) {
			neighbor.push_back(idx);
		}
		if (idx < from || idx >= to) {
			flagSame[idx] = false;
		}
	}
	
	int from_id(from);
	for (int idx(from); idx < m_sols.size(); ++idx) {
		if (!flagSame[idx] ) {
			if(from_id != idx) swapSol(from_id, idx);
			++from_id;
		}
	}
	resize(from_id);


	//std::sort(sortedId.begin(), sortedId.end(), [&](int a, int b) {
	//	if (m_fitness[a] == m_fitness[b]) {
	//		if (hashVal[a] == hashVal[b]) {

	//		}
	//		else return hashVal[a] < hashVal[b];
	//	}
	//	else return m_fitness[a] < m_fitness[b];
	//});

	

	
}


void ofec::NBN_EdgeMultiThreadDivision::divideByEdge(
	NetworkNodeSolBatch& node,
	std::vector<int>& belong,
	std::vector<double>& dis2parent,
	Random *rnd,
	const std::vector<int>& edge_seq,
	const std::vector<int>& sol_dir,
	int divDim, bool divFlag)const
{
	if (node.m_indis_ids.empty())return;
	if (node.m_indis_ids.size() >= m_divide_threadhold && divFlag) {
		std::vector<NetworkNodeSolBatch> treeDiv(m_dim);
		int edgeIdx = edge_seq[divDim];
		for (auto solId : node.m_indis_ids) {
			int edgeDi = sol_dir[solId];
			treeDiv[m_sols_edges[solId][edgeIdx][edgeDi]].m_indis_ids.push_back(solId);
		}
		{
			for (auto& node : treeDiv) {
				divideByEdge(node, belong, dis2parent, rnd, edge_seq, sol_dir,
					divDim + 1, divFlag);
			}
		}
		node.m_indis_ids.clear();
		for (auto& it : treeDiv) {
			if (it.m_best_solId != -1)
				node.m_indis_ids.push_back(it.m_best_solId);
		}
	}
	{
		int curdiv = std::max<int>(1, node.m_indis_ids.size() / m_edge_division);
		rnd->uniform.shuffle(node.m_indis_ids.begin(), node.m_indis_ids.end());
		std::vector<std::vector<int>> edgeDivision(curdiv);
		std::vector<int>  edgeDivisionBest;
		for (int nodeId(0); nodeId < node.m_indis_ids.size(); ++nodeId) {
			edgeDivision[nodeId % curdiv].push_back(node.m_indis_ids[nodeId]);
		}

		for (int divId(0); divId < curdiv; ++divId) {
			auto& space(edgeDivision[divId]);
			//auto curBest(edgeDivisionBest[divId]);
			auto curBest = udpateNeighbor(space, belong,dis2parent, rnd);
			if (curBest != -1) edgeDivisionBest.push_back(curBest);
		}

		node.m_best_solId = udpateNeighbor(edgeDivisionBest, belong, dis2parent, rnd);
	}
}

void ofec::NBN_EdgeMultiThreadDivision::updateDivisionTask(
	double seed,
	std::vector<int>& belong,
	std::vector<double>& dis2parent,
	int numLoop 
)const {
	Random randNum(seed);
	std::vector<std::vector<NetworkNodeSolBatch>> networkSols(m_network.size());
	for (int idx(0); idx < networkSols.size(); ++idx) {
		networkSols[idx].resize(m_network[idx].size());
	}
	std::vector<int> edge_seq(m_dim);
	for (int idx(0); idx < edge_seq.size(); ++idx) {
		edge_seq[idx] = idx;
	}
	std::vector<int> sol_direction(m_sols.size(),0);

	while (numLoop--) {
		clearNetworkSols(networkSols);
		randNum.uniform.shuffle(edge_seq.begin(), edge_seq.end());
		for (auto& it : sol_direction) {
			it = randNum.uniform.nextNonStd<int>(0, 2);
		}
		std::pair<int, int> nodeId = { m_divide_times,0 };
		//auto& networkLevel(networkSols[m_divide_times]);
		for (int solId(0); solId < m_sols.size(); ++solId) {
			solToIdx(solId,edge_seq,sol_direction,nodeId.second);
			networkSols[m_divide_times][nodeId.second].m_indis_ids.push_back(solId);
		}

		bool divFlag(true);
		for (int idDiv(m_divide_times); idDiv >= 0; --idDiv) {
			auto& curNetwork(m_network[idDiv]);
			for (auto& networkNode : curNetwork) {
				auto& node = networkSols[networkNode.m_network_id.front()][networkNode.m_network_id.back()];
				divideByEdge(node,belong,dis2parent,&randNum,edge_seq,sol_direction, idDiv + 1, divFlag);
				if (node.m_best_solId != -1) {
					for (auto& parentNetworkNode : networkNode.m_parents) {
						auto& parentNode=(networkSols[parentNetworkNode->m_network_id.front()][parentNetworkNode->m_network_id.back()]);
						parentNode.m_indis_ids.push_back(node.m_best_solId);
					}
				}
			}
		//	divFlag = false;
		}

		for (auto& curNetwork : m_network) {
			for (auto& networkNode : curNetwork) {
				auto& node = networkSols[networkNode.m_network_id.front()][networkNode.m_network_id.back()];
				if (node.m_best_solId == -1) continue;
				int parentId = m_belong[node.m_best_solId];
				for (auto& curSolId : node.m_indis_ids) {
					if (belong[curSolId] == curSolId) {
						if (m_fitness[curSolId] < m_fitness[parentId]) {
							double curDis = solDis(curSolId, parentId);
							belong[curSolId] = parentId;
							dis2parent[curSolId] = curDis;
						}
					}
				}
			}
		}
	}
	


}


void ofec::NBN_EdgeMultiThreadDivision::updateDivisionSubRegionTask(
	double seed,
	std::vector<int>& belong,
	std::vector<double>& dis2parent,
	int numLoop 
)const {

	ofec::Random objRnd(seed);

	std::vector<int> selectedEdge(m_dim);
	for (int idx(0); idx < selectedEdge.size(); ++idx) {
		selectedEdge[idx] = idx;
	}
	
	std::vector<int> totalSolIds;
	for (int idx(0); idx < belong.size(); ++idx) {
		totalSolIds.push_back(idx);
	}

	while (numLoop--) {
		objRnd.uniform.shuffle(selectedEdge.begin(), selectedEdge.end());
		updateDivisionTaskSubRegionDivisionByEdge(
			&objRnd, selectedEdge, 0, totalSolIds,
			belong, dis2parent
		);
	}
}


void ofec::NBN_EdgeMultiThreadDivision::updateDivisionSubRegionBandEdgeTask(
	double seed,
	std::vector<int>& belong,
	std::vector<double>& dis2parent,
	int numLoop 
)const {
	ofec::Random objRnd(seed);

	std::vector<int> selectedEdge(m_dim);
	for (int idx(0); idx < selectedEdge.size(); ++idx) {
		selectedEdge[idx] = idx;
	}

	std::vector<int> totalSolIds;
	for (int idx(0); idx < belong.size(); ++idx) {
		totalSolIds.push_back(idx);
	}

	while (numLoop--) {
		objRnd.uniform.shuffle(selectedEdge.begin(), selectedEdge.end());
		updateDivisionTaskSubRegionDivisionByEdgeBandEdge(
			&objRnd, selectedEdge, 0, totalSolIds,
			belong, dis2parent
		);
	}
}
int ofec::NBN_EdgeMultiThreadDivision::updateDivisionTaskSubRegionDivisionByEdge(
	ofec::Random* rnd,
	std::vector<int>& selectedEdgeIds,
	int from,
	const std::vector<int>& solIds,
	std::vector<int>& belong,
	std::vector<double>& dis2parent
)const {
	int edgeId = selectedEdgeIds[from];
	
	std::vector<std::vector<int>> division(m_dim);
	std::set<int> bestSolIds;
	for (auto& solId : solIds) {
		for (auto& solEdgeId : m_sols_edges[solId][edgeId]) {
			division[solEdgeId].push_back(solId);
		}
	}


	for (auto& subDivision : division) {
		if (!subDivision.empty()) {
			if (subDivision.size() <= m_numSubIndis) {
				bestSolIds.insert(
					udpateNeighbor(
						subDivision, belong, dis2parent,rnd
					));
			}
			else {

			/*	if (std::round(m_dim * (1.0 - m_err)) == from + 1) {
					int bestSolId = -1;
					double maxFit = -std::numeric_limits<double>::max();
					int curId = 0;
					int subLoop = m_numSubLoop;

					
					while (subLoop--) {
						rnd->uniform.shuffle(selectedEdgeIds.begin() + from + 1, selectedEdgeIds.end());
						curId = updateDivisionTaskSubRegionDivisionByEdge(
							rnd, selectedEdgeIds, from + 1,
							subDivision, belong, dis2parent
						);
						if (maxFit < m_fitness[curId]) {
							maxFit = m_fitness[curId];
							bestSolId = curId;
						}
					}
					bestSolIds.push_back(bestSolId);
				}
				else */{
					bestSolIds.insert(
						updateDivisionTaskSubRegionDivisionByEdge(
							rnd, selectedEdgeIds, from + 1,
							subDivision, belong, dis2parent
						));
				}

			}
		}
	}
	std::vector<int> curSolIds;
	for (auto& it : bestSolIds) {
		curSolIds.push_back(it);
	}
	return udpateNeighbor(curSolIds, belong, dis2parent, rnd);
}




int ofec::NBN_EdgeMultiThreadDivision::updateDivisionTaskSubRegionDivisionByEdgeBandEdge(
	ofec::Random* rnd,
	std::vector<int>& selectedEdgeIds,
	int from,
	const std::vector<int>& solIds,
	std::vector<int>& belong,
	std::vector<double>& dis2parent
)const {
	int edgeId = selectedEdgeIds[from];

	std::vector<std::vector<int>> division(m_dim);
	std::set<int> bestSolIds;
	std::vector<bool> feasibleEdge(m_dim,true);

	//std::vector<int> visitedSolId()
	for (auto& it: m_centerSolEdge[edgeId]) {
		feasibleEdge[it] = false;
	}
	for (auto& solId : solIds) {
		for (auto& solEdgeId : m_sols_edges[solId][edgeId]) {

			division[solEdgeId].push_back(solId);
			//if (feasibleEdge[solEdgeId]) {
			//	division[solEdgeId].push_back(solId);
			//}
		}
	}


	for (int edgeId(0); edgeId < division.size(); ++edgeId) {
		auto& subDivision = division[edgeId];
		if (feasibleEdge[edgeId]) {
			if (!subDivision.empty()) {
				if (subDivision.size() == 1) {
					bestSolIds.insert(subDivision.front());
				}
				else if (subDivision.size() <= m_numSubIndis) {
					int bestSolId = udpateNeighbor(
						subDivision, belong, dis2parent, rnd
					);
					bestSolIds.insert(bestSolId);


					if (bestSolId == -1) {
						std::cout << "error" << std::endl;
					}
				}
				else {
					int bestSolId = updateDivisionTaskSubRegionDivisionByEdgeBandEdge(
						rnd, selectedEdgeIds, from + 1,
						subDivision, belong, dis2parent
					);
					bestSolIds.insert(bestSolId);

					if (bestSolId == -1) {
						std::cout << "error" << std::endl;


						int bestSolId = updateDivisionTaskSubRegionDivisionByEdgeBandEdge(
							rnd, selectedEdgeIds, from + 1,
							subDivision, belong, dis2parent
						);
					}

				}
			}
		}
		else {
			if (!subDivision.empty()) {
		//		double maxFit = -std::numeric_limits<double>::max();
				int bestSolId(-1);
				bestSolId = subDivision.front();
				for (size_t idx(1); idx < subDivision.size(); ++idx) {
					auto& curSolId = subDivision[idx];
					if (compareSol(curSolId, bestSolId)) {
						bestSolId = curSolId;
					}
				}

				if (bestSolId == -1) {
					std::cout << "error" << std::endl;
				}

				updateSolsDis2Sol(subDivision, bestSolId,
					dis2parent, belong, rnd);
				//for (auto& curSolId : subDivision) {
				//	if (m_fitness[curSolId] > maxFit) {
				//		maxFit = m_fitness[curSolId];
				//		bestSolId = curSolId;
				//	}
				//	else if (m_fitness[curSolId] == maxFit && rnd->uniform.next() < 0.5) {
				//		bestSolId = curSolId;
				//	}
				//}
				bestSolIds.insert(bestSolId);
			}
		}
	}
	//for (auto& subDivision : division) {
	//	if (!subDivision.empty()) {
	//		if (subDivision.size() <= m_numSubIndis) {
	//			bestSolIds.insert(
	//				udpateNeighbor(
	//					subDivision, belong, dis2parent, rnd
	//				));
	//		}
	//		else {
	//			bestSolIds.insert(
	//				updateDivisionTaskSubRegionDivisionByEdge(
	//					rnd, selectedEdgeIds, from + 1,
	//					subDivision, belong, dis2parent
	//				));
	//		}
	//	}
	//}
	std::vector<int> curSolIds;
	for (auto& it : bestSolIds) {
		curSolIds.push_back(it);
	}
	return udpateNeighbor(curSolIds, belong, dis2parent, rnd);
}

int ofec::NBN_EdgeMultiThreadDivision::updateDivisionTaskSubRegionDivisionByEdgeOnCenterEdge(
	ofec::Random* rnd,
	std::vector<int>& selectedEdgeIds,
	int from,
	const std::vector<int>& sortedSolIds,
	std::vector<int>& belong,
	std::vector<double>& dis2parent
)const {
	if (from >=m_dim) {
		for (int idx(0); idx < sortedSolIds.size(); ++idx) {
			if (idx + 1 < sortedSolIds.size()) {
				int otherId = rnd->uniform.nextNonStd<int>(idx + 1, sortedSolIds.size());
				updateInfo(sortedSolIds[idx], sortedSolIds[otherId], belong, dis2parent, rnd);
			}

		}
		
		return 0;
	}


	std::vector<int> fileterSortedSolIds;

	{
		int edgeId = selectedEdgeIds[from];

		std::vector<bool> feasibleEdge(m_dim, false);
		int edgeTo = m_centerSolEdge[edgeId][rnd->uniform.nextNonStd<int>(0, 2)];

		std::vector<bool> activeSolFlag(sortedSolIds.size(), false);

		for (int idx(0); idx < activeSolFlag.size(); ++idx) {
			auto solId = sortedSolIds[idx];
			for (auto& solEdgeId : m_sols_edges[solId][edgeId]) {
				if (solEdgeId == edgeTo) {
					activeSolFlag[idx] = true;
					break;
				}
			}
		}

		int lastFilterId = -1;

		for (int idx(0); idx < activeSolFlag.size(); ++idx) {
			if (activeSolFlag[idx]) {
				fileterSortedSolIds.push_back(sortedSolIds[idx]);
				lastFilterId = idx;
			}
			else if (idx + 1 < activeSolFlag.size()) {
				int otherId = rnd->uniform.nextNonStd<int>(idx + 1, activeSolFlag.size());
				updateInfo(sortedSolIds[idx], sortedSolIds[otherId], belong, dis2parent, rnd);
			}

		}

		if (lastFilterId + 1 < activeSolFlag.size() && lastFilterId != -1) {
			int otherId = rnd->uniform.nextNonStd<int>(lastFilterId + 1, activeSolFlag.size());
			updateInfo(sortedSolIds[lastFilterId], sortedSolIds[otherId], belong, dis2parent, rnd);
		}

		//std::cout << "end line" << std::endl;
		//std::cout << "end line" << std::endl;
		//std::cout << "end line" << std::endl;

		//display(fileterSortedSolIds);

	}
	if (fileterSortedSolIds.size()<=1)return 0;
	
	return updateDivisionTaskSubRegionDivisionByEdgeOnCenterEdge(rnd, selectedEdgeIds,
		from + 1, fileterSortedSolIds, belong, dis2parent);

}


void ofec::NBN_EdgeMultiThreadDivision::updateDivisionTaskSubRegionDivisionByEdgeOnCenterEdgeTask(
	double seed,
	std::vector<int>& belong,
	std::vector<double>& dis2parent,
	int numLoop
)const {
	ofec::Random objRnd(seed);

	std::vector<int> selectedEdge(m_dim);
	for (int idx(0); idx < selectedEdge.size(); ++idx) {
		selectedEdge[idx] = idx;
	}

	//std::vector<int> totalSolIds;
	//for (int idx(0); idx < belong.size(); ++idx) {
	//	totalSolIds.push_back(idx);
	//}

	while (numLoop--) {
		objRnd.uniform.shuffle(selectedEdge.begin(), selectedEdge.end());
		updateDivisionTaskSubRegionDivisionByEdgeOnCenterEdge(
			&objRnd, selectedEdge, 0, m_sortedIds,
			belong, dis2parent
		);
	}
}

void ofec::NBN_EdgeMultiThreadDivision::updateDivision() {

	updateDivisionTask(
		m_random->uniform.next(),
		m_belong,
		m_dis2parent,
		m_num_loop
	);
	
	//clearNetworkSol();
	//m_random->uniform.shuffle(m_edge_seq.begin(), m_edge_seq.end());
	//for (int solId(0); solId < m_sols.size(); ++solId) {
	//	m_sol_direction[solId] = m_random->uniform.nextNonStd<int>(0, 2);
	//}
	//int nodeId(0);
	//auto& networkLevel(m_network[m_divide_times]);
	//for (int solId(0); solId < m_sols.size(); ++solId) {
	//	solToIdx(solId, nodeId);
	//	networkLevel[nodeId].m_indis_ids.push_back(solId);
	//}

	//bool divFlag(true);
	//for (int idDiv(m_divide_times); idDiv >= 0; --idDiv) {
	//	auto& curNetwork(m_network[idDiv]);
	//	for (auto& node : curNetwork) {
	//		divideByEdge(node, idDiv + 1, divFlag);
	//		if (node.m_best_solId != -1) {
	//			for (auto& parentNode : node.m_parents) {
	//				parentNode->m_indis_ids.push_back(node.m_best_solId);
	//			}
	//		}
	//	}
	//	divFlag = false;
	//}

	//for (auto& curNetwork : m_network) {
	//	for (auto& node : curNetwork) {
	//		if (node.m_best_solId == -1) continue;
	//		int parentId = m_belong[node.m_best_solId];
	//		for (auto& curSolId : node.m_indis_ids) {
	//			if (m_belong[curSolId] == curSolId) {
	//				if (m_fitness[curSolId] < m_fitness[parentId]) {
	//					double curDis = solDis(curSolId, parentId);
	//					m_belong[curSolId] = parentId;
	//					m_dis2parent[curSolId] = curDis;
	//				}
	//			}
	//		}
	//	}
	//}



}

void ofec::NBN_EdgeMultiThreadDivision::mergeBelongInfoThreadTask(
	double seed,
	int from, int to,
	const std::vector<std::vector<int>>& total_belong,
	const std::vector<std::vector<double>> total_dis2parent
	/*, std::vector<int>& belong,
	std::vector<double>& dis2parent*/) {
	Random randNum(seed);
	for (int idx(from); idx < to; ++idx) {

		auto& cur_dis = m_dis2parent[idx];
		auto& cur_belong = m_belong[idx];
		for (int idthrd(0); idthrd < total_belong.size(); ++idthrd) {
			if (cur_dis > total_dis2parent[idthrd][idx]) {
				cur_dis = total_dis2parent[idthrd][idx];
				cur_belong = total_belong[idthrd][idx];
			}
			else if (cur_dis == total_dis2parent[idthrd][idx]
				&& randNum.uniform.next() < 0.5) {
				//cur_dis = total_dis2parent[idthrd][idx];
				cur_belong = total_belong[idthrd][idx];
			}
		}
	}
}

void ofec::NBN_EdgeMultiThreadDivision::updateDivisionMultithread()
{
	std::vector<std::vector<int>> belong(m_numThread,m_belong);
	std::vector<std::vector<double>> dis2parent(m_numThread, m_dis2parent);

	int num_loop = m_num_loop;
	int num_task = m_numThread;

	{
		std::vector<std::thread> thrds;
		for (size_t i = 0; i < num_task; ++i) {
			double randomSeed(m_random->uniform.nextNonStd<double>(0.01, 1.0));
			thrds.push_back(std::thread(
				&NBN_EdgeMultiThreadDivision::updateDivisionTask, this,
				randomSeed,
				std::ref(belong[i]),
				std::ref(dis2parent[i]),
				num_loop));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}
	{
		std::vector<std::thread> thrds;
		std::vector<int> tasks;
		UTILITY::assignThreads(m_sols.size(), num_task, tasks);
		std::pair<int, int> from_to;
		for (size_t i = 0; i < num_task; ++i) {
			from_to.first = tasks[i];
			from_to.second = tasks[i + 1];
			double randomSeed(m_random->uniform.nextNonStd<double>(0.01, 1.0));
			thrds.push_back(std::thread(
				&NBN_EdgeMultiThreadDivision::mergeBelongInfoThreadTask,
				this, randomSeed,
				tasks[i], tasks[i + 1],
				std::cref(belong),std::cref(dis2parent)
				));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}


	
}


void ofec::NBN_EdgeMultiThreadDivision::updateDivisionSubRegionMultiThread() 
{
	std::vector<std::vector<int>> belong(m_numThread, m_belong);
	std::vector<std::vector<double>> dis2parent(m_numThread, m_dis2parent);

	int num_loop = m_num_loop;
	int num_task = m_numThread;

	{
		std::vector<std::thread> thrds;
		for (size_t i = 0; i < num_task; ++i) {
			double randomSeed(m_random->uniform.nextNonStd<double>(0.01, 1.0));
			thrds.push_back(std::thread(
				&NBN_EdgeMultiThreadDivision::updateDivisionSubRegionTask, this,
				randomSeed,
				std::ref(belong[i]),
				std::ref(dis2parent[i]),
				num_loop));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}
	{
		std::vector<std::thread> thrds;
		std::vector<int> tasks;
		// for test
	//	num_task = 1;
		UTILITY::assignThreads(m_sols.size(), num_task, tasks);
		std::pair<int, int> from_to;
		for (size_t i = 0; i < num_task; ++i) {
			from_to.first = tasks[i];
			from_to.second = tasks[i + 1];
			double randomSeed(m_random->uniform.nextNonStd<double>(0.01, 1.0));
			thrds.push_back(std::thread(
				&NBN_EdgeMultiThreadDivision::mergeBelongInfoThreadTask,
				this, randomSeed,
				tasks[i], tasks[i + 1],
				std::cref(belong), std::cref(dis2parent)
			));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}

}




void ofec::NBN_EdgeMultiThreadDivision::updateDivisionSubRegionBandEdgeMultiThread() {
	std::vector<std::vector<int>> belong(m_numThread, m_belong);
	std::vector<std::vector<double>> dis2parent(m_numThread, m_dis2parent);

	//std::vector<std::vector<int>> selectedEdge(m_numThread,std::vector<int>(m_dim));
	//for (auto& it : selectedEdge) {
	//	for (int idx(0); idx < m_dim; ++idx) {
	//		it[idx] = idx;
	//	}
	//}

	int num_loop = m_num_loop;
	int num_task = m_numThread;

	{
		std::vector<std::thread> thrds;
		// for test
	//	num_task = 1;
		//updateDivisionTaskSubRegionDivisionByEdgeBandEdge()

		for (size_t i = 0; i < num_task; ++i) {
			double randomSeed(m_random->uniform.nextNonStd<double>(0.01, 1.0));
			thrds.push_back(std::thread(
				&NBN_EdgeMultiThreadDivision::updateDivisionSubRegionBandEdgeTask, this,
			//	std::ref(selectedEdge[i]),
			//	0,
				randomSeed,
				std::ref(belong[i]),
				std::ref(dis2parent[i]),
				num_loop));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}
	{
		std::vector<std::thread> thrds;
		std::vector<int> tasks;

		UTILITY::assignThreads(m_sols.size(), num_task, tasks);
		std::pair<int, int> from_to;
		for (size_t i = 0; i < num_task; ++i) {
			from_to.first = tasks[i];
			from_to.second = tasks[i + 1];
			double randomSeed(m_random->uniform.nextNonStd<double>(0.01, 1.0));
			thrds.push_back(std::thread(
				&NBN_EdgeMultiThreadDivision::mergeBelongInfoThreadTask,
				this, randomSeed,
				tasks[i], tasks[i + 1],
				std::cref(belong), std::cref(dis2parent)
			));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}
}



void ofec::NBN_EdgeMultiThreadDivision::updateDivisionSubRegionOnCenteredEdgeMultiThread() {
	std::vector<std::vector<int>> belong(m_numThread, m_belong);
	std::vector<std::vector<double>> dis2parent(m_numThread, m_dis2parent);

	//std::vector<std::vector<int>> selectedEdge(m_numThread,std::vector<int>(m_dim));
	//for (auto& it : selectedEdge) {
	//	for (int idx(0); idx < m_dim; ++idx) {
	//		it[idx] = idx;
	//	}
	//}

	int num_loop = m_num_loop;
	int num_task = m_numThread;

	{
		std::vector<std::thread> thrds;
		// for test
		//num_task = 1;
		//updateDivisionTaskSubRegionDivisionByEdgeBandEdge()

		for (size_t i = 0; i < num_task; ++i) {
			double randomSeed(m_random->uniform.nextNonStd<double>(0.01, 1.0));
			thrds.push_back(std::thread(
				&NBN_EdgeMultiThreadDivision::updateDivisionTaskSubRegionDivisionByEdgeOnCenterEdgeTask, this,
				//	std::ref(selectedEdge[i]),
				//	0,
				randomSeed,
				std::ref(belong[i]),
				std::ref(dis2parent[i]),
				num_loop));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}
	{
		std::vector<std::thread> thrds;
		std::vector<int> tasks;

		UTILITY::assignThreads(m_sols.size(), num_task, tasks);
		std::pair<int, int> from_to;
		for (size_t i = 0; i < num_task; ++i) {
			from_to.first = tasks[i];
			from_to.second = tasks[i + 1];
			double randomSeed(m_random->uniform.nextNonStd<double>(0.01, 1.0));
			thrds.push_back(std::thread(
				&NBN_EdgeMultiThreadDivision::mergeBelongInfoThreadTask,
				this, randomSeed,
				tasks[i], tasks[i + 1],
				std::cref(belong), std::cref(dis2parent)
			));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}
}
int ofec::NBN_EdgeMultiThreadDivision::udpateNeighbor(
	std::vector<int>& indis,
	std::vector<int>& belong,
	std::vector<double>& dis2parent, 
	ofec::Random *rnd
	)const
{
	if (indis.empty())return -1;

	std::sort(indis.begin(), indis.end(), [&](int a, int b) {
		return compareSol(a, b);
	//	return m_fitness[a] > m_fitness[b];
	});

	for (int idx(0); idx < indis.size(); ++idx) {
		int solIdx(indis[idx]);
		for (int idy(0); idy < idx; ++idy) {
			int solIdy(indis[idy]);
			/*if (m_fitness[solIdy] > m_fitness[solIdx]) */{
				double curDis = solDis(solIdx, solIdy);
				if (dis2parent[solIdx] > curDis) {
					belong[solIdx] = solIdy;
					dis2parent[solIdx] = curDis;
				}
				else if (dis2parent[solIdx] == curDis && rnd->uniform.next() < 0.5) {
					belong[solIdx] = solIdy;
				}
			}
		}
	}
	return indis.front();
}

void ofec::NBN_EdgeMultiThreadDivision::updateInfo(int solIdx, int solIdy,
	std::vector<int>& belong,
	std::vector<double>& dis2parent,
	Random* rnd
)const {
	double curDis = solDis(solIdx, solIdy);
	if (dis2parent[solIdx] > curDis) {
		belong[solIdx] = solIdy;
		dis2parent[solIdx] = curDis;
	}
	else if (dis2parent[solIdx] == curDis && rnd->uniform.next() < 0.5) {
		belong[solIdx] = solIdy;
	}
}
double ofec::NBN_EdgeMultiThreadDivision::solDis(int solIdx, int solIdy) const{
	return ofec::TravellingSalesman::tspDis(m_sols_edges[solIdx], m_sols_edges[solIdy]);
	//auto& edgex(m_sols_edges[solIdx]);
	//auto& edgey(m_sols_edges[solIdy]);
	//double dis(0);
	//for (int idv(0); idv < m_dim; ++idv) {
	//	for (auto& to1 : edgex[idv]) {
	//		for (auto& to2 : edgey[idv]) {
	//			if (to1 == to2) {
	//				dis += 1;
	//				break;
	//			}
	//		}
	//	}
	//}
	//return dis / double(m_dim * edgex.front().size());
}


size_t ofec::NBN_EdgeMultiThreadDivision::size() const
{
	return m_sols.size();
	//return size_t();
}

void ofec::NBN_EdgeMultiThreadDivision::addSol(const SolutionBase& new_sol, int& belong_id, bool flag_opt, int popIter, int popSolId, int algId)
{
	auto pro = m_environment->problem();
	
	int solId(m_sols.size());
	resize(m_sols.size() + 1);
	m_sols[solId].reset(pro->createSolution(new_sol));
	ofec::TravellingSalesman::transferEdgeSol(*m_sols[solId], m_sols_edges[solId]);
	//dynamic_cast<TravellingSalesman&>(m_pro).transferEdgeSol(*m_sols[solId], m_sols_edges[solId]);

	m_fitness[solId] = m_sols[solId]->fitness();
	m_belong[solId] = solId;
	m_dis2parent[solId] = std::numeric_limits<double>::max();
	m_flagOpt[solId] = false;
	//m_sol_direction[solId] = 0;

	belong_id = -1;
}

void ofec::NBN_EdgeMultiThreadDivision::getSharedNBN(
	std::vector<std::shared_ptr<SolutionBase>>& sols,
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

void ofec::NBN_EdgeMultiThreadDivision::getSharedNBN(std::vector<std::shared_ptr<SolutionBase>>& sols, std::vector<double>& fitness, std::vector<int>& belong, std::vector<double>& dis2parent, std::vector<int> popIters, std::vector<int> popSolIds, std::vector<int> algIds) const
{

	sols = m_sols;
	fitness = m_fitness;
	belong = m_belong;
	dis2parent = m_dis2parent;
	//	flagOpt = m_flagOpt;
}
