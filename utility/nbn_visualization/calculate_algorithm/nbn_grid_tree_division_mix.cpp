#include "nbn_grid_tree_division_mix.h"
//#include "../../core/problem/continuous/continuous.h"

namespace ofec {


	void NBN_GridTreeDivisionMix::vecToIdx(const std::vector<int>& cur,
		int& idx, const std::vector<int>& dim_div) {
		idx = 0;
		for (int idDim(0); idDim < dim_div.size(); ++idDim) {
			idx *= dim_div[idDim];
			idx += cur[idDim];
		}
	}
	void NBN_GridTreeDivisionMix::solToVec(const std::vector<double>& solX,
		std::vector<int>& vec, const std::vector<int>& dim_div,
		const std::vector<double>& range_div,
		const std::vector<double>& range_from
	) {
		vec.resize(dim_div.size());

		for (int idDim(0); idDim < dim_div.size(); ++idDim) {
			vec[idDim] = (solX[idDim] - range_from[idDim]) / range_div[idDim];
			if (vec[idDim] >= dim_div[idDim]) vec[idDim] = dim_div[idDim] - 1;
		}


	}
	void NBN_GridTreeDivisionMix::vecToCenterSol(std::vector<double>& solX, const std::vector<int>& vec, const std::vector<int>& dim_div, const std::vector<double>& range_div, const std::vector<double>& range_from)const
	{
		solX.resize(m_boundary.size());
		std::fill(solX.begin(), solX.end(), 0);
		//auto& cur_sol(dynamic_cast<const Solution<>&>(sol));
		for (int idDim(0); idDim < dim_div.size(); ++idDim) {
			solX[idDim] = range_from[idDim] + (vec[idDim] + 0.5) * range_div[idDim];
		}
	}
	void NBN_GridTreeDivisionMix::idxToVec(int idx, std::vector<int>& cur,
		const std::vector<int>& dim_div) {
		for (int idDim(dim_div.size() - 1); idDim >= 0; --idDim) {
			cur[idDim] = idx % dim_div[idDim];
			idx /= dim_div[idDim];
		}
	}



	void NBN_GridTreeDivisionMix::getSolInNodeId(const std::vector<double>& solX,
		const TreeNode& node, int& belong_id) {
		std::vector<int> cur_vec;
		solToVec(solX, cur_vec, node.m_dim_div, node.m_range_div, node.m_range_from);
		vecToIdx(cur_vec, belong_id, node.m_dim_div);
	}
	void NBN_GridTreeDivisionMix::generateSols()
	{   
		// can not work in this method

		//if (!m_flag_multiThread || m_pro->hasTag(ProTag::kCSIWDN)) {
		//	SolBase* ptr_cur_sol;
		//	int id(0);
		//	int from_size = m_sols.size();
		//	assignedSols(m_maxSample - from_size);
		//	for (int idx(m_sols.size()); idx < m_maxSample; ++idx) {
		//		ptr_cur_sol = m_pro->createSolution();
		//		m_pro->initSolution(*ptr_cur_sol, m_random.get());
		//		m_eval_fun(*ptr_cur_sol, m_pro.get());
		//		std::shared_ptr<SolBase> cur_sol;
		//		cur_sol.reset(ptr_cur_sol);
		//		setSol(idx, cur_sol);
		//		//addRandomSol(*cur_sol, id);
		//	}
		//}
		//else {
		//	int id(0);
		//	int from_size = m_sols.size();
		//	int initNum = m_maxSample - from_size;
		//	std::vector<SolBase*> sols(initNum);

		//	assignedSols(initNum);
		//	UTILITY::generateRandomSolutionsMultiThreads(sols, m_random.get(), m_pro.get(), m_eval_fun);
		//	for (int idx(from_size); idx < sols.size(); ++idx) {
		//		std::shared_ptr<SolBase> cur_sol;
		//		cur_sol.reset(sols[idx]);

		//		setSol(idx, cur_sol);
		//		//addRandomSol(*cur_sol, id);
		//	}
		//	sols.clear();
		//}
	}
	void NBN_GridTreeDivisionMix::getSolId(const std::vector<double>& solX, std::vector<int>& nodePath) const{
		TreeNode* cur_node(m_root.get());
		nodePath.clear();
		int belong_idx(0);
		while (true) {
			getSolInNodeId(solX, *cur_node, belong_idx);
			nodePath.push_back(belong_idx);
			if (cur_node->m_trees[belong_idx] != nullptr) {
				cur_node = cur_node->m_trees[belong_idx].get();
			}
			else break;
		}
	}

	void NBN_GridTreeDivisionMix::assignedSols(int numSols) {
		int originSize(m_sols.size());
		m_sols.resize(m_sols.size() + numSols);
		m_solX.resize(m_sols.size());
		m_fitness.resize(m_sols.size(), 0);
		m_nbn_info.m_dis2parent.resize(m_sols.size(), std::numeric_limits<double>::max());
		m_nbn_info.m_belong.resize(m_sols.size(), -1);
		for (int idx(originSize); idx < m_sols.size(); ++idx) {
			m_nbn_info.m_belong[idx] = idx;
		//	m_solRndId[idx] = getRndId();
		}
	}
	void NBN_GridTreeDivisionMix::updateFitnessThreadTask(int from, int to) {
		for (size_t idx(from); idx < to; ++idx) {
			m_eval_fun(*m_sols[idx], m_pro.get());
			m_fitness[idx] = m_sols[idx]->fitness();
		}
	}


	void NBN_GridTreeDivisionMix::updateFitness() {
		if (!m_flag_multiThread || !m_flag_evaluate_multitask || m_pro->hasTag(ProTag::kCSIWDN)) {
			updateFitnessThreadTask(m_solsFitnessUpdatedIdxs, m_sols.size());
		}
		else {
			std::vector<std::thread> thrds;
			int num_task = std::thread::hardware_concurrency();
			int num_samples = m_sols.size() - m_solsFitnessUpdatedIdxs;
			std::vector<int> tasks;
			UTILITY::assignThreads(num_samples, num_task, tasks);
			std::pair<int, int> from_to;
			for (size_t i = 0; i < num_task; ++i) {
				from_to.first = tasks[i] + m_solsFitnessUpdatedIdxs;
				from_to.second = tasks[i + 1] + m_solsFitnessUpdatedIdxs;
				thrds.push_back(std::thread(
					&NBN_GridTreeDivisionMix::updateFitnessThreadTask, this,
					from_to.first, from_to.second));
			}
			for (auto& thrd : thrds)
				thrd.join();
		}
		m_solsFitnessUpdatedIdxs = m_sols.size();
	}

	int NBN_GridTreeDivisionMix::updateLocalNetwork(
		const std::vector<int>& solIdxs, NBN_Info& nbn_info, ofec::Random* rand)const
	{
		std::vector<int> sortedIds(solIdxs);
		std::sort(sortedIds.begin(), sortedIds.end(), [&](int a, int b) {
			return compareTwoSol(a, b);
			//return m_fitness[a] < m_fitness[b];
		});

		auto& dis2parent = nbn_info.m_dis2parent;
		auto& belong = nbn_info.m_belong;

		double cur_dis(0);
		for (int from_idx(0); from_idx < sortedIds.size(); ++from_idx) {
			auto cur_solId(sortedIds[from_idx]);
			auto& cur_sol(m_sols[cur_solId]);
			for (int to_idx(from_idx + 1); to_idx < sortedIds.size(); ++to_idx) {
				auto better_solId(sortedIds[to_idx]);

				/*if (m_fitness[better_solId] >= m_fitness[cur_solId]) */{
					auto& better_sol(m_sols[better_solId]);
					cur_dis = getSolDis(cur_solId, better_solId);
					//cur_dis = m_pro->norVariableDistance(*cur_sol, *better_sol);
					if (cur_dis < dis2parent[cur_solId] || (cur_dis == dis2parent[cur_solId] && rand->uniform.next() < 0.5)) {
						dis2parent[cur_solId] = cur_dis;
						belong[cur_solId] = better_solId;
					}
				}
			}

		}

		return sortedIds.back();
	}


	void NBN_GridTreeDivisionMix::generateNetwork() {

		m_dim = m_boundary.size();
		m_nbn_info.init(m_sols.size());

		std::vector<int> solIds(m_sols.size());
		for (int idx(0); idx < m_sols.size(); ++idx) {
			solIds[idx] = idx;
		}
		if (judgeDivisionFlag(solIds.size())) {
			m_root.reset(new TreeNode);
			const auto& boundary = m_boundary;
			
			m_root->m_nbn_info = m_nbn_info;
			DivisionInfo divisionInfo;
			divisionInfo.m_solX = m_solX;
			updateNetwork(*m_root, solIds, boundary, divisionInfo, m_random.get());

			m_nbn_info = m_root->m_nbn_info;
		}
		else {
			updateLocalNetwork(solIds, m_nbn_info,m_random.get());
		}
	}


	void NBN_GridTreeDivisionMix::mergeNetworkInfoTask(
		const std::vector<TreeNode>& totalInfo,
		NBN_Info& cur, int from, int to, ofec::Random* rnd)const
	{
		for (int idx(from); idx < to; ++idx) {

			auto& cur_dis = cur.m_dis2parent[idx];
			auto& cur_belong = cur.m_belong[idx];
			for (int idthrd(0); idthrd < totalInfo.size(); ++idthrd) {

				auto& curInfo = totalInfo[idthrd].m_nbn_info;

				if (cur_dis > curInfo.m_dis2parent[idx]) {
					cur_dis = curInfo.m_dis2parent[idx];
					cur_belong = curInfo.m_belong[idx];
				}
				else if (cur_dis == curInfo.m_dis2parent[idx]
					&& rnd->uniform.next() < 0.5) {
					//cur_dis = total_dis2parent[idthrd][idx];
					cur_belong = curInfo.m_belong[idx];
				}
			}
		}
	}

	void NBN_GridTreeDivisionMix::updateDivisionInfo(TreeNode& node, 
		std::vector<std::pair<double, double>>& boundary, DivisionInfo& divisionInfo,
		std::vector<int>& selectedDim, int numSelectedDim)const {
		node.m_nbn_info = m_nbn_info;
		boundary.resize(numSelectedDim);
		for (int idx(0); idx < boundary.size(); ++idx) {
			boundary[idx] = m_boundary[selectedDim[idx]];
		}
		divisionInfo.m_solX.resize(m_solX.size());
		for (int idS(0); idS < m_solX.size(); ++idS) {
			divisionInfo.m_solX[idS].resize(m_boundary.size());
			auto& curX = m_solX[idS];
			for (int idx(0); idx < divisionInfo.m_solX[idS].size(); ++idx) {
				divisionInfo.m_solX[idS][idx] = curX[selectedDim[idx]];
			}
		}
	}


	void NBN_GridTreeDivisionMix::generateNetworkMutiDivisionHD(
		int numSelectedDim) {

		m_dim = m_boundary.size();

		int numThread = std::thread::hardware_concurrency();
		std::vector<DivisionInfo> divisionInfo(numThread);
		std::vector<TreeNode> treeNodes(numThread);
		std::vector<std::vector<std::pair<double, double>>> boundary(numThread);
		std::vector<std::shared_ptr<ofec::Random>> rnds(numThread);
		for (auto& it : rnds) {
			it.reset(new ofec::Random(m_random->uniform.next()));
		}

		std::vector<int> solIds(m_sols.size());
		for (int idx(0); idx < m_sols.size(); ++idx) {
			solIds[idx] = idx;
		}

		// update network info
		{

			std::vector<std::vector<int>> selectedDims(numThread);
			// random num for different thread
			std::vector<int> selectedDim(m_boundary.size(), 0);
			for (int idx(0); idx < selectedDim.size(); ++idx) {
				selectedDim[idx] = idx;
			}
			for (auto& it : selectedDims) {
				it = selectedDim;
				m_random->uniform.shuffle(it.begin(), it.end());
			}
			std::vector<std::thread> thrds;
			int num_task = numThread;
			for (size_t i = 0; i < num_task; ++i) {
				thrds.push_back(std::thread(
					&NBN_GridTreeDivisionMix::updateDivisionInfo, 
					this, std::ref(treeNodes[i]), std::ref(boundary[i]),
					std::ref(divisionInfo[i]), std::ref(selectedDims[i]), numSelectedDim));
			}
			for (auto& thrd : thrds)
				thrd.join();
		}

		// calculate network by multithread
		{
			std::vector<std::thread> thrds;
			int num_task = numThread;
			for (size_t i = 0; i < num_task; ++i) {
				thrds.push_back(std::thread(
					&NBN_GridTreeDivisionMix::updateNetwork,
					this, std::ref(treeNodes[i]), std::cref(solIds),
					std::cref(boundary[i]),
					std::cref(divisionInfo[i]), rnds[i].get()));
			}
			for (auto& thrd : thrds)
				thrd.join();
		}


		// merge nbn info
		{
			std::vector<std::thread> thrds;
			int num_task = std::thread::hardware_concurrency();
			int num_samples = m_nbn_info.m_dis2parent.size();

			std::vector<int> tasks;
			UTILITY::assignThreads(num_samples, num_task, tasks);
			std::pair<int, int> from_to;
			for (size_t i = 0; i < num_task; ++i) {
				from_to.first = tasks[i] ;
				from_to.second = tasks[i + 1];
				thrds.push_back(std::thread(
					&NBN_GridTreeDivisionMix::mergeNetworkInfoTask, 
					this,
					std::cref(treeNodes), std::ref(m_nbn_info),
					from_to.first, from_to.second, rnds[i].get()));
			}
			for (auto& thrd : thrds)
				thrd.join();
		}
		
		
		
	}

	void NBN_GridTreeDivisionMix::calculateNetworkAccurate(TreeNode& curTreeNode) const{
		std::vector<int> sorted_idx(curTreeNode.m_division_nodes.size());
		for (int idx(0); idx < sorted_idx.size(); ++idx) {
			sorted_idx[idx] = idx;
		}
		auto& nbn_nodes(curTreeNode.m_division_nodes);
		std::sort(sorted_idx.begin(), sorted_idx.end(), [&](
			int a, int  b
			) {
			return compareTwoSol(nbn_nodes[a], nbn_nodes[b]);
			//return nbn_nodes[a].m_representative_solFit <
			//	nbn_nodes[b].m_representative_solFit;
		});

		curTreeNode.m_cur_visited2 = 0;
		for (auto& curNode : curTreeNode.m_division_nodes) {
			curNode.m_better_range_nodeId.clear();
			curNode.m_better_range_nodeId.push_back(curNode.m_node_id);
		}
		for (int sidx(0); sidx + 1 < sorted_idx.size(); ++sidx) {
			//std::cout << "curId\t" << sidx << std::endl;
			auto sortedId = sorted_idx[sidx];
			auto& cur_node(nbn_nodes[sortedId]);

			updateNeighborAndParent(curTreeNode, cur_node);

		}
		curTreeNode.m_best_divNodeId = sorted_idx.back();
		curTreeNode.m_best_solId = nbn_nodes[sorted_idx.back()].m_representative_solId;
	}

	//void NBN_GridTreeDivision_allSize::calSolsNetworkThreadTask(int from, int to, std::vector<int>& nodeId, const std::vector<int>& solIdxs)
	//{
	//}

	void NBN_GridTreeDivisionMix::divideSolsNetwork(TreeNode& curNode, 
		const std::vector<int>& solIdxs, const DivisionInfo& divisionInfo)const
	{

		std::vector<int> belongNodeIds(solIdxs.size(), -1);
		int belongNodeId = 0;
		for (int idx(0); idx < solIdxs.size(); ++idx) {
			auto& solId = solIdxs[idx];
			getSolInNodeId(divisionInfo.m_solX[solId], curNode, belongNodeId);
			//if (belongNodeId < 0) {
			//	getSolInNodeId(*m_sols[solId], curNode, belongNodeId);
			//}
			belongNodeIds[idx] = belongNodeId;
		}

		for (int idx(0); idx < solIdxs.size(); ++idx) {
			auto& solId = solIdxs[idx];
			auto& nodeId = belongNodeIds[idx];
			curNode.m_division_nodes[nodeId].m_sample_solIds.push_back(solId);
		}

		for (auto& it : curNode.m_division_nodes) {

			if (it.m_sample_solIds.empty()) {
				vecToCenterSol(it.m_representative_solX, it.m_vec, curNode.m_dim_div, curNode.m_range_div, curNode.m_range_from);
				it.m_representative_solFit = m_worstFit;
			//	it.m_sortedRndId = getRndId();
			}
		}

	}

	void NBN_GridTreeDivisionMix::updateNetwork(
		TreeNode& curNode,
		const std::vector<int>& solIdxs,
		const std::vector<std::pair<double, double>>& boundary,
		const DivisionInfo& divisionInfo,
		ofec::Random* rnd
		) const {

		curNode.initialize(boundary, solIdxs.size());

		divideSolsNetwork(curNode, solIdxs,divisionInfo);


		std::vector<std::pair<double, double>> curBoundary;

		for (auto& it : curNode.m_division_nodes) {

			if (!it.m_sample_solIds.empty()) {
				if (judgeDivisionFlag(it.m_sample_solIds.size())) {
					curNode.getSubSpace(it.m_node_id, curBoundary);
					curNode.m_trees[it.m_node_id].reset(new TreeNode);
					updateNetwork(
						*curNode.m_trees[it.m_node_id],
						it.m_sample_solIds,
						curBoundary, divisionInfo,rnd);

					it.m_representative_solId = curNode.m_trees[it.m_node_id]->m_best_solId;

				}
				else {
					it.m_representative_solId = updateLocalNetwork(it.m_sample_solIds,curNode.m_nbn_info,rnd);
				//	it.m_representative_solX = m_solX[it.m_representative_solId];
				//	it.m_representative_solFit = m_fitness[it.m_representative_solId];
				}

				it.m_representative_solX = divisionInfo.m_solX[it.m_representative_solId];
				it.m_representative_solFit = m_fitness[it.m_representative_solId];
				//it.m_sortedRndId = m_solRndId[it.m_representative_solId];
			}
		}
		calculateNetworkAccurate(curNode);
	}



	void NBN_GridTreeDivisionMix::TreeNode::getSubSpace(int nodeId, std::vector<std::pair<double, double>>& boundary) {
		boundary.resize(m_boundary.size());
		std::vector<int> vec(m_dim_div);
		idxToVec(nodeId, vec, m_dim_div);
		for (int idDim(0); idDim < boundary.size(); ++idDim) {
			boundary[idDim].first = m_boundary[idDim].first + vec[idDim] * m_range_div[idDim];
			boundary[idDim].second = boundary[idDim].first + m_range_div[idDim];
		}
	}



	void NBN_GridTreeDivisionMix::TreeNode::resizeNumberDivision(int numDiv) {
		m_trees.resize(numDiv, nullptr);
		m_neighbors.clear();
		//m_division_nodes.clear();
		m_division_nodes.resize(numDiv);
		for (int idx(0); idx < numDiv; ++idx) {
			m_division_nodes[idx].initialize(idx);
		}
	}

	void NBN_GridTreeDivisionMix::TreeNode::initialize(
		const std::vector<std::pair<double, double>>& boundary, int sampleSize) {
		m_boundary = boundary;
		int dim(m_boundary.size());
		int numSample_int = sampleSize;
		m_dim_div.resize(dim);
		int dim_div = std::max(2.0, exp(log(double(numSample_int)) / double(dim)));
		numSample_int = 1;
		for (auto& it : m_dim_div) {
			it = dim_div;
			numSample_int *= it;
		}
		m_numberDivision = numSample_int;

		m_range_div.resize(dim);
		m_range_from.resize(dim);
		for (int idDim(0); idDim < dim; ++idDim) {
			auto& range_var = m_boundary[idDim];
			m_range_div[idDim] = (range_var.second - range_var.first) / double(m_dim_div[idDim]);
			m_range_from[idDim] = range_var.first;
		}
		resizeNumberDivision(m_numberDivision);

		std::vector<int> curVec(dim, 0);
		int curIdx(0);
		// set solutions
		for (size_t idx(0); idx < m_numberDivision; ++idx) {
			idxToVec(idx, curVec, m_dim_div);
			vecToIdx(curVec, curIdx, m_dim_div);
			m_division_nodes[idx].m_vec = curVec;
		}
		// generate neighbor infos
		{
			std::vector<int> neiVec;
			int neiIdx(0);
			for (int idx(0); idx < m_numberDivision; ++idx) {
				auto& cur_node = m_division_nodes[idx];
				idxToVec(idx, curVec, m_dim_div);
				vecToIdx(curVec, curIdx, m_dim_div);
				neiVec = curVec;
				for (int idDim(0); idDim < dim; ++idDim) {
					--neiVec[idDim];
					if (neiVec[idDim] >= 0) {
						vecToIdx(neiVec, neiIdx, m_dim_div);
						auto& nei_node = m_division_nodes[neiIdx];
						cur_node.m_neighbor_nodeIds.push_back(neiIdx);
						nei_node.m_neighbor_nodeIds.push_back(curIdx);
						//			m_neighbors.push_back({ idx, neiIdx });
					}
					++neiVec[idDim];
				}
			}

		}

	}
}