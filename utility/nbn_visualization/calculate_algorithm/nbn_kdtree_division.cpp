#include "nbn_kdtree_division.h"

#include "../function/custom_function.h"
#include "../../utility/nondominated_sorting/fast_sort.h"
#include "../../utility/nondominated_sorting/filter_sort.h"

namespace ofec {


	/*void NBN_KDTreeDivision::vecToIdx(const std::vector<int>& cur,
		int& idx, const std::vector<int>& dim_div) {
		idx = 0;
		for (int idDim(0); idDim < dim_div.size(); ++idDim) {
			idx *= dim_div[idDim];
			idx += cur[idDim];
		}
	}
	void NBN_KDTreeDivision::solToVec(const SolutionBase& sol,
		std::vector<int>& vec, const std::vector<int>& dim_div,
		const std::vector<double>& range_div,
		const std::vector<double>& range_from
	) {
		auto& cur_sol(dynamic_cast<const Solution<>&>(sol));
		vec.resize(dim_div.size());

		for (int idDim(0); idDim < dim_div.size(); ++idDim) {
			vec[idDim] = (cur_sol.variable()[idDim] - range_from[idDim]) / range_div[idDim];
			if (vec[idDim] >= dim_div[idDim]) vec[idDim] = dim_div[idDim] - 1;
		}


	}
	void NBN_KDTreeDivision::vecToCenterSol(SolContinousType& sol, const std::vector<int>& vec, const std::vector<int>& dim_div, const std::vector<double>& range_div, const std::vector<double>& range_from)
	{
		for (int idDim(0); idDim < dim_div.size(); ++idDim) {
			sol.variable()[idDim] = range_from[idDim] + (vec[idDim] + 0.5) * range_div[idDim];
		}

	}
	void NBN_KDTreeDivision::idxToVec(int idx, std::vector<int>& cur, const std::vector<int>& dim_div) {

		for (int idDim(dim_div.size() - 1); idDim >= 0; --idDim) {
			cur[idDim] = idx % dim_div[idDim];
			idx /= dim_div[idDim];
		}
	}*/

	void NBN_KDTreeDivision::generateSols(){

		if (!m_flag_multiThread || m_pro->hasTag(ProblemTag::kCSIWDN)) {
			SolutionBase* ptr_cur_sol;
			int id(0);
			int from_size = m_sols.size();
			assignedSols(m_maxSample - from_size);
			for (int idx(m_sols.size()); idx < m_maxSample; ++idx) {
				ptr_cur_sol = m_pro->createSolution();
				m_pro->initializeSolution(*ptr_cur_sol, m_random.get());
				m_eval_fun(*ptr_cur_sol, m_pro.get());
				std::shared_ptr<SolutionBase> cur_sol;
				cur_sol.reset(ptr_cur_sol);
				setSol(idx, cur_sol);
				//addRandomSol(*cur_sol, id);
			}
		}
		else {
			int id(0);
			int from_size = m_sols.size();
			int initNum = m_maxSample - from_size;
			std::vector<SolutionBase*> sols(initNum);

			assignedSols(initNum);
			UTILITY::generateRandomSolutionsMultiThreads(sols, m_random.get(), m_pro.get(), m_eval_fun);
			for (int idx(from_size); idx < sols.size(); ++idx) {
				std::shared_ptr<SolutionBase> cur_sol;
				cur_sol.reset(sols[idx]);

				setSol(idx, cur_sol);
				//addRandomSol(*cur_sol, id);
			}
			sols.clear();
		}
	}
	void NBN_KDTreeDivision::getSolId(const SolutionBase& new_sol, std::vector<int>& nodePath) {
		TreeNode& cur_node(*m_root);
		nodePath.clear();
		int belong_idx(0);
		while (true) {
			getSolInNodeId(new_sol, cur_node, belong_idx);
			nodePath.push_back(belong_idx);
			if (cur_node.m_trees[belong_idx] != nullptr) {
				cur_node = *cur_node.m_trees[belong_idx];
			}
			else break;
		}
	}

	void NBN_KDTreeDivision::assignedSols(int numSols) {
		int originSize(m_sols.size());
		m_sols.resize(m_sols.size() + numSols);
		m_belong.resize(m_sols.size(), -1);
		m_fitness.resize(m_sols.size(), 0);
		m_dis2parent.resize(m_sols.size(), std::numeric_limits<double>::max());
		//m_flagOpt.resize(m_sols.size(), false);
		//m_popIter.resize(m_sols.size(), -1);
		//m_popSolId.resize(m_sols.size(), -1);
		//m_algId.resize(m_sols.size(), -1);

		for (int idx(originSize); idx < m_belong.size(); ++idx) {
			m_belong[idx] = idx;
		}
	}
	void NBN_KDTreeDivision::updateFitnessThreadTask(int from, int to) {
		for (size_t idx(from); idx < to; ++idx) {
			m_eval_fun(*m_sols[idx], m_pro.get());
			m_fitness[idx] = m_sols[idx]->fitness();
		}
	}


	void NBN_KDTreeDivision::updateFitness() {
		if (!m_flag_multiThread || !m_flag_evaluate_multitask || m_pro->hasTag(ProblemTag::kCSIWDN)) {
			updateFitnessThreadTask(m_solsFitnessUpdatedIdxs, m_sols.size());
		}
		else {
			if (m_pro->numberObjectives() > 1) {
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
						&NBN_KDTreeDivision::updateFitnessThreadTask, this,
						from_to.first, from_to.second));
				}
				for (auto& thrd : thrds)
					thrd.join();
			}
			else {
				for (auto& it : m_sols) {
					it->evaluate(m_pro.get());
				}

				std::vector<std::vector<Real>*> objs;
				for (size_t i = 0; i < m_sols.size(); ++i) {
					objs.emplace_back(&m_sols[i]->objective());
				}
				std::vector<int> rank;
				nd_sort::fastSortP<Real>(objs, rank, CAST_CONOP(m_pro.get())->optimizeMode());
				for (int idx(0); idx < m_sols.size(); ++idx) {
					m_fitness[idx] = -rank[idx];
					m_sols[idx]->setFitness(m_fitness[idx]);
				}
			}
		}
		m_solsFitnessUpdatedIdxs = m_sols.size();
	}

	int NBN_KDTreeDivision::updateLocalNetwork(const std::vector<int>& solIdxs)
	{
		std::vector<int> sortedIds(solIdxs);
		std::sort(sortedIds.begin(), sortedIds.end(), [&](int a, int b) {
			return m_fitness[a] < m_fitness[b];
			});

		double cur_dis(0);
		for (int from_idx(0); from_idx < sortedIds.size(); ++from_idx) {
			auto cur_solId(sortedIds[from_idx]);
			//if (cur_solId == 3235) {
			//	int stop = -1;
			//}
			auto& cur_sol(m_sols[cur_solId]);
			for (int to_idx(from_idx + 1); to_idx < sortedIds.size(); ++to_idx) {
				auto better_solId(sortedIds[to_idx]);

				if (m_fitness[better_solId] >= m_fitness[cur_solId]) {
					auto& better_sol(m_sols[better_solId]);
					cur_dis = m_pro->normalizedVariableDistance(*cur_sol, *better_sol);
					if (cur_dis < m_dis2parent[cur_solId] || (cur_dis == m_dis2parent[cur_solId] && m_random->uniform.next() < 0.5)) {
						m_dis2parent[cur_solId] = cur_dis;
						m_belong[cur_solId] = better_solId;
					}
				}
			}

		}

		return sortedIds.back();
	}

	void NBN_KDTreeDivision::generateNetwork() {

		m_dim = m_pro->numberVariables();
		std::vector<int> solIds(m_sols.size());
		for (int idx(0); idx < m_sols.size(); ++idx) {
			solIds[idx] = idx;
		}
		if (divisionFlag(solIds.size())) {
			m_root.reset(new TreeNode);
			const auto& boundary = CAST_CONOP(m_pro.get())->boundary();
			updateNetwork(*m_root, solIds, boundary);
		}
		else {
			updateLocalNetwork(solIds);
		}
	}

	void NBN_KDTreeDivision::calculateNetworkAccurate(TreeNode& curTreeNode) {
		std::vector<int> sorted_idx(curTreeNode.m_division_nodes.size());
		for (int idx(0); idx < sorted_idx.size(); ++idx) {
			sorted_idx[idx] = idx;
		}
		auto& nbn_nodes(curTreeNode.m_division_nodes);
		std::sort(sorted_idx.begin(), sorted_idx.end(), [&](
			int a, int  b
			) {
				return nbn_nodes[a].m_representative_sol->fitness() <
					nbn_nodes[b].m_representative_sol->fitness();
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

	//void NBN_KDTreeDivision::calSolsNetworkThreadTask(int from, int to, std::vector<int>& nodeId, const std::vector<int>& solIdxs)
	//{
	//}

	void NBN_KDTreeDivision::divideSolsNetwork(TreeNode& curNode, const std::vector<int>& solIdxs){

		std::vector<int> belongNodeIds(solIdxs.size(), -1);
		int belongNodeId = 0;
		for (int idx(0); idx < solIdxs.size(); ++idx) {
			auto& solId = solIdxs[idx];
			getSolInNodeId(*m_sols[solId], curNode, belongNodeId);
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
				it.m_representative_sol.reset(m_pro->createSolution());
				SolContinousType& cursol = dynamic_cast<SolContinousType&>(*it.m_representative_sol);
				solToVec(cursol, it.m_vec, curNode.m_dim_div, curNode.m_range_div, curNode.m_range_from);
				cursol.setFitness(m_worstFit);
			}
		}

	}

	void NBN_KDTreeDivision::updateNetwork(
		TreeNode& curNode,
		const std::vector<int>& solIdxs,
		const std::vector<std::pair<double, double>>& boundary) {

		curNode.initialize(boundary, solIdxs.size());

		divideSolsNetwork(curNode, solIdxs);


		std::vector<std::pair<double, double>> curBoundary;

		for (auto& it : curNode.m_division_nodes) {

			if (it.m_node_id == 36) {
				int stop = -1;
			}
			if (!it.m_sample_solIds.empty()) {
				if (judgeDivisionFlag(it.m_sample_solIds.size())) {
					curNode.getSubSpace(it.m_node_id, curBoundary);
					curNode.m_trees[it.m_node_id].reset(new TreeNode);
					updateNetwork(
						*curNode.m_trees[it.m_node_id],
						it.m_sample_solIds,
						curBoundary);

					it.m_representative_solId = curNode.m_trees[it.m_node_id]->m_best_solId;
					it.m_representative_sol = m_sols[it.m_representative_solId];

				}
				else {
					it.m_representative_solId = updateLocalNetwork(it.m_sample_solIds);
					it.m_representative_sol = m_sols[it.m_representative_solId];
				}
			}
		}
		calculateNetworkAccurate(curNode);
	}

	void NBN_KDTreeDivision::TreeNode::initialize(
		const std::vector<std::pair<double, double>>& boundary, int sampleSize) {
		m_boundary = boundary;
		int dim(m_boundary.size());
		/*int numSample_int = sampleSize;
		m_dim_div.resize(dim);
		int dim_div = std::max(2.0, exp(log(double(numSample_int)) / double(dim)));
		numSample_int = 1;
		for (auto& it : m_dim_div) {
			it = dim_div;
			numSample_int *= it;
		}*/
		m_numberDivision = sampleSize;

		/*m_range_div.resize(dim);
		m_range_from.resize(dim);
		for (int idDim(0); idDim < dim; ++idDim) {
			auto& range_var = m_boundary[idDim];
			m_range_div[idDim] = (range_var.second - range_var.first) / double(m_dim_div[idDim]);
			m_range_from[idDim] = range_var.first;
		}*/

		resizeNumberDivision(m_numberDivision);

		std::vector<int> curVec(dim, 0);
		int curIdx(0);
		// set solutions
		/*for (size_t idx(0); idx < m_numberDivision; ++idx) {
			idxToVec(idx, curVec, m_dim_div);
			vecToIdx(curVec, curIdx, m_dim_div);
			m_division_nodes[idx].m_vec = curVec;
		}*/

		// set neighbor sols
		for (int i = 0; i < m_kd_tree->size(); ++i) {
			std::list<int> neighs;
			m_kd_tree->findNeighbor(i, neighs);

			for (auto nei : neighs) {
				m_neighbors.push_back({ i, nei });
				m_division_nodes[i].m_neighbor_nodeIds.push_back(nei);
			}
			/*for (auto nei : neighs) {
				double cur_dis = m_pro->normalizedVariableDistance(*m_nbn_nodes[i].m_sol, *m_nbn_nodes[nei].m_sol);
				m_nbn_nodes[i].m_neighbor_id_dis.emplace_back(nei, cur_dis);
			}*/
		}
	}

	void NBN_KDTreeDivision::TreeNode::getSubSpace(int nodeId, std::vector<std::pair<double, double>>& boundary) {
		boundary.resize(m_boundary.size());
		std::vector<int> vec(m_dim_div);
		idxToVec(nodeId, vec, m_dim_div);
		for (int idDim(0); idDim < boundary.size(); ++idDim) {
			boundary[idDim].first = m_boundary[idDim].first + vec[idDim] * m_range_div[idDim];
			boundary[idDim].second = boundary[idDim].first + m_range_div[idDim];
		}

	}

	void NBN_KDTreeDivision::TreeNode::resizeNumberDivision(int numDiv) {
		m_trees.resize(numDiv, nullptr);
		m_neighbors.clear();
		//m_division_nodes.clear();
		m_division_nodes.resize(numDiv);
		for (int idx(0); idx < numDiv; ++idx) {
			m_division_nodes[idx].initialize(idx);
		}
	}
}
