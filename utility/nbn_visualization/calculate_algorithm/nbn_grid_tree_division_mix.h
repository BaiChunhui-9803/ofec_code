

#ifndef NBN_GRID_DIVISION_TREE_MIX_H
#define NBN_GRID_DIVISION_TREE_MIX_H

#include <vector>
#include <memory>
#include <queue>

#include "../../core/problem/solution.h"
#include "../../utility/typevar/typevar.h"
#include "nbn_division_base.h"
#include <thread>
#include "../../utility/function/custom_function.h"
#include <mutex>


//#include "core/problem/solution.h"
//#include "utility/typevar/typevar.h"



namespace ofec {
	class NBN_GridTreeDivisionMix : public NBN_DivisionBase {
	public:
		using SolutionType = SolBase;


		struct Node {
			int m_node_id = -1;
			int m_representative_solId = -1;
			std::vector<double> m_representative_solX;
		//	unsigned m_sortedRndId = 0;
			double m_representative_solFit;
			//std::shared_ptr<SolutionType> m_representative_sol = nullptr;
		
			int m_belong = -1;

			std::vector<int> m_vec;
			unsigned long long m_visited_stamp2 = 0;
			std::vector<int> m_neighbor_nodeIds;
			int m_parent_nodeId = -1;
			int m_direct_parent_nodeId = -1;
			std::vector<int> m_better_range_nodeId;
			std::vector<int> m_sample_solIds;

			void clear() {
				m_node_id = -1;
				m_representative_solId = -1;
				m_representative_solX.clear();
				m_representative_solFit = 0;
				//m_representative_sol = nullptr;
				m_vec.clear();
				m_visited_stamp2 = 0;
				m_neighbor_nodeIds.clear();
				m_parent_nodeId = -1;
				m_direct_parent_nodeId = -1;
				m_better_range_nodeId.clear();
				m_sample_solIds.clear();
			}
			void initialize(int id) {
				clear();
				m_node_id = id;
			}
		};


		struct NBN_Info {
			std::vector<int> m_belong;
			std::vector<double> m_dis2parent;


			void init(int num) {
				m_belong.resize(num);
				m_dis2parent.resize(num);
				for (int idx(0); idx < num; ++idx) {
					m_belong[idx] = idx;
					m_dis2parent[idx] = std::numeric_limits<double>::max();
				}
			}
		};


		struct DivisionInfo {
			std::vector<std::vector<double>> m_solX;
		};


		struct TreeNode {

			int m_son_nodeId = -1;
			//			std::shared_ptr<TreeNode> m_parent = nullptr;
			int m_best_solId = -1;
			int m_best_divNodeId = -1;
			int m_numberDivision = 0;
			//	unsigned long long m_cur_visited = 0;
			unsigned long long m_cur_visited2 = 0;
			//std::vector<double> m_center;
			std::vector<std::pair<double, double>> m_boundary;

			std::vector<int> m_dim_div;
			std::vector<double> m_range_div;
			std::vector<double> m_range_from;

			// neighbor info
			std::vector<std::array<int, 2>> m_neighbors;
			std::vector<Node> m_division_nodes;
			std::vector<std::shared_ptr<TreeNode>> m_trees;

			NBN_Info m_nbn_info;

			void resizeNumberDivision(int numDiv);
			void resetVisited2() {
				m_cur_visited2 = 0;
				for (auto& it : m_division_nodes) it.m_visited_stamp2 = m_cur_visited2;
			}

			void clear() {
				m_neighbors.clear();
				m_division_nodes.clear();
				m_trees.clear();
			}

			void initialize(const std::vector<std::pair<double, double>>& boundary, int sampleSize);
			void getSubSpace(int nodeId, std::vector<std::pair<double, double>>& boundary);
		};





	protected:

		std::shared_ptr<TreeNode> m_root;
		int m_division_threadhold = 1000;
		int m_min_div = 100;

		double m_worstFit = -std::numeric_limits<double>::max();


		std::vector<std::pair<double, double>> m_boundary;
		std::vector<std::vector<double>> m_solX;
		std::vector<std::shared_ptr<SolBase>> m_sols;
		std::vector<double> m_fitness;
		NBN_Info m_nbn_info;
		//int m_selected_dim_num = 5;
		int m_solsFitnessUpdatedIdxs = 0;




		bool compareTwoSol(int a, int b)const{
			if (m_fitness[a] == m_fitness[b]) {
				for (int idx(0); idx < m_dim; ++idx) {
					if (m_solX[a][idx] != m_solX[b][idx]) {
						return m_solX[a][idx] < m_solX[b][idx];
					}
				}
				return a < b;
			}
			else {
				return m_fitness[a] < m_fitness[b];
			}
		}

		bool compareTwoSol(const Node& a, const Node& b)const {


			if (a.m_representative_solId !=-1&& b.m_representative_solId!=-1) {
				return compareTwoSol(a.m_representative_solId, b.m_representative_solId);
			}
			else {

				if (a.m_representative_solFit == b.m_representative_solFit) {
					for (int idx(0); idx < m_dim; ++idx) {
						if (a.m_representative_solX[idx] != b.m_representative_solX[idx]) {
							return a.m_representative_solX[idx] < b.m_representative_solX[idx];
						}
					}
					return a.m_representative_solId < b.m_representative_solId;

				}
				else {
					return a.m_representative_solFit < b.m_representative_solFit;
				}
			}
		}

		//unsigned getRndId() {
		//	return m_random->uniform.nextNonStd<unsigned>(1e9, std::numeric_limits<unsigned>::max());
		//}

		double getSolDis(const Node& a, const Node& b)const {
			return euclideanDistance(a.m_representative_solX.begin(),a.m_representative_solX.end(), 
				b.m_representative_solX.begin());

			if (a.m_representative_solId != -1 && b.m_representative_solId != -1) {
			//	return
				//	ofec::euclideanDistance(m_solX[a.m_representative_solId], m_solX[b.m_representative_solId]);
					//m_pro->variableDistance(*m_sols[a.m_representative_solId], *m_sols[b.m_representative_solId]);
			}
			else return 0;
		}


		double getSolDis(int a, int b)const {
			return euclideanDistance(m_solX[a].begin(),m_solX[a].end(),
				m_solX[b].begin());
		}

		// nbn info
	protected:
		

		//double distanceFun(const std::vector<double>& x1, const std::vector<double>& x2) {
		//	return euclideanDistance(x1.begin(), x2.begin());
		//}

		void updateNeighborInfo(int idx, int dim, 
			std::vector<Node>& division_nodes, const std::vector<int>& m_dim_div) {
			std::vector<int> neiVec;
			int neiIdx(0);

			std::vector<int> curVec(dim, 0);
			int curIdx(0);
			/*for (int idx(0); idx < m_numberDivision; ++idx) */ {
				auto& cur_node = division_nodes[idx];
				idxToVec(idx, curVec, m_dim_div);
				vecToIdx(curVec, curIdx, m_dim_div);
				neiVec = curVec;
				for (int idDim(0); idDim < dim; ++idDim) {
					--neiVec[idDim];
					if (neiVec[idDim] >= 0) {
						vecToIdx(neiVec, neiIdx, m_dim_div);
						//	auto& nei_node = division_nodes[neiIdx];
						cur_node.m_neighbor_nodeIds.push_back(neiIdx);
						//	nei_node.m_neighbor_nodeIds.push_back(curIdx);
							//			m_neighbors.push_back({ idx, neiIdx });
					}
					++neiVec[idDim];

					++neiVec[idDim];

					if (neiVec[idDim] < dim) {
						vecToIdx(neiVec, neiIdx, m_dim_div);
						//	auto& nei_node = division_nodes[neiIdx];
						cur_node.m_neighbor_nodeIds.push_back(neiIdx);
						//nei_node.m_neighbor_nodeIds.push_back(curIdx);
						//m_neighbors.push_back({ idx, neiIdx });
					}

					--neiVec[idDim];
				}
			}
		}
		bool judgeDivisionFlag(int sampleSize)const {
			return sampleSize >= m_division_threadhold || (double(sampleSize) > pow(2, m_dim) && sampleSize <= m_min_div);
		}

		void updateFitnessThreadTask(int from, int to);
		void assignedSols(int numSols);

		void updateNeighborAndParent(TreeNode& curTreeNode, Node& cur_node)const
		{
			if (++curTreeNode.m_cur_visited2 == 0) {
				curTreeNode.resetVisited2();
				++curTreeNode.m_cur_visited2;
			}
			auto& nbn_nodes(curTreeNode.m_division_nodes);
			cur_node.m_visited_stamp2 = curTreeNode.m_cur_visited2;

			//auto& cur_solX(cur_node.m_representative_solX);

			std::priority_queue<NodeDisStruct> nei_que;
			NodeDisStruct curQueNode;
			for (auto& nei_info : cur_node.m_neighbor_nodeIds) {
				auto& nei_idx = nei_info;
				auto nei_node(&nbn_nodes[nei_idx]);
				if (nei_node->m_visited_stamp2 != curTreeNode.m_cur_visited2) {
					nei_node->m_visited_stamp2 = curTreeNode.m_cur_visited2;
					curQueNode.m_cur = nei_node;
					curQueNode.m_cur_id = nei_idx;
					curQueNode.m_cur_dis = getSolDis(cur_node, *nei_node);

					nei_que.push(curQueNode);
				}
			}


			std::vector<Node*> attraction;
			while (!nei_que.empty()) {
				curQueNode = nei_que.top();
				auto nei_node = curQueNode.m_cur;
				if (nei_node->m_parent_nodeId !=-1)
					//if (m_belong[nei_node->m_belong] != nei_node->m_representative_solId)
					/*if (m_fitness[cur_node.m_representative_solId] > m_fitness[nei_node->m_representative_solId])*/ {
					nei_que.pop();
					attraction.push_back(curQueNode.m_cur);
					for (auto& son_idx : nei_node->m_better_range_nodeId) {
						auto son_node(&nbn_nodes[son_idx]);
						if (son_node->m_visited_stamp2 != curTreeNode.m_cur_visited2) {
							son_node->m_visited_stamp2 = curTreeNode.m_cur_visited2;
							curQueNode.m_cur = son_node;
							curQueNode.m_cur_id = son_idx;
							{

								curQueNode.m_cur_dis = getSolDis(cur_node, *son_node);
								//auto& son_solX(son_node->m_representative_solX);
							//	curQueNode.m_cur_dis = distanceFun(cur_solX, son_solX);
								//cur_node.m_dis2parent = cur_sol->norVariableDistance(*par_sol, m_id_pro);
							}
							nei_que.push(curQueNode);
						}
					}
				}
				else break;
			}
			{
				//cur_node.m_better_range_nodeId = peak_ranges;
				cur_node.m_better_range_nodeId.clear();
				NodeDisStruct parentNode = nei_que.top();
				while (!nei_que.empty()) {
					curQueNode = nei_que.top();
					nei_que.pop();
					cur_node.m_better_range_nodeId.push_back(curQueNode.m_cur_id);
				}
				cur_node.m_direct_parent_nodeId = cur_node.m_parent_nodeId = parentNode.m_cur_id;
				if (cur_node.m_representative_solId != -1) {

					int curSolId = cur_node.m_representative_solId;
					int parSolId = nbn_nodes[cur_node.m_direct_parent_nodeId].m_representative_solId;
					


					auto& dis2parent = curTreeNode.m_nbn_info.m_dis2parent;
					auto& belong = curTreeNode.m_nbn_info.m_belong;

					if (dis2parent[curSolId] > parentNode.m_cur_dis) {
						dis2parent[curSolId] = parentNode.m_cur_dis;
						belong[curSolId] = parSolId;
					}

				}
			}
		}
		void calculateNetworkAccurate(TreeNode& curNode)const;
		void mergeNetworkInfoTask(
			const std::vector<TreeNode>& totalInfo,
			NBN_Info& cur, 
			int from, int to,ofec::Random* rnd)const;
		void updateDivisionInfo(TreeNode& node, std::vector<std::pair<double, double>>& boundary, DivisionInfo& divisionInfo, std::vector<int>& selectedDim, int numSelectedDim)const;

		void divideSolsNetwork(TreeNode& curNode, const std::vector<int>& solIdxs, const DivisionInfo& divisionInfo)const;
		void updateNetwork(
			TreeNode& curNode,
			const std::vector<int>& solIdxs,
			const std::vector<std::pair<double, double>>& boundary,
			const DivisionInfo& divisionInfo,
			ofec::Random* rnd
			)const;

		int  updateLocalNetwork(
			const std::vector<int>& solIdxs,
			NBN_Info& nbn_info, 
			ofec::Random* rand)const;
		void getSolId(const std::vector<double>& solX, std::vector<int>& nodeIds)const;

		virtual void addSol(const SolBase& new_sol, int& belong_id, bool flag_opt = false,
			int popIter = -1, int popSolId = -1, int algId = -1) override {
			//		insertSol(new_sol, belong_id, flag_opt, popIter, popSolId, algId);
		}
		void setSol(int solId, std::shared_ptr<SolBase>& sol) {
			m_fitness[solId] = sol->fitness();
			m_sols[solId] = sol;
		}
	public:

		virtual void initialize_(bool flag_grid_sample = true) override {
			if (m_root != nullptr) {
				m_root->m_trees.clear();
			}

			m_root = nullptr;
		}


		static void vecToIdx(const std::vector<int>& cur, int& idx, const std::vector<int>& dim_div);
		static void solToVec(const std::vector<double>& solX, std::vector<int>& vec, const std::vector<int>& dim_div,
			const std::vector<double>& range_div,
			const std::vector<double>& range_from);

		void vecToCenterSol(std::vector<double>& solX, 
			const std::vector<int>& vec,
			const std::vector<int>& dim_div,
			const std::vector<double>& range_div,
			const std::vector<double>& range_from)const;

		static void idxToVec(int idx, std::vector<int>& vec, const std::vector<int>& dim_div);
		static void getSolInNodeId(const std::vector<double>& solX, const TreeNode& node, int& belong_id);

	public:

		struct NodeDisStruct {
			int m_cur_id = 0;
			Node* m_cur = nullptr;
			double m_cur_dis = 0;
			bool operator<(const NodeDisStruct& a) const
			{
				return m_cur_dis > a.m_cur_dis;
			}
		};


		NBN_GridTreeDivisionMix() = default;
		~NBN_GridTreeDivisionMix() {
			if (m_root != nullptr) {
				m_root->m_trees.clear();
			}
			m_root = nullptr;
		}
		virtual size_t size() const override {
			return m_sols.size();
		}


		void updateBestSolSub(int numSubSols) {

			std::vector<int> solIds(m_sols.size());
			for (int idx(0); idx < solIds.size(); ++idx) {
				solIds[idx] = idx;
			}
			std::sort(solIds.begin(), solIds.end(), [&](int a, int b) {
				return compareTwoSol(a, b);
				//	return m_fitness[a] > m_fitness[b];
			});

			solIds.resize(numSubSols);

			updateLocalNetwork(
				solIds, m_nbn_info,m_random.get());
		}
		//void setNumberDivisionDim(int num) {
		//	m_selected_dim_num = num;
		//}
		void setBoundary(const std::vector<std::pair<double, double>>& boundary) {
			m_boundary = boundary;
		}
		//void setSolX(int solId, const std::vector<double>& solX) {
		//	m_solX[solId] = solX;
		//}
		//void updateGlobalNetwork();
		virtual void generateSols();
		virtual void updateFitness();
		virtual void addSol(const SolBase& new_sol,
			const std::vector<double>& solX) {
			std::shared_ptr<SolBase> cur_sol(m_pro->createSolution(new_sol));
			assignedSols(1);
			m_sols.back() = cur_sol;
			m_solX.back() = solX;
			m_fitness.back() = cur_sol->fitness();
		}
		int getSolNodeId(const std::vector<double>& solX)const {
			int nodeId(-1);
			getSolInNodeId(solX, *m_root, nodeId);
			return m_root->m_division_nodes[nodeId].m_representative_solId;
			//return nodeId;
		}
		void generateNetwork();
		void generateNetworkMutiDivisionHD(int selectedDim);
		void getNearestBetterNetworkShareMemory(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2par) {
			sols = m_sols;
			fitness = m_fitness;
			belong = m_nbn_info.m_belong;
			dis2par = m_nbn_info.m_dis2parent;
		}

		void getResult(
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			std::vector<double>& fitness
		) const{
			dis2parent = m_nbn_info.m_dis2parent;
			belong = m_nbn_info.m_belong;
			fitness = m_fitness;
		}



		void getNearestBetterNetworkShareMemory(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2par,
			std::vector<int> popIters,
			std::vector<int> popSolIds,
			std::vector<int> algIds
		) {
			sols = m_sols;
			fitness = m_fitness;
			belong = m_nbn_info.m_belong;
			dis2par = m_nbn_info.m_dis2parent;
			//	popIters = m_popIter;
			//	popSolIds = m_popSolId;
		//		algIds = m_algId;
		}

		// belong[idx]==idx  peakInfo
		virtual void getSharedNBN(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2par,
			std::vector<bool>& flagOpt
		)const override {
			sols = m_sols;
			fitness = m_fitness;
			belong = m_nbn_info.m_belong;
			dis2par = m_nbn_info.m_dis2parent;
			flagOpt.resize(m_sols.size());
			std::fill(flagOpt.begin(), flagOpt.end(), false);
			//	flagOpt = m_flagOpt;
		}

		// belong[idx]==idx  peakInfo
		virtual void getSharedNBN(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			std::vector<int> popIters,
			std::vector<int> popSolIds,
			std::vector<int> algIds
		)const override {
			sols = m_sols;
			fitness = m_fitness;
			belong = m_nbn_info.m_belong;
			dis2parent = m_nbn_info.m_dis2parent;
		}

	};

}

#endif 