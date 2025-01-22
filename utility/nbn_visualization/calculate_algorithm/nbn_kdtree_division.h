

#ifndef NBN_KDTREE_DIVISION_H
#define NBN_KDTREE_DIVISION_H

#include <vector>
#include <memory>
#include <queue>

#include "../../core/problem/solution.h"
#include "../../utility/parameter/param_map.h"
#include "nbn_division_base.h"

#include <thread>
#include "../../utility/function/custom_function.h"
#include "../../core/problem/continuous/continuous.h"
#include <mutex>

#include "../../utility/kd-tree/kdtree_space.h"
#include "../../utility/random/newran.h"

//#include "core/problem/solution.h"
//#include "utility/parameter/param_map.h"



namespace ofec {
	class NBN_KDTreeDivision : public NBN_DivisionBase {
	public:
		using SolutionType = SolutionBase;
		using SolContinousType = Solution<>;
		struct Node {

			int m_node_id = -1;
			int m_representative_solId = -1;
			std::shared_ptr<SolutionType> m_representative_sol = nullptr;
			int m_belong = -1;
			//double m_dis2parent = 0;
			std::vector<int> m_vec;
			//	unsigned long long m_visited_stamp = 0;
			unsigned long long m_visited_stamp2 = 0;
			std::vector<int> m_neighbor_nodeIds;
			//std::vector<std::pair<int, double>> m_neighbor_nodeId_dis;
			int m_parent_nodeId = -1;
			int m_direct_parent_nodeId = -1;
			std::vector<int> m_better_range_nodeId;
			std::vector<int> m_sample_solIds;

			void clear() {
				m_node_id = -1;
				m_representative_solId = -1;
				m_representative_sol = nullptr;
				//m_parent_solId = -1;
			//	m_dis2parent = 0;
				m_vec.clear();
				m_visited_stamp2 = 0;
				m_neighbor_nodeIds.clear();
				//	m_neighbor_nodeId_dis.clear();
				m_parent_nodeId = -1;
				m_direct_parent_nodeId = -1;
				m_better_range_nodeId.clear();
				m_sample_solIds.clear();
			}
			void initialize(int id) {
				clear();
				m_node_id = id;

				//	m_dis2parent = std::numeric_limits<double>::max();
			}
		};

		struct TreeNode {

			int m_son_nodeId = -1;
			//std::shared_ptr<TreeNode> m_parent = nullptr;

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


			//int m_updateNode = 0;
			//std::mutex m_updateMtx;


			void resizeNumberDivision(int numDiv);

			//void resetVisited() {
			//	m_cur_visited = 0;
			//	for (auto& it : m_division_nodes) it.m_visited_stamp = m_cur_visited;
			//}
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
			//	void mergeParent(Node* curNode);
			//	void findParent(Node* curNode, Node*& parent);

		};

	protected:
		std::shared_ptr<nanoflann::KDTreeSpace<Real>> m_kd_tree;
		//int m_dim = 2;
		std::shared_ptr<TreeNode> m_root;
		//// updateNodeInfo
		//std::vector<int> m_node_paths;
		//bool m_flag_node_division = false;
		int m_division_threadhold = 1000;
		int m_min_div = 100;
		int m_num_samples;
		int m_num_division;
		int m_num_optima;

		double m_worstFit = -std::numeric_limits<double>::max();

		//std::vector<std::shared_ptr<SolutionType>> m_added_sols;
		std::vector<std::shared_ptr<SolutionType>> m_sols;
		std::vector<double> m_fitness;
		std::vector<int> m_belong;
		std::vector<double> m_dis2parent;


		int m_solsFitnessUpdatedIdxs = 0;

		// nbn info
	public:
		NBN_KDTreeDivision() = default;
		~NBN_KDTreeDivision() {
			if (m_root != nullptr) {
				m_root->m_trees.clear();
			}
			m_root = nullptr;
		}

		// sample in each subspace
		void NBN_KDTreeDivision::generateSols(size_t num_space,size_t num_samples,Random* rnd) {
			setDivisionNumber(num_space);
			int fromSolId = m_sols.size();
			assignedSols(num_space*num_samples);
			auto bound = CAST_CONOP(m_pro.get())->boundary();
			std::vector<Real> ratios(num_space, 1. / num_space);
			m_kd_tree=std::make_shared<nanoflann::KDTreeSpace<Real>>(ratios, bound);
			m_kd_tree->buildIndex();
			for (size_t i = 0; i < m_kd_tree->size(); ++i) {
				auto box = m_kd_tree->getBox(i);
				for (size_t j = 0; j < num_samples; ++j) {
					std::vector<Real> temp_var;
					for (size_t j = 0; j < box.size(); ++j) {
						temp_var.push_back(box[j].first + rnd->uniform.next()*(box[j].second-box[j].first));
					}
					auto& cur_sol(m_sols[i*num_samples+j + fromSolId]);
					cur_sol.reset(new Solution<>(m_pro->numberObjectives(),m_pro->numberConstraints(),m_pro->numberVariables()));
					auto& curConSol(dynamic_cast<Solution<>&>(*cur_sol));
					curConSol.variable().vect() = temp_var;
				}
			}
		}

		void updateFitness();

		void updateFitnessThreadTask(int from, int to);

		void generateNetwork();

		void assignedSols(int numSols);

		void updateNeighborInfo(int idx, std::vector<Node>& division_nodes) {
			for (int idx(0); idx < m_num_division; ++idx) {
				auto& cur_node = division_nodes[idx];
				std::list<int> neighs;
				m_kd_tree->findNeighbor(idx, neighs);
				for (auto nei : neighs) {
					cur_node.m_neighbor_nodeIds.push_back(nei);
				}
			}
		}

		bool judgeDivisionFlag(int sampleSize) {
			return sampleSize >= m_division_threadhold || (double(sampleSize) > pow(2, m_dim) && sampleSize <= m_min_div);
		}

		void updateNetwork(
			TreeNode& curNode,
			const std::vector<int>& solIdxs,
			const std::vector<std::pair<double, double>>& boundary);

		int  updateLocalNetwork(const std::vector<int>& solIdxs);
		void getSolId(const SolutionBase& new_sol, std::vector<int>& nodeIds);

		//void divideSpace(std::shared_ptr<TreeNode>&, const std::vector<std::pair<double, double>>& boundary, int sampleSize);
		//void udpateSubspace(std::shared_ptr<TreeNode>&, int curIdx);
		//void divideSubSpace(std::shared_ptr<TreeNode>&, int curIdx);
		//void generateSolution(TreeNode& curNode);
		//void calculateNetwork(TreeNode& curNode);

		void updateNeighborAndParent(TreeNode& curTreeNode, Node& cur_node) {
			if (++curTreeNode.m_cur_visited2 == 0) {
				curTreeNode.resetVisited2();
				++curTreeNode.m_cur_visited2;
			}
			auto& nbn_nodes(curTreeNode.m_division_nodes);
			cur_node.m_visited_stamp2 = curTreeNode.m_cur_visited2;

			auto& cur_sol(cur_node.m_representative_sol);

			std::priority_queue<NodeDisStruct> nei_que;
			NodeDisStruct curQueNode;
			for (auto& nei_info : cur_node.m_neighbor_nodeIds) {
				auto& nei_idx = nei_info;
				auto nei_node(&nbn_nodes[nei_idx]);
				if (nei_node->m_visited_stamp2 != curTreeNode.m_cur_visited2) {
					nei_node->m_visited_stamp2 = curTreeNode.m_cur_visited2;
					curQueNode.m_cur = nei_node;
					curQueNode.m_cur_id = nei_idx;
					{
						auto& nei_sol(nei_node->m_representative_sol);
						curQueNode.m_cur_dis = m_pro->normalizedVariableDistance(*cur_sol, *nei_sol);
						//cur_node.m_dis2parent = cur_sol->normalizedVariableDistance(*par_sol, m_problem.get());
					}
					nei_que.push(curQueNode);
				}
			}


			std::vector<Node*> attraction;
			while (!nei_que.empty()) {
				curQueNode = nei_que.top();
				auto nei_node = curQueNode.m_cur;
				if (nei_node->m_direct_parent_nodeId != -1)
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

								auto& son_sol(son_node->m_representative_sol);
								curQueNode.m_cur_dis = m_pro->normalizedVariableDistance(*cur_sol, *son_sol);
								//cur_node.m_dis2parent = cur_sol->normalizedVariableDistance(*par_sol, m_problem.get());
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
					m_dis2parent[curSolId] = parentNode.m_cur_dis;
					m_belong[curSolId] = parSolId;

				}
			}
		}

		void calculateNetworkAccurate(TreeNode& curNode);

		void calSolsNetworkThreadTask(int from, int to, std::vector<int>& nodeId, const std::vector<int>& solIdxs) {}

		void divideSolsNetwork(TreeNode& curNode, const std::vector<int>& solIdxs);

		virtual void addSol(const SolutionBase& new_sol, int& belong_id, bool flag_opt = false,
			int popIter = -1, int popSolId = -1, int algId = -1) override {
			//		insertSol(new_sol, belong_id, flag_opt, popIter, popSolId, algId);
		}

		void setSol(int solId, std::shared_ptr<SolutionBase>& sol) {
			m_fitness[solId] = sol->fitness();
			m_sols[solId] = sol;
		}

		virtual void initialize_(bool flag_grid_sample = true) override {
			if (m_root != nullptr) {
				m_root->m_trees.clear();
			}

			m_root = nullptr;
		}

		static void getSolInNodeId(const SolutionBase& new_sol, const TreeNode& node, int& belong_id) {
			std::vector<int> cur_vec;
			auto& var = new_sol.variableBase();
			auto space_idx = m_kd_tree->getRegionIdx(var);
			belong_id = m_num_optima + space_idx;
			//solToVec(new_sol, cur_vec, node.m_dim_div, node.m_range_div, node.m_range_from);
			//vecToIdx(cur_vec, belong_id, node.m_dim_div);
		}

		struct NodeDisStruct {
			int m_cur_id = 0;
			Node* m_cur = nullptr;
			double m_cur_dis = 0;
			bool operator<(const NodeDisStruct& a) const
			{
				return m_cur_dis > a.m_cur_dis;
			}
		};

		virtual size_t size()  const  override {
			return m_sols.size();
		}

		//void updateGlobalNetwork();

		//virtual void generateSols();
		//virtual void updateFitness();
		virtual void addSol(const SolutionBase& new_sol) {
			std::shared_ptr<SolutionBase> cur_sol(m_pro->createSolution(new_sol));
			assignedSols(1);
			m_sols.back() = cur_sol;
		}

		virtual int inputSol(const SolutionBase& new_sol) {
			std::shared_ptr<SolutionBase> cur_sol(m_pro->createSolution(new_sol));
			int solId = m_sols.size();
			assignedSols(1);
			m_sols.back() = cur_sol;
			return solId;
		}

		virtual void addRandom(int numSample) {
			using namespace ofec;
			SolContinousType sol(m_pro->numberObjectives(),
				m_pro->numberConstraints(),
				m_pro->numberVariables());
			for (int idx(0); idx < numSample; ++idx) {
				sol.initialize(m_pro.get(), m_random.get());
				int belong_id(0);
				addSol(sol, belong_id);
			}
		}

		int getSolNodeId(const SolutionBase& sol) {
			int nodeId(-1);
			getSolInNodeId(sol, *m_root, nodeId);
			return nodeId;
		}

		void setDivisionNumber(size_t num_subspaces) {
			m_num_division = num_subspaces;
		}

		void setSampleNumber(size_t num_samples) {
			m_num_samples = num_samples;
		}

		void setOptimaNum(size_t v) {
			m_num_optima = v;
		}

		size_t getOptimaNum() {
			return m_num_optima;
		}

		bool divisionFlag(size_t num_sample) {
			return num_sample > m_min_div;
		}

		//void updateDivision();
		//void updateDivision(std::shared_ptr<TreeNode>&, const std::vector<std::pair<double, double>>& boundary, int sampleSize);
		//

		void getNearestBetterNetworkShareMemory(
			std::vector<std::shared_ptr<SolutionBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2par) {
			sols = m_sols;
			fitness = m_fitness;
			belong = m_belong;
			dis2par = m_dis2parent;
		}

		void getNearestBetterNetworkShareMemory(
			std::vector<std::shared_ptr<SolutionBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2par,
			std::vector<int> popIters,
			std::vector<int> popSolIds,
			std::vector<int> algIds
		) {
			sols = m_sols;
			fitness = m_fitness;
			belong = m_belong;
			dis2par = m_dis2parent;
			//	popIters = m_popIter;
			//	popSolIds = m_popSolId;
		//		algIds = m_algId;
		}

		// belong[idx]==idx  peakInfo
		virtual void getSharedNBN(
			std::vector<std::shared_ptr<SolutionBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2par,
			std::vector<bool>& flagOpt
		)const override {
			sols = m_sols;
			fitness = m_fitness;
			belong = m_belong;
			dis2par = m_dis2parent;
			flagOpt.resize(m_sols.size());
			std::fill(flagOpt.begin(), flagOpt.end(), false);
			//	flagOpt = m_flagOpt;
		}

		// belong[idx]==idx  peakInfo
		virtual void getSharedNBN(
			std::vector<std::shared_ptr<SolutionBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			std::vector<int> popIters,
			std::vector<int> popSolIds,
			std::vector<int> algIds
		)const override {
			sols = m_sols;
			fitness = m_fitness;
			belong = m_belong;
			dis2parent = m_dis2parent;
		}
	};
}

#endif