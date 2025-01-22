

#ifndef SAMPLE2D_GRAPH_ALGORITHM_H
#define SAMPLE2D_GRAPH_ALGORITHM_H

#include <vector>
#include <memory>

#include "../../../core/problem/solution.h"
#include "sample_nearest_better_network.h"
#include "../../../utility/parameter/param_map.h"


namespace ofec {
	struct node {
	public:
		using SolutionType = Solution<>;
		int m_node_idx = -1;
		int m_sol_id = -1;
		std::vector<int> m_vec;
	//	bool m_near_par = false;
		//std::shared_ptr<SolutionType>  m_sol;
		unsigned long long m_visited_id = 0;
		//unsigned long long m_visited_times = 0;
		std::vector<int> m_neighbor;
		std::vector<double> m_neighbor_dis;
		int m_parent = -1;
		int m_direct_parent = -1;
		double m_dis2parent = 0;
		//std::vector<node*> m_sons;
		std::vector<int> m_better_range;

		void clear() {
			m_visited_id = 0;
			m_neighbor.clear();
			m_neighbor_dis.clear();
			m_parent = -1;
			m_direct_parent = -1;
			//m_sons.clear();
			m_better_range.clear();
		}


		void init(int id) {
			m_node_idx = m_sol_id = id;
			m_direct_parent = m_parent = id;
			m_dis2parent = std::numeric_limits<double>::max();
		}
	};


	struct treeNode {
		using SolutionType = Solution<>;
		std::vector<treeNode*> m_leaves;
		std::shared_ptr<SolutionType>  m_sol;
		
	};


	class Sample2D_Graph_Algorithm  {
	public:
		using SolutionType = Solution<>;
	protected:

	
		int m_random.get();
		int m_problem.get();
		std::function<void(Solution<>& sol, Problem *pro)> m_eval_fun;

		std::vector<std::shared_ptr<SolutionType>> m_cur_sols;
		std::vector<node> m_nbn_nodes;
		unsigned long long m_cur_visited = 0;


		std::vector<int> m_belongs;

		// generated neighbor solutions
		std::vector<std::array<int, 2>> m_neighbors;
		int m_dim = 2;
		std::vector<int> m_dim_div;
		int m_maxSample = 1e5;
		//std::vector<std::vector<int>> m_idx2vec;
		std::vector<treeNode> m_tree;
		treeNode* m_head;

		// opt related solutions
		std::vector<int> m_marker_idxs;
		//std::ve
		double m_filterDis = 0.02;


	protected:
		
		void resetVisited() {
			m_cur_visited = 0;
			for (auto& it : m_nbn_nodes) it.m_visited_id = m_cur_visited;
		}
		void initNodes2D();
		void initNodesNormal();
		void setParents();

		void updateParent(node* cur_node);

		void idxToVec(int idx, std::vector<int>& cur);
		void vecToIdx(const std::vector<int>& cur, int& idx);
		bool judgeFeasible(const std::vector<int>& cur);
		
	//	void insertSearchRange();

	public:

		void setMaxSampleSize(int sampleSize) {
			m_maxSample = sampleSize;
		}

		int addRandomSol(const std::shared_ptr<SolutionType>& new_sol,int to_idx);
		int addRandomSol(const std::shared_ptr<SolutionType>& new_sol);

		~Sample2D_Graph_Algorithm() {
			for (auto& it : m_nbn_nodes) it.clear();
		}

		const std::vector<std::array<int, 2>>& neighbors()const {
			return m_neighbors;
		}
		void init(Problem *pro, Random *rnd,const std::function<void(Solution<>& sol, Problem *pro)>& eval_fun);
		void calculate();
		void udpate_network();
		void getNearestBetterNetworkShareMemory(
			std::vector<std::shared_ptr<SolutionBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong);

		void getNearestBetterNetworkShareMemory(
			std::vector<std::shared_ptr<SolutionBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2par
			);

		void getNearestBetterNetworkShareMemoryFilter(
			std::vector<std::shared_ptr<SolutionBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong);
		void outputData(SampleNearestBetterNetworkRecord& record);
		void getTree();
		void insertOpt();

	};
}

#endif