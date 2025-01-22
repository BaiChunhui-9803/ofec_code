

#ifndef NBN_GRID_SAMPLE_H
#define NBN_GRID_SAMPLE_H

#include <vector>
#include <memory>
#include "../../core/problem/solution.h"
#include "../../utility/typevar/typevar.h"
#include "nbn_division_base.h"

//#include "core/problem/solution.h"
//#include "utility/typevar/typevar.h"



namespace ofec {
	class NBN_GridDivision: public NBN_DivisionBase {
	public:
		using SolutionType = SolBase;
		using SolContinousType = Solution<>;
		struct Node {

			int m_id = -1;
			std::shared_ptr<SolutionType> m_sol;
			std::shared_ptr<SolutionType> m_representative;
			int m_popIter = -1;
			int m_popSolId = -1;
			int m_algId = -1;
			std::vector<int> m_vec;

			unsigned long long m_visited_id = 0;
			std::vector<std::pair<int,double>> m_neighbor_id_dis;
			int m_parent = -1;
			int m_direct_parent = -1;
			double m_dis2parent = 0;
			//double m_dis2center = 0;
			double m_fitness = 0;
			std::vector<int> m_better_range;
			bool m_flag_opt = false;

			void clear() {
				m_id = -1;
				m_sol = nullptr;
				m_representative = nullptr;

				m_popIter = -1;
				m_popSolId = -1;
				m_algId = -1;
				
				
				m_vec.clear();
				m_visited_id = 0;
				m_neighbor_id_dis.clear();
				m_parent = -1;
				m_direct_parent = -1;
				m_dis2parent = 0;
				m_fitness = 0;
				m_better_range.clear();
				m_flag_opt = false;
			}


			void init(int id) {
				clear();
				m_id = m_direct_parent = m_parent = id;
				m_dis2parent = std::numeric_limits<double>::max();
			}
		};

	protected:




		unsigned long long m_cur_visited = 0;
		std::vector<Node> m_nbn_nodes;

		std::array<std::vector<int>, 2> m_fitness_landscape_idxs;
		// neighbor info
		std::vector<std::array<int, 2>> m_neighbors;
		int m_dim = 2;
		std::vector<int> m_dim_div;
		std::vector<double> m_range_div;
		std::vector<double> m_range_from;

	//	std::vector<int> m_max_divs;
		double m_filterDis = 0.02;


		// nbn info
	protected:

		void resetVisited() {
			m_cur_visited = 0;
			for (auto& it : m_nbn_nodes) it.m_visited_id = m_cur_visited;
		}

		//void setParents();
		
		void mergeParent(Node* cur_node);

		void idxToVec(int idx, std::vector<int>& cur);
		void vecToIdx(const std::vector<int>& cur, int& idx);
		bool judgeFeasible(const std::vector<int>& cur);

		void solToVec(const SolBase& sol, std::vector<int>& vec);
		//	void insertSearchRange();
		
		void updateFitnessThreadTask(int from,int to);



		//void init(Problem *pro, Random *rnd, const std::function<void(SolBase& sol, Problem *pro)>& eval_fun);
		void initNodesNormal();
		void updateFitness();
		void calculate();



		virtual void initialize_(bool flag_grid_sample = true) override;
		

	public:

		NBN_GridDivision() = default;
		~NBN_GridDivision() {
			for (auto& it : m_nbn_nodes) it.clear();
		}

		virtual size_t size() const  override {
			return m_nbn_nodes.size();
		}
		void updateSol(const SolBase& new_sol, int &belong_id,bool flag_opt=false,
			int popIter = -1, int popSolId = -1,int algId=-1);
		virtual void addSol(const SolBase& new_sol, int& belong_id, bool flag_opt = false,
			int popIter = -1, int popSolId = -1, int algId = -1) override {
			updateSol(new_sol, belong_id, flag_opt, popIter, popSolId, algId);
		}


		void getNearestBetterNetworkShareMemory(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<bool>& flagOpt
			);

		void getNearestBetterNetworkShareMemory(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<std::shared_ptr<SolBase>>& representative,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2par
		);


		void getNearestBetterNetworkShareMemory(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<std::shared_ptr<SolBase>>& representative,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2par,
			std::vector<int> popIters,
			std::vector<int> popSolIds,
			std::vector<int> algIds
		);


		void getNearestBetterNetworkShareMemoryFilter(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong);


		virtual void getSharedNBN(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2parent, std::vector<bool>& flagOpt
		)const override;


		virtual void getSharedNBN(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			std::vector<int> popIters,
			std::vector<int> popSolIds,
			std::vector<int> algIds
		)const override;

	};
}

#endif