#ifndef NBN_EDGE_DIVISION_H
#define NBN_EDGE_DIVISION_H

#include <vector>
#include <memory>
#include "../../core/problem/solution.h"
#include "../../core/typevar/typevar.h"
#include "../../core/random/newran.h"
#include "nbn_division_base.h"


//#include "core/problem/solution.h"
//#include "utility/typevar/typevar.h"

// for edge tsp 

namespace ofec {
	class NBN_EdgeDivision : public NBN_DivisionBase {
		struct NetworkNode {
			std::array<int,2> m_network_id = { -1,-1 };
			std::vector<NetworkNode*> m_parents;
			std::vector<NetworkNode*> m_sons;
			int m_best_solId = -1;
			std::vector<int> m_indis_ids;
			std::vector<int> m_repre_edge;
		};

		std::vector<std::vector<NetworkNode>> m_network;

		// parameters 
		int m_divide_threadhold = 2e2;
		int m_edge_division = 50;
		int m_num_division = 1e3;
		int m_divide_times = 0;
		int m_dim = 0;
		std::vector<int> m_edge_seq;

		//std::vector<int> 
		//std::vector<int> m_selected_edges;
		//std::vector<std::vector<NetworkNode>> m_networks;
		////std::vector<std::vector<int>> ;
		//std::vector<std::vector<int>> m_edgeToEdgeId;
		//std::vector<int> m_edge_indis_ids;
		std::vector<int> m_sol_direction;
		std::vector<std::vector<std::array<int, 2>>> m_sols_edges;


		std::vector<std::shared_ptr<SolBase>> m_sols;
		std::vector<double> m_fitness;
		std::vector<int> m_belong;
		std::vector<double> m_dis2parent;
		std::vector<bool> m_flagOpt;


	protected:


		void resize(int size) {
			m_sol_direction.resize(size);
			m_sols_edges.resize(size);
			m_sols.resize(size);
			m_fitness.resize(size, 0);
			m_belong.resize(size, -1);
			m_dis2parent.resize(size, 1.0);
			m_flagOpt.resize(size, false);
		}

		void solToIdx(int solId,int & nodeId);
		void idxToVec(int idx, std::vector<int>& cur);
		void vecToIdx(const std::vector<int>& cur, int& idx);


		void Cnm(std::vector<std::vector<bool>>& selectedIdx,
			std::vector<bool>& cur, int m, int from);
		void Cnm(std::vector<unsigned>& selectedIdx, unsigned cur, int curSize,int m, int from);

		void updateNetwork();

		void clearNetworkSol();


		void generateSols();
		void generateSolsThreadTask(int from,int to, double seed);
		void generateSol(int solId, Random *rnd);

	
		void divideByEdge( NetworkNode& node, int divDim, bool divFlag);


		// return the best idxs
		int udpateNeighbor(std::vector<int>& indis);
		
		double solDis(int solIdx, int solIdy);
	public:




		void updateDivision();
		// Í¨¹ý NBN_DivisionBase ¼Ì³Ð
		virtual void initialize_(bool flag_grid_sample = true) override;

		virtual size_t size() const override;

		virtual void addSol(const SolBase& new_sol, int& belong_id, bool flag_opt = false, int popIter = -1, int popSolId = -1, int algId = -1) override;

	
		virtual void getSharedNBN(
			std::vector<std::shared_ptr<SolBase>>& sols, 
			std::vector<double>& fitness, 
			std::vector<int>& belong, 
			std::vector<double>& dis2parent, 
			std::vector<bool>& flagOpt) const override;

		virtual void getSharedNBN(std::vector<std::shared_ptr<SolBase>>& sols, std::vector<double>& fitness, std::vector<int>& belong, std::vector<double>& dis2parent, std::vector<int> popIters, std::vector<int> popSolIds, std::vector<int> algIds) const override;

	};
}

#endif