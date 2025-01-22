#ifndef NBN_EDGE_SIMPLE_DIVISION_H
#define NBN_EDGE_SIMPLE_DIVISION_H

#include <vector>
#include <memory>
#include "../../core/problem/solution.h"
#include "../../utility/typevar/typevar.h"
#include "../../utility/random/newran.h"
#include "../../utility/hash_table/hash_table.h"
#include "nbn_division_base.h"


//#include "core/problem/solution.h"
//#include "utility/typevar/typevar.h"

// for edge tsp 

namespace ofec {
	class NBN_EdgeSimpleDivison : public NBN_DivisionBase {

	public:

		struct Node {
			double m_dis2parent = 0;
			int m_belong = 0;
		};

	protected:
		std::vector<std::vector<std::array<int, 2>>> m_sols_edges;


		std::vector<std::shared_ptr<SolBase>> m_sols;
		std::vector<double> m_fitness;
		std::vector<int> m_insertId;

		std::vector<Node> m_solInfoBetter;
		std::vector<Node> m_solInfoBetterEqual;
		//std::vector<int> m_belong;
		//std::vector<double> m_dis2parent;
	//	std::vector<bool> m_flagOpt;

		int m_numLoop = 1;
		int m_numTaskUpdateNetwork = 1;
	public:

		void setNumberLoop(int numLoop) {
			m_numLoop = numLoop;
		}
		void setNumTaskUpdateNetwork(int numTask) {
			m_numTaskUpdateNetwork = numTask;
		}

		void resize(int size) {
			//m_sol_direction.resize(size);
			m_sols_edges.resize(size);
			m_sols.resize(size);
			m_fitness.resize(size, 0);
			m_insertId.resize(size, -1);
			m_solInfoBetter.resize(size);
			m_solInfoBetterEqual.resize(size);
		}

		void swapSol(int idx, int idy) {
			std::swap(m_sols_edges[idx], m_sols_edges[idy]);
			std::swap(m_sols[idx], m_sols[idy]);
			std::swap(m_fitness[idx], m_fitness[idy]);
			std::swap(m_insertId[idx], m_insertId[idy]);
			std::swap(m_solInfoBetter[idx], m_solInfoBetter[idy]);
			std::swap(m_solInfoBetterEqual[idx], m_solInfoBetterEqual[idy]);
		}
		void setSol(const std::vector<std::shared_ptr<SolBase>>& sols,
			const std::vector<int> & solId
			) {
			
			resize(sols.size());
			m_insertId = solId;
			for (int idx(0); idx < sols.size(); ++idx) {
				m_sols[idx] = sols[idx];
			}

			updateSols();
		}
		void updateSol(int solId);
		void updateSolsFromTo(int from, int to) {
			for (int idx(from); idx < to; ++idx) {
				updateSol(idx);
			}
		}
		void updateSols();

		void initNetwork(std::vector<Node>& network) {
			network.resize(m_sols.size());
			//double maxVal = std::numeric_limits<double>::max();

			for (int idx(0); idx < network.size(); ++idx) {
				network[idx].m_dis2parent = 1.0;
				network[idx].m_belong = idx;
			}
		}

		void initNBN() {
			initNetwork(m_solInfoBetter);
			initNetwork(m_solInfoBetterEqual);

		}

		void udpateNBNmultithread();

		void mergeNBNmultithread(
			int from, int to,
			const std::vector<std::vector<Node>>& networkBetter,
			const std::vector<std::vector<Node>>& networkBetterEqual
		);

		void updateNetWorkLoop(int numLoop, std::vector<Node>& networkBetter, std::vector<Node>& networkBetterEqual);
		void updateNetwork(std::vector<Node>& networkBetter, std::vector<Node>& networkBetterEqual);
		int updateDivision(const std::vector<int>& totalSol,
			std::vector<int>& sortedEdgeId,
			int sortedEdgeIdFrom,
			const std::vector<int>& solDirection,
			std::vector<Node>& networkBetter, std::vector<Node>& networkBetterEqual
		);
		//void betterSwap(int& idx, int& idy);
		int updateSolNetworkInfo(int idx, int belong, double curDis,
			std::vector<Node>& networkBetter);
		int updateRelation(int idx, int idy, std::vector<Node>& networkBetter, std::vector<Node>& networkBetterEqual);


		// Í¨¹ý NBN_DivisionBase ¼Ì³Ð
		virtual void initialize_(bool flag_grid_sample = true) override {

		}
		virtual size_t size() const override {
			return m_sols.size();
		}
		virtual void addSol(const SolBase& new_sol, int& belong_id,
			bool flag_opt = false, int popIter = -1, int popSolId = -1, int algId = -1) override {

		}
		virtual void getSharedNBN(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			std::vector<bool>& flagOpt) const override {

		}
		virtual void getSharedNBN(std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<double>& fitness, 
			std::vector<int>& belong, 
			std::vector<double>& dis2parent, 
			std::vector<int> popIters, std::vector<int> popSolIds, std::vector<int> algIds) const override {

		}


		void getResult(std::vector<int>& belong, std::vector<double>& dis2parent) {
			belong.resize(m_solInfoBetter.size());
			dis2parent.resize(m_solInfoBetter.size());
			for (int idx(0); idx < m_solInfoBetter.size(); ++idx) {
				belong[idx] = m_solInfoBetter[idx].m_belong;
				dis2parent[idx] = m_solInfoBetter[idx].m_dis2parent;
			}

			
		}


		void getResultEqualBetter(std::vector<int>& belong, std::vector<double>& dis2parent) {
			belong.resize(m_solInfoBetterEqual.size());
			dis2parent.resize(m_solInfoBetterEqual.size());
			for (int idx(0); idx < m_solInfoBetterEqual.size(); ++idx) {
				belong[idx] = m_solInfoBetterEqual[idx].m_belong;
				dis2parent[idx] = m_solInfoBetterEqual[idx].m_dis2parent;
			}


		}



	};
}


#endif